library(purrr)
library(dplyr)

match_worms <- function(tax_names) {
  name_batches <- split(tax_names, as.integer((seq_along(tax_names) - 1) / 50))
  get_matches <- insistently(function(batch) {
    res <- worrms::wm_records_names(batch, marine_only = FALSE, on_error = warning)
    if (is.null(res)) {
      # worrms returns NULL if any other error than 204
      stop("Error trying to match names")
    }
    res
  }, quiet = FALSE)
  matched_batches <- purrr::map(name_batches, function(batch) {
    res <- get_matches(batch)
    names(res) <- batch
    df <- bind_rows(res, .id = "input")
    if (nrow(df) > 0) {
      df <- df %>%
        filter(match_type == "exact" | match_type == "exact_genus" | match_type == "exact_subgenus") %>%
        filter(!is.na(scientificname)) %>%
        mutate(priority = ifelse(status == "accepted", 1, 0)) %>%
        group_by(input) %>%
        arrange(desc(priority)) %>%
        filter(row_number() == 1)
    }
    df
  })
  matches <- bind_rows(matched_batches)
  return(matches)
}
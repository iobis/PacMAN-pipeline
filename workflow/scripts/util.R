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

# # Clean set of taxon names into taxonomy as a named list
# clean_taxonomy <- function(taxa, prefixed, ranks) {
#   if (prefixed) {
#     taxa <- taxa[taxa != "" & taxa != "NA" & taxa != "nan"]
#     parts <- str_match(taxa, "([a-z]+)__(.*)")
#     recoded_ranks <- recode(parts[,2], "sk" = "superkingdom", "k" = "kingdom", "p" = "phylum", "c" = "class", "o" = "order", "f" = "family", "g" = "genus", "s" = "species")
#     taxon_names <- as.list(parts[,3])
#     names(taxon_names) <- recoded_ranks
#     exported_ranks <- intersect(recoded_ranks, ranks)
#     if (length(exported_ranks) == 0) {
#       return(setNames(list(NA), ranks[1]))
#     }
#     return(taxon_names[exported_ranks])
#   } else {
#     if (length(taxa) == 0) {
#       return(setNames(list(NA), ranks[1]))
#     }
#     taxa[taxa %in% c("", "NA", "nan", "unknown", "Unknown")] <- NA
#     taxa[grepl("uncultured", taxa, ignore.case = TRUE)] <- NA
#     taxa[grepl("sp\\.", taxa, ignore.case = TRUE)] <- NA
#     taxon_names <- setNames(as.list(taxa), ranks[1:length(taxa)])
#     return(taxon_names)
#   }
# }

# taxonomy_to_taxmat <- function(taxonomy) {
#   # Make educated guess about taxonomy format and clean taxonomies
#   # TODO: fix, this logic will only work for some reference databases

#   taxonomies <- str_split(str_replace(taxonomy, ";+$", ""), ";")

#   max_taxonomy_length <- max(sapply(taxonomies, length))
#   most_frequent_names <- names(head(sort(table(unlist(taxonomies)), decreasing = TRUE, na.last = TRUE)))
#   most_frequent_names <- most_frequent_names[most_frequent_names != ""]
#   frequent_names_prefixed <- all(str_detect(most_frequent_names, "([a-z]+)__(.*)"))

#   if (frequent_names_prefixed) {
#     prefixed <- TRUE
#   } else {
#     prefixed <- FALSE
#   }

#   if (max_taxonomy_length == 8 | "Metazoa" %in% most_frequent_names) {
#     ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
#   } else if (max_taxonomy_length == 7) {
#     ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
#   }

#   cleaned <- lapply(taxonomies, clean_taxonomy, prefixed = prefixed, ranks = ranks)

#   taxmat <- cleaned %>%
#     bind_rows() %>%
#     as.data.frame() %>%
#     select(!!!ranks)

#   return(list(taxmat = taxmat, ranks = ranks))
# }

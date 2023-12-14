# Written by Saara Suominen and Pieter Provoost (saara.suominen.work@gmail.com)
# for OBIS and PacMAN

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)
library(tidyr)
library(yaml)
library(purrr)

REF_DB_PATTERN <- "identity_filtered/\\s*(.*?)\\s*_blca_tax_table"

# Parse arguments
args <- commandArgs(trailingOnly = T)
message("args <- ", capture.output(dput(args))) # output for debugging

outpath <- args[1]
tax_file_path <- args[2]
rep_seqs_path <- args[3]
config_path <- args[4]

if (length(args) == 6) {
    basta_file_path <- args[5]
    blast_date <- args[6]
} else {
    basta_file_path <- NULL
}

########################### 1. Read input files  ################################################################################################################

tax_file <- read.csv(tax_file_path, sep = "\t", header = T)

rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

config <- read_yaml(config_path)

# If blast was performed on the unknown sequences:
if (!is.null(basta_file_path)) {
  if (file.size(basta_file_path) > 0) {
    message("0. Results of Blast annotation read")
    basta_file <- read.csv(basta_file_path, sep = "\t", header = F)
    colnames(basta_file) <- c("asv", "sum.taxonomy")
    #basta_file$identificationRemarks <- basta_file$sum.taxonomy
  }
}

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) #####################################

message("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

# Clean set of taxon names into taxonomy as a named list
clean_taxonomy <- function(taxa, prefixed, ranks) {
  if (prefixed) {
    taxa <- taxa[taxa != "" & taxa != "NA" & taxa != "nan"]
    parts <- str_match(taxa, "([a-z]+)__(.*)")
    recoded_ranks <- recode(parts[,2], "sk" = "superkingdom", "k" = "kingdom", "p" = "phylum", "c" = "class", "o" = "order", "f" = "family", "g" = "genus", "s" = "species")
    taxon_names <- as.list(parts[,3])
    names(taxon_names) <- recoded_ranks
    exported_ranks <- intersect(recoded_ranks, ranks)
    if (length(exported_ranks) == 0) {
      return(setNames(list(NA), ranks[1]))
    }
    return(taxon_names[exported_ranks])
  } else {
    if (length(taxa) == 0) {
      return(setNames(list(NA), ranks[1]))
    }
    taxa[taxa %in% c("", "NA", "nan", "unknown", "Unknown")] <- NA
    taxa[grepl("uncultured", taxa, ignore.case = TRUE)] <- NA
    taxa[grepl("sp\\.", taxa, ignore.case = TRUE)] <- NA
    taxon_names <- setNames(as.list(taxa), ranks[1:length(taxa)])
    return(taxon_names)
  }
}

taxonomy_to_taxmat <- function(taxonomy) {
  # Make educated guess about taxonomy format and clean taxonomies
  # TODO: fix, this logic will only work for some reference databases

  taxonomies <- str_split(str_replace(taxonomy, ";+$", ""), ";")

  max_taxonomy_length <- max(sapply(taxonomies, length))
  most_frequent_names <- names(head(sort(table(unlist(taxonomies)), decreasing = TRUE, na.last = TRUE)))
  most_frequent_names <- most_frequent_names[most_frequent_names != ""]
  frequent_names_prefixed <- all(str_detect(most_frequent_names, "([a-z]+)__(.*)"))

  if (frequent_names_prefixed) {
    prefixed <- TRUE
  } else {
    prefixed <- FALSE
  }

  if (max_taxonomy_length == 8 | "Metazoa" %in% most_frequent_names) {
    ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  } else if (max_taxonomy_length == 7) {
    ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  }

  cleaned <- lapply(taxonomies, clean_taxonomy, prefixed = prefixed, ranks = ranks)

  taxmat <- cleaned %>%
    bind_rows() %>%
    as.data.frame() %>%
    select(!!!ranks)

  return(list(taxmat = taxmat, ranks = ranks))
}

taxmat_result <- taxonomy_to_taxmat(tax_file$sum.taxonomy)
ranks <- taxmat_result$ranks
taxmat <- taxmat_result$taxmat %>%
    mutate(verbatimIdentification = tax_file$identificationRemarks)
row.names(taxmat) <- tax_file$asv
row.names(tax_file) <- tax_file$asv

# Clean separately basta file and add to taxmat together

if (exists("basta_file")) {
    taxmat_result_basta <- taxonomy_to_taxmat(basta_file$sum.taxonomy)
    ranks_basta <- taxmat_result_basta$ranks
    taxmat_basta <- taxmat_result_basta$taxmat %>%
      mutate(verbatimIdentification = paste("Identification based on blastn against the full nt database (downloaded on", blast_date, "), and with basta-lca with filtering on:", config$BLAST$percent_identity, "percent identity,", config$BLAST$"e-value", "e-value, and", config$BLAST$alignment_length, "alignment length"))
    row.names(taxmat_basta) <- basta_file$asv
    row.names(basta_file) <- basta_file$asv

    # Remove ASVs present in basta file from tax file
    taxmat <- taxmat[!rownames(taxmat) %in% rownames(taxmat_basta),]

    # Merge
    taxmat <- bind_rows(taxmat, taxmat_basta)
    ranks <- unique(c(ranks, ranks_basta))
}

rep_seqs_unknown <- names(rep_seqs[!names(rep_seqs) %in% row.names(taxmat),])
taxmat_unknown <- data.frame(
  otu_seq_comp_appr = rep(NA, length(rep_seqs_unknown)),
  otu_db = rep(NA, length(rep_seqs_unknown))
)
row.names(taxmat_unknown) <- rep_seqs_unknown
taxmat <- bind_rows(taxmat, taxmat_unknown)

# Match distinct names

tax_names <- taxmat %>% select(!!!ranks) %>% unlist() %>% na.omit() %>% unique() %>% sort()

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
      mutate(priority = ifelse(status == "accepted", 1, 0)) %>%
      group_by(input) %>%
      arrange(desc(priority)) %>%
      filter(row_number() == 1)
  }
  df
})
matches <- bind_rows(matched_batches) %>%
  split(f = .$input)

message("2. Matched names across all ranks")

taxmat$scientificName <- NA
taxmat$scientificNameID <- NA
taxmat$taxonRank <- NA

for (i in 1:nrow(taxmat)) {
  
  lsids <- taxmat[i, ranks] %>%
    as.character() %>%
    sapply(function(x) { matches[[x]]$lsid }) %>%
    sapply(function(x) { ifelse(is.null(x), NA, x) }) %>%
    unlist()
  if (all(is.na(lsids))) next

  most_specific_name <- taxmat[i, max(which(!is.na(lsids)))]
  best_match <- as.data.frame(matches[[most_specific_name]])
  taxmat$scientificName[i] <- best_match$scientificname
  taxmat$scientificNameID[i] <- best_match$lsid
  taxmat$taxonRank[i] <- tolower(best_match$rank)
  for (rank in ranks) {
    if (rank == "species" & best_match$rank == "Species") {
      taxmat[i, rank] <- best_match$scientificname
    } else if (rank %in% names(best_match)) {
      taxmat[i, rank] <- best_match[,rank]
    } else {
      taxmat[i, rank] <- NA
    }
  }
}

# Add Incertae LSID in case there is no last value
taxmat$scientificName[is.na(taxmat$scientificName)] <- "Incertae sedis"
taxmat$taxonRank[is.na(taxmat$scientificNameID)] <- "kingdom"
taxmat$scientificNameID[is.na(taxmat$scientificNameID)] <- "urn:lsid:marinespecies.org:taxname:12"

# Names not in WoRMS
names_not_in_worms <- names(matches)[sapply(matches, is.null)]
message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))

# Add sequence to the tax_table slot (linked to each asv)
message("Add sequences to table")
taxmat$DNA_sequence <- as.character(rep_seqs[row.names(taxmat)])

# Add identificationRemarks, containing also vsearch annotation results
#taxmat$identificationRemarks <- tax_file[row.names(taxmat), "identificationRemarks"]
message(paste("The taxonomy file contains the following fields:", colnames(tax_file), head(tax_file)))

# Write table of unknown names to make manual inspection easier:
message("Write names not recognized in worrms to separate file")
write.table(names_not_in_worms, paste0(outpath, "Taxa_not_in_worms.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# Write tax table
message("Write tax table")
write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

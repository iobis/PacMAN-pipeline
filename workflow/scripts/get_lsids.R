# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)
library(tidyr)

REF_DB_PATTERN <- "identity_filtered/\\s*(.*?)\\s*_blca_tax_table"

# Parse arguments
args <- commandArgs(trailingOnly = T)
message("args <- ", capture.output(dput(args))) # output for debugging

outpath <- args[1]
tax_file_path <- args[2]
rep_seqs_path <- args[3]
if (length(args) == 5) {
    basta_file_path <- args[4]
    blast_date <- args[5]
} else {
    basta_file_path <- NULL
}

########################### 1. Read input files  ################################################################################################################

tax_file <- read.csv(tax_file_path, sep = "\t", header = T)

rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# If blast was performed on the unknown sequences:
if (!is.null(basta_file_path)) {
  if (file.size(basta_file_path) > 0) {
    message("0. Results of Blast annotation read")
    basta_file <- read.csv(basta_file_path, sep = "\t", header = F)
    colnames(basta_file) <- c("asv", "sum.taxonomy")
    basta_file$verbatimIdentification <- basta_file$sum.taxonomy
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
    taxon_names <- setNames(as.list(taxa), ranks[1:length(taxa)])
    return(taxon_names)
  }
}

if (exists("basta_file")) {
  # Remove ASVs present in basta file from tax file
  tax_file <- tax_file[!tax_file$asv %in% basta_file$asv,]
  # Merge
  tax_file <- bind_rows(tax_file, basta_file)
}

taxonomies <- str_split(str_replace(tax_file$sum.taxonomy, ";+$", ""), ";")

# Make educated guess about taxonomy format and clean taxonomies
# TODO: fix, this logic will only work for some reference databases

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

# Create output table

taxmat <- cleaned %>%
  bind_rows() %>%
  as.data.frame() %>%
  select(!!!ranks) %>%
  mutate(verbatimIdentification = tax_file$verbatimIdentification)

row.names(taxmat) <- tax_file$asv
row.names(tax_file) <- tax_file$asv

rep_seqs_unknown <- names(rep_seqs[!names(rep_seqs) %in% row.names(taxmat),])
taxmat_unknown <- data.frame(
  otu_seq_comp_appr = rep(NA, length(rep_seqs_unknown)),
  otu_db = rep(NA, length(rep_seqs_unknown))
)
row.names(taxmat_unknown) <- rep_seqs_unknown
taxmat <- bind_rows(taxmat, taxmat_unknown)

# TODO: failing with rate limit, submit in batches with at least version 0.4.3 of worrms (https://anaconda.org/conda-forge/r-worrms)
match_name <- function(name) {
  lsid <- tryCatch({
    res <- wm_records_names(name, marine_only = FALSE)
    # TODO: fix
    Sys.sleep(1)
    matches <- res[[1]] %>%
      filter(match_type == "exact" | match_type == "exact_genus" | match_type == "exact_subgenus")
    if (nrow(matches) > 1) {
      message(paste0("Multiple matches for ", name))
    }
    return(matches[1,])
  }, error = function(cond) {
    message(cond)
    return(NULL)
  })
}

# Taxon names across all ranks
tax_names <- taxmat %>% select(!!!ranks) %>% unlist() %>% na.omit() %>% unique() %>% sort()
matches <- sapply(tax_names, match_name)

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
  scientificnameid <- matches[[most_specific_name]]$lsid

  taxmat$scientificName[i] <- matches[[most_specific_name]]$scientificname
  taxmat$scientificNameID[i] <- matches[[most_specific_name]]$lsid
  taxmat$taxonRank[i] <- tolower(matches[[most_specific_name]]$rank)
  taxmat$kingdom[i] <- matches[[most_specific_name]]$kingdom
  taxmat$phylum[i] <- matches[[most_specific_name]]$phylum
  taxmat$class[i] <- matches[[most_specific_name]]$class
  taxmat$order[i] <- matches[[most_specific_name]]$order
  taxmat$family[i] <- matches[[most_specific_name]]$family
  taxmat$genus[i] <- matches[[most_specific_name]]$genus
}

# Add Biota LSID in case there is no last value
# Kingdom is used as taxonRank so that "Biota" is also recognized correctly by GBIF
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
taxmat$identificationRemarks <- tax_file[row.names(taxmat), "identificationRemarks"]
message(paste("The taxonomy file contains the following fields:", colnames(tax_file), head(tax_file)))

# Write table of unknown names to make manual inspection easier:
message("Write names not recognized in worrms to separate file")
write.table(names_not_in_worms, paste0(outpath, "Taxa_not_in_worms.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# Write tax table
message("Write tax table")
write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

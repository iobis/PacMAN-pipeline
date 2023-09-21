# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)
library(tidyr)

RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
PIPELINE <- "bowtie2;2.4.4;ANACAPA-blca;2021"
BLAST_VERSION <- "blastn;2.12.0"
BLAST_DB <- "NCBI-nt;"
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

tax_file <- read.csv(tax_file_path, sep = "\t", header = T) %>%
  rename("verbatimIdentification" = "taxonomy_confidence")
# Add annotation pipeline and reference database
tax_file$otu_seq_comp_appr <- PIPELINE
result <- regmatches(tax_file_path, regexec(REF_DB_PATTERN, tax_file_path))
otu_db <- result[[1]][2]
tax_file$otu_db <- otu_db

rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# If blast was performed on the unknown sequences:
if (!is.null(basta_file_path)) {
  if (file.size(basta_file_path) > 0) {
    message("0. Results of Blast annotation read")
    basta_file <- read.csv(basta_file_path, sep = "\t", header = F)
    colnames(basta_file) <- c("rowname", "sum.taxonomy")
    basta_file$verbatimIdentification <- basta_file$sum.taxonomy
    # Add annotation pipeline and reference database
    basta_file$otu_seq_comp_appr <- BLAST_VERSION
    basta_file$otu_db <- paste0(BLAST_DB, blast_date)
  }
}

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) #####################################

message("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

# Clean set of taxon names into taxonomy as a named list
clean_taxonomy <- function(taxa) {
  if (all(str_detect(taxa, "([a-z]+)__(.*)_[0-9]"))) {
    taxa <- taxa[taxa != "" & taxa != "NA" & taxa != "nan"]
    if (length(taxa) == 0) {
      return(list(kingdom = NA))
    }
    parts <- str_match(taxa, "([a-z]+)__(.*)_[0-9]")
    ranks <- recode(parts[,2], "k" = "kingdom", "p" = "phylum", "c" = "class", "o" = "order", "f" = "family", "g" = "genus", "s" = "species")
    taxon_names <- as.list(parts[,3])
    names(taxon_names) <- ranks
    return(taxon_names)
  } else {
    if (length(taxa) == 0) {
      return(list(kingdom = NA))
    }
    taxa[taxa %in% c("", "NA", "nan", "unknown", "Unknown")] <- NA
    taxon_names <- setNames(as.list(taxa), RANKS[1:length(taxa)])
    return(taxon_names)
  }
}

if (exists("basta_file")) {
  # Remove ASVs present in basta file from tax file
  tax_file <- tax_file %>% filter(!rowname %in% basta_file$rowname)
  # Merge
  tax_file <- bind_rows(tax_file, basta_file)
}

taxonomies <- str_split(str_replace(tax_file$sum.taxonomy, ";+$", ""), ";")
cleaned <- lapply(taxonomies, clean_taxonomy)

taxmat <- cleaned %>%
  bind_rows() %>%
  as.data.frame() %>%
  select(!!!RANKS) %>%
  mutate(verbatimIdentification = tax_file$verbatimIdentification)%>%
  mutate(otu_seq_comp_appr = tax_file$otu_seq_comp_appr)%>%
  mutate(otu_db = tax_file$otu_db)

row.names(taxmat) <- tax_file$rowname

# Add possible remaining unknowns to the taxmat based on asvs in the rep_seqs (keep all ASVs in the final dataset)
rep_seqs_unknown <- names(rep_seqs[!names(rep_seqs)%in%row.names(taxmat),])
taxmat_unknown <- data.frame(
  otu_seq_comp_appr = rep(PIPELINE, length(rep_seqs_unknown)),
  otu_db = rep(otu_db, length(rep_seqs_unknown))
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
tax_names <- taxmat %>% select(!!!RANKS) %>% unlist() %>% na.omit() %>% unique() %>% sort()
matches <- sapply(tax_names, match_name)

taxmat$scientificName <- NA
taxmat$scientificNameID <- NA

for (i in 1:nrow(taxmat)) {

  lsids <- taxmat[i, RANKS] %>%
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

# Add Incertae sedis LSID in case there is no last value
taxmat$scientificName[is.na(taxmat$scientificName)] <- "Incertae sedis"
taxmat$scientificNameID[is.na(taxmat$scientificNameID)] <- "urn:lsid:marinespecies.org:taxname:12"

# Names not in WoRMS
names_not_in_worms <- names(matches)[sapply(matches, is.null)]
message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))

# Add sequence to the tax_table slot (linked to each asv)
taxmat$DNA_sequence <- as.character(rep_seqs[row.names(taxmat)])

# Write table of unknown names to make manual inspection easier:
write.table(names_not_in_worms, paste0(outpath, "Taxa_not_in_worms.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# Write tax table
write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

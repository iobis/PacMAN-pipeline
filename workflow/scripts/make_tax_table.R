# Written by Saara Suominen and Pieter Provoost (saara.suominen.work@gmail.com)
# for OBIS and PacMAN

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)
library(tidyr)
library(yaml)
library(purrr)
library(glue)
library(jsonlite)
source("workflow/scripts/util.R")

# Parse arguments

if (!exists("cmd_args")) {
  cmd_args <- commandArgs(trailingOnly = T)
  message("cmd_args <- ", capture.output(dput(cmd_args)))
}

outpath <- cmd_args[1]
tax_file_path <- cmd_args[2]
rep_seqs_path <- cmd_args[3]
config_path <- cmd_args[4]
config <- read_yaml(config_path)

param_rdp_confidence_threshold <- config$Rdp$cutoff
param_vsearch_identity_threshold <- config$Vsearch$pident
param_vsearch_cover_threshold <- config$Vsearch$query_cov
param_basta_identity_threshold <- ifelse(config$BLAST$include, config$BLAST$percent_identity, NA)

# Read input files

tax_file <- read.csv(tax_file_path, sep = "\t", header = T) %>%
      mutate(across(where(is.character), ~na_if(.x, "")))
rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# Combine annotations into identificationRemarks

remarks <- tax_file %>%
  group_by(asv) %>%
  summarize(identificationRemarks = as.character(
    toJSON(
      list(
        pipeline = "PacMAN pipeline",
        parameters = list(
          rdp_confidence_threshold = param_rdp_confidence_threshold,
          vsearch_identity_threshold = param_vsearch_identity_threshold,
          vsearch_cover_threshold = param_vsearch_identity_threshold,
          basta_identity_threshold = param_basta_identity_threshold
        ),
        annotations = tibble(method, scientificName, scientificNameID, accession = seqid, confidence = confidence, identity, query_cover)
      ),
      dataframe = "rows", auto_unbox = TRUE
    )
  )) %>%
  ungroup()
  
# Construct taxonomy table (in one case superkingdom not found.)

tax <- tax_file %>%
  filter(method == "RDP classifier") %>%
  select(any_of(c("asv", "scientificName", "scientificNameID", "domain", "superkingdom", "phylum", "class", "order", "family", "genus", "species"))) %>%
  left_join(remarks, by = "asv") %>%
  rowwise() %>%
  mutate(
    taxonRank = case_when(
      scientificName == domain ~ "domain",
      #scientificName == superkingdom ~ "superkingdom",
      scientificName == coalesce(cur_data()[["superkingdom"]], NA_character_) ~ "superkingdom",
      scientificName == phylum ~ "phylum",
      scientificName == class ~ "class",
      scientificName == order ~ "order",
      scientificName == family ~ "family",
      scientificName == genus ~ "genus",
      scientificName == species ~ "species",
      TRUE ~ NA
    )
  )

# Handle unidentified ASVs

unidentified_asvs <- names(rep_seqs)[which(!names(rep_seqs) %in% tax$asv)]

unidentified <- data.frame(
  asv = unidentified_asvs,
  scientificName = "Incertae sedis",
  scientificNameID = "urn:lsid:marinespecies.org:taxname:12"
)

tax <- bind_rows(tax, unidentified)

# Add sequence

message("Add sequences to table")
tax$DNA_sequence <- as.character(rep_seqs[tax$asv])

# Names not in WoRMS
# names_not_in_worms <- names(matches)[sapply(matches, is.null)]
# message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))
# write.table(names_not_in_worms, paste0(outpath, "taxa_not_in_worms.tsv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# TODO: fix
write.table(data.frame(scientificName = character(0)), paste0(outpath, "taxa_not_in_worms.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

message(paste("The taxonomy file contains the following fields:", colnames(tax), head(tax)))

# Write tax table

message("Write tax table")
write.table(tax, paste0(outpath, "full_tax_table_with_lsids.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

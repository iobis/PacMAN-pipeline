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

# Read input files

tax_file <- read.csv(tax_file_path, sep = "\t", header = T) %>%
      mutate(across(where(is.character), ~na_if(.x, "")))
rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# Combine annotations into identificationRemarks

remarks <- tax_file %>%
  group_by(asv) %>%
  summarize(identificationRemarks = as.character(toJSON(tibble(method, scientificName, scientificNameID, accession = seqid, rdp_confidence = confidence, identity, rdp_confidence_threshold, query_cover, vsearch_identity_threshold, vsearch_cover_threshold, basta_identity_threshold), dataframe = "rows", auto_unbox = TRUE))) %>%
  ungroup()
  
# Construct taxonomy table

tax <- tax_file %>%
  filter(method == "RDP classifier") %>%
  select(asv, scientificName, scientificNameID, domain, superkingdom, phylum, class, order, family, genus) %>%
  left_join(remarks, by = "asv")

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
write.table(data.frame(), paste0(outpath, "taxa_not_in_worms.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

message(paste("The taxonomy file contains the following fields:", colnames(tax), head(tax)))

# Write tax table

message("Write tax table")
write.table(tax, paste0(outpath, "full_tax_table_with_lsids.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

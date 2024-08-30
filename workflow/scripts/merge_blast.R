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
source("workflow/scripts/util.R")

REF_DB_PATTERN <- "identity_filtered/\\s*(.*?)\\s*_blca_tax_table"

# Parse arguments

if (!exists("cmd_args")) {
  cmd_args <- commandArgs(trailingOnly = T)
  message("cmd_args <- ", capture.output(dput(cmd_args)))
}

outpath <- cmd_args[1]
annotations_file_path <- cmd_args[2]
rep_seqs_path <- cmd_args[3]
config_path <- cmd_args[4]
config <- read_yaml(config_path)

basta_identity_threshold <- config$BLAST$percent_identity

if (length(cmd_args) == 6) {
    basta_file_path <- cmd_args[5]
    blast_date <- cmd_args[6]
} else {
    basta_file_path <- NULL
}

# Read input files

config <- read_yaml(config_path)
annotations_file <- read.csv(annotations_file_path, sep = "\t", header = T)
rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# Clean basta file

if (!is.null(basta_file_path)) {
  if (file.size(basta_file_path) > 0) {
    message("0. Results of Blast annotation read")
    basta_file <- read.csv(basta_file_path, sep = "\t", header = T) #%>%
      #setNames(c("asv", "sum.taxonomy"))
  }
}

if (exists("basta_file")) {

    basta_clean <- basta_file %>%
      dplyr::rename(taxonomy = Taxonomy, asv=ASV) %>%
      mutate(
        taxonomy = str_replace_all(taxonomy, "_", " ")
      ) %>%
      mutate(
        taxonomy = str_replace_all(taxonomy, "[Uu]nknown(;|$)", ";"),
      ) %>%
      separate(taxonomy, into = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
      mutate(across(where(is.character), ~na_if(.x, ""))) %>%
      rowwise() %>%
      mutate(scientificName = coalesce(species, genus, family, order, class, phylum)) %>%
      ungroup() %>%
      filter(!is.na(scientificName) & scientificName != "")

    tax_names <- basta_clean$scientificName %>% na.omit() %>% unique() %>% sort()
    matches <- match_worms(unique(tax_names)) %>%
      select(scientificName = input, scientificNameID = lsid)
    basta_clean <- basta_clean %>%
    left_join(matches, by = "scientificName") %>%
      mutate(
        method = "BASTA",
        basta_identity_threshold = basta_identity_threshold
      )

    annotations_file <- annotations_file %>%
      bind_rows(basta_clean)

}

write.table(annotations_file, paste0(outpath, "annotation_results_blast.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t", na = "")

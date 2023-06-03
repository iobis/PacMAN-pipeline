# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN
#Last update 21.6.2021

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = T)

outpath <- args[1]
tax_file_path <- args[2]
rep_seqs_path <- args[3]
if (length(args) == 5) {
    basta_file_path <- args[4]
    blast_date <- args[5]
} else {
    basta_file_path <- NULL
}

RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

########################### 1. Read input files  ################################################################################################################

tax_file <- read.csv(tax_file_path, sep = "\t", header = T)
rep_seqs <- Biostrings::readDNAStringSet(rep_seqs_path)

# If blast was performed on the unknown sequences:
if (!is.null(basta_file_path)) {
  if (file.size(basta_file_path) > 0) {
    message("0. Results of Blast annotation read")
    basta_file <- read.csv(basta_file_path, sep = "\t", header = F)
  }
}

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) #####################################

message("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

clean_taxonomy_prefixed <- function(taxa) {
  taxa <- taxa[taxa != "" & taxa != "NA" & taxa != "nan"]
  if (length(taxa) == 0) {
    return(list(kingdom = NA))
  }
  parts <- str_match(taxa, "([a-z]+)__(.*)_[0-9]")
  ranks <- recode(parts[,2], "k" = "kingdom", "p" = "phylum", "c" = "class", "o" = "order", "f" = "family", "g" = "genus", "s" = "species")
  taxon_names <- as.list(parts[,3])
  names(taxon_names) <- ranks
  return(taxon_names)
}

clean_taxonomy <- function(taxa) {
  if (length(taxa) == 0) {
    return(list(kingdom = NA))
  }
  taxa[taxa %in% c("", "NA", "nan")] <- NA
  taxon_names <- setNames(as.list(taxa), RANKS[1:length(taxa)])
  return(taxon_names)
}

taxonomies <- str_split(tax_file$sum.taxonomy, ";")

if (grepl("MIDORI_UNIQ", tax_file_path, fixed = TRUE, ignore.case = TRUE)) {
  # Move names to the appropriate columns if using the MIDORI database
  cleaned <- lapply(taxonomies, clean_taxonomy_prefixed)
} else {
  cleaned <- lapply(taxonomies, clean_taxonomy)
}

taxmat <- cleaned %>%
  bind_rows() %>%
  as.data.frame() %>%
  select(!!!RANKS) %>%
  mutate(identificationRemarks = tax_file$sum.taxonomy)

row.names(taxmat) <- tax_file$rowname

# TODO: remove
# Collect the highest known taxonomic value to the last column
# taxmat$lastvalue <- as.matrix(taxmat)[cbind(seq(1, nrow(taxmat)), max.col(!is.na(taxmat), "last"))]

# Add information on the classification approach
# TODO: remove hard coding
taxmat$otu_seq_comp_appr <- "bowtie2;2.4.4;ANACAPA-blca;2021"

# Get otu_db name from the input filename
pattern <- "identity_filtered/\\s*(.*?)\\s*_blca_tax_table"
result <- regmatches(tax_file_path, regexec(pattern, tax_file_path))
taxmat$otu_db <- result[[1]][2]

# If Blast was performed on unknown sequences and gave results after lca
if (exists("basta_file")) {
  message("1.1 Adding Blast results to taxonomic table")
  colnames(basta_file) <- c("rowname", "sum.taxonomy")
  taxmat2 <- separate(basta_file, "sum.taxonomy", into = RANKS, sep = ";")
  taxmat2[taxmat2 == ""] <- NA
  taxmat2$species <- gsub("_", " ", taxmat2$species)
  taxmat2$lastvalue <- as.matrix(taxmat2)[cbind(seq(1, nrow(taxmat2)), max.col(!is.na(taxmat2), "last"))]
  # Add also the fields that separate the otu assignment methods
  # Note: read these from the config file?
  # TODO: remove hard coding
  taxmat2$otu_seq_comp_appr <- "blastn;2.12.0"
  taxmat2$otu_db <- paste0("NCBI-nt;", blast_date)
  rownames(taxmat2) <- taxmat2$rowname
  taxmat2 <- subset(taxmat2, select = -c(rowname))

  # Combine this taxmat to the original one
  # At the moment the unknowns are there twice! So first remove duplicates before you combine them
  taxmat <- taxmat[!(rownames(taxmat) %in% rownames(taxmat2)),]
  taxmat <- rbind(taxmat, taxmat2)
}

# Add possible remaining unknowns to the taxmat based on asvs in the rep_seqs (keep all ASVs in the final dataset)
rep_seqs_unknown <- names(rep_seqs[!names(rep_seqs) %in% row.names(taxmat),])
rn <- row.names(taxmat)
taxmat[nrow(taxmat) + seq_along(rep_seqs_unknown), ] <- NA 
row.names(taxmat) <- c(rn, rep_seqs_unknown)
# Fill the original database values on the otu_seq_comp_appr and otu_db for the unknown sequences as well
taxmat[is.na(taxmat$otu_seq_comp_appr), "otu_seq_comp_appr"] <- "bowtie2;2.4.4;ANACAPA-blca;2021"
taxmat[is.na(taxmat$otu_db), "otu_db"] <- taxmat[1,"otu_db"]

# TODO: failing with rate limit, submit in batches

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
lsid_map <- lapply(tax_names, function(x) {
  list(lsid = matches[[x]]$lsid, scientificname =  matches[[x]]$scientificname, rank =  matches[[x]]$rank)
})
names(lsid_map) <- tax_names

# lookup_lsid <- function(input) {
#   lsids <- sapply(input, function(x) { lsid_map[[x]] })
#   lsids[sapply(lsids, is.null)] <- NA
#   return(unlist(lsids))
# }

# lsid_table <- taxmat %>%
#   select(!!!RANKS) %>%
#   mutate(across(everything(), lookup_lsid))
# best_column_index <- max.col(!is.na(lsid_table), "last")

taxmat$scientificName <- NA
taxmat$scientificNameID <- NA

for (i in 1:nrow(taxmat)) {

  lsids <- taxmat[i, RANKS] %>%
    as.character() %>%
    sapply(function(x) { lsid_map[[x]]$lsid }) %>%
    sapply(function(x) { ifelse(is.null(x), NA, x) }) %>%
    unlist()
  if (all(is.na(lsids))) next

  most_specific_name <- taxmat[i, max(which(!is.na(lsids)))]

  scientificname <- lsid_map[[most_specific_name]]$scientificname
  scientificnameid <- lsid_map[[most_specific_name]]$lsid
  rank <- lsid_map[[most_specific_name]]$rank

  if (!is.na(scientificnameid)) {
      taxmat$scientificName[i] <- scientificname
      taxmat$scientificNameID[i] <- scientificnameid
      taxmat$taxonRank[i] <- tolower(rank)
  }

}

# Add Biota LSID in case there is no last value
# Kingdom is used as taxonRank so that "Biota" is also recognized correctly by GBIF

taxmat$scientificName[is.na(taxmat$scientificName)] <- "Biota"
taxmat$taxonRank[is.na(taxmat$scientificNameID)] <- "kingdom"
taxmat$scientificNameID[is.na(taxmat$scientificNameID)] <- "urn:lsid:marinespecies.org:taxname:1"

# Names not in WoRMS

names_not_in_worms <- names(matches)[sapply(matches, is.null)]
message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))

# TODO: disabled for now, matching already happens at all levels
# Some of these are common contaminants (e.g. Homo sapiens, Canis lupus)
# Some of these are not marine species (e.g. insects), that could be found with terrestrial matches
# Some have a higher taxonomy that can be found in the database
# Check at genus level:
# genus_names <- unique(na.omit(taxmat$genus[is.na(taxmat$lsid)]))
# genus_matches <- sapply(genus_names, match_name)
# missing_lsids <- is.na(taxmat$lsid)
# for (i in 1:nrow(taxmat)) {
#   if (is.na(taxmat$lsid[i]) & !is.na(taxmat$genus[i]) & !is.null(genus_matches[[taxmat$genus[i]]])) {
#     taxmat$lsid[i] <- genus_matches[[taxmat$genus[i]]]$lsid
#     taxmat$taxonRank[i] <- "genus"
#   }
# }
# This way unknowns decreases. However, these values should be checked for possible addition to WoRMS?
# Will require manual inspection
# genera_not_in_worms <- names(genus_matches)[sapply(genus_matches, is.null)]
# message("Number of genus names not recognized in WORMS: ", length(genera_not_in_worms))
# message("These taxa can be found in the table: ", outpath, "Taxa_not_in_worms.csv")
# not_in_worms <- taxmat[is.na(taxmat$lsid),]

# Add sequence to the tax_table slot (linked to each asv)
taxmat$DNA_sequence <- as.character(rep_seqs[row.names(taxmat)])

# Write table of unknown names to make manual inspection easier:

write.table(names_not_in_worms, paste0(outpath, "Taxa_not_in_worms.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# Write tax table
#taxmat <- apply(taxmat,2,as.character)
write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

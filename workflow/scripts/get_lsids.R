# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN
#Last update 21.6.2021

library(worrms)
library(stringr)
library(Biostrings)
library(dplyr)

args <- commandArgs(trailingOnly = T)

outpath <- args[1]

########################### 1. Read input files  ################################################################################################################

tax_file <- read.csv(args[2], sep = "\t", header = T)
rep_seqs <- Biostrings::readDNAStringSet(args[3])

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) ######################################

message("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
taxmat <- str_split(tax_file$sum.taxonomy, ";", simplify = T)
colnames(taxmat) <- levels
rownames(taxmat) <- tax_file$rowname

# First remove extra characters from strings

taxmat <- gsub(".*[__]([^.]+)[_].*", "\\1", taxmat)

# Replace all empty strings with na

taxmat[taxmat == ""] <- NA
taxmat[taxmat == "NA"] <- NA

# Collect the highest known taxonomic value to the last column

taxmat <- as.data.frame(taxmat)
taxmat$lastvalue <- taxmat[cbind(seq(1, nrow(taxmat)), max.col(!is.na(taxmat), "last"))]

# Because WORMS doesn't recognize Eukaryota, change those that have this in the lastvalue to Biota:

taxmat$lastvalue <- recode(taxmat$lastvalue, "Eukaryota" = "Biota")

# Find worrms value with worrms
# Function from Pieter, but I removed the exact match requirement so be careful with this.
# I removed it, because I was getting 'like' matches even for exact matches, but I am not sure why.

tax_names <- unique(na.omit(taxmat$lastvalue))

match_name <- function(name) {
  lsid <- tryCatch({
    res <- wm_records_names(name, marine_only = FALSE)
    matches <- res[[1]] %>%
      filter(match_type == "exact" | match_type == "exact_genus")
    if (nrow(matches) > 1) {
      message(paste0("Multiple matches for ", name))
    }
    return(matches$lsid[1])
  }, error = function(cond) {
    return(NA)
  })
}

matches <- sapply(tax_names, match_name)
taxmat$lsid <- matches[taxmat$lastvalue]

# Add Biota LSID in case there is no last value

taxmat$lsid[is.na(taxmat$lastvalue)] <- "urn:lsid:marinespecies.org:taxname:1"

# Names not in WoRMS

names_not_in_worms <- names(matches)[is.na(matches)]
message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))

# Some of these are common contaminants (e.g. Homo sapiens, Canis lupus)
# Some of these are not marine species (e.g. insects), that could be found with terrestrial matches
# Some have a higher taxonomy that can be found in the database
# Check at genus level:

genus_names <- unique(na.omit(taxmat$genus[is.na(taxmat$lsid)]))
genus_matches <- sapply(genus_names, match_name)
missing_lsids <- is.na(taxmat$lsid)
taxmat$lsid[missing_lsids] <- matches[taxmat$genus[missing_lsids]]

# This way unknowns decreases. However, these values should be checked for possible addition to WoRMS?
# Will require manual inspection

genera_not_in_worms <- names(genus_matches)[is.na(genus_matches)]
message("Number of genus names not recognized in WORMS: ", length(genera_not_in_worms))
message("These taxa can be found in the table: ", outpath, "Taxa_not_in_worms.csv")

not_in_worms <- taxmat[is.na(taxmat$lsid),]

# We have the dilemma of keeping sequences that are completely unknown (as 'biota'),
# while removing sequences that have a known non-marine origin?

# Add 'Biota' as the name for the unknown sequences

taxmat$lsid[is.na(taxmat$lsid)] <- "urn:lsid:marinespecies.org:taxname:1"

#Add sequence to the tax_table slot (linked to each asv)

taxmat$DNA_sequence <- as.character(rep_seqs[rownames(taxmat)])

# Write table of unknown names to make manual inspection easier:

write.table(not_in_worms, paste0(outpath, "Taxa_not_in_worms.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# Write tax table

write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

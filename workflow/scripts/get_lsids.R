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

########################### 1. Read input files  ################################################################################################################

tax_file <- read.csv(args[2], sep = "\t", header = T)
rep_seqs <- Biostrings::readDNAStringSet(args[3])

#If blast was performed on the unknown sequences:
if (length(args) == 5 & file.size(args[4]) > 0) {
  message("0. Results of Blast annotation read")
  basta_file <- read.csv(args[4], sep = "\t", header = F)
}

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) #####################################

message("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

# Move names to the appropriate columns if using the MIDORI database

if (grepl("MIDORI_UNIQ", args[2], fixed = TRUE, ignore.case = TRUE)) {

  taxonomies <- str_split(tax_file$sum.taxonomy, ";")

  clean_taxonomy <- function(taxa) {
    taxa <- taxa[taxa != "" & taxa != "NA"]
    if (length(taxa) == 0) return(list(kingdom = NA))
    parts <- str_match(taxa, "([a-z]+)__(.*)_[0-9]")
    ranks <- recode(parts[,2], "k" = "kingdom", "p" = "phylum", "c" = "class", "o" = "order", "f" = "family", "g" = "genus", "s" = "species")
    taxon_names <- as.list(parts[,3])
    names(taxon_names) <- ranks
    return(taxon_names)
  }

  cleaned <- lapply(taxonomies, clean_taxonomy)
  taxmat <- as.data.frame(bind_rows(cleaned))
  row.names(taxmat) <- tax_file$rowname

} else {

  taxmat <- separate(tax_file, "sum.taxonomy", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")
  rownames(taxmat) <- taxmat$rowname
  taxmat <- subset(taxmat, select = -c(rowname))

}

taxmat <- taxmat %>%
  select(kingdom, phylum, class, order, family, genus, species)

# Collect the highest known taxonomic value to the last column
taxmat$lastvalue <- as.matrix(taxmat)[cbind(seq(1, nrow(taxmat)), max.col(!is.na(taxmat), "last"))]

# Add information on the classification approach
taxmat$otu_seq_comp_appr <- "bowtie2;2.4.4;ANACAPA-blca;2021"

# Get otu_db name from the input filename
pattern <- "identity_filtered/\\s*(.*?)\\s*_blca_tax_table"
result <- regmatches(args[2], regexec(pattern, args[2]))
taxmat$otu_db <- result[[1]][2]

# If Blast was performed on unknown sequences and gave results after lca
if (length(args) == 5 & file.size(args[4]) > 0) {
  message("1.1 Adding Blast results to taxonomic table")
  colnames(basta_file) <- c("rowname", "sum.taxonomy")
  taxmat2 <- separate(basta_file, "sum.taxonomy", into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")
  taxmat2[taxmat2 == ""] <- NA
  taxmat2$species <- gsub("_", " ", taxmat2$species)
  taxmat2$lastvalue <- as.matrix(taxmat2)[cbind(seq(1, nrow(taxmat2)), max.col(!is.na(taxmat2), "last"))]
  # Add also the fields that separate the otu assignment methods
  # Note: read these from the config file?
  taxmat2$otu_seq_comp_appr <- "blastn;2.12.0"
  taxmat2$otu_db <- paste0("NCBI-nt;", args[5])
  rownames(taxmat2) <- taxmat2$rowname
  taxmat2 <- subset(taxmat2, select = -c(rowname))
  
  # Combine this taxmat to the original one
  # At the moment the unknowns are there twice! So first remove duplicates before you combine them
  taxmat <- taxmat[!(rownames(taxmat) %in% rownames(taxmat2)),]
  taxmat <- rbind(taxmat, taxmat2)

}

# Add possible remaining unknowns to the taxmat based on asvs in the rep_seqs (keep all ASVs in the final dataset)
rep_seqs_unknown <- names(rep_seqs[!names(rep_seqs)%in%row.names(taxmat),])
rn <- row.names(taxmat)
taxmat[nrow(taxmat) + seq_along(rep_seqs_unknown), ] <- NA 
row.names(taxmat) <- c(rn, rep_seqs_unknown)
#Fill the original database values on the otu_seq_comp_appr and otu_db for the unknown sequences as well
taxmat[is.na(taxmat$otu_seq_comp_appr), "otu_seq_comp_appr"] <- "bowtie2;2.4.4;ANACAPA-blca;2021"
taxmat[is.na(taxmat$otu_db), "otu_db"] <- taxmat[1,"otu_db"]

# Because WORMS doesn't recognize Eukaryota, change those that have this in the lastvalue to Biota:
taxmat$lastvalue <- recode(taxmat$lastvalue, "Eukaryota" = "Biota")

# Find WoRMS LSID with worrms

tax_names <- unique(na.omit(taxmat$lastvalue))

match_name <- function(name) {
  lsid <- tryCatch({
    res <- wm_records_names(name, marine_only = FALSE)
    matches <- res[[1]] %>%
      filter(match_type == "exact" | match_type == "exact_genus" | match_type == "exact_subgenus")
    if (nrow(matches) > 1) {
      message(paste0("Multiple matches for ", name))
    }
    return(matches[1,])
  }, error = function(cond) {
    return(NULL)
  })
}

matches <- sapply(tax_names, match_name)
taxmat$lsid <- sapply(taxmat$lastvalue, function(name) { ifelse(!is.null(matches[[name]]), matches[[name]]$lsid, NA) })
taxmat$taxonRank <- sapply(taxmat$lastvalue, function(name) { ifelse(!is.null(matches[[name]]), tolower(matches[[name]]$rank), NA) })
taxmat$specificEpithet <- sapply(taxmat$lastvalue, function(name) {
  if (!is.null(matches[[name]])) {
    if (!is.na(matches[[name]]$rank) & matches[[name]]$rank == "Species") {
      return(sub(".*\\s", "", matches[[name]]$scientificname))
    }
  }
  return(NA)
})

# Add Biota LSID in case there is no last value
# Kingdom is used as taxonRank so that "Biota" is also recognized correctly by GBIF

taxmat$lsid[is.na(taxmat$lastvalue)] <- "urn:lsid:marinespecies.org:taxname:1"
taxmat$taxonRank[is.na(taxmat$lastvalue)] <- "kingdom"

# Names not in WoRMS

names_not_in_worms <- names(matches)[sapply(matches, is.null)]
message("Number of species names not recognized in WORMS: ", length(names_not_in_worms))

# Some of these are common contaminants (e.g. Homo sapiens, Canis lupus)
# Some of these are not marine species (e.g. insects), that could be found with terrestrial matches
# Some have a higher taxonomy that can be found in the database
# Check at genus level:

genus_names <- unique(na.omit(taxmat$genus[is.na(taxmat$lsid)]))
genus_matches <- sapply(genus_names, match_name)
missing_lsids <- is.na(taxmat$lsid)

for (i in 1:nrow(taxmat)) {
  if (is.na(taxmat$lsid[i]) & !is.na(taxmat$genus[i]) & !is.null(genus_matches[[taxmat$genus[i]]])) {
    taxmat$lsid[i] <- genus_matches[[taxmat$genus[i]]]$lsid
    taxmat$taxonRank[i] <- "genus"
  }
}

# This way unknowns decreases. However, these values should be checked for possible addition to WoRMS?
# Will require manual inspection

genera_not_in_worms <- names(genus_matches)[sapply(genus_matches, is.null)]
message("Number of genus names not recognized in WORMS: ", length(genera_not_in_worms))
message("These taxa can be found in the table: ", outpath, "Taxa_not_in_worms.csv")

not_in_worms <- taxmat[is.na(taxmat$lsid),]

# Add 'Biota' as the name for the unknown sequences

taxmat$lastvalue[is.na(taxmat$lastvalue)] <- "Biota"

# Add sequence to the tax_table slot (linked to each asv)
taxmat$DNA_sequence <- as.character(rep_seqs[row.names(taxmat)])

# Write table of unknown names to make manual inspection easier:

write.table(not_in_worms, paste0(outpath, "Taxa_not_in_worms.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

# Write tax table
#taxmat <- apply(taxmat,2,as.character)
write.table(taxmat, paste0(outpath, "Full_tax_table_with_lsids.csv"), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE, na = "")

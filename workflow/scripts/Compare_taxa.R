library(dplyr)
library(tidyr)
library(yaml)
library(stringr)
library(glue)

#setwd("Dropbox (IPOfI)/Saara/github/Mock communities/Hleap2021/New_pacman/")

args <- commandArgs(trailingOnly = T)
config <- read_yaml(args[2])

outpath <- args[1]
rdp_file_path <- args[3]
vsearch_file_path <- args[4]
vsearch_lca_file_path <- args[5]

rdp <- read.csv(rdp_file_path, sep = "\t", header = F)
vsearch <- read.csv(vsearch_file_path, sep = "\t", header = F)
vsearch_lca <- read.csv(vsearch_lca_file_path, sep = "\t", header = F)

############# Modify Vsearch results ########################

colnames(vsearch) <- c("asv", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#Simplify the results table
#vsearch=vsearch %>% tidytable::separate_wider_delim(sseqid, ";", names = c("seqID", "taxonomy"), too_few = "align_start", too_many = "debug")

# Split sequence ID and taxonomy
vsearch <- vsearch %>%
  separate(sseqid, sep=";", into = c("vsearch_seqid", "vsearch_taxonomy")) %>%
  mutate(pident = as.numeric(pident)) %>%
  rename(vsearch_pident = pident)

# Remove rank prefixes from taxonomy
vsearch_lca <- vsearch_lca %>%
  mutate(across(2, gsub, pattern = "[a-z]:", replacement = "")) %>% 
  mutate(across(2, gsub, pattern = ",", replacement = ";")) %>%
  setNames(c("asv", "vsearch_consensus"))

#In each case we have multiple matches for the asvs. There is an lca output which already gives a consolidated taxonomy based on these results
# i.e. this is the common taxonomy based on all of the alignments. Use this output for the taxonomy but keep the seqids and percent identities from the original

# Aggregate vsearch results by asv, remove asvs without results
vsearch_mod <- vsearch %>%
  rowwise() %>%
  mutate(vsearch_taxonomy = str_replace(tail(strsplit(vsearch_taxonomy, ",")[[1]], 1), ".*:", "")) %>%
  group_by(asv) %>%
  summarize_all(~paste(na.omit(.), collapse = ",")) %>% 
  select(asv, vsearch_seqid, vsearch_pident, vsearch_taxonomy) %>%
  filter(vsearch_seqid != "*")

vsearch_combined <- vsearch_mod %>%
  inner_join(vsearch_lca, by = "asv")

################### Modify rdp results ####################################

confidence_threshold <- config$Rdp$cutoff

process_rdp_taxonomy <- function(row) {
  taxa <- row[seq(6, ncol(row), 3)]
  ranks <- row[seq(7, ncol(row), 3)]
  confidences <- row[seq(8, length(values), 3)]
  retained <- which(confidences >= confidence_threshold)
  rdp_taxonomy <- paste0(taxa, collapse=";")
  rdp_confidences <- paste0(confidences, collapse=";")
  rdp_consensus <- paste0(taxa[retained], collapse=";")
  return(data.frame(rdp_taxonomy, rdp_confidences, rdp_consensus))
}

rdp_mod <- rdp %>%
  rename(asv = 1) %>%
  rowwise() %>%
  mutate(process_rdp_taxonomy(across(everything()))) %>%
  ungroup() %>%
  select(asv, rdp_taxonomy, rdp_confidences, rdp_consensus)

######### Combine info from vsearch to rdp ###############

#vsearch_combi_mod = vsearch_combined %>% separate(taxonomy, sep=";", names = c("kingdom", "phylum", "class", "order","family", "genus", "species"), too_few = "align_start")

#Build output of asv, sum.taxonomy, and identificationRemarks
#sum.taxonomy is the taxonomic information separated by ;
#Identification remarks, should contain both the taxonomy and the confidence from rdp,
#as well as the results of the vsearch high confidence search. 
#Return also those that are fully unknown from rdp for a possible blast step. 
               
results <- rdp_mod %>%
  left_join(vsearch_combined, by = "asv") %>%
  mutate(
    rdp = glue("Identification based on the RDP classifier at the confidence level {config$Rdp$cutoff}: taxonomy {rdp_taxonomy}, confidences {rdp_confidences}."),
    vsearch = glue("Confirmation with VSEARCH against the {config$Vsearch$vsearch_db_name} database at 97% similarity: hits {vsearch_seqid}, identities {vsearch_pident}, taxonomy {vsearch_taxonomy}, consensus {vsearch_consensus}.")
  ) %>%
  mutate(
    identificationRemarks = ifelse(!is.na(vsearch_pident), glue("{rdp} {vsearch}"), glue("{rdp} No VSEARCH hits."))
  ) %>%
  select(asv, sum.taxonomy = rdp_consensus, identificationRemarks)

unknowns <- rdp_mod$asv[rdp_mod$rdp_consensus=="" | rdp_mod$rdp_consensus=="Eukaryota" | rdp_mod$rdp_consensus=="sk__Eukaryota"]
  
write.table(results, paste0(outpath, "04-taxonomy/rdp_vsearch_results.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep="\t")
writeLines(unknowns, paste0(outpath, "04-taxonomy/unknown_asvs.txt"))

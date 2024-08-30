library(dplyr)
library(tidyr)
library(Biostrings)

#This script removes species annotations where the top blast hit is less than 99% (and coverage is less than 90%) and genus annotations were the top blast hit is less than 97%
#LCA already calculates a consensus taxonomy based on filtering parameters and hits, but when we don't have good hits, we sometimes end up with species 
#annotations with low identity to the database. Therefore this filtering step is necessary so we don't have clearly incorrect annotations left in the dataset. 


if (!exists("cmd_args")) {
  cmd_args <- commandArgs(trailingOnly = T)
  message("cmd_args <- ", capture.output(dput(cmd_args)))
}

blast_hits <- read.csv(cmd_args[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
blast_lca <- read.csv(cmd_args[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#fasta_file <- cmd_args[3]
sequences <- readDNAStringSet(cmd_args[3])
outpath <- cmd_args[4]

#blast_hits <- read.csv("./runs/Bangladesh_12SMifish/04-taxonomy/blast/unclassified_blast_results.tab", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blast_hits) <- c("ASV", "SubjectID", "PID", "AlignmentLength", "Mismatches", "GapOpens", "QStart", "QEnd", "SStart", "SEnd", "EValue", "BitScore")

#blast_lca <- read.csv("./runs/Bangladesh_12SMifish/04-taxonomy/blast/unclassified_blast_results_lca.tab", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blast_lca)<-c("ASV", "Taxonomy")

blast_separated <- blast_lca %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right", extra = "drop")
head(blast_separated)

#Select the top hits (first PID, then alignmentLength, then EValue)
#Collect the first 5 NCBIids and these top values to the table
top_hits_full <- blast_hits %>%
  group_by(ASV) %>%
  filter(PID == max(PID)) %>%
  filter(AlignmentLength == max(AlignmentLength)) %>%
  filter(EValue == min(EValue)) %>%
  summarise(
    SubjectIDs = paste(unique(SubjectID)[1:5], collapse = ";"),
    PID = dplyr::first(PID),
    AlignmentLength = dplyr::first(AlignmentLength),
    EValue = dplyr::first(EValue)
  )%>%  
  ungroup()

#Remove NAs from the IDs
top_hits_full$SubjectIDs = gsub(";NA","", top_hits_full$SubjectIDs)
 

blast_lca_comb <- blast_separated %>% left_join(top_hits_full, by="ASV")

#Get the length of the sequences from the fasta file
# Extract the lengths of the sequences
sequence_lengths <- width(sequences)

# Create a data frame with ASV IDs and sequence lengths
asv_lengths <- data.frame(
    ASV = names(sequences),
    Length = sequence_lengths
)

blast_lca_comb <- blast_lca_comb %>% left_join(asv_lengths, by="ASV")

blast_lca_comb$coverage <- blast_lca_comb$AlignmentLength/blast_lca_comb$Length


blast_lca_comb <- blast_lca_comb %>%
    mutate(
        Species = ifelse(PID < 99, NA, Species),
        Species = ifelse(coverage < 0.95, NA, Species),
        Genus = ifelse(PID < 97, NA, Genus)
    )

blast_lca_mod <- blast_lca_comb %>%
    unite("Taxonomy", Kingdom:Species, sep = ";", remove = TRUE, na.rm = TRUE)
head(blast_lca_mod)

write.table(blast_lca_mod, outpath, sep="\t", row.names=F, quote=F)

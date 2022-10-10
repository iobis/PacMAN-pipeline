# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN
#Last update 21.6.2021

library("stringr")
library("phyloseq")
library("dplyr")
library("xml2")

args <- commandArgs(trailingOnly = T)

outpath <- args[1]

# LOAD data (this will be done from input later on)

tax_file <- read.csv(args[3], sep = "\t", header = T, row.names = 1)
otu_file <- read.csv(args[2], sep = "\t", header = T, row.names = 1)
#Rep_seqs=Biostrings::readDNAStringSet(args[4])
sample_file <- read.csv(args[4], sep = ";", header = T, row.names = 1)

# Add checks!
# 1. Make sure that the values are given
# 2. Make sure that the sample names in the otu-file (the original sample table), and the sample_data match!
# 3. Make sure the required fields are not empty
#       BasisOfRecord, eventdate, occurrencestatus (what about latitude and longitude?)

########################### 2. Add user provided fields to sample data ############################################################################################

# First change empty values to NA
args[args[5:15] == "None"] <- NA

sample_file$target_gene <- args[5]
sample_file$subfragment <- args[6]
sample_file$pcr_primer_forward <- args[7]
sample_file$pcr_primer_reverse <- args[8]
sample_file$pcr_primer_name_forw <- args[9]
sample_file$pcr_primer_name_reverse <- args[10]
sample_file$pcr_primer_reference <- args[11]
sample_file$lib_layout <- args[12]
sample_file$seq_meth <- args[13]
sample_file$sop <- args[14]
sample_file$votu_db <- args[15]

# Addition of possible extra fields:
args_name_value <- data.frame(command = 1, value = 1)

if (args[16] != "None") {
  extra_args <- str_split(args[16], ",", simplify = T)
  for (i in 1:length(extra_args)) {
    args_name_value[i,] <- str_split(extra_args[i], ":", simplify = T)
  }
  for(i in 1:nrow(args_name_value)) {
    sample_file[,args_name_value[i, 1]] <- args_name_value[i, 2]
  }
}

########################### 3. Collect all values to phyloseq object ######################################

tax_table <- phyloseq::tax_table(as(tax_file, "matrix"))
otu_table <- phyloseq::otu_table(otu_file, taxa_are_rows = T)
sample_data <- phyloseq::sample_data(sample_file)

# Here I make a phyloseq object with the three files
phydata <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
# Print the amount of information stored:
phydata

# Save the phyloseq rdata object for easier access in the future
print("Saving the phyloseq -table to an Rdata object, for ease of access for data analysis later")
print("The object can be loaded with readRDS, while the phyloseq package and library is loaded")
print(paste("The saved object can be found here: ", outpath, "phyloseq_object.rds", sep = ""))
saveRDS(phydata, paste0(outpath, "phyloseq_object.rds"))

# Remove the OTUs that are found in the control samples (occurrenceStatus==absent)
if ("absent" %in% sample_data$occurrenceStatus) {
  control_taxa <- taxa_names(filter_taxa(subset_samples(phydata, occurrenceStatus == "absent"), function(x) sum(x) > 0, TRUE))
  good_taxa <- taxa_names(phydata)[!(taxa_names(phydata) %in% control_taxa)]
  phydata_no_control <- prune_taxa(good_taxa, phydata)
}

# Here add the total read counts in each sample to the sample_data table:
# Here we should add a check that the samples are in the right order:
if ("absent" %in% sample_data$occurrenceStatus) {
  sample_data(phydata_no_control)$sampleSizeValue <- sample_sums(phydata)
  sample_data(phydata_no_control)$organismQuantityType <- "DNA Sequence reads"
  sample_data(phydata_no_control)$sampleSizeUnit <- "DNA Sequence reads"
  phydf <- psmelt(phydata_no_control)
} else {
  sample_data(phydata)$sampleSizeValue <- sample_sums(phydata)
  sample_data(phydata)$organismQuantityType <- "DNA Sequence reads"
  sample_data(phydata)$sampleSizeUnit <- "DNA Sequence reads"
  phydf <- psmelt(phydata)
}

phydf$occurrenceID <- paste(phydf$OTU, phydf$Sample, sep = "_")

# Change names where necessary
#names(phydf[names(phydf)=="lastvalue"])="ScientificName"
phydf <- phydf %>%
          rename(
              scientificName = lastvalue,
              organismQuantity = Abundance,
              scientificNameID = lsid
            )
            # %>%
            # add_column(organismQuantityType = "DNA Sequence reads") %>%
            # add_column(sampleSizeUnit = "DNA Sequence reads")

# Remove those occurrences which did not return an lsid:
phydf <- phydf[!is.na(phydf$scientificNameID),]

# Remove 0 Abundance data (not valuable for us)
phydf_present <- phydf[phydf$organismQuantity > 0,]

# Write tables with all the fields found in the current tables:

get_dwc_fields <- function(spec_url) {
  doc <- read_xml(spec_url)
  doc %>%
    xml_ns_strip() %>%
    xml_find_all("//property") %>%
    xml_attr(attr = "name")
}

spec_occurrence <- "https://rs.gbif.org/core/dwc_occurrence_2020-07-15.xml"
occurrence_table_fields <- get_dwc_fields(spec_occurrence)

spec_dna <- "https://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2021-07-05.xml"
DNA_extension_fields <- get_dwc_fields(spec_dna)

occurrence_table <- phydf_present[,colnames(phydf_present) %in% occurrence_table_fields]
DNA_derived_data_extension <- phydf_present[,colnames(phydf_present) %in% c('occurrenceID', DNA_extension_fields)]

write.table(occurrence_table, paste0(outpath, "Occurence_table.csv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")
write.table(DNA_derived_data_extension, paste0(outpath, "DNA_extension_table.csv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "")

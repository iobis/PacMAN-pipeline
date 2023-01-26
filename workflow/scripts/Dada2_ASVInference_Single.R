#####
# Dada2 Step 2: ASV inference

# @author: Saara Suominen
# last major modification on the 8.12.2021

#input=snakemake@input[[1]] Not working with the shell command, files are read as {input}
library(plyr)
library(dplyr)
library(yaml)
library(dada2)
library(Biostrings)
library(ggplot2)
library(tidyr)
library(tibble)

args <- commandArgs(trailingOnly = T)
config <- read_yaml(args[2])

# Set the different paths for all the supplied libraries
# NOTICE: only forward files given as input
filtFs <- args[3:length(args)]
filtFs_single <- gsub("_1P", "_1U", filtFs)

filtRs <- gsub("_1P", "_2P", filtFs)
filtRs_single <- gsub("_2P", "_2U", filtRs)

outpath <- args[1]

# The parameters nti/ntj, need the chosen nucleotides in a vector format:
nti <- strsplit(config$DADA2$plotERRORS$nti, "")[[1]]
if (length(nti) == 0) {
  message("Default value given to nti, all nucleotide transitions will be shown in the error plot")
  nti <- c("A", "C", "G", "T")
}

ntj <- strsplit(config$DADA2$plotERRORS$ntj, "")[[1]]
if (length(nti) == 0) {
  message("Default value given to ntj, all nucleotide transitions will be shown in the error plot")
  ntj <- c("A", "C", "G", "T")
}

# Get sample names
sample.names <- gsub("_1P.fastq.gz", "", basename(filtFs))
message(paste0("Sample ", sample.names, " will be analyzed", collapse = "\n"))

# Assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names
names(filtFs_single) <- sample.names
names(filtRs_single) <- sample.names

allfiles <- list(filtFs, filtRs, filtFs_single, filtRs_single)
names(allfiles) <- c("filtFs", "filtRs", "filtFs_single", "filtRs_single")
errs <- list()
dereps <- list()
dadas <- list()
seqtab <- list()

# Loop through all file types (forward, reverse, unpaired forward, unpaired reverse) for learning errors and dereplicating
for (i in 1:4) {

  if (any(file.exists(allfiles[[i]]))) {

    message(paste("learning error rates of files:", i, ": (1) forward paired (2) reverse paired (3) forward single and (4) reverse single " , sep=" "))

    errs[[i]] <- learnErrors(allfiles[[i]][file.exists(allfiles[[i]])],
                            multithread = config$DADA2$learnERRORS$multithread,
                            nbases = as.numeric(config$DADA2$learnERRORS$nbases),
                            randomize = config$DADA2$learnERRORS$randomize,
                            MAX_CONSIST = as.numeric(config$DADA2$learnERRORS$MAX_CONSIST),
                            OMEGA_C = as.numeric(config$DADA2$learnERRORS$OMEGA_C),
                            verbose = config$DADA2$learnERRORS$verbose)

    message("Making error estimation plots of reads")
    png(filename = paste0(outpath, "06-report/dada2/error_profile_", names(allfiles)[i], ".png"))
      p_ERR <- plotErrors(
        errs[[i]],
        nti = nti,
        ntj = ntj,
        obs = config$DADA2$plotERRORS$obs,
        err_out = config$DADA2$plotERRORS$err_out,
        err_in = config$DADA2$plotERRORS$err_in,
        nominalQ = config$DADA2$plotERRORS$nominalQ)
      print(p_ERR)
    dev.off()

    message("Running dereplication of reads")
    dereps[[i]] <- derepFastq(allfiles[[i]][file.exists(allfiles[[i]])],
                              n = as.numeric(config$DADA2$derepFastq$num),
                              config$DADA2$learnERRORS$verbose)

    message("Running dada on reads")
    dadas[[i]] <- dada(dereps[[i]],
                      errs[[i]],
                      selfConsist = config$DADA2$dada$selfConsist,
                      pool = config$DADA2$dada$pool,
                      priors = config$DADA2$dada$priors,
                      multithread = config$DADA2$learnERRORS$multithread,
                      verbose = config$DADA2$learnERRORS$verbose)

    # Convert to list in case there's only one file
    if (is(dadas[[i]], "dada")) {
      dadas[[i]] <- list(dadas[[i]])
      names(dadas[[i]]) <- names(allfiles[[i]][file.exists(allfiles[[i]])])
    }
    if (is(dereps[[i]], "derep")) {
      dereps[[i]] <- list(dereps[[i]])
      names(dereps[[i]]) <- names(allfiles[[i]][file.exists(allfiles[[i]])])
    }

  # If no files found for the paired reads (should not be the case!)
  } else if (!grepl("single", names(allfiles)[i])) {

    message("Error: no paired reads to process")
    stop()

  # If no files found for the unpaired reads, continue with the workflow
  # It has to be made sure in the snakefile that this step is run despite not requiring output files
  } else {

    message("No further unpaired reads to process")

  }
}

# Merge forward and reverse paired reads
message("Attempting merge")

# Define function for merging and formatting seqtabs: next steps require an integer matrix
merge_format_seqtab <- function(seqtab1, seqtab2) {
  df1 <- reshape2::melt(seqtab1, varnames = c("sample", "sequence"))
  df2 <- reshape2::melt(seqtab2, varnames = c("sample", "sequence"))
  df <- bind_rows(df1, df2) %>%
    group_by(sample, sequence) %>%
    summarize(value = sum(value)) %>%
    ungroup()
  m <- reshape2::acast(df, sample ~ sequence, value.var = "value")
  mode(m) <- "integer"
  return(m)
}

if (config$DADA2$mergePairs$include) {
  message("merging pairs")
  mergers <- mergePairs(dadas[[1]],
                        dereps[[1]],
                        dadas[[2]],
                        dereps[[2]],
                        minOverlap = as.numeric(config$DADA2$mergePairs$minOverlap),
                        maxMismatch = as.numeric(config$DADA2$mergePairs$maxMismatch),
                        returnRejects = config$DADA2$mergePairs$returnRejects,
                        propagateCol = config$DADA2$mergePairs$propagateCol,
                        justConcatenate = config$DADA2$mergePairs$justConcatenate,
                        trimOverhang = config$DADA2$mergePairs$trimOverhang,
                        verbose = config$DADA2$learnERRORS$verbose)
  if (is.data.frame(mergers)) {
    mergers <- list(mergers) %>% setNames(sample.names)
  }

  # Create read/ASV mapping
  mergers_all <- mergePairs(dadas[[1]], dereps[[1]], dadas[[2]], dereps[[2]],
    minOverlap = as.numeric(config$DADA2$mergePairs$minOverlap), maxMismatch = as.numeric(config$DADA2$mergePairs$maxMismatch),
    returnRejects = TRUE, propagateCol = config$DADA2$mergePairs$propagateCol,
    justConcatenate = config$DADA2$mergePairs$justConcatenate, trimOverhang = config$DADA2$mergePairs$trimOverhang,
    verbose = config$DADA2$learnERRORS$verbose)
  if (is.data.frame(mergers_all)) {
    mergers_all <- list(mergers_all) %>% setNames(sample.names)
  }

  mapping <- lapply(names(mergers_all), function(name) {
    merger_to_dada <- lapply(mergers_all[[name]]$forward, function(x) {
      which(dadas[[1]][[name]]$map == x)
    })
    merger_to_read <- lapply(merger_to_dada, function(x) {
      which(dereps[[1]][[name]]$map %in% x)
    })
    names(merger_to_read) <- mergers_all[[name]]$sequence
    merger_to_read
  })
  names(mapping) <- names(mergers_all)
  for (name in names(mapping)) {
    fq <- microseq::readFastq(filtFs[[name]])
    headers <- sub(" .*", "", fq$Header)
    mapping[[name]] <- mapping[[name]] %>%
      enframe(name = "sequence", value = "read") %>%
      unnest() %>%
      filter(sequence != "") %>%
      mutate(read = headers[read]) %>%
      arrange(sequence, read)
  }

  seqtab <- makeSequenceTable(mergers)

  # When merging is done with returnRejects=TRUE, the abundance of the rejected merges is returned, but not the sequence
  # We want to collect also these single sequences and add them to the seqtab (to avoid loosing ANY data)
  if (config$DADA2$mergePairs$returnRejects == TRUE) {

    unmerged_f <- list()
    unmerged_r <- list()
    concatenated <- list()

    for (i in 1:length(sample.names)) {
      unmerged_f[[i]] <- dadas[[1]][[sample.names[i]]]$sequence[mergers[[sample.names[i]]]$forward[!mergers[[sample.names[i]]]$accept]]
      unmerged_r[[i]] <- dadas[[2]][[sample.names[i]]]$sequence[mergers[[sample.names[i]]]$reverse[!mergers[[sample.names[i]]]$accept]]
      # Here for the rejected reads (!merger$sample$accept) the indices are collected (merger$sample$forward, merger$sample$reverse)
      # The sequences are sourced from the original dada-file (dadaF$sample$denoised, dadaR$sample$denoised)
      # It seems that concatenating these reads and keeping them for further analyses can result in better taxonomic coverage (Dacey et al. 2021 https://doi.org/10.1186/s12859-021-04410-2)
      # Abundances for these reads is taken from the merged abundances.
      # reverse complement reverse reads so that the following taxonomic assignment will work optimally.
      unmerged_r[[i]] <- sapply(sapply(sapply(unmerged_r[[i]], DNAString), reverseComplement), toString)
      sequence <- paste0(unmerged_f[[i]], unmerged_r[[i]])
      abundance <- mergers[[sample.names[i]]]$abundance[!mergers[[sample.names[i]]]$accept]
      concatenated[[i]] <- tibble(sequence, abundance)
    }
    names(concatenated) <- sample.names
    #names(unmerged_r)=sample.names

    # Make sequence table
    seqtab_unmerged <- makeSequenceTable(concatenated)

    # The merged returnrejects = T seqtab also contains a column with an empty header,
    # This is all rejected (non-merged) abundances combined.
    # This column messes with future steps, so we want to remove it.
    # We have instead collected the abundances of the unmerged reads to add to the table with sequences.
    seqtab <- seqtab[,-which(colnames(seqtab) == "")]

    # Merge with paired reads and format for the next steps (integer matrix)
    seqtab <- merge_format_seqtab(seqtab, seqtab_unmerged)

  }

} else {

  message("no merging of paired reads")
  # Combine forward and reverse sequences to one table
  seqtab1 <- makeSequenceTable(dadas[[1]])
  seqtab2 <- makeSequenceTable(dadas[[2]])
  # Reverse complement reverse reads so that the following taxonomic assignment will work optimally.
  colnames(seqtab2) <- sapply(sapply(sapply(colnames(seqtab2), DNAString), reverseComplement), toString)
  seqtab <- cbind(seqtab1, seqtab2)

}

# Add ASVs from single reads to full table, and format table to the right format to continue with the pipeline
# It has to be an integer matrix with samples as rownames and sequences as column names
if (any(file.exists(allfiles[[3]]))) {
  message("Adding ASVs from unpaired forward reads to ASV-table")
  seqtab3 <- makeSequenceTable(dadas[[3]])
  seqtab <- merge_format_seqtab(seqtab, seqtab3)
}

if (any(file.exists(allfiles[[4]]))) {
  message("Adding ASVs from unpaired reverse reads to ASV-table")
  seqtab4 <- makeSequenceTable(dadas[[4]])
  # Here also the sequences from the reverse reads are reverse complemented before they are added to the sequence table
  colnames(seqtab4) <- sapply(sapply(sapply(colnames(seqtab4), DNAString), reverseComplement), toString)
  seqtab <- merge_format_seqtab(seqtab, seqtab4)
}

# Chimeras removed from the full combined table as per: https://github.com/benjjneb/dada2/issues/1235:
# We think the best way (in most cases, using current common techs -- such is the challenge of recommendations) is to combine the tables from multiple runs and then remove chimeras on that combined table.
# Look into details of how chimera removal is done to understand if this is smart
# In their case they are looking at equal length reads, possibly the single reads and paired reads should not be combined for chimera removal?

message("removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = config$DADA2$removeBimeraDenovo$method,
                                    multithread = config$DADA2$learnERRORS$multithread,
                                    verbose = config$DADA2$learnERRORS$verbose)

print(dim(seqtab.nochim))
# ASVs denominated by the actual sequence, we want to simplify the names.
new.names <- c(paste("asv.", 1:length(colnames(seqtab.nochim)), sep = ""))
message(head(new.names))

# Save fasta, before changing names in the seqtab table
message("Making fasta table")
uniquesToFasta(seqtab.nochim, fout = paste0(outpath, "03-dada2/rep-seqs.fna"), ids = new.names)

# Export read/asv mapping
if (exists("mapping")) {
  message("Exporting read/asv mapping")
  asv_sequences <- data.frame(asv = new.names, sequence = colnames(seqtab.nochim))
  for (name in names(mapping)) {
    mapping[[name]] <- mapping[[name]] %>%
      left_join(asv_sequences, by = "sequence") %>%
      select(read, asv) %>%
      filter(!is.na(asv))
    write.table(mapping[[name]], paste0(outpath, "03-dada2/mapping/", name, "/", name, "_mapping.txt"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

# Finally, change names in otu table
message("Changing sequence names")
colnames(seqtab.nochim) <- new.names

# This show sequence length distributions (see if you should include this)
#seq_hist <- table(nchar(getSequences(seqtab)))
#fname_seqh <- paste(args[6],"seq_hist.txt",sep="")
#write.table(seq_hist, file = fname_seqh  , sep = "\t", quote=FALSE, col.names = FALSE)

# Collect results of how many reads are available at each step in a table:
getN <- function(x) sum(getUniques(x))

# Make a table with all information on the reads retained from the run, if paired reads were merged:
message("Making summary table")
if (config$DADA2$mergePairs$include) {
  track <- cbind(sapply(dadas[[1]], getN), sapply(dadas[[2]], getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim != 0))
  colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim", "ASVs")
# If paired reads were not merged:
} else {
  track <- cbind(sapply(dadas[[1]], getN), sapply(dadas[[2]], getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim != 0))
  colnames(track) <- c("denoisedF", "denoisedR", "nonchim", "ASVs")
}

# Add also information on the reads that came from possibly evaluated single reads
if (any(file.exists(allfiles[[3]]))) {
  message("Adding info from unpaired forward reads to the summary table")
  denoisedF_single <- sapply(dadas[[3]], getN)
  track <- merge(track, data.frame(denoisedF_single), by = 0, all = TRUE)
  rownames(track) <- track$Row.names
  track <- subset(track, select = -Row.names)
  # Reorder columns:
  track <- track %>% relocate(denoisedF_single, .before = nonchim)
}

if (any(file.exists(allfiles[[4]]))) {
  message("Adding info from unpaired reverse reads to the summary table")
  denoisedR_single <- sapply(dadas[[4]], getN)
  track <- merge(track, data.frame(denoisedR_single), by = 0, all = TRUE)
  rownames(track) <- track$Row.names
  track <- subset(track, select = -Row.names)
  track <- track %>% relocate(denoisedR_single, .before = nonchim)
}

rownames(track)
rownames(track) <- sample.names
message(track)

# Read results of filtering step and append the results of ASV step:
out <- read.table(paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), header = TRUE)
track <- cbind(out, track)

# Write output tables of reads and the otu_table:
write.table(track, paste0(outpath, "06-report/dada2/dada2_stats.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(t(seqtab.nochim), paste0(outpath,"03-dada2/seqtab-nochim.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Fix for https://github.com/tidyverse/ggplot2/issues/2787
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

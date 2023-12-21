#####
# Dada2 Step 2: ASV inference
#
# @author: J. Engelmann
# @author: A. ABdala
# modified by Saara Suominen on the 15.6.2021

#input=snakemake@input[[1]] Not working with the shell command, files are read as {input}
library(yaml)
library(dada2)
library(Biostrings)
library(ggplot2)

args <- commandArgs(trailingOnly = T)
message("args <- ", capture.output(dput(args))) # output for debugging
config <- read_yaml(args[2])

# Set the different paths for all the supplied libraries
#NOTICE: only forward files given as input
filtFs <- args[3:length(args)]
filtRs <- gsub("_1P", "_2P", filtFs)

outpath <- args[1]

# Get sample names
sample.names <- gsub("_1P.fastq.gz", "", basename(filtFs))
message(paste0("Sample ", sample.names, " will be analyzed", collapse = "\n"))

# assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print("learning error rates")
errF <- learnErrors(filtFs, multithread = config$DADA2$learnERRORS$multithread,
  nbases = as.numeric(config$DADA2$learnERRORS$nbases), randomize = config$DADA2$learnERRORS$randomize,
  MAX_CONSIST = as.numeric(config$DADA2$learnERRORS$MAX_CONSIST), OMEGA_C = as.numeric(config$DADA2$learnERRORS$OMEGA_C),
  verbose = config$DADA2$learnERRORS$verbose)

errR <- learnErrors(filtRs, multithread = config$DADA2$learnERRORS$multithread,
  nbases = as.numeric(config$DADA2$learnERRORS$nbases), randomize = config$DADA2$learnERRORS$randomize,
  MAX_CONSIST = as.numeric(config$DADA2$learnERRORS$MAX_CONSIST), OMEGA_C = as.numeric(config$DADA2$learnERRORS$OMEGA_C),
  verbose = config$DADA2$learnERRORS$verbose)

# Plot estimated erros as a sanity check:

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

print("Making error estimation plots")
png(filename = paste0(outpath, "06-report/dada2/error_profile_forward.png"))
plotErrors(errF, nti = nti, ntj = ntj,
  obs = config$DADA2$plotERRORS$obs, err_out = config$DADA2$plotERRORS$err_out,
  err_in = config$DADA2$plotERRORS$err_in, nominalQ = config$DADA2$plotERRORS$nominalQ)
dev.off()

png(filename = paste0(outpath, "06-report/dada2/error_profile_reverse.png"))
plotErrors(errR, nti = nti, ntj = ntj,
  obs = config$DADA2$plotERRORS$obs, err_out = config$DADA2$plotERRORS$err_out,
  err_in = config$DADA2$plotERRORS$err_in, nominalQ = config$DADA2$plotERRORS$nominalQ)
dev.off()

message("Running dereplication")
derepFs <- derepFastq(filtFs, n = as.numeric(config$DADA2$derepFastq$num), config$DADA2$learnERRORS$verbose)
derepRs <- derepFastq(filtRs, n = as.numeric(config$DADA2$derepFastq$num), config$DADA2$learnERRORS$verbose)

message("Running dada")
dadaFs <- dada(derepFs, errF, selfConsist = config$DADA2$dada$selfConsist, pool = config$DADA2$dada$pool,
  priors = config$DADA2$dada$priors, multithread = config$DADA2$learnERRORS$multithread,
  verbose = config$DADA2$learnERRORS$verbose)
dadaRs <- dada(derepRs, errR, selfConsist = config$DADA2$dada$selfConsist, pool = config$DADA2$dada$pool,
  priors = config$DADA2$dada$priors, multithread = config$DADA2$learnERRORS$multithread,
  verbose = config$DADA2$learnERRORS$verbose)
print(names(dadaFs)); print(names(dadaRs))

message("Attempting merge")
if (tolower(config$meta$sequencing$lib_layout) == "paired") {
    message("Merging pairs")
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
      minOverlap = as.numeric(config$DADA2$mergePairs$minOverlap), maxMismatch = as.numeric(config$DADA2$mergePairs$maxMismatch),
      returnRejects = config$DADA2$mergePairs$returnRejects, propagateCol = config$DADA2$mergePairs$propagateCol,
      justConcatenate = config$DADA2$mergePairs$justConcatenate, trimOverhang = config$DADA2$mergePairs$trimOverhang,
      verbose = config$DADA2$learnERRORS$verbose)

    seqtab <- makeSequenceTable(mergers)

} else {
  message("No merging")
  seqtab1 <- makeSequenceTable(dadaFs)
  seqtab2 <- makeSequenceTable(dadaRs)
  seqtab <- cbind(seqtab1, seqtab2)
}

message("Removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = config$DADA2$removeBimeraDenovo$method,
  multithread = config$DADA2$learnERRORS$multithread, verbose = config$DADA2$learnERRORS$verbose)

print(dim(seqtab.nochim))
# ASVs denominated by the actual sequence, we want to simplify the names.
new.names <- c(paste("asv.", 1:length(colnames(seqtab.nochim)), sep = ""))
message(head(new.names))

# Save fasta, before changing names in the seqtab table
message("Making fasta table")
uniquesToFasta(seqtab.nochim, fout = paste0(outpath, "03-dada2/rep-seqs.fna"), ids = new.names)

# Finally, change names in otu table
colnames(seqtab.nochim) <- new.names

#this show sequence length distributions (see if you should include this)
#seq_hist <- table(nchar(getSequences(seqtab)))
#fname_seqh <- paste(args[6],"seq_hist.txt",sep="")
#write.table(seq_hist, file = fname_seqh  , sep = "\t", quote=FALSE, col.names = FALSE)

# Collect results of how many reads are available at each step in a table:
getN <- function(x) sum(getUniques(x))

#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim!=0))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","ASVs")
#rownames(track) <- sample.names
if (config$meta$sequencing$lib_layout=="Paired"|config$meta$sequencing$lib_layout=="paired"|config$meta$sequencing$lib_layout=="PAIRED") {
  track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim != 0))
  colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim", "ASVs")
} else {
  track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim != 0))
  colnames(track) <- c("denoisedF", "denoisedR", "nonchim", "ASVs")
}

rownames(track) <- sample.names
message(track)
# Read results of filtering step and append the results of ASV step:
out <- read.table(paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), header = TRUE)
print(out)
track <- cbind(out, track)

# Write output tables of reads and the otu_table:
write.table(track, paste0(outpath, "06-report/dada2/dada2_stats.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(t(seqtab.nochim), paste0(outpath,"03-dada2/seqtab-nochim.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Fix for https://github.com/tidyverse/ggplot2/issues/2787
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

#####
# Dada2 Step 2: ASV inference
#
# @author: J. Engelmann
# @author: A. ABdala
# modified by Saara Suominen on the 15.6.2021

#input=snakemake@input[[1]] Not working with the shell command, files are read as {input}

library(dada2)
library(Biostrings)
library(ggplot2)

args <- commandArgs(trailingOnly = T)

# Add all arguments for dada2 parameters from configfile!
#args[1]... = outpath for dada2 files (project/run/dada2/)
#learnERRORS:
multithread <-     as.logical(args[2])
nbases <-          as.numeric(args[3])
randomize <-       as.logical(args[4])
MAX_CONSIST <-     as.numeric(args[5])
OMEGA_C <-         as.numeric(args[6])
verbose <-         as.logical(args[7])
# PlotERRORS
obs <-             as.logical(args[10])
err_out <-         as.logical(args[11])
err_in <-          as.logical(args[12])
nominalQ <-        as.logical(args[13])
#derepFastq
n <-               as.numeric(args[14])
#dada
selfConsist <-     as.logical(args[15])
pool <-            as.logical(args[16])
priors <-          if (args[17] == "None") "" else args[17]
mergePairs <-      as.logical(args[18])
minOverlap <-      as.numeric(args[19])
maxMismatch <-     as.numeric(args[20])
returnRejects <-   as.logical(args[21])
propagateCol <-    if (args[22] == "None") "" else args[22]
justConcatenate <- as.logical(args[23])
trimOverhang <-    as.logical(args[24])
#removeBimeraDenovo
method <-          args[25]
#args[25:n]...= input.files}

# The parameters nti/ntj, need the chosen nucleotides in a vector format:
nti <- strsplit(args[8], "")[[1]]
if (length(nti) == 0) {
  message("Default value given to nti, all nucleotide transitions will be shown in the error plot")
  nti <- c("A", "C", "G", "T")
}

ntj <- strsplit(args[9], "")[[1]]
if (length(nti) == 0) {
  message("Default value given to ntj, all nucleotide transitions will be shown in the error plot")
  ntj <- c("A", "C", "G", "T")
}

# Set the different paths for all the supplied libraries
#NOTICE: only forward files given as input
filtFs <- args[26:length(args)]
#print(filtFs)
filtRs <- gsub("_1P", "_2P", filtFs)
#print(filtRs)

outpath <- args[1]
#print(outpath)

#List files
#filtFs <- sort(list.files(paths, pattern="_1P.fastq.gz", full.names = TRUE))
#filtRs <- sort(list.files(paths, pattern="_2P.fastq.gz", full.names = TRUE))
#filtFs <- gsub("//", "/", filtFs)
#filtRs <- gsub("//", "/", filtRs)


# Get sample names
sample.names <- gsub("_1P.fastq.gz", "", basename(filtFs))
message(paste0("Sample ", sample.names, " will be analyzed", collapse = "\n"))

#print(args[1:26])

# assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# learn error rates
errF <- learnErrors(filtFs, multithread = multithread, nbases = nbases, randomize = randomize,
  MAX_CONSIST = MAX_CONSIST, OMEGA_C = OMEGA_C, verbose = verbose)
errR <- learnErrors(filtRs, multithread = multithread, nbases = nbases, randomize = randomize,
  MAX_CONSIST = MAX_CONSIST, OMEGA_C = OMEGA_C, verbose = verbose)

# Plot estimated erros as a sanity check:
png(filename = paste0(outpath, "06-report/dada2/error_profile_forward.png"))
plotErrors(errF, nti = nti, ntj = ntj, obs = obs, err_out = err_out,
  err_in = err_in, nominalQ = nominalQ)
dev.off()

png(filename = paste0(outpath, "06-report/dada2/error_profile_reverse.png"))
plotErrors(errR, nti = nti, ntj = ntj, obs = obs, err_out = err_out,
  err_in = err_in, nominalQ = nominalQ)
dev.off()

derepFs <- derepFastq(filtFs, n = n, verbose = verbose)
derepRs <- derepFastq(filtRs, n = n, verbose = verbose)

dadaFs <- dada(derepFs, err = errF, selfConsist = selfConsist, pool = pool, priors = priors,
  multithread = multithread, verbose = verbose)
dadaRs <- dada(derepRs, err = errR, selfConsist = selfConsist, pool = pool, priors = priors,
  multithread = multithread, verbose = verbose)

if (mergePairs) {

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
  minOverlap = minOverlap, maxMismatch = maxMismatch, returnRejects = returnRejects,
  propagateCol = propagateCol, justConcatenate = justConcatenate, trimOverhang = trimOverhang, verbose = verbose)

seqtab <- makeSequenceTable(mergers)

} else {

seqtab <- makeSequenceTable(dadaFs)

}

seqtab.nochim <- removeBimeraDenovo(seqtab, method = method, multithread = multithread, verbose = verbose)

# ASVs denominated by the actual sequence, we want to simplify the names.
new.names <- c(paste("asv.", 1:length(colnames(seqtab.nochim)), sep = ""))
message(head(new.names))

# Save fasta, before changing names in the seqtab table
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
if (mergePairs) {
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


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

#Add all arguments for dada2 parameters from configfile!
#args[1]... = outpath for dada2 files (project/run/dada2/)
#learnERRORS:
multithread =     as.logical(args[2])
nbases =          as.numeric(args[3])
randomize =       as.logical(args[4])
MAX_CONSIST =     as.numeric(args[5])
OMEGA_C =         as.numeric(args[6])
verbose =         as.logical(args[7])
#PlotERRORS
#nti=args[8] See later: built with if statement
#ntj=args[9]
obs =             as.logical(args[10])
err_out =         as.logical(args[11])
err_in =          as.logical(args[12])
nominalQ =        as.logical(args[13])
#derepFastq
n =               as.numeric(args[14])
#dada
selfConsist =     as.logical(args[15])
pool =            as.logical(args[16])
if(args[17]=="None"){ priors="" } else { propagateCol=args[17] }
#mergePairs
minOverlap =      as.numeric(args[18])
maxMismatch =     as.numeric(args[19])
returnRejects =   as.logical(args[20])
if(args[21]=="None"){ propagateCol="" } else { propagateCol=args[21] }
justConcatenate = as.logical(args[22])
trimOverhang =    as.logical(args[23])
#removeBimeraDenovo
method =          args[24]
#args[25:n]...= input.files}

#The parameters nti/ntj, need the chosen nucleotides in a vector format:
nti=c()
if (grepl("A", args[8])){
  nti=c(nti, "A")}
if (grepl("C", args[8])){
  nti=c(nti, "C")}
if (grepl("T", args[8])){
  nti=c(nti, "T")}
if (grepl("G", args[8])){
  nti=c(nti, "G")}
if (length(nti)==0) {
  print("Default value given to nti, all nucleotide transitions will be shown in the error plot")
  nti=c("A", "C", "G", "T")}


ntj=c()
if (grepl("A", args[9])){
  ntj=c(ntj, "A")}
if (grepl("C", args[9])){
  ntj=c(ntj, "C")}
if (grepl("T", args[9])){
  ntj=c(ntj, "T")}
if (grepl("G", args[9])){
  ntj=c(ntj, "G")}
if (length(ntj)==0) {
  print("Default value given to ntj, all nucleotide transitions will be shown in the error plot")
  ntj=c("A", "C", "G", "T")}


#Set the different paths for all the supplied libraries
#NOTICE: only forward files given as input
filtFs=args[25:length(args)]
#print(filtFs)
filtRs=gsub('_1P','_2P',filtFs)
#print(filtRs)

outpath <- args[1]
#print(outpath)

#List files
#filtFs <- sort(list.files(paths, pattern="_1P.fastq.gz", full.names = TRUE))
#filtRs <- sort(list.files(paths, pattern="_2P.fastq.gz", full.names = TRUE))
#filtFs <- gsub("//", "/", filtFs)
#filtRs <- gsub("//", "/", filtRs)


#Get sample names
sample.names <- gsub('_1P.fastq.gz', '', basename(filtFs))
print(paste("Sample", sample.names, "will be analysed"))
#print(args[1:26])

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# learn error rates
errF <- learnErrors(filtFs, multithread=multithread, nbases=nbases, randomize=randomize,
  MAX_CONSIST=MAX_CONSIST, OMEGA_C=OMEGA_C, verbose=verbose)
errR <- learnErrors(filtRs, multithread=multithread, nbases=nbases, randomize=randomize,
  MAX_CONSIST=MAX_CONSIST, OMEGA_C=OMEGA_C, verbose=verbose)

#Plot estimated erros as a sanity check:
png(filename=paste0(outpath, "06-report/dada2/error_profile_forward.png"))
plotErrors(errF, nti=nti, ntj=ntj, obs=obs, err_out=err_out,
  err_in=err_in, nominalQ=nominalQ)
dev.off()

png(filename=paste0(outpath, "06-report/dada2/error_profile_reverse.png"))
plotErrors(errR, nti=nti, ntj=ntj, obs=obs, err_out=err_out,
  err_in=err_in, nominalQ=nominalQ)
dev.off()

derepFs <- derepFastq(filtFs, n=n, verbose=verbose)
derepRs <- derepFastq(filtRs, n=n, verbose=verbose)

dadaFs <- dada(derepFs, err=errF, selfConsist=selfConsist, pool=pool, priors=priors,
  multithread=multithread, verbose=verbose)
dadaRs <- dada(derepRs, err=errR, selfConsist=selfConsist, pool=pool, priors=priors,
  multithread=multithread, verbose=verbose)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
  minOverlap=minOverlap, maxMismatch=maxMismatch, returnRejects=returnRejects,
  propagateCol=propagateCol, justConcatenate=justConcatenate, trimOverhang=trimOverhang, verbose=verbose)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method=method, multithread=multithread, verbose=verbose)

#ASVs denominated by the actual sequence, we want to simplify the names.
new.names=c(paste("asv.",1:length(colnames(seqtab.nochim)),sep=""))
print(head(new.names))

#Save fasta, before changing names in the seqtab table
uniquesToFasta(seqtab.nochim, fout=paste0(outpath, '03-dada2/rep-seqs.fna'), ids=new.names)

#Finally, change names in otu table
colnames(seqtab.nochim)=new.names

#this show sequence length distributions (see if you should include this)
#seq_hist <- table(nchar(getSequences(seqtab)))
#fname_seqh <- paste(args[6],"seq_hist.txt",sep="")
#write.table(seq_hist, file = fname_seqh  , sep = "\t", quote=FALSE, col.names = FALSE)

#Collect results of how many reads are available at each step in a table:
getN <- function(x) sum(getUniques(x))

#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim!=0))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","ASVs")
#rownames(track) <- sample.names

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim!=0))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim","ASVs")
rownames(track) <- sample.names
print(track)
#Read results of filtering step and append the results of ASV step:
out=read.table(paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), header=TRUE)
print(out)
track=cbind(out, track)

#Write output tables of reads and the otu_table:
write.table(track, paste0(outpath, "06-report/dada2/dada2_stats.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(t(seqtab.nochim), paste0(outpath,"03-dada2/seqtab-nochim.txt"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

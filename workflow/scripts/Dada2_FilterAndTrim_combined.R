#####
# Dada2 Step 1: filter and trim
#
# @author: J. Engelmann
# @author: A. ABdala
# modified by Saara Suominen on the 15.6.2021

#input=snakemake@input[[1]] Not working with the shell command, files are read as {input}
library(yaml)
library(dada2)
#library(Biostrings)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# Add all arguments for dada2 parameters from configfile!
config <- read_yaml(args[2])

# Set the different paths for all the supplied libraries
paths <- args[3:length(args)]
#print(paths)

outpath <- args[1]
#print(outpath)

# Set the different paths for all the supplied libraries
# NOTICE: only forward files given as input
filesForw <- paths
filesRev <- gsub("_1P", "_2P", filesForw)

filesForw_single <- gsub("_1P", "_1U", filesForw)
filesRev_single <- gsub("_1P", "_2U", filesForw)

# Get sample names
#sample.names <- gsub("_1P.fastq.gz", "", basename(filesForw))
#message(paste("Sample", sample.names, "will be analysed", collapse = "\n"))
sample.names <- list()

# Make folder for quality plots, and report
dir.create(paste0(outpath, "03-dada2/quality"), recursive = TRUE)
dir.create(paste0(outpath, "06-report/dada2/"))

allfiles <- list(filesForw, filesRev, filesForw_single, filesRev_single)
names(allfiles) <- c("filesForw", "filesRev", "filesForw_single", "filesRev_single")
files_exist <- list()
filts <- list()
quals <- list()

# loop through 1. forward paired, 2. reverse paired 3. forward single 4. reverse single reads
# First check that files are not empty (cutadapt returns empty files for those that don't pass any filters.)
# Filter and trim does not make empty files that don't pass filters, so make these separately.
message("Analyse different files group separately")

for (i in 1:4) {

  message("Check that files are not empty")

  files_loop <- c()

  for (j in 1:length(allfiles[[i]])) {
    info <- file.info(allfiles[[i]][j])
    if (!is.na(info$size) & info$size > 36) { # 36 is the size of an empty .gz file
      files_loop <- c(files_loop, allfiles[[i]][j])
    }
  }

  if (is.null(files_loop)) {
    files_exist[i] <- list(NULL)
  } else {
    files_exist[[i]] <- files_loop
  }

  if (length(files_exist[[i]] != 0)) {
    message(paste("making quality plots of raw reads:", i, ": (1) forward paired (2) reverse paired (3) forward single and (4) reverse single " , sep = " "))

    sample.names[[i]] <- gsub("_1P.fastq.gz", "", basename(files_exist[[i]]))
    sample.names[[i]] <- gsub("/checked", "", sample.names[[i]])
    sample.names[[i]] <- gsub("_2P.fastq.gz", "", sample.names[[i]])
    sample.names[[i]] <- gsub("_1U.fastq.gz", "", sample.names[[i]])
    sample.names[[i]] <- gsub("_2U.fastq.gz", "", sample.names[[i]])

    message(paste(names(allfiles)[i], "reads of sample", sample.names[[i]], "will be analysed", collapse = "\n"))

    quals[[i]] <- gsub("02-cutadapt/checked/", "03-dada2/quality/", files_exist[[i]])
    quals[[i]] <- gsub(".fastq.gz", ".png", quals[[i]])

    plot_list <- list()
    message("Making quality plots")
    for (j in 1:length(quals[[i]])) {
      dir.create(dirname(quals[[i]][j]), showWarnings = FALSE)
      p <- plotQualityProfile(files_exist[[i]][j])
      ggsave(quals[[i]][j], plot = p, dpi = 150, width = 10, height = 10, units = "cm")
    }

    #Make aggregate plot of all reads
    print(plotQualityProfile(files_exist[[i]], aggregate = T))
    ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_", names(allfiles)[i], ".png"), dpi = 300, width = 10, height = 10, units = "cm")

    #Create path and file names for filtered samples"
    filts[[i]] <- gsub("02-cutadapt/checked", "03-dada2/filtered/", files_exist[[i]])

    #assign names to files
    names(filts[[i]]) <- sample.names[[i]]

  } else {

    message(paste("No reads to process for", names(allfiles)[i], sep=" "))

  }
}

# Run filtering on reads:
# Paired reads together and single reads separately
if (config$meta$sequencing$lib_layout=="Paired") {
  message("Filtering and Trimming paired reads based on parameter set in the config file")

  out <- filterAndTrim(files_exist[[1]], filts[[1]], files_exist[[2]], filts[[2]],
    truncLen = c(config$DADA2$filterAndTrim$Trunc_len_f,config$DADA2$filterAndTrim$Trunc_len_r),
    truncQ = config$DADA2$filterAndTrim$TruncQ,
    trimRight = config$DADA2$filterAndTrim$Trim_right,
    trimLeft = config$DADA2$filterAndTrim$Trim_left,
    maxLen = config$DADA2$filterAndTrim$maxLen,
    minLen = config$DADA2$filterAndTrim$minLen,
    maxN = config$DADA2$filterAndTrim$maxN,
    minQ = config$DADA2$filterAndTrim$minQ,
    maxEE = config$DADA2$filterAndTrim$MaxEE,
    rm.phix = config$DADA2$filterAndTrim$Rm.phix,
    orient.fwd = config$DADA2$filterAndTrim$orient.fwd,
    matchIDs = config$DADA2$filterAndTrim$matchIDs,
    id.sep = config$DADA2$filterAndTrim$id.sep,
    id.field = config$DADA2$filterAndTrim$id.field,
    compress = config$DADA2$filterAndTrim$compress,
    multithread = config$DADA2$filterAndTrim$multithread,
    n = config$DADA2$filterAndTrim$num,
    OMP = config$DADA2$filterAndTrim$OMP,
    verbose = config$DADA2$filterAndTrim$verbose
  )

# Write out to save the effect of filtering on the reads:
rownames(out) <- sample.names[[1]]
out <- as.data.frame(out)
colnames(out)[2] <- "reads.out.paired"
stats_reads <- out
#write.table(out, paste0(outpath, "06-report/dada2/dada2_filtering_stats_paired_reads.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)

qualfiltsFs <- gsub(".png", "_filtered.png", quals[[1]])
qualfiltsRs <- gsub(".png", "_filtered.png", quals[[2]])

# Make quality profile plots.
message("Making quality plots of the filtered reads")
filts_passedF <- c()
filts_passedR <- c()
plot_list <- list()

for (j in 1:length(qualfiltsFs)) {
  if (file.exists(filts[[1]][j])) {
      filts_passedF <- c(filts_passedF, filts[[1]][j])
      dir.create(dirname(qualfiltsFs[j]), showWarnings = FALSE)
      p <- plotQualityProfile(filts[[1]][j])
      ggsave(qualfiltsFs[j], plot = p, dpi = 150, width = 10, height = 10, units = "cm")
  }
  if (file.exists(filts[[2]][j])) {
    filts_passedR <- c(filts_passedR, filts[[2]][j])
    dir.create(dirname(qualfiltsRs[j]), showWarnings = FALSE)
    q <- plotQualityProfile(filts[[2]][j])
    ggsave(qualfiltsRs[j], plot = q, dpi = 150, width = 10, height = 10, units = "cm")
  }
}

print(plotQualityProfile(filts_passedF, aggregate = T))
ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_paired_filtered_forward.png"), dpi = 300, width = 10, height = 10, units = "cm")
print(plotQualityProfile(filts_passedR, aggregate = T))
ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_paired_filtered_reverse.png"), dpi = 300, width = 10, height = 10, units = "cm")

} else {

message("Paired read files will be analysed in single-end mode")

#Append files_exist 1P and files_exist 2P to the single reads so that they are analysed in the same workflow
#Forward reads:
if (length(files_exist[[3]] != 0)) { 
files_exist[[3]]=c(files_exist[[1]], files_exist[[3]])
filts[[3]]=c(filts[[1]], filts[[3]])
quals[[3]]=c(quals[[1]], quals[[3]])
sample.names[[3]]=c(sample.names[[1]], sample.names[[3]])
} else {
files_exist[[3]]=files_exist[[1]]
filts[[3]]=filts[[1]]
quals[[3]]=quals[[1]]
sample.names[[3]]=sample.names[[1]]
}

#Reverse reads:
if (length(files_exist[[4]] != 0)) { 
files_exist[[4]]=c(files_exist[[2]], files_exist[[4]])
filts[[4]]=c(filts[[2]], filts[[4]])
quals[[4]]=c(quals[[2]], quals[[4]])
sample.names[[4]]=c(sample.names[[2]], sample.names[[4]])
} else {
files_exist[[4]]=files_exist[[2]]
filts[[4]]=filts[[2]]
quals[[4]]=quals[[2]]
sample.names[[4]]=sample.names[[2]]
}

} 

# Same for single reads:

message("Filtering and Trimming unpaired forward reads based on parameter set in the config file")

if (length(files_exist[[3]]) != 0) {

  #filts[[3]] <- filts[[3]][file.exists(filts[[3]])]
  out <- filterAndTrim(files_exist[[3]], filts[[3]],
    truncLen = config$DADA2$filterAndTrim$Trunc_len_f,
    truncQ = config$DADA2$filterAndTrim$TruncQ,
    trimRight = config$DADA2$filterAndTrim$Trim_right,
    trimLeft = config$DADA2$filterAndTrim$Trim_left,
    maxLen = config$DADA2$filterAndTrim$maxLen,
    minLen = config$DADA2$filterAndTrim$minLen,
    maxN = config$DADA2$filterAndTrim$maxN,
    minQ = config$DADA2$filterAndTrim$minQ,
    maxEE = config$DADA2$filterAndTrim$MaxEE,
    rm.phix = config$DADA2$filterAndTrim$Rm.phix,
    orient.fwd = config$DADA2$filterAndTrim$orient.fwd,
    matchIDs = config$DADA2$filterAndTrim$matchIDs,
    id.sep = config$DADA2$filterAndTrim$id.sep,
    id.field = config$DADA2$filterAndTrim$id.field,
    compress = config$DADA2$filterAndTrim$compress,
    multithread = config$DADA2$filterAndTrim$multithread,
    n = config$DADA2$filterAndTrim$num,
    OMP = config$DADA2$filterAndTrim$OMP,
    verbose = config$DADA2$filterAndTrim$verbose
  )

  # Write out to save the effect of filtering on the reads:
  rownames(out) <- sample.names[[3]]
  out <- as.data.frame(out)
  colnames(out) <- c("reads.in.forward.single", "reads.out.forward.single")
  
  if (config$meta$sequencing$lib_layout=="Paired") {
  stats_reads <- cbind(stats_reads, out[match(rownames(stats_reads), rownames(out)),])
  #stats_reads$reads.out.forward.single = out$reads.out.forward.single[match(rownames(stats_reads), rownames(out))]
  #write.table(out, paste0(outpath, "06-report/dada2/dada2_filtering_stats_unpaired_forward_reads.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
  } else {
    stats_reads <- out
  }


  qualfiltsFs_single <- gsub(".png", "_filtered.png", quals[[3]])

  # Make quality profile plots.
  message("Making quality plots of the filtered reads")
  filts_passedF <- c()
  plot_list <- list()

  for (j in 1:length(qualfiltsFs_single)) {
    if (file.exists(filts[[3]][j])) {
      filts_passedF <- c(filts_passedF, filts[[3]][j])
      dir.create(dirname(qualfiltsFs_single[j]), showWarnings = FALSE)
      p <- plotQualityProfile(filts[[3]][j])
      ggsave(qualfiltsFs_single[j], plot = p, dpi = 150, width = 10, height = 10, units = "cm")
    }
  }

if (length(filts_passedF)!=0){ 
    print(plotQualityProfile(filts_passedF, aggregate = T))
    ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_filtered_unpaired_forward.png"), dpi = 300, width = 10, height = 10, units = "cm")
  }

}

message("Filtering and Trimming unpaired reverse reads based on parameter set in the config file")

if (length(files_exist[[4]]) != 0) {

  out <- filterAndTrim(files_exist[[4]], filts[[4]],
    truncLen = config$DADA2$filterAndTrim$Trunc_len_r,
    truncQ = config$DADA2$filterAndTrim$TruncQ,
    trimRight = config$DADA2$filterAndTrim$Trim_right,
    trimLeft = config$DADA2$filterAndTrim$Trim_left,
    maxLen = config$DADA2$filterAndTrim$maxLen,
    minLen = config$DADA2$filterAndTrim$minLen,
    maxN = config$DADA2$filterAndTrim$maxN,
    minQ = config$DADA2$filterAndTrim$minQ,
    maxEE = config$DADA2$filterAndTrim$MaxEE,
    rm.phix = config$DADA2$filterAndTrim$Rm.phix,
    orient.fwd = config$DADA2$filterAndTrim$orient.fwd,
    matchIDs = config$DADA2$filterAndTrim$matchIDs,
    id.sep = config$DADA2$filterAndTrim$id.sep,
    id.field = config$DADA2$filterAndTrim$id.field,
    compress = config$DADA2$filterAndTrim$compress,
    multithread = config$DADA2$filterAndTrim$multithread,
    n = config$DADA2$filterAndTrim$num,
    OMP = config$DADA2$filterAndTrim$OMP,
    verbose = config$DADA2$filterAndTrim$verbose
  )

  # Write out to save the effect of filtering on the reads:
  rownames(out) <- sample.names[[4]]
  out <- as.data.frame(out)
  colnames(out) <- c("reads.in.reverse.single", "reads.out.reverse.single")
  stats_reads <- cbind(stats_reads, out[match(rownames(stats_reads), rownames(out)),])
  #stats_reads$reads.out.reverse.single = out$reads.out.reverse.single[match(rownames(stats_reads), rownames(out))]
  #write.table(out, paste0(outpath, "06-report/dada2/dada2_filtering_stats_unpaired_reverse_reads.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)

  qualfiltsRs_single <- gsub(".png", "_filtered.png", quals[[4]])

  # Make quality profile plots.
  message("Making quality plots of the filtered reads")
  filts_passedR <- c()
  plot_list <- list()


  for (j in 1:length(qualfiltsRs_single)) {
    if (file.exists(filts[[4]][j])){
      filts_passedR <- c(filts_passedR, filts[[4]][j])
      dir.create(dirname(qualfiltsRs_single[j]), showWarnings = FALSE)
      p <- plotQualityProfile(filts[[4]][j])
      ggsave(qualfiltsRs_single[j], plot = p, dpi = 150, width = 10, height = 10, units = "cm")
    }
  }

  if (length(filts_passedR)!=0){ 
    print(plotQualityProfile(filts_passedR, aggregate = T))
    ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_filtered_unpaired_reverse.png"), dpi = 300, width = 10, height = 10, units = "cm")
  }
}

write.table(stats_reads, paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)

# Fix for https://github.com/tidyverse/ggplot2/issues/2787
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

# Make empty files for the those samples that did not pass filtering (and no files therefore created)

## There is no error if all reads have been eliminated, but no files
  ## are written in this case. Check for the output files, and if they
  ## don't exist, create empty ones.
 
 for (i in 1:length(filts)) {
  for(fn in filts[[i]]){
    if(!file.exists(fn)){
      cat(gettextf('creating empty file %s\n', fn))
      gzf = gzfile(fn)
      cat('', file=gzf, fill=FALSE)
      close(gzf)
    }}}
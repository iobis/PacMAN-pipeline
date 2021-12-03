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
#Add all arguments for dada2 parameters from configfile!
config=read_yaml(args[2])

#Set the different paths for all the supplied libraries
paths <- args[3:length(args)]
#print(paths)

outpath <- args[1]
#print(outpath)

#Set the different paths for all the supplied libraries
#NOTICE: only forward files given as input
filesForw <- paths
filesRev <- gsub("_1P", "_2P", filesForw)

#Get sample names
sample.names <- gsub("_1P.fastq.gz", "", basename(filesForw))
message(paste("Sample", sample.names, "will be analysed", collapse = "\n"))

qualFs <- gsub("02-cutadapt/", "03-dada2/quality/", filesForw)
qualFs <- gsub(".fastq.gz", ".png", qualFs)
qualRs <- gsub("02-cutadapt/", "03-dada2/quality/", filesRev)
qualRs <- gsub(".fastq.gz", ".png", qualRs)

#print(qualFs)
#print(qualRs)

# Make quality profile plots.

dir.create(paste0(outpath, "03-dada2/quality"), recursive = TRUE)

plot_list_F <- list()
plot_list_R <- list()

for (i in 1:length(qualFs)) {
  dir.create(dirname(qualFs[i]))
  p <- plotQualityProfile(filesForw[i])
  q <- plotQualityProfile(filesRev[i])
  plot_list_F[[i]] <- p
  plot_list_R[[i]] <- q
}

# Save plots to png. Makes a separate file for each plot.
for (i in 1:length(qualFs)) {

  #print(qualFs[i])
  ggsave(qualFs[i], plot = plot_list_F[[i]], dpi = 150, width = 10, height = 10, units = "cm")
  #dev.off()

  #print(qualRs[i])
  ggsave(qualRs[i], plot = plot_list_R[[i]], dpi = 150, width = 10, height = 10, units = "cm")
  #dev.off()
}


dir.create(paste0(outpath,"06-report/dada2/"))

# Make also an aggregate quality plot which can be added to the report
print(plotQualityProfile(filesForw, aggregate = T))
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_forward.png"), dpi = 300, width = 10, height = 10, units = "cm")

print(plotQualityProfile(filesRev, aggregate = T))
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_reverse.png"), dpi = 300, width = 10, height = 10, units = "cm")

#Create path and file names for filtered samples"
filtFs <- gsub("02-cutadapt/", "03-dada2/filtered/", filesForw)
filtRs <- gsub("02-cutadapt/", "03-dada2/filtered/", filesRev)

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(filesForw, filtFs, filesRev, filtRs,
  truncLen = c(config$DADA2$filterAndTrim$Trunc_len_f,config$DADA2$filterAndTrim$Trunc_len_r),
  truncQ=config$DADA2$filterAndTrim$TruncQ,
  trimRight = config$DADA2$filterAndTrim$Trim_right, trimLeft = config$DADA2$filterAndTrim$Trim_left,
  maxLen = config$DADA2$filterAndTrim$maxLen, minLen = config$DADA2$filterAndTrim$minLen,
  maxN = config$DADA2$filterAndTrim$maxN, minQ = config$DADA2$filterAndTrim$minQ,
  maxEE = config$DADA2$filterAndTrim$MaxEE,
  rm.phix = config$DADA2$filterAndTrim$Rm.phix, orient.fwd = config$DADA2$filterAndTrim$orient.fwd,
  matchIDs = config$DADA2$filterAndTrim$matchIDs, id.sep = config$DADA2$filterAndTrim$id.sep,
  id.field = config$DADA2$filterAndTrim$id.field,
  compress = config$DADA2$filterAndTrim$compress, multithread = config$DADA2$filterAndTrim$multithread,
  n = config$DADA2$filterAndTrim$num, OMP = config$DADA2$filterAndTrim$OMP,
  verbose = config$DADA2$filterAndTrim$verbose
)

#Write out to save the effect of filtering on the reads:
rownames(out) <- sample.names
#print(paste0(outpath, "dada2_stats.txt"))
write.table(out, paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


qualfiltFs <- gsub(".png", "_filtered.png", qualFs)
qualfiltRs <- gsub(".png", "_filtered.png", qualRs)

# Make quality profile plots.
plot_list_F <- list()
plot_list_R <- list()

for (i in 1:length(qualfiltFs)) {
  p <- plotQualityProfile(filtFs[i])
  q <- plotQualityProfile(filtRs[i])
  plot_list_F[[i]] <- p
  plot_list_R[[i]] <- q
}

# Save plots to png. Makes a separate file for each plot.
for (i in 1:length(qualfiltFs)) {
  png(qualfiltFs[i])
  print(plot_list_F[[i]])
  dev.off()

  png(qualfiltRs[i])
  print(plot_list_R[[i]])
  dev.off()
}

print(plotQualityProfile(filtFs, aggregate=T))
ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_filtered_forward.png"), dpi = 300, width = 10, height = 10, units = "cm")

print(plotQualityProfile(filtRs, aggregate=T))
ggsave(paste0(outpath, "06-report/dada2/aggregate_quality_profiles_filtered_reverse.png"), dpi = 300, width = 10, height = 10, units = "cm")

# Fix for https://github.com/tidyverse/ggplot2/issues/2787
if (file.exists("Rplots.pdf")) {
  file.remove("Rplots.pdf")
}

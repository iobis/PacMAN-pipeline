
#####
# Dada2 Step 1: filter and trim
#
# @author: J. Engelmann
# @author: A. ABdala
# modified by Saara Suominen on the 15.6.2021

#input=snakemake@input[[1]] Not working with the shell command, files are read as {input}

library(dada2)
#library(Biostrings)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

#Add all arguments for dada2 parameters from configfile!
#args[1]... = outpath for dada2 files (project/run/dada2/)
truncLen_F =  as.numeric(args[2])
truncLen_R =  as.numeric(args[3])
truncQ =      as.numeric(args[4])
trimRight=    as.numeric(args[5])
trimLeft=     as.numeric(args[6])
maxLen=       as.numeric(args[7])
minLen=       as.numeric(args[8])
maxN=         as.numeric(args[9])
minQ=         as.numeric(args[10])
maxEE=        as.numeric(args[11])
rm.phix=      as.logical(args[12])
if(args[13]=="None"){ orient.fwd=NULL } else { orient.fwd=args[13] }
matchIDs=     as.logical(args[14])
id.sep=args[15]
if(args[16]=="None"){ id.field=NULL } else { id.field=args[16] }
compress=     as.logical(args[17])
#args[18]...= multithread]} \
n=            as.numeric(args[19])
OMP=          as.logical(args[20])
verbose=      as.logical(args[21])
#args[22:n]...= input.files}

#Multithread can either be logical or a numeric: probably will be False and True..
if (args[18]=="FALSE"|args[18]=="TRUE"|args[18]=="T"|args[18]=="F"|args[18]=="False"|args[18]=="True"){
  multithread=as.logical(args[18])
} else {
  multithread=as.numeric(args[18])
}

#print(args[1:22])

#Set the different paths for all the supplied libraries
paths=args[22:length(args)]
#print(paths)

outpath <- args[1]
#print(outpath)

#List files
#filesForw <- sort(list.files(paths, pattern="_1P.fastq.gz", full.names = TRUE))
#filesRev <- sort(list.files(paths, pattern="_2P.fastq.gz", full.names = TRUE))
#filesForw <- gsub("//", "/", filesForw)
#filesRev <- gsub("//", "/", filesRev)
#print(filesForw)
#print(filesRev)

#Set the different paths for all the supplied libraries
#NOTICE: only forward files given as input
filesForw=paths
#print(filtFs)
filesRev=gsub('_1P','_2P',filesForw)
#print(filtRs)

#Get sample names
sample.names <- gsub('_1P.fastq.gz', '', basename(filesForw))
print(paste("Sample", sample.names, "will be analysed"))

qualFs <- gsub('02-cutadapt/','03-dada2/quality/',filesForw)
qualFs <- gsub('.fastq.gz','.png',qualFs)
qualRs <- gsub('02-cutadapt/','03-dada2/quality/',filesRev)
qualRs <- gsub('.fastq.gz','.png',qualRs)

#print(qualFs)
#print(qualRs)

# Make quality profile plots.

dir.create(paste0(outpath, "03-dada2/quality"))
plot_list_F = list()
plot_list_R = list()
for (i in 1:length(qualFs)) {
    dir.create(dirname(qualFs[i]))
    p = plotQualityProfile(filesForw[i])
    q = plotQualityProfile(filesRev[i])
    plot_list_F[[i]] = p
    plot_list_R[[i]] = q
}

# Save plots to png. Makes a separate file for each plot.
for (i in 1:length(qualFs)) {

    #print(qualFs[i])
    ggsave(qualFs[i], plot=plot_list_F[[i]], dpi=150, width=10, height=10, units="cm")
    #dev.off()

    #print(qualRs[i])
    ggsave(qualRs[i], plot=plot_list_R[[i]], dpi=150, width=10, height=10, units="cm")
    #dev.off()
}


dir.create(paste0(outpath,"06-report/dada2/"))

# Make also an aggregate quality plot which can be added to the report
print(plotQualityProfile(filesForw, aggregate=T))
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_forward.png"), dpi=300, width=10, height=10, units="cm")

print(plotQualityProfile(filesRev, aggregate=T))
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_reverse.png"), dpi=300, width=10, height=10, units="cm")



#Create path and file names for filtered samples"
filtFs <- gsub('02-cutadapt/','03-dada2/filtered/',filesForw)
filtRs <- gsub('02-cutadapt/','03-dada2/filtered/',filesRev)

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(filesForw, filtFs, filesRev, filtRs,
  truncLen=c(truncLen_F,truncLen_R),truncQ=truncQ,
  trimRight=trimRight, trimLeft=trimLeft,
  maxLen=maxLen, minLen=minLen,
  maxN=maxN, minQ=minQ, maxEE=maxEE,
  rm.phix=rm.phix, orient.fwd=orient.fwd,
  matchIDs=matchIDs, id.sep=id.sep, id.field=id.field,
  compress=compress, multithread=multithread,
  n=n, OMP=OMP,
  verbose=verbose)

#Write out to save the effect of filtering on the reads:
rownames(out) <- sample.names
#print(paste0(outpath, "dada2_stats.txt"))
write.table(out, paste0(outpath, "06-report/dada2/dada2_filtering_stats.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE)


qualfiltFs <- gsub('.png','_filtered.png',qualFs)
qualfiltRs <- gsub('.png','_filtered.png',qualRs)

# Make quality profile plots.
plot_list_F = list()
plot_list_R = list()

for (i in 1:length(qualfiltFs)) {
    p = plotQualityProfile(filtFs[i])
    q = plotQualityProfile(filtRs[i])
    plot_list_F[[i]] = p
    plot_list_R[[i]] = q
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
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_filtered_forward.png"), dpi=300, width=10, height=10, units="cm")

print(plotQualityProfile(filtRs, aggregate=T))
ggsave(paste0(outpath,"06-report/dada2/aggregate_quality_profiles_filtered_reverse.png"), dpi=300, width=10, height=10, units="cm")

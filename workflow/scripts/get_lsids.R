# Written by Saara Suominen (saara.suominen.work@gmail.com)
# for OBIS and PacMAN
#Last update 21.6.2021

library('worrms')
library('stringr')
library('Biostrings')

args <- commandArgs(trailingOnly = T)

outpath=args[1]

########################### 1. Read input files  ################################################################################################################

Tax_file=read.csv(args[2], sep='\t', header=T)
Rep_seqs=Biostrings::readDNAStringSet(args[3])

########################### 2. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package) ######################################

print("1. Modify tax table for taxonomic ranks, and find all possible worms ids (using worrms package)")

levels = c("kingdom","phylum", "class", "order", "family", "genus", "species")
taxmat=str_split(Tax_file$sum.taxonomy, ";", simplify =T)
colnames(taxmat)=levels
rownames(taxmat)=Tax_file$rowname

#There are some missing values in between which means that taxa are not in the right category
#The following loops corrects this based on the first letter of the resuts:

#It would be better to write a function here, or correct this in the previous step (why are there empty values in between?):
#match_level=function(taxmat){}

for(i in 1:dim(taxmat)[1]) {
    if ((!grepl("p__", taxmat[i,"phylum"]))&taxmat[i,"phylum"]!="") {
        for(j in dim(taxmat)[2]:2) {
              taxmat[i,j]<-taxmat[i,j-1]
      }
      taxmat[i,"phylum"]=NA
   }
   if ((!grepl("c__", taxmat[i,"class"]))&taxmat[i,"class"]!="") {
        for(j in dim(taxmat)[2]:3) {
             taxmat[i,j]<-taxmat[i,j-1]
        }
      taxmat[i,"class"]=NA
   }
   if ((!grepl("o__", taxmat[i,"order"]))&taxmat[i,"order"]!="") {
        for(j in dim(taxmat)[2]:4) {
             taxmat[i,j]<-taxmat[i,j-1]
        }
      taxmat[i,"order"]=NA
   }
   if ((!grepl("f__", taxmat[i,"family"]))&taxmat[i,"family"]!="") {
        for(j in dim(taxmat)[2]:5) {
             taxmat[i,j]<-taxmat[i,j-1]
        }
      taxmat[i,"family"]=NA
   }
    if ((!grepl("g__", taxmat[i,"genus"]))&taxmat[i,"genus"]!="") {
        for(j in dim(taxmat)[2]:6) {
             taxmat[i,j]<-taxmat[i,j-1]
        }
      taxmat[i,"genus"]=NA
   }
  }

#First remove extra characters from strings
taxmat=gsub(".*[__]([^.]+)[_].*", "\\1", taxmat)
#Replace all empty strings with na
taxmat[taxmat==""]=NA
taxmat[taxmat=="NA"]=NA

#Collect the highest known taxonomic value to the last column
#Requires some wrangling, as the lastvalue is returned as a list:
lastValue <- function(x)  { tail(x[!is.na(x)], 1)}
lastvalue=apply(taxmat, 1, lastValue)
taxmat=as.data.frame(taxmat)
taxmat$lastvalue=NA

for(i in 1:length(lastvalue)) {
  if(length(unlist(lastvalue[i]))>0) {
    taxmat$lastvalue[i]=unlist(lastvalue[i])
  } else {
    taxmat$lastvalue[i]=NA
  }}

#Because WORMS doesn't recognize Eukaryota, change those that have this in the lastvalue to Biota:
taxmat[taxmat$lastvalue=='Eukaryota'&!is.na(taxmat$lastvalue), 'lastvalue']='Biota'

#Find worrms value with worrms
#Function from Pieter, but I removed the exact match requirement so be careful with this.
#I removed it, because I was getting 'like' matches even for exact matches, but I am not sure why.
match_name <- function(name) {
  lsid <- tryCatch({
    #res <- wm_records_names(name) #Note, this one not returning all, i.e. Peronospora, what is the difference to the other command?
    res <- wm_records_name(name)
    matches <- res[[1]] #%>%
      #filter(match_type == "exact" | match_type == "exact_genus") # Even for the exact species name I get 'like' matches, what is the issue?
    if (nrow(matches) > 1) {
      message(paste0("Multiple matches for ", name))
    }
    return(matches$lsid[1])
  }, error = function(cond) {
    return(NA)
  })
}

#Use loop to fill in table with the matches. Easy way to skip the NA's (wm_records_names, takes NAs as characters)
#For each NA add instead the ID for 'biota', so that all of these are kept for the data submission
for(i in 1:length(lastvalue)){
  if(!is.na(taxmat$lastvalue[i])) {
    taxmat$lsid[i]=match_name(taxmat$lastvalue[i])
  } else {
    taxmat$lsid[i]='urn:lsid:marinespecies.org:taxname:1'
  }}

#The loop returns NAs for those names that are not found in WORMS
names_not_in_worms=taxmat[is.na(taxmat$lsid),'lastvalue']
print(paste("Number of species names not recognized in WORMS: ", nrow(unique(names_not_in_worms))))

#Some of these are common contaminants (e.g. Homo sapiens, Canis lupus)
#Some of these are not marine species (e.g. insects), that could be found with terrestrial matches
#Some have a higher taxonomy that can be found in the database
#Check at genus level:
for(i in 1:length(lastvalue)){
  if(is.na(taxmat$lsid[i])) {
    if(!is.na(taxmat$genus[i])) {
      taxmat$lsid[i]=match_name(taxmat$genus[i])
    }}}

#This way unknowns decreases. However, these values should be checked for possible addition to WoRMS?
#Will require manual inspection
not_in_worms=taxmat[is.na(taxmat$lsid),]
print(paste("Number of genus names not recognized in WORMS: ", nrow(unique(not_in_worms))))
print(paste("These taxa can be found in the table: ", paste0(outpath,"Taxa_not_in_worms.csv")))

#We have the dilemma of keeping sequences that are completely unknown (as 'biota'),
#while removing sequences that have a known non-marine origin?

#Add 'Biota' as the name for the unknown sequences
taxmat[is.na(taxmat$lastvalue), 'lastvalue']="Biota"

#Add sequence to the tax_table slot (linked to each asv)
DNA_sequence=as.character(Rep_seqs[rownames(taxmat)])
taxmat=cbind(taxmat, DNA_sequence)

#Return table of unknown names to make manual inspection easier:
write.table(not_in_worms, paste0(outpath,"Taxa_not_in_worms.csv"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
#Return tax table
write.table(taxmat, paste0(outpath,"Full_tax_table_with_lsids.csv"), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

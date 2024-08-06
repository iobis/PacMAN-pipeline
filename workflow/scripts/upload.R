library(dplyr)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
data_path <- args
stopifnot(length(data_path) > 0)

occurrence_files <- list.files(path = data_path, pattern = "Occurrence_table.tsv", recursive = TRUE, full.names = TRUE)
dna_files <- list.files(path = data_path, pattern = "DNA_extension_table.tsv", recursive = TRUE, full.names = TRUE)

occurrence <- purrr::map(occurrence_files, function(file) {
    read.table(file, sep = "\t", header = TRUE, colClasses = "character")
}) %>% bind_rows()

dna <- purrr::map(dna_files, function(file) {
    read.table(file, sep = "\t", header = TRUE, colClasses = "character")
}) %>% bind_rows()

stopifnot(!any(duplicated(occurrence$occurrenceID)))
stopifnot(!any(duplicated(dna$occurrenceID)))

tmp <- tempdir()
occurrence_path <- file.path(tmp, "occurrence.txt")
dna_path <- file.path(tmp, "dna.txt")
write.table(occurrence, file = occurrence_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(dna, file = dna_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

system(glue("aws s3 cp {occurrence_path} s3://pacman-dataset/edna/occurrence.txt"))
system(glue("aws s3 cp {dna_path} s3://pacman-dataset/edna/dna.txt"))

library(dplyr)
library(tidyr)
library(yaml)
library(stringr)
library(glue)
library(Biostrings)
source("workflow/scripts/util.R")

max_vsearch_results <- 5

if (!exists("cmd_args")) {
  cmd_args <- commandArgs(trailingOnly = T)
  message("cmd_args <- ", capture.output(dput(cmd_args)))
}

config <- read_yaml(cmd_args[2])
outpath <- cmd_args[1]
rdp_file_path <- cmd_args[3]
vsearch_file_path <- cmd_args[4]
vsearch_lca_file_path <- cmd_args[5]
rep_seqs_path <- cmd_args[6]

rdp_confidence_threshold <- config$Rdp$cutoff
vsearch_identity_threshold <- config$Vsearch$pident
vsearch_cover_threshold <- config$Vsearch$query_cov

rdp <- read.csv(rdp_file_path, sep = "\t", header = F)
vsearch <- read.csv(vsearch_file_path, sep = "\t", header = F)
vsearch_lca <- read.csv(vsearch_lca_file_path, sep = "\t", header = F)

dna <- readDNAStringSet(rep_seqs_path)
seq_lengths <- data.frame(asv = names(dna), seq_length = width(dna))

# process vsearch results

colnames(vsearch) <- c("asv", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

vsearch <- vsearch %>%
  separate(sseqid, sep=";", into = c("vsearch_seqid", "vsearch_taxonomy")) %>%
  mutate(pident = as.numeric(pident) / 100) %>%
  left_join(seq_lengths, by = "asv") %>%
  mutate(
    query_cover = (qend - qstart) / seq_length
  ) %>%
  select(asv, taxonomy = vsearch_taxonomy, seqid = vsearch_seqid, identity = pident, query_cover) %>%
  filter(!is.na(taxonomy))

split_taxonomy <- function(taxonomy_string) {
  pairs <- str_split(taxonomy_string, ",")[[1]]
  pairs_list <- str_split(pairs, ":")
  keys <- sapply(pairs_list, `[`, 1)
  values <- sapply(pairs_list, `[`, 2)
  names(values) <- keys
  return(values)
}

vsearch_clean <- vsearch %>%
  mutate(
    taxonomy = str_replace(taxonomy, "tax=", ""),
    taxonomy = str_replace(taxonomy, "_", " "),
  ) %>%
  rowwise() %>%
  mutate(parsed = list(split_taxonomy(taxonomy))) %>%
  unnest_wider(parsed, names_sep = "_") %>%
  ungroup() %>%
  dplyr::rename(
    domain = parsed_d,
    phylum = parsed_p,
    class = parsed_c,
    order = parsed_o,
    family = parsed_f,
    genus = parsed_g,
    species = parsed_s
  ) %>%
  select(-taxonomy) %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  group_by(asv) %>%
  slice_max(identity, n = 5, with_ties = FALSE) %>%
  rowwise() %>%
  mutate(scientificName = coalesce(species, genus, family, order, class, phylum)) %>%
  ungroup()

tax_names <- vsearch_clean$scientificName %>% na.omit() %>% unique() %>% sort()

message("Matching names")

matches <- match_worms(unique(tax_names)) %>%
  select(scientificName = input, scientificNameID = lsid)
vsearch_clean <- vsearch_clean %>%
  left_join(matches, by = "scientificName") %>%
  mutate(
    method = "VSEARCH",
    vsearch_identity_threshold = vsearch_identity_threshold,
    vsearch_cover_threshold = vsearch_cover_threshold
  )

# process vsearch consensus results?

# process rdp results
# TODO: check terrimporter taxonomy format

process_rdp_taxonomy <- function(row) {
  taxa <- as.character(row[seq(6, ncol(row), 3)])
  ranks <- as.character(row[seq(7, ncol(row), 3)])
  confidences <- as.numeric(row[seq(8, ncol(row), 3)])

  df <- data.frame(scientificName = taxa, rank = ranks, confidence = confidences) %>%
    arrange(desc(confidence)) %>%
    filter(confidence >= rdp_confidence_threshold) %>%
    mutate(across(where(is.character), ~str_replace_all(.x, "^([a-z]_)+.*", ""))) %>%
    mutate(across(where(is.character), ~str_replace_all(.x, "_", " "))) %>%
    mutate(across(where(is.character), ~na_if(.x, ""))) %>%
    filter(!is.na(scientificName))

  last_confidence <- tail(df, 1)$confidence
  last_scientificName <- tail(df, 1)$scientificName

  df <- df %>%
    select(-confidence) %>%
    mutate(rank = tolower(rank)) %>%
    pivot_wider(names_from = rank, values_from = scientificName) %>%
    mutate(
      asv = as.character(row[,1]),
      scientificName = last_scientificName,
      confidence = last_confidence
    )

  return(df)
}

rdp_clean <- rdp %>%
  rowwise() %>%
  group_split() %>%
  map(process_rdp_taxonomy) %>%
  bind_rows()

tax_names <- rdp_clean$scientificName %>% na.omit() %>% unique() %>% sort()
matches <- match_worms(unique(tax_names)) %>%
  select(scientificName = input, scientificNameID = lsid)
rdp_clean <- rdp_clean %>%
  left_join(matches, by = "scientificName") %>%
  mutate(
    method = "RDP classifier",
    rdp_confidence_threshold = rdp_confidence_threshold
  )

# combine vsearch and rdp results

combined <- bind_rows(rdp_clean, vsearch_clean) %>%
  arrange(asv, method) %>%
  relocate(asv, method, scientificName, confidence, identity)

unknowns <- combined %>% filter(method == "RDP classifier" & is.na(phylum)) %>% pull(asv)

write.table(combined, paste0(outpath, "04-taxonomy/annotation_results.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t", na = "")
writeLines(unknowns, paste0(outpath, "04-taxonomy/unknown_asvs.txt"))

######### Combine info from vsearch to rdp ###############

#Build output of asv, sum.taxonomy, and identificationRemarks
#sum.taxonomy is the taxonomic information separated by ;
#Identification remarks, should contain both the taxonomy and the confidence from rdp,
#as well as the results of the vsearch high confidence search. 
#Return also those that are fully unknown from rdp for a possible blast step. 
               
# results <- rdp_mod %>%
#   left_join(vsearch_combined, by = "asv") %>%
#   mutate(
#     rdp = glue("Identification based on the RDP classifier at the confidence level {config$Rdp$cutoff}: taxonomy {rdp_taxonomy}, confidences {rdp_confidences}."),
#     vsearch = glue("Confirmation with VSEARCH against the {config$Vsearch$vsearch_db_name} database at {config$Vsearch$pident} similarity: hits {vsearch_seqid}, identities {vsearch_pident}, taxonomy {vsearch_taxonomy}, consensus {vsearch_consensus}.")
#   ) %>%
#   mutate(
#     identificationRemarks = ifelse(!is.na(vsearch_pident), glue("{rdp} {vsearch}"), glue("{rdp} No VSEARCH hits at {config$Vsearch$pident} identity."))
#   ) %>%
#   select(asv, sum.taxonomy = rdp_consensus, identificationRemarks)

# unknowns <- rdp_mod$asv[rdp_mod$rdp_consensus=="" | rdp_mod$rdp_consensus=="Eukaryota" | rdp_mod$rdp_consensus=="sk__Eukaryota"]
  
# write.table(results, paste0(outpath, "04-taxonomy/rdp_vsearch_results.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep="\t")
# writeLines(unknowns, paste0(outpath, "04-taxonomy/unknown_asvs.txt"))

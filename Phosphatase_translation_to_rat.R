library(tidyverse)
library(readxl)
library(biomaRt)

# Load DEPOD interaction data
depod_data <- readxl::read_excel("PPase_protSubtrates_201903.xls")

# ---------------------------------------------------------------------------- #

# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# 
# devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")

library(MotrpacRatTraining6moData)
# data("RAT_TO_HUMAN_GENE")

head(RAT_TO_HUMAN_GENE)

RAT_TO_HUMAN_GENE %>%
  summarise(
    rat_genes = sum(!is.na(RAT_SYMBOL)),
    human_genes = sum(!is.na(HUMAN_ORTHOLOG_SYMBOL)),
    distinct_rat_genes = n_distinct(RAT_SYMBOL),
    distinct_human_genes = n_distinct(HUMAN_ORTHOLOG_SYMBOL)
  )

depod_data %>%
  summarise(
    distinct_human_genes = n_distinct(`Substrate entry names`),
    distinct_human_phosphatases = n_distinct(`Phosphatase entry names`)
  )

dat_RAT_TO_HUMAN_GENE <- RAT_TO_HUMAN_GENE %>%
  dplyr::select(HUMAN_ORTHOLOG_SYMBOL, RAT_SYMBOL)

depod_with_genes <- depod_data %>%
  left_join(
    dat_RAT_TO_HUMAN_GENE, 
    by = c("Substrate entry names" = "HUMAN_ORTHOLOG_SYMBOL"),
    relationship = "many-to-many"
    )

dat_RAT_TO_HUMAN_GENE <- RAT_TO_HUMAN_GENE %>%
  dplyr::select(HUMAN_ORTHOLOG_SYMBOL, RAT_SYMBOL) %>%
  rename(RAT_SYMBOL_phosp = RAT_SYMBOL)

depod_with_genes <- depod_with_genes %>%
  left_join(
    dat_RAT_TO_HUMAN_GENE, 
    by = c("Phosphatase entry names" = "HUMAN_ORTHOLOG_SYMBOL"),
    relationship = "many-to-many"
  )

# # Connect to Ensembl database
# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
# 
# # Get human-to-rat homologs
# human_genes <- c("TP53", "BRCA1", "EGFR") # Replace with your HGNC symbols
# homologs <- getLDS(
#   attributes = c("hgnc_symbol"),
#   filters = "hgnc_symbol",
#   values = human_genes,
#   mart = human_mart,
#   attributesL = c("external_gene_name"),
#   martL = rat_mart
# )
# 
# # View results
# print(homologs)

# ---------------------------------------------------------------------------- #

# library(gprofiler2)
# 
# human_genes <- c("TP53", "BRCA1", "EGFR") # Replace with your genes
# 
# result <- gconvert(
#   query = human_genes,
#   organism = "hsapiens", # Source organism
#   target = "rnorvegicus" # Target organism
# )
# 
# validated_genes <- gconvert(
#   query = human_genes,
#   organism = "hsapiens",
#   target = NULL
# )
# 
# print(validated_genes)

# ---------------------------------------------------------------------------- #

# # Map protein accessions to gene names using biomaRt
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# accession_to_gene <- getBM(
#   attributes = c("uniprotswissprot", "hgnc_symbol"),
#   filters = "uniprotswissprot",
#   values = unique(c(depod_data$`Phosphatase accession numbers`, depod_data$`Substrate accession numbers`)),
#   mart = mart
# ) %>%
#   tibble::as_tibble()
# 
# # Join dataframes
# depod_with_genes <- depod_data %>%
#   # left_join(
#   #   accession_to_gene,
#   #   by = c("Substrate accession numbers" = "uniprotswissprot"),
#   #   relationship = "many-to-many"
#   # ) %>%
#   dplyr::select(`Phosphatase entry names`, `Substrate entry names`, Dephosphosites)

# ---------------------------------------------------------------------------- #

readr::write_csv(depod_with_genes, "depod_with_gene_names.csv")

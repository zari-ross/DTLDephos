library(tidyverse)
library(readxl)
library(biomaRt)

# Load DEPOD interaction data
depod_data <- readxl::read_excel("PPase_protSubtrates_201903.xls")

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

# Map protein accessions to gene names using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

accession_to_gene <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = unique(c(depod_data$`Phosphatase accession numbers`, depod_data$`Substrate accession numbers`)),
  mart = mart
) %>%
  tibble::as_tibble()

# Join dataframes
depod_with_genes <- depod_data %>%
  # left_join(
  #   accession_to_gene,
  #   by = c("Substrate accession numbers" = "uniprotswissprot"),
  #   relationship = "many-to-many"
  # ) %>%
  dplyr::select(`Phosphatase entry names`, `Substrate entry names`, Dephosphosites)

readr::write_csv(depod_with_genes, "depod_with_gene_names.csv")

library(tidyverse)
library(readxl)
library(biomaRt)

# Load DEPOD interaction data
depod_data <- readxl::read_excel("PPase_protSubtrates_201903.xls")

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

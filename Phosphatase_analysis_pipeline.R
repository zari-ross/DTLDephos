library(readr)
library(tidyverse)

path_to_files <- "C:/Users/NITru/OneDrive/Documents/PhD_work/GitHub/Meta-Proteomics/arsenite_phosph/"

dat_peptides <- read_csv(paste0(path_to_files, "lfq.peptides.csv"))
dat_protein_peptides <- read_csv(paste0(path_to_files, "lfq.protein-peptides.csv"))
dat_proteins <- read_csv(paste0(path_to_files, "lfq.proteins.csv"))

# ---------------------------------------------------------------------------- #

dat_proc <- dat_protein_peptides %>%
  dplyr::filter(USED == "Y") %>%
  dplyr::select(Start, End, Peptide, `Protein Accession`) %>%
  dplyr::inner_join(dat_proteins %>% dplyr::select(Accession, Description), by = c(`Protein Accession` = "Accession")) %>%
  dplyr::filter(grepl("57", Peptide)) %>%
  dplyr::mutate(
    Phosphopeptide = Peptide,
    Phosphopeptide = gsub("\\(\\+57\\.02\\)", "", Phosphopeptide),
    Phosphopeptide = gsub("\\(\\+15\\.99\\)", "", Phosphopeptide),
    Phosphopeptide = gsub("\\(\\+42\\.01\\)", "", Phosphopeptide),
    Gene = sub(" .*", "", sub(".*GN=", "", Description)),
    Gene = gsub("_", "-", Gene),
    Full_description = Description,
    Description = sub("OS.*", "", Description)
  ) %>%
  dplyr::mutate(`Protein Accession` = sub("\\|.*", "", `Protein Accession`))

# ---------------------------------------------------------------------------- #

# install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)

dat_proteome <- Biostrings::readAAStringSet("C:/Users/NITru/OneDrive/Documents/PhD_work/Projects/NT_projects/ProteinDatabaseCleanup/RatProteomeCleanup/20231112_rat_proteome_dedupl_aCas.fasta")

dat_proteome <- data.frame(
  Accession = sub("^[^|]*\\|([^|]*)\\|.*", "\\1", names(dat_proteome)),
  Sequence = as.character(dat_proteome)
)

dat_proc <- dat_proc %>%
  dplyr::left_join(dat_proteome, by = c(`Protein Accession` = "Accession"))

# ---------------------------------------------------------------------------- #

dat_proc <- dat_proc %>%
  dplyr::mutate(
    Phosphosites = stringr::str_locate_all(Phosphopeptide, "\\(\\+79\\.97\\)"),
    Base_peptide = gsub("\\(\\+79\\.97\\)", "", Phosphopeptide)
  ) %>%
  tidyr::unnest_longer(Phosphosites) %>%
  dplyr::mutate(
    Phosphopeptide = purrr::map2_chr(
      Base_peptide, Phosphosites[, 1],
      ~ paste0(substr(.x, 1, .y - 1), "(+79.97)", substr(.x, .y, nchar(.x)))
    ),
    Position = Start + Phosphosites[, 1] - 1,
    Phosphomark = paste0(substr(Phosphopeptide, Phosphosites[, 1] - 1, Phosphosites[, 1] - 1), "", Position)
  ) %>%
  dplyr::select(-Phosphosites, -Base_peptide)

# ---------------------------------------------------------------------------- #

dat_proc <- dat_proc %>%
  dplyr::mutate(
    Residue_window = purrr::map2_chr(
      Sequence, Position,
      ~ substr(.x, max(.y - 15, 1), min(.y + 15, nchar(.x)))
    )
  )

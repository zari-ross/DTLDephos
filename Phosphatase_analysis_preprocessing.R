library(readr)
library(readxl)
library(tidyverse)
library(Biostrings)

source("config.R")

window_size <- 31
half_window <- (window_size - 1) / 2

# ---------------------------------------------------------------------------- #

dat_orig <- readxl::read_excel(path_to_file)

dat_preproc <- dat_orig %>%
  dplyr::filter(Change == "Decreased") %>%
  dplyr::mutate(Accession = sub("\\|.*", "", Accession))

# ---------------------------------------------------------------------------- #

dat_proteome <- Biostrings::readAAStringSet(path_to_proteome)

dat_proteome <- data.frame(
  Accession = sub("^[^|]*\\|([^|]*)\\|.*", "\\1", names(dat_proteome)),
  `Fasta name` = sub(" .*", "", names(dat_proteome)),
  check.names = FALSE,
  Sequence = as.character(dat_proteome)
)

dat_proc <- dat_preproc %>%
  dplyr::left_join(dat_proteome, by = "Accession") %>%
  dplyr::mutate(Position = as.numeric(sub("^[A-Z]", "", Phosphosite)))

# ---------------------------------------------------------------------------- #

dat_proc <- dat_proc %>%
  dplyr::mutate(
    Residue_window = purrr::map2_chr(
      Sequence, Position,
      ~ substr(.x, max(.y - half_window, 1), min(.y + half_window, nchar(.x)))
    )
  ) %>%
  dplyr::filter(nchar(Residue_window) == window_size) %>%
  dplyr::mutate(
    Start = Position - half_window,
    End = Position + half_window,
    `Fasta name` = paste0(`Fasta name`, "%", Start, "%", End)
  )

# ---------------------------------------------------------------------------- #

fasta_lines_y <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "Y")) %>%
  dplyr::mutate(
    Fasta_line = paste0(">", `Fasta name`, "\n", Residue_window)
  ) %>%
  dplyr::pull(Fasta_line)

fasta_lines_st <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "[ST]")) %>%
  dplyr::mutate(
    Fasta_line = paste0(">", `Fasta name`, "\n", Residue_window)
  ) %>%
  dplyr::pull(Fasta_line)

writeLines(fasta_lines_y, "output_y.fasta")
writeLines(fasta_lines_st, "output_st.fasta")

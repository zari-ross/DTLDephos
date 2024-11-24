library(readr)
library(readxl)
library(tidyverse)

source("config.R")

prob_cut_off <- 0.5

# ---------------------------------------------------------------------------- #

# Load preprocessed dataset
dat_proc <- readr::read_csv(paste0("dat_preproc/", condition, "_full_table.csv"))

# Load model outputs
output_y <- readr::read_csv(paste0("model_output/", condition, "_Y.csv")) %>%
  dplyr::mutate(`Sequence ID` = paste0(">", `Sequence ID`)) %>%
  dplyr::distinct(`Sequence ID`, .keep_all = TRUE)
output_st <- readr::read_csv(paste0("model_output/", condition, "_ST.csv")) %>%
  dplyr::mutate(`Sequence ID` = paste0(">", `Sequence ID`)) %>%
  dplyr::distinct(`Sequence ID`, .keep_all = TRUE)  #dplyr::distinct(.keep_all = TRUE)

# Sometimes, the score is a bit different
# your_data <- readr::read_csv(paste0("model_output/", condition, "_ST.csv")) %>%
#   dplyr::mutate(`Sequence ID` = paste0(">", `Sequence ID`)) 
# distinct_sequence_id <- your_data %>%
#   distinct(`Sequence ID`, .keep_all = TRUE)
# distinct_all <- your_data %>%
#   distinct(.keep_all = TRUE)
# difference_rows <- anti_join(distinct_all, distinct_sequence_id)
# your_data %>%
#   filter(stringr::str_detect(`Sequence ID`, "TAU_RAT")) %>% View()

# Add dephosphorylation probabilities to the dataset
dat_proc_y <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "Y")) %>%
  dplyr::left_join(output_y, by = c("Fasta name" = "Sequence ID"))

dat_proc_st <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "[ST]")) %>%
  dplyr::left_join(output_st, by = c("Fasta name" = "Sequence ID"))

# Combine results into one table
dat_proc_combined <- dplyr::bind_rows(dat_proc_y, dat_proc_st) %>%
  dplyr::mutate(Dephosphorylation = ifelse(`Dephosphorylation Probability` > prob_cut_off, "yes", "no"))

# Save combined table
readr::write_csv(dat_proc_combined, paste0("dat_preproc/", condition, "_processed_with_predictions.csv"))

# ---------------------------------------------------------------------------- #

dat_orig <- readxl::read_excel(path_to_file) %>%
  dplyr::left_join(dat_proc_combined %>% select(Phosphosite, Dephosphorylation, Change), 
                   by = c("Change", "Phosphosite"))

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

cols <- c(
  "Decreased" = col_control,
  "NS" = "gray",
  "Increased" = col_exp
)

norm_opt <- "notNORM"
pValue_cutoff <- 0.05
Sign_cutoff <- -10*log10(pValue_cutoff)
FC_cutoff <- 2

dat_plot <- dat_orig

if (norm_opt != "norm") {
  names(dat_plot)[names(dat_plot) == "Significance"] <- "Significance_new"
}

dat_plot <- dat_orig %>%
  dplyr::mutate(
    logFC = as.numeric(logFC),
    Significance_new = as.numeric(Significance),
    Significance_new = ifelse(Significance_new >= 60, Inf, Significance_new),
    logFC = ifelse(logFC >= 6, Inf, ifelse(logFC <= -6, -Inf, logFC))
  )

plot_volcano <- ggplot(dat_plot, aes(x = logFC, y = Significance_new, colour = Change,
                                     text = paste("Gene:", Gene, "\n", "Description:", Description))) +
  geom_point(size = 0.7, alpha = 0.4) +
  scale_colour_manual(values = cols) +
  
  geom_hline(yintercept = Sign_cutoff, colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = log2(FC_cutoff), colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = -log2(FC_cutoff), colour = "gray", linetype = "dashed") +
  
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(color = "black", size = 10),
    axis.ticks = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  geom_point(
    data = dat_plot %>% filter(stringr::str_starts(Phosphosite, "[Y]") & Dephosphorylation == "yes"),
    aes(x = logFC, y = Significance_new, fill = "Y residues"),
    size = 1.2,
    shape = 21,
    alpha = 0.8
  ) +
  geom_point(
    data = dat_plot %>% filter(stringr::str_starts(Phosphosite, "[ST]") & Dephosphorylation == "yes"),
    aes(x = logFC, y = Significance_new, fill = "ST residues"),
    size = 1.2,
    shape = 21,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c("Y residues" = "blue", "ST residues" = "red"),
    name = paste("Residue Type\nDephosphorylation cut-off", prob_cut_off)
  )

plot_volcano

ggsave(plot_volcano, filename = paste("dat_proc/plot_volcano_", norm_opt, ".png", sep=""), width = 8, height = 6)


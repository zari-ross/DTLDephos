library(readr)
library(readxl)
library(tidyverse)

source("config.R")

# ---------------------------------------------------------------------------- #

# Load preprocessed dataset
dat_proc <- readr::read_csv(paste0("dat_preproc/", condition, "_full_table.csv"))

# Load model outputs
output_y <- readr::read_csv(paste0("model_output/", condition, "_Y.csv")) %>%
  dplyr::mutate(`Sequence ID` = paste0(">", `Sequence ID`)) %>%
  dplyr::distinct(`Sequence ID`, .keep_all = TRUE)
output_st <- readr::read_csv(paste0("model_output/", condition, "_ST.csv")) %>%
  dplyr::mutate(`Sequence ID` = paste0(">", `Sequence ID`)) %>%
  dplyr::distinct(`Sequence ID`, .keep_all = TRUE)

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

# ---------------------------------------------------------------------------- #

# Calculate fractions of dephosphorylated sites
# Add relative fractions and total counts to the data
fraction_dephos <- bind_rows(
  output_y %>% 
    mutate(Residue = "Y"),
  output_st %>% 
    mutate(Residue = "ST")
) %>%
  mutate(
    Dephosphorylated = ifelse(`Dephosphorylation Probability` >= prob_cut_off, "Yes", "No")
  ) %>%
  group_by(Residue, Dephosphorylated) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Residue) %>%
  mutate(
    Total = sum(Count),
    Fraction = Count / Total
  ) %>%
  ungroup()

# Update plot to include relative fractions and total counts in axis labels
plot_fractions <- ggplot(fraction_dephos, aes(x = Residue, y = Fraction, fill = Dephosphorylated)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.5) +
  scale_fill_manual(values = c("Yes" = "purple", "No" = "gray"), name = "Dephosphorylated") +
  labs(
    title = paste("Relative fractions of dephosphorylated sites\nModel cut-off =", prob_cut_off),
    x = "Residue Type (Total Count)",
    y = "Fraction"
  ) +
  geom_text(
    aes(label = Count),
    position = position_fill(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  scale_x_discrete(labels = function(x) {
    totals <- fraction_dephos %>% group_by(Residue) %>% summarise(Total = unique(Total))
    sapply(x, function(residue) paste0(residue, " (", totals$Total[totals$Residue == residue], ")"))
  }) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5)
    , axis.text = element_text(color = "black")
    , legend.position = "top"
  )

plot_fractions

ggsave(plot_fractions, filename = paste0("dat_proc/plot_fractions_", condition, "_", prob_cut_off, ".png"), width = 4, height = 4)

# ---------------------------------------------------------------------------- #

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
readr::write_csv(dat_proc_combined, paste0("dat_proc/", condition, "_proc_with_predictions.csv"))

# ---------------------------------------------------------------------------- #

dat_proc_combined <- dat_proc_combined %>% 
  dplyr::select(Gene, Phosphosite, Dephosphorylation, Change) %>%
  filter(Change == "Decreased" & Dephosphorylation == "yes") %>%
  dplyr::distinct(.keep_all = TRUE)

dat_orig <- readxl::read_excel(path_to_file) %>%
  dplyr::left_join(dat_proc_combined, 
                   by = c("Gene", "Change", "Phosphosite"))

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

y_limit <- 60
x_limit <- 6

dat_plot <- dat_orig %>%
  dplyr::mutate(
    logFC = as.numeric(logFC),
    Significance_new = as.numeric(Significance),
    Significance_new = ifelse(Significance_new >= y_limit, Inf, Significance_new),
    logFC = ifelse(logFC >= x_limit, Inf, ifelse(logFC <= -x_limit, -Inf, logFC))
  )

plot_volcano <- ggplot(dat_plot, aes(x = logFC, y = Significance_new, colour = Change,
                                     text = paste("Gene:", Gene, "\n", "Description:", Description))) +
  geom_point(size = 0.7, alpha = 0.4) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(limits = c(-x_limit, x_limit)) +
  scale_y_continuous(limits = c(0, y_limit)) +
  
  geom_hline(yintercept = Sign_cutoff, colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = log2(FC_cutoff), colour = "gray", linetype = "dashed") + 
  geom_vline(xintercept = -log2(FC_cutoff), colour = "gray", linetype = "dashed") +
  
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(color = "black", size = 10),
    axis.ticks = element_blank(),
    axis.text = element_text(color = "black"),
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
    name = paste("Residue type\nModel cut-off", prob_cut_off)
  )

plot_volcano

ggsave(plot_volcano, filename = paste("dat_proc/plot_volcano_", norm_opt, "_", condition, "_", prob_cut_off, ".png", sep=""), width = 8, height = 6)

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

depod_with_genes <- readr::read_csv("depod_with_gene_names.csv")

dat_plot_decreased <- dat_plot %>%
  filter(Change == "Decreased" & Dephosphorylation == "yes") %>%
  dplyr::distinct(.keep_all = TRUE) %>%
  mutate(Gene = toupper(Gene)) %>%
  left_join(
    depod_with_genes,
    by = c("Gene" = "Substrate entry names"),
      relationship = "many-to-many"
  ) %>%
  filter(!is.na(`Phosphatase entry names`))

# ---------------------------------------------------------------------------- #

phosphatase_count <- dat_plot_decreased %>%
  mutate(`Phosphatase entry names` = as.factor(`Phosphatase entry names`)) %>%
  count(`Phosphatase entry names`) %>%
  arrange(desc(n))

# Create a plot to visualize phosphatase mapping
plot_phosphatase_mapping <- ggplot(phosphatase_count, aes(x = reorder(`Phosphatase entry names`, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, max(phosphatase_count$n) + 1, 2)) +
  labs(
    title = "Phosphatase mapping",
    x = "Phosphatase entry names",
    y = "Count"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

plot_phosphatase_mapping

ggsave(plot_phosphatase_mapping, filename = paste0("dat_proc/plot_phosphatase_mapping_", condition, "_", prob_cut_off, ".png"), width = 6, height = 4)


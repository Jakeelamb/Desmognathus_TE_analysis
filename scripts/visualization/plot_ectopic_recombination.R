#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(script_flag) > 0) sub("^--file=", "", script_flag[[1]]) else "scripts/visualization/plot_ectopic_recombination.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
helper_path <- file.path(dirname(script_dir), "R", "path_config_utils.R")
source(helper_path)

project_root <- find_project_root()
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

input_path <- file.path(project_root, "results", "data", "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv")
figure_dir <- file.path(project_root, "results", "figures", "ectopic_recombination")
data_dir <- file.path(project_root, "results", "data")
summary_path <- file.path(data_dir, "ectopic_recombination_species_summary.csv")
stats_path <- file.path(data_dir, "ectopic_recombination_species_tests.txt")
figure_path <- file.path(figure_dir, "ectopic_ratio_violin_log10.png")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Filtered ectopic recombination table not found: ", input_path, call. = FALSE)
}

ectopic <- readr::read_tsv(input_path, show_col_types = FALSE) %>%
  mutate(
    ratio_terminal_internal = suppressWarnings(as.numeric(ratio_terminal_internal)),
    species = as.character(species)
  ) %>%
  filter(!is.na(ratio_terminal_internal), ratio_terminal_internal > 0)

if (nrow(ectopic) == 0) {
  stop("No positive ectopic ratios available after filtering.", call. = FALSE)
}

species_summary <- ectopic %>%
  group_by(species) %>%
  summarise(
    n_elements = n(),
    mean_ratio = mean(ratio_terminal_internal),
    median_ratio = median(ratio_terminal_internal),
    sd_ratio = sd(ratio_terminal_internal),
    mean_log10_ratio = mean(log10(ratio_terminal_internal)),
    complete_fraction_known = mean(Complete == "yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_ratio))

readr::write_csv(species_summary, summary_path)

ordered_species <- species_summary$species
ectopic$species <- factor(ectopic$species, levels = ordered_species)

plot_data <- ectopic %>%
  mutate(log10_ratio = log10(ratio_terminal_internal))

p <- ggplot(plot_data, aes(x = species, y = ratio_terminal_internal, fill = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#b2182b", linewidth = 0.7) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.85, color = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "LTR terminal to internal depth ratio by species",
    subtitle = "Filtered to elements >= 3000 bp with five or more annotated domains",
    x = "Species",
    y = "Terminal : internal depth ratio (log10 scale)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(figure_path, p, width = 11, height = 8, dpi = 300)

kruskal_result <- kruskal.test(ratio_terminal_internal ~ species, data = ectopic)
pairwise_result <- suppressWarnings(pairwise.wilcox.test(
  x = ectopic$ratio_terminal_internal,
  g = ectopic$species,
  p.adjust.method = "BH"
))

stats_lines <- c(
  "Ectopic recombination species comparison",
  paste("Input table:", input_path),
  paste("Rows after filtering:", nrow(ectopic)),
  paste("Species represented:", length(unique(ectopic$species))),
  "",
  "Kruskal-Wallis test on raw terminal:internal ratios:",
  capture.output(print(kruskal_result)),
  "",
  "Pairwise Wilcoxon tests with BH adjustment:",
  capture.output(print(pairwise_result))
)

writeLines(stats_lines, stats_path)

cat("Ectopic recombination plot saved to:", figure_path, "\n")
cat("Species summary saved to:", summary_path, "\n")
cat("Nonparametric tests saved to:", stats_path, "\n")

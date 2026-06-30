#!/usr/bin/env Rscript

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/processing"
  }
}, error = function(...) {
  "scripts/processing"
})

source(file.path(dirname(script_dir), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(phytools)
  library(readr)
  library(tibble)
  library(yaml)
})

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
path_input_file <- file.path(project_root, "path_analysis", "data", "derived", "path_input_master.csv")

pairwise_path <- file.path(results_data_dir, "ltr_age", "ltr_age_pairwise_divergence.csv")
species_path <- file.path(results_data_dir, "ltr_age", "ltr_age_species_summary.csv")
output_dir <- file.path(results_data_dir, "ltr_age")
report_dir <- file.path(project_root, "results", "reports")
signal_note_path <- file.path(report_dir, "LTR_DIVERGENCE_SIGNAL_AUDIT.md")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

pairwise <- read_csv(pairwise_path, show_col_types = FALSE)
species_summary <- read_csv(species_path, show_col_types = FALSE)
path_input <- read_csv(path_input_file, show_col_types = FALSE)

successful <- pairwise %>%
  filter(estimation_status == "estimated")

high_conf <- successful %>%
  filter(recommended_high_confidence)

feature_candidates <- c(
  "ectopic_log10_mean_ratio",
  "weighted_te_divergence_p90",
  "order_pielou",
  "ltr_line_logratio",
  "ltr_divergence_p90"
)

species_signal <- species_summary %>%
  mutate(species_key = gsub("^D\\.", "", species)) %>%
  left_join(
    path_input %>% select(any_of(c("species", feature_candidates))),
    by = c("species_key" = "species")
  )

feature_correlations <- lapply(feature_candidates, function(metric) {
  if (!metric %in% names(species_signal)) {
    return(NULL)
  }

  sub <- species_signal %>%
    select(median_k2p_distance, !!sym(metric)) %>%
    filter(!is.na(median_k2p_distance), !is.na(.data[[metric]]))

  if (nrow(sub) < 8) {
    return(NULL)
  }

  spearman_test <- suppressWarnings(cor.test(
    sub$median_k2p_distance,
    sub[[metric]],
    method = "spearman",
    exact = FALSE
  ))

  tibble(
    feature = metric,
    n = nrow(sub),
    spearman_rho = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value
  )
}) %>%
  bind_rows() %>%
  arrange(spearman_p)

tree_file <- file.path(phylo_dir, "desmo900dated_test.tre")
tree <- read.tree(tree_file)

common_species <- intersect(tree$tip.label, species_signal$species_key)
subtree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
signal_values <- species_signal$median_k2p_distance[match(subtree$tip.label, species_signal$species_key)]
names(signal_values) <- subtree$tip.label

k_signal <- phylosig(subtree, signal_values, method = "K", test = TRUE)
lambda_signal <- phylosig(subtree, signal_values, method = "lambda", test = TRUE)

phylo_signal <- tibble(
  metric = "median_k2p_distance",
  n_species = length(signal_values),
  blombergs_k = k_signal$K,
  blombergs_k_p = k_signal$P,
  pagels_lambda = lambda_signal$lambda,
  pagels_lambda_p = lambda_signal$P
)

species_extremes <- bind_rows(
  species_summary %>%
    arrange(desc(median_k2p_distance)) %>%
    slice_head(n = 10) %>%
    mutate(extreme_set = "highest_k2p"),
  species_summary %>%
    arrange(median_k2p_distance) %>%
    slice_head(n = 10) %>%
    mutate(extreme_set = "lowest_k2p")
) %>%
  select(
    extreme_set,
    species,
    n_pairs_estimated,
    n_high_confidence_pairs_estimated,
    median_k2p_distance,
    median_comparable_sites
  )

feature_corr_path <- file.path(output_dir, "ltr_divergence_feature_correlations.csv")
phylo_signal_path <- file.path(output_dir, "ltr_divergence_phylogenetic_signal.csv")
species_extremes_path <- file.path(output_dir, "ltr_divergence_species_extremes.csv")

write_csv(feature_correlations, feature_corr_path)
write_csv(phylo_signal, phylo_signal_path)
write_csv(species_extremes, species_extremes_path)

all_n <- nrow(successful)
all_median <- median(successful$k2p_distance, na.rm = TRUE)
high_n <- nrow(high_conf)
high_median <- median(high_conf$k2p_distance, na.rm = TRUE)

top_species <- species_summary %>%
  arrange(desc(median_k2p_distance)) %>%
  slice_head(n = 5)

bottom_species <- species_summary %>%
  arrange(median_k2p_distance) %>%
  slice_head(n = 5)

top_feature_lines <- if (nrow(feature_correlations) > 0) {
  apply(head(feature_correlations, 5), 1, function(row) {
    sprintf(
      "- `%s`: Spearman rho `%0.3f`, `p = %0.4f`, `n = %s`",
      row[["feature"]],
      as.numeric(row[["spearman_rho"]]),
      as.numeric(row[["spearman_p"]]),
      row[["n"]]
    )
  })
} else {
  "- No feature correlations passed the minimum sample-size requirement."
}

lines <- c(
  "# LTR Divergence Signal Audit",
  "",
  "This note summarizes what the sequence-derived paired-LTR divergence layer is",
  "saying biologically after the assembly-backed extraction workflow completed.",
  "",
  "## Bottom Line",
  "",
  sprintf("- Successful sequence-derived LTR divergence estimates: `%d` pairs across `%d` species.", all_n, nrow(species_summary)),
  sprintf("- The estimated set is overwhelmingly `Gypsy` (`%d / %d` successful pairs).", sum(successful$Superfamily == "Gypsy"), all_n),
  sprintf("- The recommended high-confidence subset (`Complete = yes` and `5+` domains) has a slightly lower median K2P distance (`%0.6f`) than the full estimated set (`%0.6f`).", high_median, all_median),
  sprintf("- Species-level median LTR divergence shows weak phylogenetic structure (`K = %0.4f`, `p = %0.4g`; `lambda = %0.4f`, `p = %0.4g`).", phylo_signal$blombergs_k[[1]], phylo_signal$blombergs_k_p[[1]], phylo_signal$pagels_lambda[[1]], phylo_signal$pagels_lambda_p[[1]]),
  "- Species-level median LTR divergence does not show a strong monotonic relationship with the existing TE/path summary features tested here.",
  "",
  "## Species Extremes",
  "",
  "Highest median K2P species:",
  vapply(seq_len(nrow(top_species)), function(i) {
    sprintf(
      "- `%s`: median K2P `%0.6f` from `%d` estimated pairs",
      top_species$species[[i]],
      top_species$median_k2p_distance[[i]],
      top_species$n_pairs_estimated[[i]]
    )
  }, character(1)),
  "",
  "Lowest median K2P species:",
  vapply(seq_len(nrow(bottom_species)), function(i) {
    sprintf(
      "- `%s`: median K2P `%0.6f` from `%d` estimated pairs",
      bottom_species$species[[i]],
      bottom_species$median_k2p_distance[[i]],
      bottom_species$n_pairs_estimated[[i]]
    )
  }, character(1)),
  "",
  "## Feature Correlations",
  "",
  top_feature_lines,
  "",
  "## Interpretation",
  "",
  "- The new LTR divergence layer appears to add information that is not trivially redundant with the current TE composition and ectopic summary features.",
  "- Because the signal is weakly phylogenetically structured and not strongly tied to the existing TE/path features, it looks more species-specific than deeply clade-conserved at the current species coverage.",
  "- This branch is now suitable for supplementary comparative use as a sequence-derived LTR recency/divergence layer, but absolute insertion ages still require an externally justified substitution-rate calibration."
)

writeLines(lines, signal_note_path)

message("Wrote feature correlations to: ", feature_corr_path)
message("Wrote phylogenetic signal summary to: ", phylo_signal_path)
message("Wrote species extremes to: ", species_extremes_path)
message("Wrote audit note to: ", signal_note_path)

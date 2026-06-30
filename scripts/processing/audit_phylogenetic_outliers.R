#!/usr/bin/env Rscript

prefer_active_conda_r_library <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
  if (!nzchar(conda_prefix)) {
    return(invisible(NULL))
  }

  conda_lib <- normalizePath(
    file.path(conda_prefix, "lib", "R", "library"),
    mustWork = FALSE
  )
  if (dir.exists(conda_lib)) {
    .libPaths(conda_lib)
  }

  invisible(NULL)
}

prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(purrr)
  library(readr)
  library(phytools)
})

script_file_from_args <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match) == 0) {
    return(NULL)
  }
  sub("^--file=", "", match[[1]])
}

standardize_species <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^D\\.\\s*", "", x)
  x
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

script_path <- script_file_from_args()
project_root <- if (!is.null(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

tree_path <- file.path(project_root, "input_data", "phylogeny", "desmo900dated_test.tre")
master_path <- file.path(project_root, "path_analysis", "data", "derived", "path_input_master.csv")
spectrum_path <- file.path(project_root, "results", "data", "te_age_spectra", "te_age_spectrum_metrics.csv")
clade_path <- file.path(project_root, "results", "data", "permanova", "species_clade_assignments.csv")
output_dir <- file.path(project_root, "results", "data", "phylo_residuals")
summary_path <- file.path(output_dir, "phylo_residual_summary.csv")
outlier_path <- file.path(output_dir, "phylo_outlier_species.csv")
report_dir <- file.path(project_root, "results", "reports")
note_path <- file.path(report_dir, "PHYLOGENETIC_OUTLIER_AUDIT.md")

trait_spec <- tibble::tribble(
  ~trait, ~label, ~transform_label,
  "ltr_line_logratio", "LTR-to-LINE balance", "identity",
  "order_pielou", "Order-level TE evenness", "identity",
  "weighted_te_divergence_p90", "Weighted TE divergence p90", "identity",
  "ectopic_log10_mean_ratio", "Ectopic ratio log10 mean", "identity",
  "ltr_history_median_k2p_distance", "Paired-LTR median K2P", "identity",
  "ltr_history_n_pairs_estimated", "Paired-LTR recoverable count", "log1p",
  "all_recent_mass_frac_0_5", "Landscape recent-mass fraction", "identity",
  "all_old_tail_frac_20plus", "Landscape old-tail fraction", "identity"
)

transform_trait <- function(values, transform_label) {
  if (transform_label == "log1p") {
    return(log1p(values))
  }
  values
}

load_tree <- function() {
  tree <- ape::read.tree(tree_path)
  tree$tip.label <- standardize_species(tree$tip.label)
  tree
}

analyse_trait <- function(tree, data_df, trait_name, label, transform_label) {
  sub <- data_df %>%
    select(species, clade, all_of(trait_name)) %>%
    rename(value = !!trait_name) %>%
    filter(is.finite(value))

  sub$value_transformed <- transform_trait(sub$value, transform_label)
  sub <- sub %>% filter(is.finite(value_transformed))

  keep <- intersect(tree$tip.label, sub$species)
  if (length(keep) < 10) {
    return(NULL)
  }

  tree_sub <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  sub <- sub %>%
    filter(species %in% keep) %>%
    distinct(species, .keep_all = TRUE)
  trait_vec <- setNames(sub$value_transformed, sub$species)[tree_sub$tip.label]
  tip_tbl <- tibble(
    species = tree_sub$tip.label,
    value_transformed = as.numeric(trait_vec)
  ) %>%
    left_join(sub, by = c("species", "value_transformed"))

  anc <- phytools::fastAnc(tree_sub, trait_vec)
  edge_tbl <- tibble(
    parent_node = tree_sub$edge[, 1],
    child_index = tree_sub$edge[, 2],
    branch_length = tree_sub$edge.length
  ) %>%
    filter(child_index <= length(tree_sub$tip.label)) %>%
    mutate(species = tree_sub$tip.label[child_index])

  min_positive_branch <- min(edge_tbl$branch_length[edge_tbl$branch_length > 0], na.rm = TRUE)
  branch_denom <- sqrt(pmax(edge_tbl$branch_length, min_positive_branch / 2))

  outliers <- edge_tbl %>%
    mutate(
      expected_parent_state = anc[as.character(parent_node)],
      observed_value = trait_vec[species],
      raw_deviation = observed_value - expected_parent_state,
      branch_scaled_deviation = raw_deviation / branch_denom,
      terminal_shift_z = safe_scale(branch_scaled_deviation)
    ) %>%
    left_join(sub %>% select(species, clade, value), by = "species") %>%
    mutate(
      trait = trait_name,
      trait_label = label,
      transform = transform_label,
      abs_terminal_shift_z = abs(terminal_shift_z),
      direction = ifelse(terminal_shift_z >= 0, "positive", "negative"),
      outlier_flag_abs_z_ge_2 = abs_terminal_shift_z >= 2
    ) %>%
    arrange(desc(abs_terminal_shift_z))

  pos_rank <- rank(-outliers$terminal_shift_z, ties.method = "min")
  neg_rank <- rank(outliers$terminal_shift_z, ties.method = "min")
  outliers$positive_rank <- ifelse(outliers$terminal_shift_z > 0, pos_rank, NA_integer_)
  outliers$negative_rank <- ifelse(outliers$terminal_shift_z < 0, neg_rank, NA_integer_)
  outliers$abs_rank <- rank(-outliers$abs_terminal_shift_z, ties.method = "min")

  k_signal <- tryCatch(
    phytools::phylosig(tree_sub, trait_vec, method = "K", test = TRUE, nsim = 999),
    error = function(e) NULL
  )
  lambda_signal <- tryCatch(
    phytools::phylosig(tree_sub, trait_vec, method = "lambda", test = TRUE),
    error = function(e) NULL
  )

  strongest_positive <- outliers %>% filter(terminal_shift_z == max(terminal_shift_z, na.rm = TRUE)) %>% slice(1)
  strongest_negative <- outliers %>% filter(terminal_shift_z == min(terminal_shift_z, na.rm = TRUE)) %>% slice(1)

  summary_row <- tibble(
    trait = trait_name,
    trait_label = label,
    transform = transform_label,
    n_species = nrow(outliers),
    blombergs_k = if (!is.null(k_signal)) unname(k_signal$K) else NA_real_,
    blombergs_k_p = if (!is.null(k_signal)) unname(k_signal$P) else NA_real_,
    pagels_lambda = if (!is.null(lambda_signal)) unname(lambda_signal$lambda) else NA_real_,
    pagels_lambda_p = if (!is.null(lambda_signal)) unname(lambda_signal$P) else NA_real_,
    residual_sd = sd(outliers$branch_scaled_deviation, na.rm = TRUE),
    n_abs_z_ge_2 = sum(outliers$outlier_flag_abs_z_ge_2, na.rm = TRUE),
    strongest_positive_species = strongest_positive$species[[1]],
    strongest_positive_z = strongest_positive$terminal_shift_z[[1]],
    strongest_negative_species = strongest_negative$species[[1]],
    strongest_negative_z = strongest_negative$terminal_shift_z[[1]]
  )

  list(summary = summary_row, outliers = outliers)
}

write_note <- function(summary_df, outlier_df) {
  recurrent <- outlier_df %>%
    filter(outlier_flag_abs_z_ge_2) %>%
    count(species, sort = TRUE, name = "n_trait_outliers")

  recurrent_lines <- if (nrow(recurrent) == 0) {
    "- No species exceeded |z| >= 2 in more than one trait."
  } else {
    paste0("- `", recurrent$species, "`: `", recurrent$n_trait_outliers, "` traits") %>%
      paste(collapse = "\n")
  }

  strongest_lines <- summary_df %>%
    transmute(
      line = sprintf(
        "- `%s`: positive `%s` (z = %.3f), negative `%s` (z = %.3f), lambda = %.3f",
        trait,
        strongest_positive_species,
        strongest_positive_z,
        strongest_negative_species,
        strongest_negative_z,
        pagels_lambda
      )
    ) %>%
    pull(line) %>%
    paste(collapse = "\n")

  note <- paste(
    "# Phylogenetic Outlier Audit",
    "",
    "## Purpose",
    "",
    "This note identifies species that sit unusually high or low for key TE and",
    "LTR-history variables relative to nearby phylogenetic expectation.",
    "",
    "## Upstream inputs",
    "",
    "- `input_data/phylogeny/desmo900dated_test.tre`",
    "- `path_analysis/data/derived/path_input_master.csv`",
    "- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`",
    "- `results/data/permanova/species_clade_assignments.csv`",
    "",
    "## Outputs",
    "",
    "- `results/data/phylo_residuals/phylo_residual_summary.csv`",
    "- `results/data/phylo_residuals/phylo_outlier_species.csv`",
    "",
    "## Method",
    "",
    "For each trait, ancestral states were reconstructed on the time-calibrated",
    "tree with `phytools::fastAnc()`. Species-level outlier scores are terminal",
    "branch deviations from the reconstructed parent state, scaled by branch",
    "length and then z-scored across species within trait.",
    "",
    "This is a branch-local outlier screen, not a claim about full causal",
    "direction or adaptive explanation.",
    "",
    "## Strongest branch-level deviations",
    "",
    strongest_lines,
    "",
    "## Recurrent outlier species",
    "",
    recurrent_lines,
    "",
    "## Bottom line",
    "",
    "The TE landscape and LTR-history layers are not only clade-level patterns.",
    "They also contain species-level terminal deviations that can anchor concrete",
    "biological follow-up and interpretation.",
    sep = "\n"
  )

  writeLines(note, note_path)
}

main <- function() {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  tree <- load_tree()
  master <- readr::read_csv(master_path, show_col_types = FALSE) %>%
    mutate(species = standardize_species(species))
  spectrum <- readr::read_csv(spectrum_path, show_col_types = FALSE) %>%
    mutate(species = standardize_species(species))
  clades <- readr::read_csv(clade_path, show_col_types = FALSE) %>%
    transmute(species = standardize_species(species), clade = as.character(clade))

  data_df <- master %>%
    left_join(
      spectrum %>% select(species, all_recent_mass_frac_0_5, all_old_tail_frac_20plus),
      by = "species"
    ) %>%
    left_join(clades, by = "species")

  results <- purrr::pmap(
    trait_spec,
    function(trait, label, transform_label) {
      analyse_trait(tree, data_df, trait, label, transform_label)
    }
  )
  results <- results[!vapply(results, is.null, logical(1))]

  summary_df <- bind_rows(purrr::map(results, "summary")) %>%
    arrange(desc(n_abs_z_ge_2), desc(abs(strongest_positive_z)), desc(abs(strongest_negative_z)))
  outlier_df <- bind_rows(purrr::map(results, "outliers")) %>%
    arrange(trait, abs_rank, desc(abs_terminal_shift_z))

  readr::write_csv(summary_df, summary_path)
  readr::write_csv(outlier_df, outlier_path)
  write_note(summary_df, outlier_df)

  message("Wrote ", summary_path)
  message("Wrote ", outlier_path)
  message("Wrote ", note_path)
}

main()

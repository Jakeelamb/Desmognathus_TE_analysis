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
  library(ggplot2)
  library(phytools)
  library(readr)
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
  x <- sub("^Desmognathus\\s+", "", x)
  x <- sub("^D\\.\\s*", "", x)
  x
}

display_species <- function(x) {
  x <- standardize_species(x)
  ifelse(nzchar(x), paste0("D. ", x), "")
}

safe_z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

safe_pct <- function(numer, denom) {
  out <- (numer / denom - 1) * 100
  out[!is.finite(out)] <- NA_real_
  out
}

format_num <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")
}

format_pct <- function(x, digits = 1) {
  ifelse(is.finite(x), paste0(formatC(x, format = "f", digits = digits), "%"), "NA")
}

script_path <- script_file_from_args()
project_root <- if (!is.null(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

tree_path <- file.path(project_root, "input_data", "phylogeny", "desmo900dated_test.tre")
balanced_dir <- file.path(
  project_root,
  "path_analysis",
  "data",
  "external",
  "derived",
  "balanced_genome_iod_sensitivity"
)
estimate_path <- file.path(balanced_dir, "balanced_genome_iod_species_estimates.csv")
comparison_path <- file.path(balanced_dir, "balanced_genome_iod_comparison.csv")
output_dir <- file.path(balanced_dir, "phylogenetic_genome_audit")
figure_dir <- file.path(output_dir, "figures")

species_audit_path <- file.path(output_dir, "genome_phylo_species_audit.csv")
signal_summary_path <- file.path(output_dir, "genome_phylo_signal_summary.csv")
branch_jump_path <- file.path(output_dir, "genome_phylo_branch_jumps.csv")
report_path <- file.path(output_dir, "GENOME_PHYLOGENETIC_AUDIT.md")
html_path <- file.path(output_dir, "index.html")

required_panels <- c("current_top50_primary", "balanced_qc_only")
primary_panel <- "balanced_qc_only"
reference_species_key <- "fuscus"
reference_species_display <- "D. fuscus"

viewer_rel_path_for_panel <- function(panel_name) {
  if (panel_name == "balanced_qc_curated") {
    return("../viewer_balanced_qc_curated/index.html")
  }
  "../viewer_balanced_qc_only/index.html"
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path, call. = FALSE)
  }
}

load_tree <- function(path) {
  tree <- ape::read.tree(path)
  tree$tip.label <- standardize_species(tree$tip.label)
  tree
}

nearest_relative_table <- function(tree_sub, species_keys, genome_pg) {
  dist_mat <- ape::cophenetic.phylo(tree_sub)
  rows <- lapply(species_keys, function(species_key) {
    distances <- dist_mat[species_key, , drop = TRUE]
    distances <- distances[names(distances) != species_key]
    distances <- distances[is.finite(distances)]
    nearest_key <- names(which.min(distances))
    nearest_distance <- unname(distances[[nearest_key]])
    current_pg <- unname(genome_pg[[species_key]])
    nearest_pg <- unname(genome_pg[[nearest_key]])
    tibble(
      species_key = species_key,
      nearest_species_key = nearest_key,
      nearest_species = display_species(nearest_key),
      nearest_patristic_distance_mya = nearest_distance,
      nearest_genome_pg = nearest_pg,
      nearest_relative_pct_diff = safe_pct(current_pg, nearest_pg),
      abs_nearest_relative_pct_diff = abs(safe_pct(current_pg, nearest_pg))
    )
  })
  bind_rows(rows)
}

branch_descendant_label <- function(tree_sub, node_id) {
  if (node_id <= length(tree_sub$tip.label)) {
    return(tree_sub$tip.label[[node_id]])
  }
  tips <- ape::extract.clade(tree_sub, node = node_id)$tip.label
  if (length(tips) <= 4) {
    paste(tips, collapse = ";")
  } else {
    paste0(length(tips), "_tip_clade:", paste(head(tips, 4), collapse = ";"), ";...")
  }
}

node_state <- function(tree_sub, trait_vec, anc, node_id) {
  if (node_id <= length(tree_sub$tip.label)) {
    return(unname(trait_vec[[tree_sub$tip.label[[node_id]]]]))
  }
  unname(anc[[as.character(node_id)]])
}

classify_support_penalty <- function(support_label) {
  dplyr::case_when(
    support_label == "low" ~ 1,
    support_label == "limited" ~ 0.5,
    TRUE ~ 0
  )
}

add_audit_triggers <- function(df) {
  df %>%
    mutate(
      trigger_fuscus_delta_ge_50pct = abs_pct_vs_fuscus >= 50,
      trigger_terminal_phylo_z_ge_2 = abs_terminal_shift_z >= 2,
      trigger_nearest_relative_jump_ge_25pct = abs_nearest_relative_pct_diff >= 25,
      trigger_method_shift_ge_40pct = abs_balanced_vs_current_pct >= 40,
      trigger_low_independent_support = support_label == "low",
      audit_trigger_count =
        as.integer(trigger_fuscus_delta_ge_50pct) +
        as.integer(trigger_terminal_phylo_z_ge_2) +
        as.integer(trigger_nearest_relative_jump_ge_25pct) +
        as.integer(trigger_method_shift_ge_40pct) +
        as.integer(trigger_low_independent_support),
      audit_trigger_summary = apply(
        cbind(
          ifelse(trigger_fuscus_delta_ge_50pct, "fuscus_delta_ge_50pct", ""),
          ifelse(trigger_terminal_phylo_z_ge_2, "terminal_phylo_z_ge_2", ""),
          ifelse(trigger_nearest_relative_jump_ge_25pct, "nearest_relative_jump_ge_25pct", ""),
          ifelse(trigger_method_shift_ge_40pct, "method_shift_ge_40pct", ""),
          ifelse(trigger_low_independent_support, "low_independent_support", "")
        ),
        1,
        function(x) paste(x[nzchar(x)], collapse = "; ")
      ),
      audit_priority_score =
        pmin(abs_pct_vs_fuscus / 50, 2) +
        pmin(abs_terminal_shift_z, 3) +
        pmin(abs_nearest_relative_pct_diff / 25, 3) +
        pmin(abs_balanced_vs_current_pct / 40, 2) +
        classify_support_penalty(support_label),
      audit_priority_tier = case_when(
        audit_trigger_count >= 3 ~ "high",
        trigger_terminal_phylo_z_ge_2 | trigger_nearest_relative_jump_ge_25pct ~ "high",
        audit_trigger_count >= 2 ~ "medium",
        audit_trigger_count == 1 ~ "watch",
        TRUE ~ "baseline"
      )
    )
}

analyse_panel <- function(tree, estimates, comparison, panel_name) {
  panel_df <- estimates %>%
    filter(panel == panel_name, is.finite(estimated_genome_pg), estimated_genome_pg > 0) %>%
    mutate(
      species_key = standardize_species(species),
      species = display_species(species_key)
    ) %>%
    distinct(species_key, .keep_all = TRUE)

  keep <- intersect(tree$tip.label, panel_df$species_key)
  if (length(keep) < 6) {
    stop("Panel ", panel_name, " has fewer than 6 species matching the tree.", call. = FALSE)
  }

  tree_sub <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  panel_df <- panel_df %>% filter(species_key %in% tree_sub$tip.label)
  genome_pg <- setNames(panel_df$estimated_genome_pg, panel_df$species_key)[tree_sub$tip.label]
  trait_vec <- log(genome_pg)

  anc <- phytools::fastAnc(tree_sub, trait_vec)
  edge_tbl <- tibble(
    panel = panel_name,
    parent_node = tree_sub$edge[, 1],
    child_node = tree_sub$edge[, 2],
    branch_length_mya = tree_sub$edge.length
  )
  min_positive_branch <- min(edge_tbl$branch_length_mya[edge_tbl$branch_length_mya > 0], na.rm = TRUE)
  if (!is.finite(min_positive_branch)) {
    min_positive_branch <- 1
  }

  branch_tbl <- edge_tbl %>%
    rowwise() %>%
    mutate(
      child_is_tip = child_node <= length(tree_sub$tip.label),
      child_species_key = ifelse(child_is_tip, tree_sub$tip.label[[child_node]], NA_character_),
      child_label = branch_descendant_label(tree_sub, child_node),
      parent_log_genome = node_state(tree_sub, trait_vec, anc, parent_node),
      child_log_genome = node_state(tree_sub, trait_vec, anc, child_node),
      parent_genome_pg = exp(parent_log_genome),
      child_genome_pg = exp(child_log_genome),
      branch_log_shift = child_log_genome - parent_log_genome,
      branch_pct_shift = (exp(branch_log_shift) - 1) * 100,
      branch_scaled_shift = branch_log_shift / sqrt(max(branch_length_mya, min_positive_branch / 2))
    ) %>%
    ungroup() %>%
    mutate(
      branch_shift_z = safe_z(branch_scaled_shift),
      abs_branch_shift_z = abs(branch_shift_z),
      abs_branch_pct_shift = abs(branch_pct_shift),
      abs_branch_shift_rank = rank(-abs_branch_shift_z, ties.method = "min")
    ) %>%
    arrange(abs_branch_shift_rank, desc(abs_branch_pct_shift))

  terminal_tbl <- branch_tbl %>%
    filter(child_is_tip) %>%
    transmute(
      panel,
      species_key = child_species_key,
      expected_parent_genome_pg = parent_genome_pg,
      terminal_log_shift = branch_log_shift,
      terminal_pct_shift = branch_pct_shift,
      terminal_branch_length_mya = branch_length_mya,
      terminal_branch_scaled_shift = branch_scaled_shift,
      terminal_shift_z = branch_shift_z,
      abs_terminal_shift_z = abs_branch_shift_z
    )

  dist_to_fuscus <- rep(NA_real_, length(tree_sub$tip.label))
  names(dist_to_fuscus) <- tree_sub$tip.label
  if (reference_species_key %in% tree_sub$tip.label) {
    dist_mat <- ape::cophenetic.phylo(tree_sub)
    dist_to_fuscus <- dist_mat[, reference_species_key]
  }
  fuscus_pg <- if (reference_species_key %in% names(genome_pg)) {
    unname(genome_pg[[reference_species_key]])
  } else {
    NA_real_
  }
  nearest_tbl <- nearest_relative_table(tree_sub, tree_sub$tip.label, genome_pg)

  species_tbl <- panel_df %>%
    select(
      panel,
      species,
      species_key,
      estimated_genome_pg,
      estimated_genome_pg_q1,
      estimated_genome_pg_q3,
      iod_ratio_to_reference,
      n_selected_pairs,
      n_selected_images,
      n_selected_specimens,
      effective_n,
      support_label,
      support_warnings
    ) %>%
    left_join(terminal_tbl, by = c("panel", "species_key")) %>%
    left_join(nearest_tbl, by = "species_key") %>%
    mutate(
      fuscus_genome_pg = fuscus_pg,
      delta_vs_fuscus_pg = estimated_genome_pg - fuscus_pg,
      pct_vs_fuscus = safe_pct(estimated_genome_pg, fuscus_pg),
      abs_pct_vs_fuscus = abs(pct_vs_fuscus),
      patristic_distance_to_fuscus_mya = dist_to_fuscus[species_key],
      viewer_species_url = paste0(viewer_rel_path_for_panel(panel_name), "?species=", URLencode(species, reserved = TRUE))
    )

  current_keyed <- estimates %>%
    filter(panel == "current_top50_primary") %>%
    mutate(species_key = standardize_species(species)) %>%
    select(
      species_key,
      support_label_current_top50_primary = support_label,
      support_warnings_current_top50_primary = support_warnings,
      estimated_genome_pg_current_top50_primary = estimated_genome_pg
    )
  original_keyed <- estimates %>%
    filter(panel == "balanced_qc_only") %>%
    mutate(species_key = standardize_species(species)) %>%
    select(
      species_key,
      estimated_genome_pg_balanced_qc_only = estimated_genome_pg,
      support_label_balanced_qc_only = support_label,
      support_warnings_balanced_qc_only = support_warnings
    )

  species_tbl <- species_tbl %>%
    left_join(current_keyed, by = "species_key") %>%
    left_join(original_keyed, by = "species_key") %>%
    mutate(
      panel_minus_current_pg = estimated_genome_pg - estimated_genome_pg_current_top50_primary,
      panel_vs_current_pct = safe_pct(estimated_genome_pg, estimated_genome_pg_current_top50_primary),
      panel_minus_original_balanced_pg = estimated_genome_pg - estimated_genome_pg_balanced_qc_only,
      panel_vs_original_balanced_pct = safe_pct(estimated_genome_pg, estimated_genome_pg_balanced_qc_only),
      balanced_minus_current_pg = panel_minus_current_pg,
      balanced_vs_current_pct = panel_vs_current_pct,
      abs_balanced_vs_current_pct = abs(panel_vs_current_pct)
    ) %>%
    add_audit_triggers() %>%
    mutate(
      fuscus_deviation_rank = rank(-abs_pct_vs_fuscus, ties.method = "min"),
      terminal_phylo_rank = rank(-abs_terminal_shift_z, ties.method = "min"),
      nearest_relative_jump_rank = rank(-abs_nearest_relative_pct_diff, ties.method = "min"),
      audit_priority_rank = rank(-audit_priority_score, ties.method = "min")
    ) %>%
    arrange(audit_priority_rank, fuscus_deviation_rank, terminal_phylo_rank)

  set.seed(20260708 + sum(utf8ToInt(panel_name)))
  k_signal <- tryCatch(
    phytools::phylosig(tree_sub, trait_vec, method = "K", test = TRUE, nsim = 999),
    error = function(e) NULL
  )
  lambda_signal <- tryCatch(
    phytools::phylosig(tree_sub, trait_vec, method = "lambda", test = TRUE),
    error = function(e) NULL
  )

  signal_tbl <- tibble(
    panel = panel_name,
    n_species = nrow(species_tbl),
    reference_species = reference_species_display,
    reference_genome_pg = fuscus_pg,
    blombergs_k = if (!is.null(k_signal)) unname(k_signal$K) else NA_real_,
    blombergs_k_p = if (!is.null(k_signal)) unname(k_signal$P) else NA_real_,
    pagels_lambda = if (!is.null(lambda_signal)) unname(lambda_signal$lambda) else NA_real_,
    pagels_lambda_p = if (!is.null(lambda_signal)) unname(lambda_signal$P) else NA_real_,
    spearman_genome_vs_fuscus_distance = suppressWarnings(
      cor(
        species_tbl$estimated_genome_pg,
        species_tbl$patristic_distance_to_fuscus_mya,
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    ),
    n_terminal_abs_z_ge_2 = sum(species_tbl$abs_terminal_shift_z >= 2, na.rm = TRUE),
    n_fuscus_delta_ge_50pct = sum(species_tbl$abs_pct_vs_fuscus >= 50, na.rm = TRUE),
    n_nearest_relative_jump_ge_25pct = sum(species_tbl$abs_nearest_relative_pct_diff >= 25, na.rm = TRUE),
    strongest_fuscus_deviation_species = species_tbl$species[which.max(species_tbl$abs_pct_vs_fuscus)],
    strongest_fuscus_deviation_pct = max(species_tbl$abs_pct_vs_fuscus, na.rm = TRUE),
    strongest_terminal_outlier_species = species_tbl$species[which.max(species_tbl$abs_terminal_shift_z)],
    strongest_terminal_outlier_z = max(species_tbl$abs_terminal_shift_z, na.rm = TRUE),
    strongest_nearest_jump_species = species_tbl$species[which.max(species_tbl$abs_nearest_relative_pct_diff)],
    strongest_nearest_jump_pct = max(species_tbl$abs_nearest_relative_pct_diff, na.rm = TRUE)
  )

  list(
    species = species_tbl,
    signal = signal_tbl,
    branches = branch_tbl,
    tree = tree_sub
  )
}

write_ranked_plot <- function(df, y_col, output_path, title, y_label, top_n = 21) {
  plot_df <- df %>%
    arrange(desc(abs(.data[[y_col]]))) %>%
    slice_head(n = top_n) %>%
    mutate(species = reorder(species, .data[[y_col]]))

  p <- ggplot(plot_df, aes(x = species, y = .data[[y_col]], fill = audit_priority_tier)) +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.15) +
    coord_flip() +
    labs(x = NULL, y = y_label, title = title, fill = "Audit tier") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )
  ggsave(output_path, p, width = 8.5, height = 6.5, dpi = 180)
}

write_fuscus_distance_plot <- function(df, output_path) {
  label_df <- df %>%
    arrange(audit_priority_rank) %>%
    slice_head(n = 8)

  p <- ggplot(
    df,
    aes(
      x = patristic_distance_to_fuscus_mya,
      y = estimated_genome_pg,
      color = terminal_shift_z,
      shape = support_label
    )
  ) +
    geom_hline(aes(yintercept = fuscus_genome_pg), color = "grey45", linetype = "dashed", linewidth = 0.35) +
    geom_point(size = 3.4, alpha = 0.92) +
    geom_text(
      data = label_df,
      aes(label = species),
      check_overlap = TRUE,
      nudge_y = 0.45,
      size = 3,
      show.legend = FALSE
    ) +
    scale_color_gradient2(low = "#3666a6", mid = "grey75", high = "#a8382d", midpoint = 0) +
    labs(
      x = "Patristic distance to D. fuscus on dated tree",
      y = "Estimated genome size (pg)",
      color = "Terminal\nphylo z",
      shape = "Support",
      title = "Genome size relative to D. fuscus and phylogenetic distance"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
  ggsave(output_path, p, width = 9, height = 6.3, dpi = 180)
}

write_tree_plot <- function(tree_sub, df, output_path) {
  ordered <- df %>% distinct(species_key, .keep_all = TRUE)
  genome_pg <- setNames(ordered$estimated_genome_pg, ordered$species_key)
  z_score <- setNames(ordered$terminal_shift_z, ordered$species_key)
  tiers <- setNames(ordered$audit_priority_tier, ordered$species_key)
  tip_pg <- genome_pg[tree_sub$tip.label]
  rng <- range(tip_pg, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    scaled <- rep(0.5, length(tip_pg))
  } else {
    scaled <- (tip_pg - rng[[1]]) / diff(rng)
  }
  palette <- colorRampPalette(c("#2e5e9e", "#f5f1df", "#b23a2f"))(101)
  tip_cols <- palette[pmax(1, pmin(101, round(scaled * 100) + 1))]
  top_tips <- which(tree_sub$tip.label %in% ordered$species_key[ordered$audit_priority_rank <= 8])

  png(output_path, width = 1500, height = 1050, res = 160)
  par(mar = c(3, 1, 3, 9), xpd = NA)
  plot.phylo(
    tree_sub,
    cex = 0.72,
    label.offset = 0.15,
    main = "Balanced genome-IOD estimates on dated phylogeny"
  )
  tiplabels(pch = 21, bg = tip_cols, col = "grey20", cex = 1.15)
  if (length(top_tips) > 0) {
    tiplabels(tip = top_tips, pch = 8, col = "black", cex = 0.85)
  }
  legend(
    "right",
    inset = c(-0.22, 0),
    legend = c(
      paste0("Genome pg: ", format_num(rng[[1]], 1)),
      paste0("Genome pg: ", format_num(mean(rng), 1)),
      paste0("Genome pg: ", format_num(rng[[2]], 1)),
      "Star: top audit priority"
    ),
    pt.bg = c(palette[[1]], palette[[51]], palette[[101]], NA),
    pch = c(21, 21, 21, 8),
    col = c("grey20", "grey20", "grey20", "black"),
    bty = "n",
    cex = 0.8
  )
  invisible(z_score)
  invisible(tiers)
  dev.off()
}

markdown_table <- function(df, columns, n = 10) {
  show_df <- df %>%
    select(all_of(columns)) %>%
    slice_head(n = n)
  header <- paste(columns, collapse = " | ")
  sep <- paste(rep("---", length(columns)), collapse = " | ")
  rows <- apply(show_df, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, sep, rows), collapse = "\n")
}

write_report <- function(species_df, signal_df) {
  primary <- species_df %>%
    filter(panel == primary_panel) %>%
    arrange(audit_priority_rank)
  top_fuscus <- primary %>%
    arrange(fuscus_deviation_rank) %>%
    transmute(
      species,
      genome_pg = format_num(estimated_genome_pg, 2),
      vs_fuscus = format_pct(pct_vs_fuscus, 1),
      terminal_z = format_num(terminal_shift_z, 2),
      nearest = nearest_species,
      nearest_jump = format_pct(nearest_relative_pct_diff, 1),
      triggers = ifelse(nzchar(audit_trigger_summary), audit_trigger_summary, "none")
    )
  top_terminal <- primary %>%
    arrange(terminal_phylo_rank) %>%
    transmute(
      species,
      genome_pg = format_num(estimated_genome_pg, 2),
      terminal_z = format_num(terminal_shift_z, 2),
      parent_expected_pg = format_num(expected_parent_genome_pg, 2),
      terminal_shift = format_pct(terminal_pct_shift, 1),
      triggers = ifelse(nzchar(audit_trigger_summary), audit_trigger_summary, "none")
    )
  top_nearest <- primary %>%
    arrange(nearest_relative_jump_rank) %>%
    transmute(
      species,
      genome_pg = format_num(estimated_genome_pg, 2),
      nearest = nearest_species,
      nearest_genome_pg = format_num(nearest_genome_pg, 2),
      nearest_jump = format_pct(nearest_relative_pct_diff, 1),
      distance_mya = format_num(nearest_patristic_distance_mya, 2)
    )
  top_method <- primary %>%
    arrange(desc(abs_balanced_vs_current_pct)) %>%
    transmute(
      species,
      current_pg = format_num(estimated_genome_pg_current_top50_primary, 2),
      primary_panel_pg = format_num(estimated_genome_pg, 2),
      shift_vs_current = format_pct(panel_vs_current_pct, 1),
      shift_vs_original_balanced = format_pct(panel_vs_original_balanced_pct, 1),
      support = support_label,
      triggers = ifelse(nzchar(audit_trigger_summary), audit_trigger_summary, "none")
    )

  signal_primary <- signal_df %>% filter(panel == primary_panel) %>% slice(1)
  note <- paste(
    "# Genome Phylogenetic Deviation Audit",
    "",
    "## Purpose",
    "",
    "This audit ranks species that most deserve image-level review because their",
    "estimated genome sizes are far from `D. fuscus`, far from local",
    "phylogenetic expectation, or sharply different from their nearest measured",
    "relative on the time-calibrated tree.",
    "",
    "## Inputs",
    "",
    paste0("- Tree: `", tree_path, "`"),
    paste0("- Genome estimates: `", estimate_path, "`"),
    paste0("- Current-vs-balanced comparison: `", comparison_path, "`"),
    "",
    "## Outputs",
    "",
    paste0("- Species audit table: `", species_audit_path, "`"),
    paste0("- Phylogenetic signal summary: `", signal_summary_path, "`"),
    paste0("- Branch jump table: `", branch_jump_path, "`"),
    paste0("- Browser dashboard: `", html_path, "`"),
    "",
    "## Method",
    "",
    "Genome size was analyzed on the log scale. For each estimate panel,",
    "`phytools::fastAnc()` reconstructed internal states on the dated tree.",
    "Terminal outlier scores are observed tip values minus reconstructed parent",
    "states, scaled by terminal branch length and then z-scored within panel.",
    "Nearest-relative jumps compare each species to the closest measured species",
    "by patristic distance. This is an audit-prioritization screen, not a final",
    "claim that any species is biologically wrong.",
    "",
    paste0("## ", primary_panel, " signal snapshot"),
    "",
    paste0("- Species: `", signal_primary$n_species, "`"),
    paste0("- Blomberg K: `", format_num(signal_primary$blombergs_k, 3), "`; p = `", format_num(signal_primary$blombergs_k_p, 3), "`"),
    paste0("- Pagel lambda: `", format_num(signal_primary$pagels_lambda, 3), "`; p = `", format_num(signal_primary$pagels_lambda_p, 3), "`"),
    paste0("- Spearman genome vs distance from D. fuscus: `", format_num(signal_primary$spearman_genome_vs_fuscus_distance, 3), "`"),
    paste0("- Terminal |z| >= 2 species: `", signal_primary$n_terminal_abs_z_ge_2, "`"),
    paste0("- Species >= 50% from D. fuscus: `", signal_primary$n_fuscus_delta_ge_50pct, "`"),
    paste0("- Species >= 25% from nearest measured relative: `", signal_primary$n_nearest_relative_jump_ge_25pct, "`"),
    "",
    "## Largest deviations from D. fuscus",
    "",
    markdown_table(top_fuscus, names(top_fuscus), n = 10),
    "",
    "## Largest terminal phylogenetic residuals",
    "",
    markdown_table(top_terminal, names(top_terminal), n = 10),
    "",
    "## Largest nearest-relative jumps",
    "",
    markdown_table(top_nearest, names(top_nearest), n = 10),
    "",
    "## Largest current-to-primary-panel method shifts",
    "",
    markdown_table(top_method, names(top_method), n = 10),
    "",
    "## Figures",
    "",
    "- `figures/genome_vs_fuscus_distance.png`",
    "- `figures/terminal_phylo_outlier_rank.png`",
    "- `figures/nearest_relative_jump_rank.png`",
    "- `figures/tree_genome_audit.png`",
    "",
    "## Audit rule",
    "",
    "Prioritize species with multiple triggers, especially if the same species is",
    "far from `D. fuscus`, has a high terminal phylogenetic residual, and also",
    "shows a large current-to-balanced method shift. Those species should be",
    "checked in the mask viewer before their genome-size values are treated as",
    "biological signal.",
    sep = "\n"
  )
  writeLines(note, report_path)
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

write_html_dashboard <- function(species_df, signal_df) {
  primary <- species_df %>%
    filter(panel == primary_panel) %>%
    arrange(audit_priority_rank)
  top_rows <- primary %>% slice_head(n = 21)

  row_html <- apply(top_rows, 1, function(row) {
    triggers <- row[["audit_trigger_summary"]]
    if (!nzchar(triggers)) {
      triggers <- "none"
    }
    paste0(
      "<tr>",
      "<td>", html_escape(row[["audit_priority_rank"]]), "</td>",
      "<td><a href=\"", html_escape(row[["viewer_species_url"]]), "\">", html_escape(row[["species"]]), "</a></td>",
      "<td>", html_escape(format_num(as.numeric(row[["estimated_genome_pg"]]), 2)), "</td>",
      "<td>", html_escape(format_pct(as.numeric(row[["pct_vs_fuscus"]]), 1)), "</td>",
      "<td>", html_escape(format_num(as.numeric(row[["terminal_shift_z"]]), 2)), "</td>",
      "<td>", html_escape(row[["nearest_species"]]), "</td>",
      "<td>", html_escape(format_pct(as.numeric(row[["nearest_relative_pct_diff"]]), 1)), "</td>",
      "<td>", html_escape(format_pct(as.numeric(row[["panel_vs_current_pct"]]), 1)), "</td>",
      "<td>", html_escape(row[["support_label"]]), "</td>",
      "<td>", html_escape(triggers), "</td>",
      "</tr>"
    )
  })

  signal_primary <- signal_df %>% filter(panel == primary_panel) %>% slice(1)
  html <- paste0(
    "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>Genome Phylogenetic Audit</title>\n",
    "<style>",
    "body{font-family:Inter,system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;background:#f6f6f2;color:#20282b;}",
    "header{position:sticky;top:0;background:#f6f6f2;border-bottom:1px solid #d8ddd9;padding:14px 18px;z-index:5;}",
    "h1{font-size:20px;margin:0 0 4px;}p{margin:4px 0;color:#59656a;}main{padding:18px;max-width:1240px;margin:0 auto;}",
    ".stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin-bottom:18px;}",
    ".stat{background:white;border:1px solid #d8ddd9;border-radius:8px;padding:10px;}.stat b{display:block;font-size:20px;}",
    ".figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:14px;margin:16px 0;}",
    ".fig{background:white;border:1px solid #d8ddd9;border-radius:8px;padding:10px;}.fig img{width:100%;height:auto;display:block;}",
    "table{width:100%;border-collapse:collapse;background:white;border:1px solid #d8ddd9;font-size:13px;}th,td{border-bottom:1px solid #e2e6e3;padding:8px;text-align:left;vertical-align:top;}th{background:#eef1ee;position:sticky;top:61px;}",
    "a{color:#1d5f8f;text-decoration:none;}a:hover{text-decoration:underline;}",
    "</style>\n</head>\n<body>\n",
    "<header><h1>Genome Phylogenetic Audit</h1>",
    "<p>", html_escape(primary_panel), " genome estimates ranked by fuscus-relative deviation, terminal phylogenetic residuals, nearest-relative jumps, and method sensitivity. Species links open the mask viewer filtered to that species.</p></header>\n",
    "<main>\n<section class=\"stats\">",
    "<div class=\"stat\"><span>Species</span><b>", html_escape(signal_primary$n_species), "</b></div>",
    "<div class=\"stat\"><span>Blomberg K</span><b>", html_escape(format_num(signal_primary$blombergs_k, 3)), "</b></div>",
    "<div class=\"stat\"><span>Pagel lambda</span><b>", html_escape(format_num(signal_primary$pagels_lambda, 3)), "</b></div>",
    "<div class=\"stat\"><span>|terminal z| >= 2</span><b>", html_escape(signal_primary$n_terminal_abs_z_ge_2), "</b></div>",
    "</section>\n",
    "<section class=\"figs\">",
    "<div class=\"fig\"><img src=\"figures/genome_vs_fuscus_distance.png\" alt=\"Genome size versus distance to fuscus\"></div>",
    "<div class=\"fig\"><img src=\"figures/terminal_phylo_outlier_rank.png\" alt=\"Terminal phylogenetic residual ranks\"></div>",
    "<div class=\"fig\"><img src=\"figures/nearest_relative_jump_rank.png\" alt=\"Nearest relative jump ranks\"></div>",
    "<div class=\"fig\"><img src=\"figures/tree_genome_audit.png\" alt=\"Genome estimates on phylogeny\"></div>",
    "</section>\n",
    "<table><thead><tr>",
    "<th>Rank</th><th>Species</th><th>Genome pg</th><th>vs fuscus</th><th>Terminal z</th><th>Nearest</th><th>Nearest jump</th><th>Panel vs current</th><th>Support</th><th>Triggers</th>",
    "</tr></thead><tbody>",
    paste(row_html, collapse = "\n"),
    "</tbody></table>\n",
    "<p>Tables: <a href=\"genome_phylo_species_audit.csv\">species audit CSV</a> | <a href=\"genome_phylo_signal_summary.csv\">signal summary CSV</a> | <a href=\"genome_phylo_branch_jumps.csv\">branch jumps CSV</a> | <a href=\"GENOME_PHYLOGENETIC_AUDIT.md\">markdown report</a></p>",
    "</main>\n</body>\n</html>\n"
  )
  writeLines(html, html_path)
}

main <- function() {
  require_file(tree_path)
  require_file(estimate_path)
  require_file(comparison_path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  tree <- load_tree(tree_path)
  estimates <- readr::read_csv(estimate_path, show_col_types = FALSE)
  comparison <- readr::read_csv(comparison_path, show_col_types = FALSE)

  missing_panels <- setdiff(required_panels, unique(estimates$panel))
  if (length(missing_panels) > 0) {
    stop("Missing required estimate panels: ", paste(missing_panels, collapse = ", "), call. = FALSE)
  }
  primary_panel <<- if ("balanced_qc_curated" %in% unique(estimates$panel)) {
    "balanced_qc_curated"
  } else {
    "balanced_qc_only"
  }

  panels <- unique(estimates$panel)
  analyses <- lapply(panels, function(panel_name) analyse_panel(tree, estimates, comparison, panel_name))
  names(analyses) <- panels

  species_df <- bind_rows(lapply(analyses, `[[`, "species")) %>%
    arrange(panel, audit_priority_rank, species)
  signal_df <- bind_rows(lapply(analyses, `[[`, "signal")) %>%
    arrange(panel)
  branch_df <- bind_rows(lapply(analyses, `[[`, "branches")) %>%
    arrange(panel, abs_branch_shift_rank)

  readr::write_csv(species_df, species_audit_path)
  readr::write_csv(signal_df, signal_summary_path)
  readr::write_csv(branch_df, branch_jump_path)

  primary_df <- species_df %>% filter(panel == primary_panel)
  write_fuscus_distance_plot(primary_df, file.path(figure_dir, "genome_vs_fuscus_distance.png"))
  write_ranked_plot(
    primary_df,
    "terminal_shift_z",
    file.path(figure_dir, "terminal_phylo_outlier_rank.png"),
    paste(primary_panel, "terminal phylogenetic residuals"),
    "Terminal branch residual z"
  )
  write_ranked_plot(
    primary_df,
    "nearest_relative_pct_diff",
    file.path(figure_dir, "nearest_relative_jump_rank.png"),
    paste(primary_panel, "nearest measured relative jumps"),
    "Percent difference from nearest measured relative"
  )
  write_tree_plot(
    analyses[[primary_panel]]$tree,
    primary_df,
    file.path(figure_dir, "tree_genome_audit.png")
  )

  write_report(species_df, signal_df)
  write_html_dashboard(species_df, signal_df)

  message("Wrote ", species_audit_path)
  message("Wrote ", signal_summary_path)
  message("Wrote ", branch_jump_path)
  message("Wrote ", report_path)
  message("Wrote ", html_path)
}

main()

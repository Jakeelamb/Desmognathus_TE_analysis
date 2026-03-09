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
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

if (!requireNamespace("phylopath", quietly = TRUE)) {
  stop("Package 'phylopath' is required for audit_model_stability.R")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

script_file_from_args <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match) == 0) {
    return(NULL)
  }
  sub("^--file=", "", match[[1]])
}

script_path <- script_file_from_args()
project_root <- if (!is.null(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
workspace_root <- file.path(project_root, "path_analysis")

source(file.path(workspace_root, "scripts", "path_model_scaffold.R"))

derived_dir <- file.path(workspace_root, "data", "derived")
panel_dir <- file.path(derived_dir, "panels")
results_dir <- file.path(workspace_root, "results")
output_summary <- file.path(derived_dir, "model_stability_summary.csv")
output_detail <- file.path(derived_dir, "model_stability_detail.csv")
note_path <- file.path(workspace_root, "MODEL_STABILITY_AUDIT.md")
clade_path <- file.path(project_root, "results", "data", "permanova", "species_clade_assignments.csv")

targets <- tibble::tribble(
  ~family, ~panel,
  "te_genome", "te_genome_primary_mediumplus",
  "te_genome_organismal", "te_genome_organismal_primary_mediumplus",
  "te_genome_ectopic_organismal", "te_genome_ectopic_organismal_primary_mediumplus",
  "te_genome_ltr_history", "te_genome_ltr_history_primary_mediumplus"
)

tracked_edges <- tibble::tribble(
  ~edge_id, ~parent, ~child,
  "ltr_balance_to_te_evenness", "ltr_balance", "te_evenness",
  "te_evenness_to_gs", "te_evenness", "gs",
  "body_size_to_gs", "body_size", "gs",
  "ltr_history_to_gs", "ltr_history", "gs",
  "ectopic_index_to_gs", "ectopic_index", "gs"
)

sign_code <- function(x) {
  ifelse(
    is.na(x), "absent",
    ifelse(x > 0, "positive", ifelse(x < 0, "negative", "zero"))
  )
}

extract_edge_value <- function(edge_tbl, parent, child) {
  value <- edge_tbl %>%
    filter(.data$parent == !!parent, .data$child == !!child) %>%
    pull(coefficient)
  if (length(value) == 0) {
    return(NA_real_)
  }
  as.numeric(value[[1]])
}

load_baseline <- function(panel) {
  ranking <- readr::read_csv(
    file.path(results_dir, paste0(panel, "_model_ranking.csv")),
    show_col_types = FALSE
  )
  edges <- readr::read_csv(
    file.path(results_dir, paste0(panel, "_best_model_edges.csv")),
    show_col_types = FALSE
  )

  baseline <- list(
    winner_model = ranking$model[[1]],
    winner_margin_cicc = if (nrow(ranking) >= 2) as.numeric(ranking$delta_CICc[[2]]) else NA_real_,
    top_weight = as.numeric(ranking$w[[1]]),
    edge_tbl = edges
  )

  for (i in seq_len(nrow(tracked_edges))) {
    edge <- tracked_edges[i, ]
    value <- extract_edge_value(edges, edge$parent[[1]], edge$child[[1]])
    baseline[[paste0(edge$edge_id[[1]], "_coef")]] <- value
    baseline[[paste0(edge$edge_id[[1]], "_sign")]] <- sign_code(value)
  }

  baseline
}

build_run_table <- function(panel_df, clade_df) {
  species_runs <- tibble(
    run_type = "leave_one_species_out",
    run_label = paste0("drop_", panel_df$species),
    removed_species = panel_df$species,
    removed_clade = NA_character_
  )

  named_clades <- clade_df %>%
    inner_join(panel_df %>% select(species), by = "species") %>%
    filter(!is.na(clade), clade != "other") %>%
    count(clade, name = "panel_clade_size") %>%
    filter(panel_clade_size >= 3)

  clade_runs <- clade_df %>%
    inner_join(named_clades, by = "clade") %>%
    inner_join(panel_df %>% select(species), by = "species") %>%
    group_by(clade) %>%
    summarise(
      removed_species = paste(sort(species), collapse = ";"),
      .groups = "drop"
    ) %>%
    transmute(
      run_type = "clade_jackknife",
      run_label = paste0("drop_clade_", clade),
      removed_species,
      removed_clade = clade
    )

  bind_rows(species_runs, clade_runs)
}

run_family_case <- function(family, panel, panel_df, run_row, baseline) {
  removed_species <- unlist(strsplit(run_row$removed_species[[1]], ";", fixed = TRUE))
  test_df <- panel_df %>%
    filter(!species %in% removed_species)

  if (nrow(test_df) < 10) {
    return(tibble(
      family = family,
      panel = panel,
      run_type = run_row$run_type[[1]],
      run_label = run_row$run_label[[1]],
      removed_species = run_row$removed_species[[1]],
      removed_clade = run_row$removed_clade[[1]],
      n_species = nrow(test_df),
      status = "skipped_too_small"
    ))
  }

  warning_messages <- character()
  fit_result <- tryCatch(
    withCallingHandlers(
      {
        analysis_df <- prepare_analysis_input(test_df, family, panel)
        tree <- load_tree(project_root, analysis_df$species)
        model_set <- model_set_for_family(family, asNamespace("phylopath"))
        analysis_mat <- tibble::column_to_rownames(analysis_df, var = "species")
        tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(analysis_mat)))
        analysis_mat <- analysis_mat[tree$tip.label, , drop = FALSE]
        fit <- phylopath::phylo_path(model_set, data = analysis_mat, tree = tree)
        ranking <- as.data.frame(summary(fit))
        best_fit <- phylopath::best(fit)
        best_edges <- matrix_to_edges(best_fit$coef, best_fit$se)
        list(
          analysis_df = analysis_df,
          ranking = ranking,
          best_edges = best_edges
        )
      },
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit_result, "error")) {
    return(tibble(
      family = family,
      panel = panel,
      run_type = run_row$run_type[[1]],
      run_label = run_row$run_label[[1]],
      removed_species = run_row$removed_species[[1]],
      removed_clade = run_row$removed_clade[[1]],
      n_species = nrow(test_df),
      status = "failed",
      error_message = conditionMessage(fit_result)
    ))
  }

  ranking <- fit_result$ranking
  best_edges <- fit_result$best_edges

  row <- tibble(
    family = family,
    panel = panel,
    run_type = run_row$run_type[[1]],
    run_label = run_row$run_label[[1]],
    removed_species = run_row$removed_species[[1]],
    removed_clade = run_row$removed_clade[[1]],
    n_species = nrow(fit_result$analysis_df),
    status = "ok",
    warning_count = length(unique(warning_messages)),
    winner_model = ranking$model[[1]],
    winner_margin_cicc = if (nrow(ranking) >= 2) as.numeric(ranking$delta_CICc[[2]]) else NA_real_,
    winner_weight = as.numeric(ranking$w[[1]]),
    baseline_winner_model = baseline$winner_model,
    baseline_winner_margin_cicc = baseline$winner_margin_cicc,
    same_winner = ranking$model[[1]] == baseline$winner_model
  )

  for (i in seq_len(nrow(tracked_edges))) {
    edge <- tracked_edges[i, ]
    edge_id <- edge$edge_id[[1]]
    value <- extract_edge_value(best_edges, edge$parent[[1]], edge$child[[1]])
    current_sign <- sign_code(value)
    baseline_sign <- baseline[[paste0(edge_id, "_sign")]]

    row[[paste0(edge_id, "_coef")]] <- value
    row[[paste0(edge_id, "_sign")]] <- current_sign
    row[[paste0(edge_id, "_same_sign_as_baseline")]] <- identical(current_sign, baseline_sign)
    row[[paste0(edge_id, "_present")]] <- !is.na(value)
  }

  row
}

summarise_family_runs <- function(detail_df, baseline) {
  ok_runs <- detail_df %>% filter(status == "ok")

  summary_row <- ok_runs %>%
    summarise(
      family = first(family),
      panel = first(panel),
      run_type = "all_runs",
      n_runs = n(),
      same_winner_n = sum(same_winner, na.rm = TRUE),
      same_winner_frac = mean(same_winner, na.rm = TRUE),
      unique_winner_models = n_distinct(winner_model),
      min_winner_margin_cicc = min(winner_margin_cicc, na.rm = TRUE),
      median_winner_margin_cicc = median(winner_margin_cicc, na.rm = TRUE),
      failed_or_skipped_runs = sum(detail_df$status != "ok"),
      baseline_winner_model = first(baseline$winner_model)
    )

  split_rows <- ok_runs %>%
    group_by(family, panel, run_type) %>%
    summarise(
      n_runs = n(),
      same_winner_n = sum(same_winner, na.rm = TRUE),
      same_winner_frac = mean(same_winner, na.rm = TRUE),
      unique_winner_models = n_distinct(winner_model),
      min_winner_margin_cicc = min(winner_margin_cicc, na.rm = TRUE),
      median_winner_margin_cicc = median(winner_margin_cicc, na.rm = TRUE),
      failed_or_skipped_runs = sum(detail_df$status != "ok" & detail_df$run_type == first(run_type)),
      baseline_winner_model = baseline$winner_model,
      .groups = "drop"
    )

  for (i in seq_len(nrow(tracked_edges))) {
    edge_id <- tracked_edges$edge_id[[i]]
    presence_col <- paste0(edge_id, "_present")
    sign_col <- paste0(edge_id, "_same_sign_as_baseline")
    baseline_sign <- baseline[[paste0(edge_id, "_sign")]]
    baseline_present <- baseline_sign != "absent"

    summary_row[[paste0(edge_id, "_baseline_sign")]] <- baseline_sign
    summary_row[[paste0(edge_id, "_presence_frac")]] <- if (baseline_present) {
      mean(ok_runs[[presence_col]], na.rm = TRUE)
    } else {
      NA_real_
    }
    summary_row[[paste0(edge_id, "_same_sign_frac")]] <- if (baseline_present) {
      mean(ok_runs[[sign_col]], na.rm = TRUE)
    } else {
      NA_real_
    }

    split_rows[[paste0(edge_id, "_baseline_sign")]] <- baseline_sign
    split_rows[[paste0(edge_id, "_presence_frac")]] <- if (baseline_present) {
      ok_runs %>%
        group_by(run_type) %>%
        summarise(value = mean(.data[[presence_col]], na.rm = TRUE), .groups = "drop") %>%
        arrange(match(run_type, split_rows$run_type)) %>%
        pull(value)
    } else {
      NA_real_
    }
    split_rows[[paste0(edge_id, "_same_sign_frac")]] <- if (baseline_present) {
      ok_runs %>%
        group_by(run_type) %>%
        summarise(value = mean(.data[[sign_col]], na.rm = TRUE), .groups = "drop") %>%
        arrange(match(run_type, split_rows$run_type)) %>%
        pull(value)
    } else {
      NA_real_
    }
  }

  bind_rows(summary_row, split_rows)
}

write_note <- function(summary_df, detail_df) {
  family_sections <- list()
  for (i in seq_len(nrow(targets))) {
    family <- targets$family[[i]]
    panel <- targets$panel[[i]]
    baseline <- load_baseline(panel)
    all_runs <- summary_df %>%
      filter(family == !!family, panel == !!panel, run_type == "all_runs")
    loo <- summary_df %>%
      filter(family == !!family, panel == !!panel, run_type == "leave_one_species_out")
    clade <- summary_df %>%
      filter(family == !!family, panel == !!panel, run_type == "clade_jackknife")
    if (nrow(all_runs) == 0) {
      next
    }

    baseline_edges <- tracked_edges %>%
      mutate(
        baseline_sign = map_chr(edge_id, ~ baseline[[paste0(.x, "_sign")]])
      ) %>%
      filter(baseline_sign != "absent")

    edge_lines <- vapply(
      seq_len(nrow(baseline_edges)),
      function(idx) {
        edge <- baseline_edges[idx, ]
        same_frac <- all_runs[[paste0(edge$edge_id[[1]], "_same_sign_frac")]]
        present_frac <- all_runs[[paste0(edge$edge_id[[1]], "_presence_frac")]]
        sprintf(
          "- `%s -> %s`: baseline `%s`, present in %.3f of runs, same sign in %.3f of runs",
          edge$parent[[1]],
          edge$child[[1]],
          edge$baseline_sign[[1]],
          present_frac,
          same_frac
        )
      },
      character(1)
    )

    family_sections[[length(family_sections) + 1]] <- paste(
      sprintf("## %s", panel),
      "",
      sprintf("- Baseline winner: `%s`", baseline$winner_model),
      sprintf("- Winner retained across all runs: `%d / %d` (%.3f)",
              all_runs$same_winner_n[[1]], all_runs$n_runs[[1]], all_runs$same_winner_frac[[1]]),
      sprintf("- Leave-one-species-out retention: `%d / %d` (%.3f)",
              loo$same_winner_n[[1]], loo$n_runs[[1]], loo$same_winner_frac[[1]]),
      sprintf("- Clade-jackknife retention: `%d / %d` (%.3f)",
              clade$same_winner_n[[1]], clade$n_runs[[1]], clade$same_winner_frac[[1]]),
      sprintf("- Alternative winners observed: `%d`", all_runs$unique_winner_models[[1]]),
      sprintf("- Median winner margin (delta CICc to runner-up): `%.3f`", all_runs$median_winner_margin_cicc[[1]]),
      "",
      "Tracked edge stability:",
      paste(edge_lines, collapse = "\n"),
      sep = "\n"
    )
  }

  note <- paste(
    "# Model Stability Audit",
    "",
    "## Purpose",
    "",
    "This note audits whether the current primary medium-plus path-model winners",
    "are robust to single-species removal and named-clade jackknife exclusion.",
    "",
    "## Upstream inputs",
    "",
    "- `path_analysis/data/derived/panels/*.csv`",
    "- `path_analysis/results/*primary_mediumplus_model_ranking.csv`",
    "- `path_analysis/results/*primary_mediumplus_best_model_edges.csv`",
    "- `results/data/permanova/species_clade_assignments.csv`",
    "- `input_data/phylogeny/desmo900dated_test.tre`",
    "",
    "## Outputs",
    "",
    "- `path_analysis/data/derived/model_stability_detail.csv`",
    "- `path_analysis/data/derived/model_stability_summary.csv`",
    "",
    paste(family_sections, collapse = "\n\n"),
    "",
    "## Bottom line",
    "",
    "Treat families as primary only if winner identity is broadly retained and",
    "their baseline edge signs remain stable under both leave-one-out and named",
    "clade jackknife perturbations.",
    sep = "\n"
  )

  writeLines(note, note_path)
}

main <- function() {
  clade_df <- readr::read_csv(clade_path, show_col_types = FALSE) %>%
    transmute(
      species = standardize_species(species),
      clade = as.character(clade)
    )

  detail_rows <- list()
  summary_rows <- list()

  for (i in seq_len(nrow(targets))) {
    family <- targets$family[[i]]
    panel <- targets$panel[[i]]
    message("Auditing stability for ", panel)
    panel_df <- readr::read_csv(file.path(panel_dir, paste0(panel, ".csv")), show_col_types = FALSE) %>%
      mutate(species = standardize_species(species))
    baseline <- load_baseline(panel)
    runs <- build_run_table(panel_df, clade_df)

    family_detail <- purrr::map_dfr(
      seq_len(nrow(runs)),
      ~ run_family_case(family, panel, panel_df, runs[.x, , drop = FALSE], baseline)
    )
    detail_rows[[length(detail_rows) + 1]] <- family_detail
    summary_rows[[length(summary_rows) + 1]] <- summarise_family_runs(family_detail, baseline)
  }

  detail_df <- bind_rows(detail_rows) %>%
    arrange(family, run_type, run_label)
  summary_df <- bind_rows(summary_rows) %>%
    arrange(family, factor(run_type, levels = c("all_runs", "leave_one_species_out", "clade_jackknife")))

  readr::write_csv(detail_df, output_detail)
  readr::write_csv(summary_df, output_summary)
  write_note(summary_df, detail_df)

  message("Wrote ", output_detail)
  message("Wrote ", output_summary)
  message("Wrote ", note_path)
}

main()

#!/usr/bin/env Rscript

# Path-analysis scaffold for the Desmognathus path-analysis workspace.
# This script is intentionally conservative: it prepares transformed inputs,
# defines a small candidate model set, and only runs phylopath if the package
# is available locally.

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
  library(readr)
  library(tibble)
})

print_help <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript path_analysis/scripts/path_model_scaffold.R [--family FAMILY] [--panel PANEL_NAME] [--derived-dir DIR] [--summary-only]",
      "",
      "Families:",
      "  te_genome",
      "  te_genome_organismal",
      "  te_genome_ltr_history",
      "  te_genome_ectopic",
      "  te_genome_ectopic_organismal",
      "  genome_morphology",
      "  te_genome_morphology",
      "",
      "Panels:",
      "  te_genome_all",
      "  te_genome_primary",
      "  te_genome_primary_mediumplus",
      "  te_genome_primary_phylofill",
      "  te_genome_primary_strict_body",
      "  te_genome_organismal_primary_mediumplus",
      "  te_genome_organismal_primary_phylofill",
      "  te_genome_organismal_primary_strict_body",
      "  te_genome_ltr_history_primary_mediumplus",
      "  te_genome_ltr_history_primary_strict_body",
      "  te_genome_ectopic_all",
      "  te_genome_ectopic_primary_mediumplus",
      "  te_genome_ectopic_organismal_primary_mediumplus",
      "  te_genome_ectopic_organismal_primary_phylofill",
      "  te_genome_ectopic_organismal_primary_strict_body",
      "  te_genome_ectopic_primary_phylofill",
      "  te_genome_ectopic_primary_strict_body",
      "  te_genome_morphology_all",
      "  te_genome_morphology_primary_mediumplus",
      "  te_genome_morphology_primary_phylofill",
      "  te_genome_morphology_primary_strict_body",
      "",
      "Examples:",
      "  Rscript path_analysis/scripts/path_model_scaffold.R --summary-only",
      "  Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --summary-only",
      "  Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_strict_body",
      sep = "\n"
    )
  )
}

parse_args <- function(args) {
  opts <- list(
    family = "te_genome",
    panel = NULL,
    derived_dir = NULL,
    summary_only = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      print_help()
      quit(save = "no", status = 0)
    } else if (arg == "--family") {
      i <- i + 1
      opts$family <- args[[i]]
    } else if (arg == "--panel") {
      i <- i + 1
      opts$panel <- args[[i]]
    } else if (arg == "--derived-dir") {
      i <- i + 1
      opts$derived_dir <- args[[i]]
    } else if (arg == "--summary-only") {
      opts$summary_only <- TRUE
    } else {
      stop("Unknown argument: ", arg)
    }
    i <- i + 1
  }

  opts
}

find_project_root <- function() {
  current <- normalizePath(getwd(), mustWork = TRUE)
  for (i in 1:10) {
    if (file.exists(file.path(current, "paths.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }

  script_dir <- normalizePath(dirname(sys.frame(1)$ofile %||% "path_analysis/scripts"), mustWork = FALSE)
  current <- script_dir
  for (i in 1:10) {
    if (file.exists(file.path(current, "paths.yaml"))) {
      return(normalizePath(current))
    }
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }

  stop("Could not locate project root containing paths.yaml")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

standardize_species <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^D\\.\\s*", "", x)
  x
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

build_transform_audit <- function(analysis_df) {
  rows <- list()
  predictor_cols <- setdiff(colnames(analysis_df), "species")

  for (col in predictor_cols) {
    values <- as.numeric(analysis_df[[col]])
    bad_idx <- which(!is.finite(values))
    if (length(bad_idx)) {
      for (idx in bad_idx) {
        rows[[length(rows) + 1]] <- tibble::tibble(
          audit_type = "non_finite_transform",
          species = analysis_df$species[[idx]],
          column = col,
          details = paste0("Value after transformation: ", as.character(analysis_df[[col]][[idx]]))
        )
      }
    }

    finite_values <- values[is.finite(values)]
    if (length(finite_values) > 1 && stats::sd(finite_values) == 0) {
      rows[[length(rows) + 1]] <- tibble::tibble(
        audit_type = "constant_predictor",
        species = NA_character_,
        column = col,
        details = "Predictor has zero variance in the prepared analysis input"
      )
    }
  }

  if (!length(rows)) {
    return(tibble::tibble(
      audit_type = character(),
      species = character(),
      column = character(),
      details = character()
    ))
  }

  dplyr::bind_rows(rows)
}

clr_transform <- function(mat) {
  mat <- as.matrix(mat)
  positive <- mat[mat > 0]
  pseudocount <- if (length(positive)) min(positive) / 2 else 1e-6
  mat <- mat + pseudocount
  gm <- exp(rowMeans(log(mat)))
  sweep(log(mat), 1, log(gm), "-")
}

dataset_file_for_family <- function(family) {
  switch(
    family,
    te_genome = "dataset_te_genome.csv",
    te_genome_organismal = file.path("panels", "te_genome_organismal_primary_mediumplus.csv"),
    te_genome_ltr_history = file.path("panels", "te_genome_ltr_history_primary_mediumplus.csv"),
    te_genome_ectopic = "dataset_te_genome_ectopic.csv",
    te_genome_ectopic_organismal = file.path("panels", "te_genome_ectopic_organismal_primary_mediumplus.csv"),
    genome_morphology = "dataset_genome_morphology.csv",
    te_genome_morphology = "dataset_te_genome_morphology.csv",
    stop("Unknown family: ", family)
  )
}

input_spec_for_run <- function(family, panel, derived_dir) {
  if (is.null(panel)) {
    return(
      list(
        input_file = file.path(derived_dir, dataset_file_for_family(family)),
        output_tag = family
      )
    )
  }

  compatible_prefixes <- c(family)
  if (family == "te_genome_organismal") {
    compatible_prefixes <- c(compatible_prefixes, "te_genome")
  }
  if (family == "te_genome_ltr_history") {
    compatible_prefixes <- c(compatible_prefixes, "te_genome")
  }
  if (family == "te_genome_ectopic_organismal") {
    compatible_prefixes <- c(compatible_prefixes, "te_genome_ectopic")
  }

  if (!any(startsWith(panel, paste0(compatible_prefixes, "_")))) {
    stop("Panel '", panel, "' does not match family '", family, "'.")
  }

  list(
    input_file = file.path(derived_dir, "panels", paste0(panel, ".csv")),
    output_tag = panel
  )
}

load_tree <- function(project_root, species) {
  tree_path <- file.path(project_root, "input_data", "phylogeny", "desmo900dated_test.tre")
  tree <- ape::read.tree(tree_path)
  tree$tip.label <- standardize_species(tree$tip.label)
  keep <- intersect(tree$tip.label, species)
  if (length(keep) < 4) {
    stop("Fewer than 4 matching species between tree and dataset")
  }
  ape::drop.tip(tree, setdiff(tree$tip.label, keep))
}

prepare_analysis_input <- function(df, family, panel = NULL) {
  df$species <- standardize_species(df$species)
  use_phylofill <- !is.null(panel) && grepl("phylofill", panel)

  needs_te_features <- family %in% c(
    "te_genome",
    "te_genome_organismal",
    "te_genome_ltr_history",
    "te_genome_ectopic",
    "te_genome_ectopic_organismal",
    "te_genome_morphology"
  )
  if (needs_te_features) {
    order_cols <- c(
      "order_dirs", "order_helitron", "order_line", "order_ltr",
      "order_maverick", "order_sine", "order_tir", "order_yr"
    )
    available_order_cols <- intersect(order_cols, colnames(df))
    complete_te_rows <- stats::complete.cases(df[, available_order_cols, drop = FALSE])

    if (length(available_order_cols) >= 3 && any(complete_te_rows)) {
      clr <- clr_transform(df[complete_te_rows, available_order_cols, drop = FALSE])
      pca <- stats::prcomp(clr, center = FALSE, scale. = FALSE)
      df$te_pc1 <- NA_real_
      df$te_pc2 <- NA_real_
      df$te_pc1[complete_te_rows] <- as.numeric(pca$x[, 1])
      df$te_pc2[complete_te_rows] <- as.numeric(pca$x[, 2])
    } else {
      df$te_pc1 <- NA_real_
      df$te_pc2 <- NA_real_
    }

    if ("ltr_line_logratio" %in% colnames(df)) {
      df$ltr_balance <- as.numeric(df$ltr_line_logratio)
    } else {
      ltr <- df$order_ltr
      line <- df$order_line
      ltr_line_positive <- c(ltr[ltr > 0], line[line > 0])
      pseudocount <- if (length(ltr_line_positive)) min(ltr_line_positive) / 2 else 1e-6
      df$ltr_balance <- log((ltr + pseudocount) / (line + pseudocount))
    }
  }

  df$gs <- zscore(log10(df$genome_size_pg))
  if ("order_pielou" %in% colnames(df)) {
    df$te_evenness <- zscore(df$order_pielou)
  }
  if ("ectopic_mean_ratio" %in% colnames(df)) {
    df$ectopic_index <- zscore(log10(df$ectopic_mean_ratio))
  } else if ("ectopic_log10_mean_ratio" %in% colnames(df)) {
    df$ectopic_index <- zscore(df$ectopic_log10_mean_ratio)
  }
  if ("morph_nucleus_area_um2" %in% colnames(df)) {
    df$ns <- zscore(log10(df$morph_nucleus_area_um2))
  }
  if ("morph_cell_area_um2" %in% colnames(df)) {
    df$cs <- zscore(log10(df$morph_cell_area_um2))
  }
  if ("body_size_proxy_mm" %in% colnames(df)) {
    body_size_col <- if (use_phylofill && "body_size_proxy_mm_phylofill" %in% colnames(df)) {
      "body_size_proxy_mm_phylofill"
    } else {
      "body_size_proxy_mm"
    }
    df$body_size <- zscore(log10(df[[body_size_col]]))
  }
  if ("aquaticity_index" %in% colnames(df)) {
    aquaticity_col <- if (use_phylofill && "aquaticity_index_phylofill" %in% colnames(df)) {
      "aquaticity_index_phylofill"
    } else {
      "aquaticity_index"
    }
    df$aquaticity <- zscore(as.numeric(df[[aquaticity_col]]))
  }
  if ("ltr_history_median_k2p_distance" %in% colnames(df)) {
    df$ltr_history <- zscore(as.numeric(df$ltr_history_median_k2p_distance))
  }

  analysis_df <- switch(
    family,
    te_genome = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness),
    te_genome_organismal = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness, body_size, aquaticity),
    te_genome_ltr_history = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness, ltr_history),
    te_genome_ectopic = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness, ectopic_index),
    te_genome_ectopic_organismal = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness, ectopic_index, body_size),
    genome_morphology = df %>%
      transmute(species, gs, ns, cs),
    te_genome_morphology = df %>%
      transmute(species, gs, ltr_balance = zscore(ltr_balance), te_evenness, ns, cs)
  )

  transform_audit <- build_transform_audit(analysis_df)

  analysis_df <- analysis_df %>%
    filter(if_all(-species, ~ is.finite(.x))) %>%
    distinct(species, .keep_all = TRUE)

  analysis_df <- as.data.frame(analysis_df)
  attr(analysis_df, "transform_audit") <- transform_audit
  analysis_df
}

model_set_for_family <- function(family, phylopath_ns) {
  switch(
    family,
    te_genome = phylopath_ns$define_model_set(
      genome_null = c(),
      ltr_balance_only = c(gs ~ ltr_balance),
      evenness_only = c(gs ~ te_evenness),
      additive_load_evenness = c(gs ~ ltr_balance + te_evenness),
      mediated_evenness = c(te_evenness ~ ltr_balance, gs ~ te_evenness)
    ),
    te_genome_organismal = phylopath_ns$define_model_set(
      te_baseline = c(te_evenness ~ ltr_balance, gs ~ te_evenness),
      body_size_additive = c(te_evenness ~ ltr_balance, gs ~ te_evenness + body_size),
      aquaticity_additive = c(te_evenness ~ ltr_balance, gs ~ te_evenness + aquaticity),
      organismal_additive = c(te_evenness ~ ltr_balance, gs ~ te_evenness + body_size + aquaticity),
      aquaticity_confounds_body_and_te = c(
        body_size ~ aquaticity,
        te_evenness ~ ltr_balance + aquaticity,
        gs ~ te_evenness + body_size
      )
    ),
    te_genome_ltr_history = phylopath_ns$define_model_set(
      te_baseline = c(te_evenness ~ ltr_balance, gs ~ te_evenness),
      history_direct = c(gs ~ ltr_history),
      history_additive = c(te_evenness ~ ltr_balance, gs ~ te_evenness + ltr_history),
      history_to_balance = c(ltr_balance ~ ltr_history, te_evenness ~ ltr_balance, gs ~ te_evenness),
      history_to_evenness = c(te_evenness ~ ltr_balance + ltr_history, gs ~ te_evenness)
    ),
    te_genome_ectopic = phylopath_ns$define_model_set(
      ectopic_only = c(gs ~ ectopic_index),
      ltr_to_ectopic = c(ectopic_index ~ ltr_balance, gs ~ ectopic_index),
      evenness_and_ectopic = c(ectopic_index ~ ltr_balance, gs ~ te_evenness + ectopic_index),
      full_mechanism = c(ectopic_index ~ ltr_balance + te_evenness, gs ~ ltr_balance + te_evenness + ectopic_index)
    ),
    te_genome_ectopic_organismal = phylopath_ns$define_model_set(
      ectopic_baseline = c(gs ~ ectopic_index),
      ectopic_body_size_additive = c(gs ~ ectopic_index + body_size),
      te_body_size_baseline = c(te_evenness ~ ltr_balance, gs ~ te_evenness + body_size),
      te_ectopic_body_size = c(te_evenness ~ ltr_balance, gs ~ te_evenness + ectopic_index + body_size),
      ltr_to_ectopic_body_size = c(ectopic_index ~ ltr_balance, gs ~ ectopic_index + body_size)
    ),
    genome_morphology = phylopath_ns$define_model_set(
      morphology_null = c(),
      genome_to_nucleus = c(ns ~ gs),
      genome_to_cell_direct = c(cs ~ gs),
      mediated_cell_size = c(ns ~ gs, cs ~ ns)
    ),
    te_genome_morphology = phylopath_ns$define_model_set(
      te_to_genome_to_nucleus_to_cell = c(gs ~ ltr_balance + te_evenness, ns ~ gs, cs ~ ns),
      te_to_genome_partial_cell = c(gs ~ ltr_balance + te_evenness, ns ~ gs, cs ~ gs + ns),
      te_direct_to_genome_nucleus = c(gs ~ ltr_balance, ns ~ gs, cs ~ ns),
      te_evenness_path = c(te_evenness ~ ltr_balance, gs ~ te_evenness, ns ~ gs, cs ~ ns)
    ),
    stop("Unknown family: ", family)
  )
}

matrix_to_edges <- function(coef_mat, se_mat = NULL) {
  coef_df <- as.data.frame(as.table(coef_mat), stringsAsFactors = FALSE)
  colnames(coef_df) <- c("parent", "child", "coefficient")
  coef_df <- coef_df %>% filter(coefficient != 0)

  if (!is.null(se_mat)) {
    se_df <- as.data.frame(as.table(se_mat), stringsAsFactors = FALSE)
    colnames(se_df) <- c("parent", "child", "std_error")
    coef_df <- coef_df %>% left_join(se_df, by = c("parent", "child"))
  }

  coef_df
}

run_summary_only <- function(family, analysis_df) {
  cat("Family:", family, "\n")
  cat("Species in prepared input:", nrow(analysis_df), "\n")
  cat("Columns:\n")
  cat(paste0("  - ", colnames(analysis_df)), sep = "\n")
  transform_audit <- attr(analysis_df, "transform_audit")
  if (!is.null(transform_audit) && nrow(transform_audit) > 0) {
    cat("\nTransform audit warnings:\n")
    print(transform_audit)
  }
  cat("\nphylopath is not required for --summary-only mode.\n")
}

write_model_warnings <- function(fit, output_path) {
  if (length(fit$warnings) == 0) {
    if (file.exists(output_path)) {
      unlink(output_path)
    }
    return(invisible(NULL))
  }

  warning_text <- paste(unlist(fit$warnings), collapse = "\n")
  writeLines(warning_text, output_path)
  invisible(NULL)
}

run_phylopath_family <- function(family, analysis_df, tree, results_dir, output_tag = family) {
  if (!requireNamespace("phylopath", quietly = TRUE)) {
    stop(
      paste(
        "Package 'phylopath' is not installed.",
        "Install it with install.packages('phylopath') and rerun this script."
      )
    )
  }

  phylopath_ns <- asNamespace("phylopath")
  model_set <- model_set_for_family(family, phylopath_ns)
  analysis_df <- tibble::column_to_rownames(analysis_df, var = "species")
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(analysis_df)))
  analysis_df <- analysis_df[tree$tip.label, , drop = FALSE]

  fit <- phylopath_ns$phylo_path(model_set, data = analysis_df, tree = tree)
  summary_tbl <- as.data.frame(summary(fit))
  readr::write_csv(summary_tbl, file.path(results_dir, paste0(output_tag, "_model_ranking.csv")))
  write_model_warnings(fit, file.path(results_dir, paste0(output_tag, "_model_warnings.txt")))

  best_fit <- phylopath_ns$best(fit)
  best_edges <- matrix_to_edges(best_fit$coef, best_fit$se)
  readr::write_csv(best_edges, file.path(results_dir, paste0(output_tag, "_best_model_edges.csv")))

  avg_fit <- try(phylopath_ns$average(fit, cut_off = 2), silent = TRUE)
  if (!inherits(avg_fit, "try-error")) {
    avg_edges <- matrix_to_edges(avg_fit$coef, avg_fit$se)
    readr::write_csv(avg_edges, file.path(results_dir, paste0(output_tag, "_average_model_edges.csv")))
  }

  pdf(file.path(results_dir, paste0(output_tag, "_model_set.pdf")), width = 10, height = 7)
  phylopath_ns$plot_model_set(model_set)
  dev.off()

  pdf(file.path(results_dir, paste0(output_tag, "_model_summary.pdf")), width = 10, height = 7)
  plot(summary(fit))
  dev.off()

  pdf(file.path(results_dir, paste0(output_tag, "_best_model.pdf")), width = 10, height = 7)
  plot(best_fit)
  dev.off()

  if (!inherits(avg_fit, "try-error")) {
    pdf(file.path(results_dir, paste0(output_tag, "_average_model.pdf")), width = 10, height = 7)
    plot(avg_fit)
    dev.off()
  }

  invisible(fit)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  project_root <- find_project_root()
  workspace_root <- file.path(project_root, "path_analysis")
  derived_dir <- opts$derived_dir %||% file.path(workspace_root, "data", "derived")
  results_dir <- file.path(workspace_root, "results")

  input_spec <- input_spec_for_run(opts$family, opts$panel, derived_dir)
  if (!file.exists(input_spec$input_file)) {
    stop("Derived dataset not found: ", input_spec$input_file, "\nRun the dataset/panel builders first.")
  }

  df <- readr::read_csv(input_spec$input_file, show_col_types = FALSE)
  analysis_df <- prepare_analysis_input(df, opts$family, opts$panel)

  if (opts$family %in% c("genome_morphology", "te_genome_morphology")) {
    warning(
      "Morphology families are scaffolded for planning. ",
      "Current genome_size_pg is still provisional for formal genome->nucleus interpretation."
    )
  }

  if (opts$summary_only) {
    run_summary_only(opts$family, analysis_df)
    return(invisible(NULL))
  }

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(analysis_df, file.path(results_dir, paste0(input_spec$output_tag, "_analysis_input.csv")))
  transform_audit <- attr(analysis_df, "transform_audit")
  if (!is.null(transform_audit) && nrow(transform_audit) > 0) {
    readr::write_csv(transform_audit, file.path(results_dir, paste0(input_spec$output_tag, "_transform_audit.csv")))
  }

  tree <- load_tree(project_root, analysis_df$species)
  run_phylopath_family(opts$family, analysis_df, tree, results_dir, input_spec$output_tag)
}

if (sys.nframe() == 0) {
  main()
}

#!/usr/bin/env Rscript

# Canonical TE compositional PCA workflow.
#
# This script is the single supported entrypoint for non-phylogenetic TE PCA.
# It intentionally replaces the older raw-proportion and clustering-heavy
# variants with a compositional workflow that:
# 1. uses frozen repo-local TE breakdown tables
# 2. filters features explicitly
# 3. applies deterministic zero replacement
# 4. runs CLR PCA with no post-CLR variance scaling
# 5. writes the exact transformed matrices used downstream

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  } else {
    "scripts/processing"
  }
}, error = function(...) {
  "scripts/processing"
})

source(file.path(script_dir, "pca_utils.R"))


parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list(
    analysis = NULL,
    list_analyses = FALSE,
    no_plots = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--analysis") {
      i <- i + 1L
      params$analysis <- args[[i]]
    } else if (arg == "--list-analyses") {
      params$list_analyses <- TRUE
    } else if (arg == "--no-plots") {
      params$no_plots <- TRUE
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  params
}


print_available_analyses <- function() {
  specs <- build_te_pca_specs()
  cat("Available TE PCA analyses:\n")
  for (spec in specs) {
    cat("- ", spec$id, ": ", spec$description, "\n", sep = "")
  }
}


relative_to_root <- function(path, project_root) {
  root_norm <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  sub(paste0(root_norm, "/"), "", path_norm, fixed = TRUE)
}


main <- function() {
  params <- parse_args()
  if (params$list_analyses) {
    print_available_analyses()
    return(invisible(NULL))
  }

  required_pca_packages("ggplot2")

  paths <- get_pca_paths()
  ensure_dir(paths$pca_table_dir)
  ensure_dir(paths$pca_figure_dir)

  specs <- select_pca_specs(parse_csv_arg(params$analysis))
  metric_context <- compute_te_context_metrics(paths)

  manifest_rows <- list()
  correlation_rows <- list()
  stability_summary_rows <- list()
  stability_detail_rows <- list()

  for (spec in specs) {
    message("\n--- Running canonical TE PCA: ", spec$id, " ---")

    raw_matrix <- load_breakdown_matrix(paths, spec$level)
    prepared <- prepare_compositional_dataset(raw_matrix, spec)
    result <- run_standard_pca(prepared$clr_matrix)
    stability <- compute_pca_stability(prepared$clr_matrix, result$pca, spec$id)
    correlations <- compute_pc_metric_correlations(result$scores, metric_context, spec$id)

    prefix <- file.path(paths$pca_table_dir, spec$id)
    write_matrix_with_species(prepared$zero_replaced_matrix, paste0(prefix, "_zero_replaced_composition.csv"))
    write_matrix_with_species(prepared$clr_matrix, paste0(prefix, "_clr_matrix.csv"))
    write_csv_if_nonempty(prepared$feature_manifest, paste0(prefix, "_feature_manifest.csv"))
    write_csv_if_nonempty(result$scores, paste0(prefix, "_scores.csv"))
    write_csv_if_nonempty(result$loadings, paste0(prefix, "_loadings.csv"))
    write_csv_if_nonempty(result$variance, paste0(prefix, "_variance.csv"))

    if (!params$no_plots) {
      plot_pca_scree(
        result$variance,
        title = paste(spec$id, "CLR PCA scree plot"),
        output_path = file.path(paths$pca_figure_dir, paste0(spec$id, "_scree_plot.png"))
      )
      plot_pca_scores(
        result$scores,
        result$variance,
        result$loadings,
        title = paste(spec$id, "CLR PCA scores"),
        output_path = file.path(paths$pca_figure_dir, paste0(spec$id, "_scores_pc1_pc2.png"))
      )
    }

    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      analysis_id = spec$id,
      level = spec$level,
      role = spec$role,
      feature_set = spec$feature_set,
      min_presence = spec$min_presence,
      n_species = nrow(prepared$clr_matrix),
      n_features = ncol(prepared$clr_matrix),
      zero_fraction = prepared$zero_fraction,
      zero_replacement = if (prepared$pseudocount > 0) "half_min_positive_global" else "not_needed",
      min_positive_fraction = prepared$min_positive,
      pseudocount_fraction = prepared$pseudocount,
      pc1_variance = result$variance$proportion_variance[result$variance$pc == "PC1"],
      pc2_variance = if ("PC2" %in% result$variance$pc) {
        result$variance$proportion_variance[result$variance$pc == "PC2"]
      } else {
        NA_real_
      },
      source_path = relative_to_root(prepared$source_path, paths$project_root),
      description = spec$description,
      stringsAsFactors = FALSE
    )

    if (nrow(correlations) > 0) {
      correlation_rows[[length(correlation_rows) + 1L]] <- correlations
    }
    if (nrow(stability$summary) > 0) {
      stability_summary_rows[[length(stability_summary_rows) + 1L]] <- stability$summary
      stability_detail_rows[[length(stability_detail_rows) + 1L]] <- stability$detail
    }
  }

  manifest_df <- do.call(rbind, manifest_rows)
  write_csv_if_nonempty(manifest_df, file.path(paths$pca_table_dir, "te_pca_analysis_manifest.csv"))

  if (length(correlation_rows) > 0) {
    write_csv_if_nonempty(
      do.call(rbind, correlation_rows),
      file.path(paths$pca_table_dir, "te_pca_pc_metric_correlations.csv")
    )
  }

  if (length(stability_summary_rows) > 0) {
    write_csv_if_nonempty(
      do.call(rbind, stability_summary_rows),
      file.path(paths$pca_table_dir, "te_pca_stability_summary.csv")
    )
    write_csv_if_nonempty(
      do.call(rbind, stability_detail_rows),
      file.path(paths$pca_table_dir, "te_pca_stability_detail.csv")
    )
  }

  message("\nCanonical TE PCA complete.")
  message("Tables: ", relative_to_root(paths$pca_table_dir, paths$project_root))
  message("Figures: ", relative_to_root(paths$pca_figure_dir, paths$project_root))
}


main()

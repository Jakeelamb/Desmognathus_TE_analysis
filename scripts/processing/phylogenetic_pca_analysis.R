#!/usr/bin/env Rscript

# Canonical supplementary phylogenetic PCA workflow.
#
# This script consumes the exact CLR matrices written by scripts/processing/pca.R
# so the standard PCA and phylogenetic PCA share the same transformed input.

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
    tree = NULL,
    list_analyses = FALSE,
    no_plots = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--analysis") {
      i <- i + 1L
      params$analysis <- args[[i]]
    } else if (arg == "--tree") {
      i <- i + 1L
      params$tree <- args[[i]]
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
  cat("Available TE phylogenetic PCA analyses:\n")
  for (spec in specs) {
    cat("- ", spec$id, " (reads canonical CLR matrix from standard PCA)\n", sep = "")
  }
}


rename_component_columns <- function(df, prefix) {
  component_names <- paste0(prefix, seq_len(ncol(df)))
  names(df) <- component_names
  df
}


extract_ppca_scores <- function(ppca_result, ordered_species) {
  for (candidate in c("S", "scores", "x")) {
    values <- ppca_result[[candidate]]
    if (!is.null(values)) {
      score_df <- as.data.frame(values, check.names = FALSE)
      if (nrow(score_df) == length(ordered_species)) {
        rownames(score_df) <- ordered_species
        score_df <- rename_component_columns(score_df, "pPC")
        score_df$species <- ordered_species
        return(score_df[, c("species", setdiff(names(score_df), "species")), drop = FALSE])
      }
    }
  }
  stop("Unable to extract phylogenetic PCA scores from the phyl.pca result.", call. = FALSE)
}


extract_ppca_loadings <- function(ppca_result) {
  for (candidate in c("L", "Evec", "rotation")) {
    values <- ppca_result[[candidate]]
    if (!is.null(values)) {
      loading_df <- as.data.frame(values, check.names = FALSE)
      loading_df <- rename_component_columns(loading_df, "pPC")
      loading_df$feature <- rownames(loading_df)
      return(loading_df[, c("feature", setdiff(names(loading_df), "feature")), drop = FALSE])
    }
  }
  stop("Unable to extract phylogenetic PCA loadings from the phyl.pca result.", call. = FALSE)
}


extract_ppca_variance <- function(ppca_result) {
  eigenvalues <- NULL
  if (!is.null(ppca_result$Eval)) {
    eigenvalues <- ppca_result$Eval
    if (is.matrix(eigenvalues)) {
      eigenvalues <- diag(eigenvalues)
    }
  } else if (!is.null(ppca_result$StDev)) {
    eigenvalues <- ppca_result$StDev ^ 2
  }

  if (is.null(eigenvalues) || length(eigenvalues) == 0) {
    stop("Unable to extract phylogenetic PCA eigenvalues.", call. = FALSE)
  }

  eigenvalues <- as.numeric(eigenvalues)
  eigenvalues[eigenvalues < 0] <- 0
  proportions <- eigenvalues / sum(eigenvalues)
  data.frame(
    component = paste0("pPC", seq_along(eigenvalues)),
    eigenvalue = eigenvalues,
    proportion_variance = proportions,
    cumulative_variance = cumsum(proportions),
    stringsAsFactors = FALSE
  )
}


plot_ppca_scores <- function(scores_df, variance_df, loadings_df, title, output_path) {
  required_pca_packages("ggplot2")

  if (!all(c("pPC1", "pPC2") %in% names(scores_df))) {
    return(invisible(NULL))
  }

  top_pc1 <- paste(extract_top_features(loadings_df, "pPC1"), collapse = ", ")
  top_pc2 <- paste(extract_top_features(loadings_df, "pPC2"), collapse = ", ")
  pc1_var <- variance_df$proportion_variance[variance_df$component == "pPC1"] * 100
  pc2_var <- variance_df$proportion_variance[variance_df$component == "pPC2"] * 100

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = pPC1, y = pPC2, label = species)) +
    ggplot2::geom_point(color = "firebrick", alpha = 0.85, size = 2.4)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
  } else {
    p <- p + ggplot2::geom_text(size = 2.6, check_overlap = TRUE, vjust = -0.4)
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = sprintf("pPC1 (%.1f%%) | top loadings: %s", pc1_var, top_pc1),
      y = sprintf("pPC2 (%.1f%%) | top loadings: %s", pc2_var, top_pc2)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::coord_equal()

  ggplot2::ggsave(output_path, p, width = 9.5, height = 7.5, dpi = 300)
  invisible(output_path)
}


plot_phylomorphospace <- function(tree, scores_df, variance_df, title, output_path) {
  if (!all(c("pPC1", "pPC2") %in% names(scores_df))) {
    return(invisible(NULL))
  }

  xy <- as.matrix(scores_df[, c("pPC1", "pPC2"), drop = FALSE])
  rownames(xy) <- scores_df$species
  pc1_var <- variance_df$proportion_variance[variance_df$component == "pPC1"] * 100
  pc2_var <- variance_df$proportion_variance[variance_df$component == "pPC2"] * 100

  grDevices::png(output_path, width = 2400, height = 1900, res = 300)
  phytools::phylomorphospace(
    tree,
    xy,
    label = "horizontal",
    node.size = c(0, 0),
    xlab = sprintf("pPC1 (%.1f%%)", pc1_var),
    ylab = sprintf("pPC2 (%.1f%%)", pc2_var),
    main = title
  )
  graphics::text(
    xy[, 1],
    xy[, 2],
    labels = rownames(xy),
    pos = 4,
    cex = 0.55
  )
  grDevices::dev.off()
  invisible(output_path)
}


run_phylogenetic_pca <- function(clr_matrix, tree_file) {
  tree <- ape::read.tree(tree_file)
  tree <- standardize_tree_labels(tree)

  common <- intersect(tree$tip.label, rownames(clr_matrix))
  if (length(common) < 4) {
    stop("Need at least four shared taxa to run phylogenetic PCA.", call. = FALSE)
  }

  pruned_tree <- ape::keep.tip(tree, common)
  ordered_species <- pruned_tree$tip.label
  ordered_matrix <- clr_matrix[ordered_species, , drop = FALSE]

  mode_used <- "cov"
  ppca <- tryCatch(
    phytools::phyl.pca(pruned_tree, ordered_matrix, method = "BM", mode = "cov"),
    error = function(...) {
      mode_used <<- "corr"
      phytools::phyl.pca(pruned_tree, ordered_matrix, method = "BM", mode = "corr")
    }
  )

  list(
    tree = pruned_tree,
    ordered_species = ordered_species,
    mode_used = mode_used,
    result = ppca,
    scores = extract_ppca_scores(ppca, ordered_species),
    loadings = extract_ppca_loadings(ppca),
    variance = extract_ppca_variance(ppca)
  )
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

  required_pca_packages(c("ape", "phytools"))

  paths <- get_pca_paths()
  ensure_dir(paths$pca_table_dir)
  ensure_dir(paths$pca_figure_dir)
  tree_file <- resolve_tree_file(paths, params$tree)
  specs <- select_pca_specs(parse_csv_arg(params$analysis))

  manifest_rows <- list()

  for (spec in specs) {
    message("\n--- Running canonical TE phylogenetic PCA: ", spec$id, " ---")

    clr_matrix <- load_saved_clr_matrix(paths, spec$id)
    ppca <- run_phylogenetic_pca(clr_matrix, tree_file)
    prefix <- file.path(paths$pca_table_dir, spec$id)

    write_csv_if_nonempty(ppca$scores, paste0(prefix, "_ppca_scores.csv"))
    write_csv_if_nonempty(ppca$loadings, paste0(prefix, "_ppca_loadings.csv"))
    write_csv_if_nonempty(ppca$variance, paste0(prefix, "_ppca_variance.csv"))

    if (!params$no_plots) {
      plot_ppca_scores(
        ppca$scores,
        ppca$variance,
        ppca$loadings,
        title = paste(spec$id, "phylogenetic CLR PCA scores"),
        output_path = file.path(paths$pca_figure_dir, paste0(spec$id, "_ppca_scores_pc1_pc2.png"))
      )
      plot_phylomorphospace(
        ppca$tree,
        ppca$scores,
        ppca$variance,
        title = paste(spec$id, "phylogenetic CLR PCA phylomorphospace"),
        output_path = file.path(paths$pca_figure_dir, paste0(spec$id, "_ppca_phylomorphospace.png"))
      )
    }

    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      analysis_id = spec$id,
      n_tree_tips = length(ppca$tree$tip.label),
      n_features = ncol(clr_matrix),
      mode_used = ppca$mode_used,
      tree_file = relative_to_root(tree_file, paths$project_root),
      source_clr_matrix = file.path("results", "tables", "pca", paste0(spec$id, "_clr_matrix.csv")),
      ppc1_variance = ppca$variance$proportion_variance[ppca$variance$component == "pPC1"],
      ppc2_variance = if ("pPC2" %in% ppca$variance$component) {
        ppca$variance$proportion_variance[ppca$variance$component == "pPC2"]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }

  write_csv_if_nonempty(
    do.call(rbind, manifest_rows),
    file.path(paths$pca_table_dir, "te_ppca_analysis_manifest.csv")
  )

  message("\nCanonical TE phylogenetic PCA complete.")
  message("Tree: ", relative_to_root(tree_file, paths$project_root))
  message("Tables: ", relative_to_root(paths$pca_table_dir, paths$project_root))
  message("Figures: ", relative_to_root(paths$pca_figure_dir, paths$project_root))
}


main()

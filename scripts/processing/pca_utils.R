#' Shared utilities for the canonical TE compositional PCA workflow.
#'
#' These helpers intentionally avoid the older raw-proportion PCA assumptions.
#' The canonical workflow is:
#' 1. load frozen TE breakdown tables
#' 2. apply explicit feature filtering
#' 3. replace zeros deterministically
#' 4. close compositions and apply CLR
#' 5. run standard PCA or phylogenetic PCA on the same transformed matrix

required_pca_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}


configure_conda_r_library <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (!nzchar(conda_prefix)) {
    return(invisible(NULL))
  }

  env_library <- file.path(conda_prefix, "lib", "R", "library")
  if (!dir.exists(env_library)) {
    return(invisible(NULL))
  }

  current_paths <- .libPaths()
  if (!(env_library %in% current_paths)) {
    .libPaths(c(env_library, current_paths))
  }
  invisible(env_library)
}


configure_conda_r_library()


find_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  start_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
  } else {
    getwd()
  }

  current <- normalizePath(start_dir, mustWork = FALSE)
  for (i in seq_len(10)) {
    if (
      file.exists(file.path(current, "config", "paths.yaml")) ||
      dir.exists(file.path(current, ".git"))
    ) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }

  stop("Could not determine project root from the current script location.", call. = FALSE)
}


get_pca_paths <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- find_project_root()
  }

  list(
    project_root = project_root,
    order_csv = file.path(project_root, "results", "data", "dnaPipeTE_order_breakdown.csv"),
    superfamily_csv = file.path(project_root, "results", "data", "dnaPipeTE_superfamily_breakdown.csv"),
    order_diversity_csv = file.path(project_root, "results", "data", "diversity_order_stats.csv"),
    pca_table_dir = file.path(project_root, "results", "tables", "pca"),
    pca_figure_dir = file.path(project_root, "results", "figures", "pca"),
    tree_candidates = c(
      file.path(project_root, "results", "data", "desmo900dated_test_cleaned_phylo.tre"),
      file.path(project_root, "input_data", "phylogeny", "desmo900dated_test.tre"),
      file.path(project_root, "results", "phylogeny", "processed_phylogeny.nwk")
    )
  )
}


ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}


clean_species_ids <- function(x) {
  cleaned <- trimws(as.character(x))
  cleaned <- sub("^D\\.\\s*", "", cleaned)
  cleaned
}


load_breakdown_matrix <- function(paths, level = c("order", "superfamily")) {
  level <- match.arg(level)
  csv_path <- switch(
    level,
    order = paths$order_csv,
    superfamily = paths$superfamily_csv
  )

  if (!file.exists(csv_path)) {
    stop("Missing TE breakdown file: ", csv_path, call. = FALSE)
  }

  df <- utils::read.csv(csv_path, check.names = FALSE)
  if (ncol(df) < 2) {
    stop("Breakdown table must contain a species column and at least one feature column: ", csv_path, call. = FALSE)
  }

  species <- clean_species_ids(df[[1]])
  if (anyDuplicated(species)) {
    duplicates <- unique(species[duplicated(species)])
    stop(
      "Species IDs are duplicated after cleaning in ",
      basename(csv_path),
      ": ",
      paste(duplicates, collapse = ", "),
      call. = FALSE
    )
  }

  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- species
  list(
    level = level,
    source_path = csv_path,
    species = species,
    matrix = mat
  )
}


build_te_pca_specs <- function() {
  list(
    list(
      id = "superfamily_primary",
      level = "superfamily",
      role = "primary",
      min_presence = 3L,
      feature_set = "all",
      feature_names = NULL,
      description = "Primary supplementary ordination: superfamily CLR PCA using features present in at least 3 species."
    ),
    list(
      id = "superfamily_presence5",
      level = "superfamily",
      role = "sensitivity",
      min_presence = 5L,
      feature_set = "all",
      feature_names = NULL,
      description = "Sensitivity ordination: superfamily CLR PCA using features present in at least 5 species."
    ),
    list(
      id = "order_major",
      level = "order",
      role = "sensitivity",
      min_presence = 1L,
      feature_set = "major_orders",
      feature_names = c("LTR", "LINE", "TIR", "DIRS", "SINE"),
      description = "Sensitivity ordination: CLR PCA on major TE orders only."
    ),
    list(
      id = "order_all",
      level = "order",
      role = "sensitivity",
      min_presence = 1L,
      feature_set = "all",
      feature_names = NULL,
      description = "Sensitivity ordination: CLR PCA on all order-level TE classes."
    )
  )
}


select_pca_specs <- function(spec_ids = NULL) {
  specs <- build_te_pca_specs()
  if (is.null(spec_ids) || length(spec_ids) == 0) {
    return(specs)
  }

  lookup <- setNames(specs, vapply(specs, `[[`, character(1), "id"))
  missing <- setdiff(spec_ids, names(lookup))
  if (length(missing) > 0) {
    stop("Unknown PCA analysis IDs: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  unname(lookup[spec_ids])
}


parse_csv_arg <- function(value) {
  if (is.null(value) || identical(value, "")) {
    return(character(0))
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}


prepare_compositional_dataset <- function(raw_matrix, spec) {
  mat <- raw_matrix$matrix

  if (!is.null(spec$feature_names)) {
    keep <- intersect(spec$feature_names, colnames(mat))
    if (length(keep) < 2) {
      stop("Analysis ", spec$id, " does not retain enough requested features.", call. = FALSE)
    }
    mat <- mat[, keep, drop = FALSE]
  }

  feature_presence <- colSums(mat > 0, na.rm = TRUE)
  mat <- mat[, feature_presence >= spec$min_presence, drop = FALSE]
  feature_presence <- feature_presence[colnames(mat)]

  if (ncol(mat) < 2) {
    stop("Analysis ", spec$id, " retains fewer than two features after filtering.", call. = FALSE)
  }
  if (any(mat < 0, na.rm = TRUE)) {
    stop("Negative TE proportions found in analysis ", spec$id, ".", call. = FALSE)
  }

  row_totals <- rowSums(mat, na.rm = TRUE)
  if (any(!is.finite(row_totals)) || any(row_totals <= 0)) {
    stop("Non-positive or non-finite row totals found in analysis ", spec$id, ".", call. = FALSE)
  }

  closed <- sweep(mat, 1, row_totals, "/")
  positive <- closed[closed > 0 & is.finite(closed)]
  if (length(positive) == 0) {
    stop("No positive values remain for analysis ", spec$id, ".", call. = FALSE)
  }

  min_positive <- min(positive)
  has_zero <- any(closed == 0, na.rm = TRUE)
  pseudocount <- if (has_zero) min_positive / 2 else 0

  zero_replaced <- closed
  if (has_zero) {
    zero_replaced[zero_replaced == 0] <- pseudocount
    zero_replaced <- sweep(zero_replaced, 1, rowSums(zero_replaced), "/")
  }

  clr <- log(zero_replaced)
  clr <- clr - rowMeans(clr)

  feature_manifest <- data.frame(
    analysis_id = spec$id,
    feature = colnames(mat),
    presence_n_species = as.integer(feature_presence),
    mean_percent = colMeans(mat),
    stringsAsFactors = FALSE
  )

  list(
    spec = spec,
    species = rownames(mat),
    original_matrix = mat,
    closed_matrix = closed,
    zero_replaced_matrix = zero_replaced,
    clr_matrix = clr,
    feature_manifest = feature_manifest,
    min_positive = min_positive,
    pseudocount = pseudocount,
    zero_fraction = mean(closed == 0),
    source_path = raw_matrix$source_path
  )
}


run_standard_pca <- function(clr_matrix) {
  pca_result <- stats::prcomp(clr_matrix, center = TRUE, scale. = FALSE)
  variance <- (pca_result$sdev ^ 2)
  variance <- variance / sum(variance)
  variance_df <- data.frame(
    pc = paste0("PC", seq_along(variance)),
    eigenvalue = pca_result$sdev ^ 2,
    proportion_variance = variance,
    cumulative_variance = cumsum(variance),
    stringsAsFactors = FALSE
  )

  scores <- as.data.frame(pca_result$x, check.names = FALSE)
  scores$species <- rownames(clr_matrix)
  scores <- scores[, c("species", setdiff(names(scores), "species")), drop = FALSE]

  loadings <- as.data.frame(pca_result$rotation, check.names = FALSE)
  loadings$feature <- rownames(pca_result$rotation)
  loadings <- loadings[, c("feature", setdiff(names(loadings), "feature")), drop = FALSE]

  list(
    pca = pca_result,
    scores = scores,
    loadings = loadings,
    variance = variance_df
  )
}


extract_top_features <- function(loadings_df, pc = "PC1", n = 3L) {
  if (!pc %in% names(loadings_df)) {
    return(character(0))
  }
  ranked <- order(abs(loadings_df[[pc]]), decreasing = TRUE)
  loadings_df$feature[ranked][seq_len(min(n, nrow(loadings_df)))]
}


plot_pca_scree <- function(variance_df, title, output_path) {
  required_pca_packages("ggplot2")

  plot_df <- variance_df
  plot_df$pc <- factor(plot_df$pc, levels = plot_df$pc)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pc)) +
    ggplot2::geom_col(ggplot2::aes(y = proportion_variance * 100), fill = "steelblue", alpha = 0.85) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = cumulative_variance * 100,
        group = 1
      ),
      color = "firebrick",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = cumulative_variance * 100),
      color = "firebrick",
      size = 1.8
    ) +
    ggplot2::geom_hline(yintercept = 80, linetype = "dashed", color = "gray45") +
    ggplot2::labs(
      title = title,
      x = "Principal component",
      y = "Variance explained (%)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(output_path, p, width = 9, height = 5.5, dpi = 300)
  invisible(output_path)
}


plot_pca_scores <- function(scores_df, variance_df, loadings_df, title, output_path) {
  required_pca_packages("ggplot2")

  if (!all(c("PC1", "PC2") %in% names(scores_df))) {
    return(invisible(NULL))
  }

  top_pc1 <- paste(extract_top_features(loadings_df, "PC1"), collapse = ", ")
  top_pc2 <- paste(extract_top_features(loadings_df, "PC2"), collapse = ", ")
  pc1_var <- variance_df$proportion_variance[variance_df$pc == "PC1"] * 100
  pc2_var <- variance_df$proportion_variance[variance_df$pc == "PC2"] * 100

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = PC1, y = PC2, label = species)) +
    ggplot2::geom_point(color = "steelblue", alpha = 0.85, size = 2.4)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
  } else {
    p <- p + ggplot2::geom_text(size = 2.6, check_overlap = TRUE, vjust = -0.4)
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = sprintf("PC1 (%.1f%%) | top loadings: %s", pc1_var, top_pc1),
      y = sprintf("PC2 (%.1f%%) | top loadings: %s", pc2_var, top_pc2)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::coord_equal()

  ggplot2::ggsave(output_path, p, width = 9.5, height = 7.5, dpi = 300)
  invisible(output_path)
}


compute_te_context_metrics <- function(paths) {
  order_df <- utils::read.csv(paths$order_csv, check.names = FALSE)
  diversity_df <- utils::read.csv(paths$order_diversity_csv, check.names = FALSE)

  species <- clean_species_ids(order_df[[1]])
  order_mat <- as.matrix(order_df[, -1, drop = FALSE])
  storage.mode(order_mat) <- "numeric"
  rownames(order_mat) <- species
  closed <- sweep(order_mat, 1, rowSums(order_mat), "/")

  positive <- closed[closed > 0 & is.finite(closed)]
  pseudocount <- min(positive) / 2

  ltr <- closed[, "LTR"]
  line <- closed[, "LINE"]
  tir <- closed[, "TIR"]
  retro <- closed[, c("LTR", "LINE", "SINE", "PLE"), drop = FALSE]
  dna <- closed[, c("TIR", "Maverick", "Helitron"), drop = FALSE]

  metrics <- data.frame(
    species = rownames(order_mat),
    ltr_line_logratio = log((ltr + pseudocount) / (line + pseudocount)),
    retro_dna_logratio = log((rowSums(retro) + pseudocount) / (rowSums(dna) + pseudocount)),
    stringsAsFactors = FALSE
  )

  diversity_species <- clean_species_ids(diversity_df[[1]])
  diversity_lookup <- data.frame(
    species = diversity_species,
    order_pielou = diversity_df$Pielou_Evenness,
    stringsAsFactors = FALSE
  )

  merge(metrics, diversity_lookup, by = "species", all.x = TRUE, sort = FALSE)
}


compute_pc_metric_correlations <- function(scores_df, metrics_df, analysis_id) {
  merged <- merge(scores_df, metrics_df, by = "species", all.x = FALSE, all.y = FALSE, sort = FALSE)
  rows <- list()

  for (pc in intersect(c("PC1", "PC2"), names(merged))) {
    for (metric in c("ltr_line_logratio", "retro_dna_logratio", "order_pielou")) {
      if (!metric %in% names(merged)) {
        next
      }
      keep <- stats::complete.cases(merged[, c(pc, metric)])
      if (sum(keep) < 4) {
        next
      }
      pearson <- suppressWarnings(stats::cor(merged[[pc]][keep], merged[[metric]][keep], method = "pearson"))
      spearman <- suppressWarnings(stats::cor(merged[[pc]][keep], merged[[metric]][keep], method = "spearman"))
      rows[[length(rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        pc = pc,
        metric = metric,
        n_species = sum(keep),
        pearson_r = pearson,
        spearman_rho = spearman,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}


compute_pca_stability <- function(clr_matrix, full_pca, analysis_id, max_pc = 2L) {
  if (nrow(clr_matrix) < 6) {
    return(list(summary = data.frame(), detail = data.frame()))
  }

  detail_rows <- list()
  max_pc <- min(max_pc, ncol(full_pca$rotation), ncol(full_pca$x))

  for (i in seq_len(nrow(clr_matrix))) {
    species_out <- rownames(clr_matrix)[i]
    loo_matrix <- clr_matrix[-i, , drop = FALSE]
    loo_pca <- stats::prcomp(loo_matrix, center = TRUE, scale. = FALSE)

    for (pc_index in seq_len(max_pc)) {
      pc_name <- paste0("PC", pc_index)
      full_loading <- full_pca$rotation[, pc_index]
      loo_loading <- loo_pca$rotation[names(full_loading), pc_index]
      loading_cor_raw <- suppressWarnings(stats::cor(full_loading, loo_loading))
      sign_flip <- if (is.na(loading_cor_raw) || loading_cor_raw >= 0) 1 else -1
      loo_loading <- loo_loading * sign_flip
      aligned_loading_cor <- suppressWarnings(stats::cor(full_loading, loo_loading))

      full_scores <- full_pca$x[rownames(loo_matrix), pc_index]
      loo_scores <- loo_pca$x[, pc_index] * sign_flip
      score_cor <- suppressWarnings(stats::cor(full_scores, loo_scores))

      detail_rows[[length(detail_rows) + 1L]] <- data.frame(
        analysis_id = analysis_id,
        species_left_out = species_out,
        pc = pc_name,
        loading_correlation = aligned_loading_cor,
        score_correlation = score_cor,
        stringsAsFactors = FALSE
      )
    }
  }

  detail_df <- do.call(rbind, detail_rows)
  summary_df <- stats::aggregate(
    detail_df[, c("loading_correlation", "score_correlation")],
    by = list(analysis_id = detail_df$analysis_id, pc = detail_df$pc),
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), min = min(x, na.rm = TRUE))
  )

  summary_df <- do.call(data.frame, summary_df)
  names(summary_df) <- c(
    "analysis_id",
    "pc",
    "loading_correlation_mean",
    "loading_correlation_min",
    "score_correlation_mean",
    "score_correlation_min"
  )

  list(summary = summary_df, detail = detail_df)
}


write_csv_if_nonempty <- function(df, path) {
  if (is.null(df)) {
    return(invisible(NULL))
  }
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}


write_matrix_with_species <- function(mat, path) {
  df <- data.frame(species = rownames(mat), mat, check.names = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}


load_saved_clr_matrix <- function(paths, analysis_id) {
  csv_path <- file.path(paths$pca_table_dir, paste0(analysis_id, "_clr_matrix.csv"))
  if (!file.exists(csv_path)) {
    stop(
      "Missing canonical CLR matrix for analysis '",
      analysis_id,
      "'. Run scripts/processing/pca.R first.",
      call. = FALSE
    )
  }

  df <- utils::read.csv(csv_path, check.names = FALSE)
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- clean_species_ids(df[[1]])
  mat
}


resolve_tree_file <- function(paths, explicit_tree = NULL) {
  candidates <- c(explicit_tree, paths$tree_candidates)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("No phylogeny file found for phylogenetic PCA.", call. = FALSE)
  }
  existing[[1]]
}


standardize_tree_labels <- function(tree) {
  tree$tip.label <- clean_species_ids(tree$tip.label)
  tree
}

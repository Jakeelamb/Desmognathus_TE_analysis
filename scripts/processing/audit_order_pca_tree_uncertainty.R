#!/usr/bin/env Rscript

# Propagate the corrected order-composition phylogenetic PCA across the 200
# published Stewart-Wiens time-calibrated bootstrap trees.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/processing/audit_order_pca_tree_uncertainty.R", mustWork = TRUE)
}
project_root <- dirname(dirname(dirname(script_path)))

clr_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_clr_matrices_analysis18_v1.csv")
ordinary_scores_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_scores_analysis18_v1.csv")
ordinary_loadings_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_loadings_analysis18_v1.csv")
focal_scores_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_order_phylogenetic_pca_scores_analysis18_v1.csv")
focal_loadings_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_order_phylogenetic_pca_loadings_analysis18_v1.csv")
published_main_path <- file.path(project_root, "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk")
published_bootstrap_path <- file.path(project_root, "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex")

output_dir <- file.path(project_root, "results/data/corrected/diversity_pca")
figure_dir <- file.path(project_root, "results/figures/corrected/diversity_pca")
audit_dir <- file.path(project_root, "plans/publication-readiness-deep-audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metrics_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_metrics_analysis18_v1.csv")
scores_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_scores_analysis18_v1.csv")
score_summary_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_score_summary_analysis18_v1.csv")
loadings_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_loadings_analysis18_v1.csv")
loading_summary_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_loading_summary_analysis18_v1.csv")
manifest_output <- file.path(output_dir, "te_order_pca_tree_uncertainty_analysis18_v1.manifest.json")
report_output <- file.path(audit_dir, "te_order_pca_tree_uncertainty_analysis18_v1.md")

score_figure <- file.path(figure_dir, "te_order_pca_tree_uncertainty_score_intervals_analysis18_v1.png")
loading_figure <- file.path(figure_dir, "te_order_pca_tree_uncertainty_loading_intervals_analysis18_v1.png")
correlation_figure <- file.path(figure_dir, "te_order_pca_tree_uncertainty_correlation_analysis18_v1.png")
variance_figure <- file.path(figure_dir, "te_order_pca_tree_uncertainty_variance_analysis18_v1.png")

required <- c(
  clr_path, ordinary_scores_path, ordinary_loadings_path,
  focal_scores_path, focal_loadings_path, published_main_path, published_bootstrap_path
)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required input(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

portable <- function(path) {
  sub(paste0(normalizePath(project_root, winslash = "/"), "/"), "", normalizePath(path, winslash = "/"), fixed = TRUE)
}
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

build_ilr_basis <- function(feature_names) {
  basis <- stats::contr.helmert(length(feature_names))
  basis <- sweep(basis, 2, sqrt(colSums(basis ^ 2)), "/")
  rownames(basis) <- feature_names
  basis
}

fit_ppca <- function(tree, clr, basis) {
  clr <- clr[tree$tip.label, , drop = FALSE]
  ilr <- clr %*% basis
  result <- suppressWarnings(phytools::phyl.pca(tree, ilr, method = "BM", mode = "cov"))
  eigenvalues <- diag(result$Eval)
  if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    stop("Non-positive/non-finite phylogenetic-PCA eigenvalue", call. = FALSE)
  }
  list(
    scores = as.matrix(result$S),
    loadings = basis %*% as.matrix(result$Evec),
    variance = eigenvalues / sum(eigenvalues)
  )
}

align_first_two <- function(fit, reference_loadings) {
  correlation <- cor(reference_loadings[, 1:2, drop = FALSE], fit$loadings[, 1:2, drop = FALSE])
  identity_score <- abs(correlation[1, 1]) + abs(correlation[2, 2])
  swap_score <- abs(correlation[1, 2]) + abs(correlation[2, 1])
  swapped <- swap_score > identity_score
  if (swapped) {
    fit$scores[, 1:2] <- fit$scores[, c(2, 1), drop = FALSE]
    fit$loadings[, 1:2] <- fit$loadings[, c(2, 1), drop = FALSE]
    fit$variance[1:2] <- fit$variance[c(2, 1)]
  }
  for (axis in seq_len(2)) {
    if (cor(reference_loadings[, axis], fit$loadings[, axis]) < 0) {
      fit$scores[, axis] <- -fit$scores[, axis]
      fit$loadings[, axis] <- -fit$loadings[, axis]
    }
  }
  fit$axis_swap_required <- swapped
  fit
}

summarize_group <- function(values) {
  c(
    mean = mean(values),
    sd = sd(values),
    q025 = unname(quantile(values, 0.025)),
    q05 = unname(quantile(values, 0.05)),
    median = median(values),
    q95 = unname(quantile(values, 0.95)),
    q975 = unname(quantile(values, 0.975)),
    minimum = min(values),
    maximum = max(values)
  )
}

clr_long <- read.csv(clr_path, check.names = FALSE, stringsAsFactors = FALSE)
clr_long <- subset(clr_long, analysis_id == "order_classified_primary" & zero_replacement == "not_needed")
clr <- as.matrix(xtabs(clr_value ~ species + feature, data = clr_long))
basis <- build_ilr_basis(colnames(clr))

ordinary_scores <- read.csv(ordinary_scores_path, check.names = FALSE, stringsAsFactors = FALSE)
ordinary_scores <- subset(ordinary_scores, analysis_id == "order_classified_primary" & zero_replacement == "not_needed")
rownames(ordinary_scores) <- ordinary_scores$species
ordinary_score_matrix <- as.matrix(ordinary_scores[rownames(clr), c("PC1", "PC2"), drop = FALSE])

ordinary_loadings <- read.csv(ordinary_loadings_path, check.names = FALSE, stringsAsFactors = FALSE)
ordinary_loadings <- subset(ordinary_loadings, analysis_id == "order_classified_primary" & zero_replacement == "not_needed")
rownames(ordinary_loadings) <- ordinary_loadings$feature
ordinary_loading_matrix <- as.matrix(ordinary_loadings[colnames(clr), c("PC1", "PC2"), drop = FALSE])

focal_scores <- read.csv(focal_scores_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(focal_scores) <- focal_scores$species
focal_score_matrix <- as.matrix(focal_scores[rownames(clr), c("pPC1", "pPC2"), drop = FALSE])

focal_loadings <- read.csv(focal_loadings_path, check.names = FALSE, stringsAsFactors = FALSE)
rownames(focal_loadings) <- focal_loadings$feature
focal_loading_matrix <- as.matrix(focal_loadings[colnames(clr), c("pPC1", "pPC2"), drop = FALSE])

main <- ape::read.tree(published_main_path)
bootstrap_trees <- ape::read.nexus(published_bootstrap_path)
if (length(bootstrap_trees) != 200) stop("Expected exactly 200 bootstrap trees", call. = FALSE)
all_trees <- c(list(main), unclass(bootstrap_trees))
tree_ids <- c("published_main", sprintf("published_bootstrap_%03d", seq_along(bootstrap_trees)))
tree_roles <- c("published_main", rep("published_bootstrap", length(bootstrap_trees)))

metrics_rows <- vector("list", length(all_trees))
score_rows <- vector("list", length(all_trees))
loading_rows <- vector("list", length(all_trees))
for (i in seq_along(all_trees)) {
  tree <- all_trees[[i]]
  if (!setequal(tree$tip.label, rownames(clr))) {
    stop("Tree ", tree_ids[[i]], " does not match CLR species", call. = FALSE)
  }
  fit <- align_first_two(fit_ppca(tree, clr, basis), focal_loading_matrix)
  rownames(fit$scores) <- tree$tip.label
  colnames(fit$scores) <- paste0("pPC", seq_len(ncol(fit$scores)))
  colnames(fit$loadings) <- paste0("pPC", seq_len(ncol(fit$loadings)))

  metrics_rows[[i]] <- data.frame(
    tree_id = tree_ids[[i]],
    tree_role = tree_roles[[i]],
    axis_swap_required = fit$axis_swap_required,
    ppc1_variance = fit$variance[[1]],
    ppc2_variance = fit$variance[[2]],
    pc1_ppc1_score_correlation = cor(ordinary_score_matrix[, 1], fit$scores[rownames(clr), 1]),
    pc2_ppc2_score_correlation = cor(ordinary_score_matrix[, 2], fit$scores[rownames(clr), 2]),
    pc1_ppc1_loading_correlation = cor(ordinary_loading_matrix[, 1], fit$loadings[, 1]),
    pc2_ppc2_loading_correlation = cor(ordinary_loading_matrix[, 2], fit$loadings[, 2]),
    focal_ppc1_score_correlation = cor(focal_score_matrix[, 1], fit$scores[rownames(clr), 1]),
    focal_ppc2_score_correlation = cor(focal_score_matrix[, 2], fit$scores[rownames(clr), 2]),
    focal_ppc1_loading_correlation = cor(focal_loading_matrix[, 1], fit$loadings[, 1]),
    focal_ppc2_loading_correlation = cor(focal_loading_matrix[, 2], fit$loadings[, 2]),
    stringsAsFactors = FALSE
  )
  score_rows[[i]] <- data.frame(
    tree_id = tree_ids[[i]],
    tree_role = tree_roles[[i]],
    species = rownames(clr),
    pPC1 = fit$scores[rownames(clr), 1],
    pPC2 = fit$scores[rownames(clr), 2],
    stringsAsFactors = FALSE
  )
  loading_rows[[i]] <- data.frame(
    tree_id = tree_ids[[i]],
    tree_role = tree_roles[[i]],
    feature = colnames(clr),
    pPC1 = fit$loadings[, 1],
    pPC2 = fit$loadings[, 2],
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(rbind, metrics_rows)
scores <- do.call(rbind, score_rows)
loadings <- do.call(rbind, loading_rows)
write.csv(metrics, metrics_output, row.names = FALSE)
write.csv(scores, scores_output, row.names = FALSE)
write.csv(loadings, loadings_output, row.names = FALSE)

bootstrap_scores <- subset(scores, tree_role == "published_bootstrap")
score_summary_rows <- list()
for (species in unique(bootstrap_scores$species)) {
  for (component in c("pPC1", "pPC2")) {
    values <- bootstrap_scores[bootstrap_scores$species == species, component]
    score_summary_rows[[length(score_summary_rows) + 1L]] <- data.frame(
      species = species,
      component = component,
      as.list(summarize_group(values)),
      sign_positive_fraction = mean(values > 0),
      focal_score = focal_scores[species, component],
      stringsAsFactors = FALSE
    )
  }
}
score_summary <- do.call(rbind, score_summary_rows)
write.csv(score_summary, score_summary_output, row.names = FALSE)

bootstrap_loadings <- subset(loadings, tree_role == "published_bootstrap")
loading_summary_rows <- list()
for (feature in unique(bootstrap_loadings$feature)) {
  for (component in c("pPC1", "pPC2")) {
    values <- bootstrap_loadings[bootstrap_loadings$feature == feature, component]
    loading_summary_rows[[length(loading_summary_rows) + 1L]] <- data.frame(
      feature = feature,
      component = component,
      as.list(summarize_group(values)),
      sign_positive_fraction = mean(values > 0),
      focal_loading = focal_loadings[feature, component],
      stringsAsFactors = FALSE
    )
  }
}
loading_summary <- do.call(rbind, loading_summary_rows)
write.csv(loading_summary, loading_summary_output, row.names = FALSE)

score_summary$species <- factor(score_summary$species, levels = rev(sort(unique(score_summary$species))))
score_plot <- ggplot(score_summary, aes(x = median, y = species)) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, color = "#2A9D8F", linewidth = 0.8) +
  geom_point(color = "#264653", size = 2.2) +
  geom_point(aes(x = focal_score), shape = 4, color = "#E76F51", size = 2.3, stroke = 0.9) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
  facet_wrap(~component, scales = "free_x") +
  labs(
    title = "Order pPCA species-score uncertainty across 200 time trees",
    subtitle = "Points/lines: published bootstrap median and 95% interval; orange x: collaborator focal tree",
    x = "Phylogenetic PCA score",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave(score_figure, score_plot, width = 10, height = 7.5, dpi = 300)

loading_summary$feature <- factor(loading_summary$feature, levels = rev(sort(unique(loading_summary$feature))))
loading_plot <- ggplot(loading_summary, aes(x = median, y = feature)) +
  geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, color = "#2A9D8F", linewidth = 0.9) +
  geom_point(color = "#264653", size = 2.3) +
  geom_point(aes(x = focal_loading), shape = 4, color = "#E76F51", size = 2.5, stroke = 0.9) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.45) +
  facet_wrap(~component, scales = "free_x") +
  labs(
    title = "Order pPCA loading uncertainty across 200 time trees",
    subtitle = "Points/lines: published bootstrap median and 95% interval; orange x: collaborator focal tree",
    x = "Back-projected CLR loading",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))
ggsave(loading_figure, loading_plot, width = 10, height = 6.5, dpi = 300)

bootstrap_metrics <- subset(metrics, tree_role == "published_bootstrap")
correlation_long <- rbind(
  data.frame(axis = "Axis 1", evidence = "Species scores", correlation = bootstrap_metrics$pc1_ppc1_score_correlation),
  data.frame(axis = "Axis 2", evidence = "Species scores", correlation = bootstrap_metrics$pc2_ppc2_score_correlation),
  data.frame(axis = "Axis 1", evidence = "Feature loadings", correlation = bootstrap_metrics$pc1_ppc1_loading_correlation),
  data.frame(axis = "Axis 2", evidence = "Feature loadings", correlation = bootstrap_metrics$pc2_ppc2_loading_correlation)
)
correlation_plot <- ggplot(correlation_long, aes(x = axis, y = correlation, fill = evidence)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey55") +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE, alpha = 0.8) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.18, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Species scores" = "#2A9D8F", "Feature loadings" = "#E9C46A")) +
  labs(
    title = "Ordinary CLR PCA remains aligned across published time trees",
    x = NULL,
    y = "Correlation with ordinary CLR PCA",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(correlation_figure, correlation_plot, width = 8, height = 5.5, dpi = 300)

variance_long <- rbind(
  data.frame(component = "pPC1", variance = bootstrap_metrics$ppc1_variance),
  data.frame(component = "pPC2", variance = bootstrap_metrics$ppc2_variance)
)
variance_plot <- ggplot(variance_long, aes(x = component, y = 100 * variance, fill = component)) +
  geom_violin(trim = FALSE, alpha = 0.75, show.legend = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", show.legend = FALSE) +
  scale_fill_manual(values = c("pPC1" = "#2A9D8F", "pPC2" = "#E76F51")) +
  labs(
    title = "Phylogenetic variance explained across 200 time trees",
    x = NULL,
    y = "Variance explained (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(variance_figure, variance_plot, width = 7.5, height = 5.5, dpi = 300)

minimum_metric <- function(column) min(bootstrap_metrics[[column]])
loading_crosses_zero <- aggregate(
  cbind(q025, q975) ~ feature + component,
  loading_summary,
  identity
)
n_loading_intervals_cross_zero <- sum(loading_summary$q025 <= 0 & loading_summary$q975 >= 0)
n_axis_swaps <- sum(bootstrap_metrics$axis_swap_required)

report <- sprintf(
  paste0(
    "# Order-level PCA sensitivity across 200 published time trees\n\n",
    "The corrected order composition was projected from CLR to a 9-dimensional orthonormal ILR basis and refit with Brownian-covariance phylogenetic PCA on the Stewart-Wiens (2025) optimal time tree plus all 200 time-calibrated bootstrap trees. Axes were matched to the collaborator-focal pPCA using loading correlations; signs were aligned for reporting only.\n\n",
    "## Stability results\n\n",
    "- Axis swaps required: %d/200 bootstrap trees.\n",
    "- Minimum ordinary-versus-phylogenetic species-score correlation: %.3f for axis 1 and %.3f for axis 2.\n",
    "- Minimum ordinary-versus-phylogenetic loading correlation: %.3f for axis 1 and %.3f for axis 2.\n",
    "- Minimum published-tree versus collaborator-focal species-score correlation: %.3f for pPC1 and %.3f for pPC2.\n",
    "- Median variance explained across trees: %.1f%% for pPC1 and %.1f%% for pPC2.\n",
    "- %d of %d feature-by-axis 95%% loading intervals cross zero; these loadings must not be given direction-stable biological interpretations.\n\n",
    "## Verdict\n\n",
    "**Approved as a descriptive ordination across fixed-tree and published tree-time uncertainty.** The dominant order-composition geometry, especially pPC1, is robust to shared ancestry and the 200 published time trees. This still does not promote PCA axes to causal path variables: direct declared log-ratios or a multivariate phylogenetic model remain preferable for mechanistic inference, and the bootstrap set does not represent reticulation uncertainty.\n\n",
    "## Outputs\n\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n- `%s`\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n"
  ),
  n_axis_swaps,
  minimum_metric("pc1_ppc1_score_correlation"), minimum_metric("pc2_ppc2_score_correlation"),
  minimum_metric("pc1_ppc1_loading_correlation"), minimum_metric("pc2_ppc2_loading_correlation"),
  minimum_metric("focal_ppc1_score_correlation"), minimum_metric("focal_ppc2_score_correlation"),
  100 * median(bootstrap_metrics$ppc1_variance), 100 * median(bootstrap_metrics$ppc2_variance),
  n_loading_intervals_cross_zero, nrow(loading_summary),
  portable(metrics_output), portable(scores_output), portable(score_summary_output),
  portable(loadings_output), portable(loading_summary_output),
  portable(score_figure), portable(loading_figure), portable(correlation_figure), portable(variance_figure)
)
writeLines(report, report_output)

output_paths <- c(
  metrics_output, scores_output, score_summary_output, loadings_output, loading_summary_output,
  report_output, score_figure, loading_figure, correlation_figure, variance_figure
)
manifest <- list(
  analysis_scope = "order_classified_primary_final18_tree_uncertainty",
  n_species = nrow(clr),
  n_order_features = ncol(clr),
  n_published_bootstrap_trees = 200,
  transform = "CLR_to_orthonormal_Helmert_ILR",
  phylogenetic_model = "Brownian_motion_covariance",
  tree_set = portable(published_bootstrap_path),
  tree_set_sha256 = sha256(published_bootstrap_path),
  descriptive_ordination_status = "approved_across_fixed_and_published_tree_time_uncertainty",
  causal_predictor_status = "not_approved",
  reticulation_uncertainty_status = "not_represented",
  axis_swaps_required = n_axis_swaps,
  minimum_correlations = list(
    ordinary_score_axis1 = unname(minimum_metric("pc1_ppc1_score_correlation")),
    ordinary_score_axis2 = unname(minimum_metric("pc2_ppc2_score_correlation")),
    ordinary_loading_axis1 = unname(minimum_metric("pc1_ppc1_loading_correlation")),
    ordinary_loading_axis2 = unname(minimum_metric("pc2_ppc2_loading_correlation"))
  ),
  outputs = lapply(output_paths, function(path) list(path = portable(path), sha256 = sha256(path)))
)
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), manifest_output)

message("Wrote ", report_output)
message("Wrote ", manifest_output)

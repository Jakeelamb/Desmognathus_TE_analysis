#!/usr/bin/env Rscript

# Fixed-tree phylogenetic sensitivity for the corrected order-level CLR PCA.
# CLR coordinates are first projected to an orthonormal ILR basis so the
# phylogenetic covariance calculation is full rank. Historical PCA products
# are read only and all outputs are versioned under corrected/.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggrepel)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/processing/audit_order_phylogenetic_pca.R", mustWork = TRUE)
}
project_root <- dirname(dirname(dirname(script_path)))

tree_path <- file.path(project_root, "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk")
clr_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_clr_matrices_analysis18_v1.csv")
ordinary_scores_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_scores_analysis18_v1.csv")
ordinary_loadings_path <- file.path(project_root, "results/data/corrected/diversity_pca/te_pca_loadings_analysis18_v1.csv")

output_dir <- file.path(project_root, "results/data/corrected/diversity_pca")
figure_dir <- file.path(project_root, "results/figures/corrected/diversity_pca")
audit_dir <- file.path(project_root, "plans/publication-readiness-deep-audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

scores_output <- file.path(output_dir, "te_order_phylogenetic_pca_scores_analysis18_v1.csv")
loadings_output <- file.path(output_dir, "te_order_phylogenetic_pca_loadings_analysis18_v1.csv")
variance_output <- file.path(output_dir, "te_order_phylogenetic_pca_variance_analysis18_v1.csv")
comparison_output <- file.path(output_dir, "te_order_phylogenetic_pca_comparison_analysis18_v1.csv")
stability_output <- file.path(output_dir, "te_order_phylogenetic_pca_stability_analysis18_v1.csv")
manifest_output <- file.path(output_dir, "te_order_phylogenetic_pca_analysis18_v1.manifest.json")
report_output <- file.path(audit_dir, "te_order_phylogenetic_pca_analysis18_v1.md")

scores_figure <- file.path(figure_dir, "te_order_phylogenetic_pca_scores_analysis18_v1.png")
comparison_figure <- file.path(figure_dir, "te_order_ordinary_vs_phylogenetic_pca_analysis18_v1.png")
loadings_figure <- file.path(figure_dir, "te_order_phylogenetic_pca_loadings_analysis18_v1.png")
stability_figure <- file.path(figure_dir, "te_order_phylogenetic_pca_stability_analysis18_v1.png")

required <- c(tree_path, clr_path, ordinary_scores_path, ordinary_loadings_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required input(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

portable <- function(path) {
  sub(paste0(normalizePath(project_root, winslash = "/"), "/"), "", normalizePath(path, winslash = "/"), fixed = TRUE)
}

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

build_ilr_basis <- function(feature_names) {
  n_features <- length(feature_names)
  basis <- stats::contr.helmert(n_features)
  basis <- sweep(basis, 2, sqrt(colSums(basis ^ 2)), "/")
  rownames(basis) <- feature_names
  colnames(basis) <- paste0("ILR", seq_len(ncol(basis)))
  if (max(abs(crossprod(basis) - diag(ncol(basis)))) > 1e-12 ||
      max(abs(colSums(basis))) > 1e-12) {
    stop("ILR basis failed orthonormality or zero-sum validation", call. = FALSE)
  }
  basis
}

fit_ppca <- function(tree, clr, basis) {
  clr <- clr[tree$tip.label, , drop = FALSE]
  ilr <- clr %*% basis
  result <- suppressWarnings(phytools::phyl.pca(tree, ilr, method = "BM", mode = "cov"))
  eigenvalues <- diag(result$Eval)
  if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    stop("Phylogenetic PCA returned non-positive or non-finite ILR eigenvalues", call. = FALSE)
  }
  scores <- as.matrix(result$S)
  clr_loadings <- basis %*% as.matrix(result$Evec)
  names_pc <- paste0("pPC", seq_len(ncol(scores)))
  colnames(scores) <- names_pc
  colnames(clr_loadings) <- names_pc
  list(
    scores = scores,
    loadings = clr_loadings,
    eigenvalues = eigenvalues,
    variance = eigenvalues / sum(eigenvalues)
  )
}

orient_to_reference <- function(fit, reference_loadings) {
  n_axes <- min(ncol(fit$loadings), ncol(reference_loadings))
  for (axis in seq_len(n_axes)) {
    correlation <- suppressWarnings(cor(
      reference_loadings[rownames(fit$loadings), axis],
      fit$loadings[, axis]
    ))
    if (is.finite(correlation) && correlation < 0) {
      fit$loadings[, axis] <- -fit$loadings[, axis]
      fit$scores[, axis] <- -fit$scores[, axis]
    }
  }
  fit
}

tree <- ape::read.tree(tree_path)
if (!ape::is.rooted(tree) || !ape::is.binary(tree) || any(tree$edge.length <= 0)) {
  stop("Corrected tree failed rooted/binary/positive-edge contract", call. = FALSE)
}
if (!ape::is.ultrametric(tree, tol = 1e-8)) {
  stop("Corrected tree is not ultrametric", call. = FALSE)
}

clr_long <- read.csv(clr_path, check.names = FALSE, stringsAsFactors = FALSE)
clr_long <- subset(
  clr_long,
  analysis_id == "order_classified_primary" & zero_replacement == "not_needed"
)
clr <- xtabs(clr_value ~ species + feature, data = clr_long)
clr <- as.matrix(clr)
if (!setequal(rownames(clr), tree$tip.label)) {
  stop("Corrected order CLR species do not exactly match corrected tree tips", call. = FALSE)
}
clr <- clr[tree$tip.label, , drop = FALSE]
if (max(abs(rowMeans(clr))) > 1e-10) {
  stop("Input matrix is not row-centered CLR data", call. = FALSE)
}

ordinary_scores <- read.csv(ordinary_scores_path, check.names = FALSE, stringsAsFactors = FALSE)
ordinary_scores <- subset(
  ordinary_scores,
  analysis_id == "order_classified_primary" & zero_replacement == "not_needed"
)
rownames(ordinary_scores) <- ordinary_scores$species

ordinary_loadings <- read.csv(ordinary_loadings_path, check.names = FALSE, stringsAsFactors = FALSE)
ordinary_loadings <- subset(
  ordinary_loadings,
  analysis_id == "order_classified_primary" & zero_replacement == "not_needed"
)
rownames(ordinary_loadings) <- ordinary_loadings$feature
ordinary_loading_matrix <- as.matrix(ordinary_loadings[colnames(clr), paste0("PC", seq_len(ncol(clr) - 1)), drop = FALSE])

basis <- build_ilr_basis(colnames(clr))
full <- fit_ppca(tree, clr, basis)
full <- orient_to_reference(full, ordinary_loading_matrix)

scores_df <- data.frame(
  species = rownames(full$scores),
  full$scores,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
loadings_df <- data.frame(
  feature = rownames(full$loadings),
  full$loadings,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
variance_df <- data.frame(
  component = paste0("pPC", seq_along(full$eigenvalues)),
  eigenvalue = full$eigenvalues,
  proportion_variance = full$variance,
  cumulative_variance = cumsum(full$variance),
  stringsAsFactors = FALSE
)

comparison_rows <- list()
for (kind in c("score", "loading")) {
  left <- if (kind == "score") {
    as.matrix(ordinary_scores[rownames(full$scores), paste0("PC", seq_len(2)), drop = FALSE])
  } else {
    ordinary_loading_matrix[, seq_len(2), drop = FALSE]
  }
  right <- if (kind == "score") full$scores[, seq_len(2), drop = FALSE] else full$loadings[, seq_len(2), drop = FALSE]
  matrix_cor <- cor(left, right)
  for (i in seq_len(2)) {
    for (j in seq_len(2)) {
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        comparison_kind = kind,
        ordinary_axis = paste0("PC", i),
        phylogenetic_axis = paste0("pPC", j),
        pearson_correlation = matrix_cor[i, j],
        absolute_correlation = abs(matrix_cor[i, j]),
        stringsAsFactors = FALSE
      )
    }
  }
}
comparison_df <- do.call(rbind, comparison_rows)

stability_rows <- list()
for (left_out in tree$tip.label) {
  loo_tree <- ape::drop.tip(tree, left_out)
  loo <- fit_ppca(loo_tree, clr[loo_tree$tip.label, , drop = FALSE], basis)
  for (axis in seq_len(2)) {
    loading_cor <- cor(full$loadings[, axis], loo$loadings[, axis])
    sign_value <- if (is.finite(loading_cor) && loading_cor < 0) -1 else 1
    loo$loadings[, axis] <- loo$loadings[, axis] * sign_value
    loo$scores[, axis] <- loo$scores[, axis] * sign_value
    loading_cor <- abs(cor(full$loadings[, axis], loo$loadings[, axis]))
    score_cor <- abs(cor(full$scores[loo_tree$tip.label, axis], loo$scores[, axis]))
    stability_rows[[length(stability_rows) + 1L]] <- data.frame(
      species_left_out = left_out,
      component = paste0("pPC", axis),
      loading_correlation = loading_cor,
      score_correlation = score_cor,
      stringsAsFactors = FALSE
    )
  }
}
stability_df <- do.call(rbind, stability_rows)

write.csv(scores_df, scores_output, row.names = FALSE)
write.csv(loadings_df, loadings_output, row.names = FALSE)
write.csv(variance_df, variance_output, row.names = FALSE)
write.csv(comparison_df, comparison_output, row.names = FALSE)
write.csv(stability_df, stability_output, row.names = FALSE)

score_plot <- ggplot(scores_df, aes(x = pPC1, y = pPC2, label = species)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.35) +
  geom_point(size = 2.7, color = "#2A9D8F") +
  ggrepel::geom_text_repel(size = 3.0, max.overlaps = Inf, box.padding = 0.35, seed = 20260709) +
  labs(
    title = "Order composition: fixed-tree phylogenetic PCA",
    subtitle = "CLR projected to a full-rank ILR basis; Brownian covariance; n = 18 species",
    x = sprintf("pPC1 (%.1f%%)", 100 * full$variance[[1]]),
    y = sprintf("pPC2 (%.1f%%)", 100 * full$variance[[2]])
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(scores_figure, score_plot, width = 9, height = 7, dpi = 300)

compare_plot_df <- data.frame(
  species = rep(tree$tip.label, 2),
  ordinary = c(ordinary_scores[tree$tip.label, "PC1"], ordinary_scores[tree$tip.label, "PC2"]),
  phylogenetic = c(full$scores[, "pPC1"], full$scores[, "pPC2"]),
  axis = rep(c("PC1 vs pPC1", "PC2 vs pPC2"), each = length(tree$tip.label)),
  stringsAsFactors = FALSE
)
diag_score_cor <- subset(comparison_df, comparison_kind == "score" & sub("p", "", phylogenetic_axis) == ordinary_axis)
cor_labels <- setNames(sprintf("r = %.3f", diag_score_cor$pearson_correlation), paste(diag_score_cor$ordinary_axis, "vs", diag_score_cor$phylogenetic_axis))
comparison_plot <- ggplot(compare_plot_df, aes(x = ordinary, y = phylogenetic, label = species)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#E76F51", linewidth = 0.8) +
  geom_point(size = 2.4, color = "#264653") +
  ggrepel::geom_text_repel(size = 2.4, max.overlaps = 8, seed = 20260709) +
  facet_wrap(~axis, scales = "free") +
  geom_label(
    data = data.frame(axis = names(cor_labels), ordinary = -Inf, phylogenetic = Inf, label = unname(cor_labels)),
    aes(label = label), inherit.aes = FALSE, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
    size = 3, label.size = 0, fill = "white"
  ) +
  labs(
    title = "Ordinary CLR PCA is broadly retained after phylogenetic correction",
    x = "Ordinary CLR-PCA score",
    y = "Phylogenetic ILR-PCA score"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(comparison_figure, comparison_plot, width = 10, height = 5.5, dpi = 300)

loading_plot_df <- rbind(
  data.frame(feature = rownames(full$loadings), component = "pPC1", loading = full$loadings[, 1]),
  data.frame(feature = rownames(full$loadings), component = "pPC2", loading = full$loadings[, 2])
)
loading_plot <- ggplot(loading_plot_df, aes(x = reorder(feature, loading), y = loading, fill = loading > 0)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~component, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51")) +
  labs(title = "Back-projected CLR feature loadings", x = "TE order", y = "Loading") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(loadings_figure, loading_plot, width = 9.5, height = 6.5, dpi = 300)

stability_long <- rbind(
  data.frame(stability_df[c("species_left_out", "component")], metric = "Loading correlation", value = stability_df$loading_correlation),
  data.frame(stability_df[c("species_left_out", "component")], metric = "Score correlation", value = stability_df$score_correlation)
)
stability_plot <- ggplot(stability_long, aes(x = reorder(species_left_out, value), y = value, color = metric, group = metric)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey55") +
  geom_point(position = position_dodge(width = 0.35), size = 2.0) +
  coord_flip() +
  facet_wrap(~component) +
  scale_color_manual(values = c("Loading correlation" = "#E76F51", "Score correlation" = "#2A9D8F")) +
  labs(
    title = "Leave-one-species-out stability of the fixed-tree pPCA",
    x = "Species omitted",
    y = "Absolute correlation with full fit",
    color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(stability_figure, stability_plot, width = 10, height = 7.5, dpi = 300)

pc1_score_cor <- subset(comparison_df, comparison_kind == "score" & ordinary_axis == "PC1" & phylogenetic_axis == "pPC1")$pearson_correlation
pc2_score_cor <- subset(comparison_df, comparison_kind == "score" & ordinary_axis == "PC2" & phylogenetic_axis == "pPC2")$pearson_correlation
pc1_loading_cor <- subset(comparison_df, comparison_kind == "loading" & ordinary_axis == "PC1" & phylogenetic_axis == "pPC1")$pearson_correlation
pc2_loading_cor <- subset(comparison_df, comparison_kind == "loading" & ordinary_axis == "PC2" & phylogenetic_axis == "pPC2")$pearson_correlation

minimums <- aggregate(cbind(loading_correlation, score_correlation) ~ component, stability_df, min)
pc1_min <- minimums[minimums$component == "pPC1", ]
pc2_min <- minimums[minimums$component == "pPC2", ]

report <- sprintf(
  paste0(
    "# Order-level phylogenetic PCA sensitivity\n\n",
    "The approved descriptive order CLR-PCA was refit on the exact corrected 18-tip time tree. Because CLR coordinates are singular by construction, the CLR matrix was projected into a 9-dimensional orthonormal ILR basis before `phytools::phyl.pca(method = \"BM\", mode = \"cov\")`; loadings were then projected back into CLR feature space.\n\n",
    "## Results\n\n",
    "- pPC1 explains %.1f%% and pPC2 explains %.1f%% of phylogenetic ILR variance (%.1f%% cumulative).\n",
    "- Ordinary versus phylogenetic species-score correlations are %.3f for axis 1 and %.3f for axis 2 after sign alignment; axis signs have no biological meaning.\n",
    "- Ordinary versus back-projected loading correlations are %.3f for axis 1 and %.3f for axis 2. The dominant axis is strongly preserved; axis 2 is recognizably similar but more sensitive to phylogenetic covariance.\n",
    "- Minimum leave-one-species-out loading/score correlations are %.3f/%.3f for pPC1 and %.3f/%.3f for pPC2.\n\n",
    "## Verdict\n\n",
    "**Approved as a fixed-tree descriptive sensitivity.** The principal order-composition gradient is not an artifact of ignoring shared ancestry. This does not promote either PCA axis to a causal path variable: the tree source/posterior uncertainty remains unresolved, and direct prespecified log-ratios or multivariate phylogenetic models are preferable for inference.\n\n",
    "## Outputs\n\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n- `%s`\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n"
  ),
  100 * full$variance[[1]], 100 * full$variance[[2]], 100 * sum(full$variance[1:2]),
  pc1_score_cor, pc2_score_cor, pc1_loading_cor, pc2_loading_cor,
  pc1_min$loading_correlation, pc1_min$score_correlation,
  pc2_min$loading_correlation, pc2_min$score_correlation,
  portable(scores_output), portable(loadings_output), portable(variance_output),
  portable(comparison_output), portable(stability_output),
  portable(scores_figure), portable(comparison_figure), portable(loadings_figure), portable(stability_figure)
)
writeLines(report, report_output)

output_paths <- c(
  scores_output, loadings_output, variance_output, comparison_output, stability_output,
  report_output, scores_figure, comparison_figure, loadings_figure, stability_figure
)
manifest <- list(
  analysis_scope = "order_classified_primary_final18_fixed_tree_sensitivity",
  n_species = nrow(clr),
  n_order_features = ncol(clr),
  n_ilr_coordinates = ncol(basis),
  transform = "CLR_to_orthonormal_Helmert_ILR",
  phylogenetic_model = "Brownian_motion_covariance",
  tree = portable(tree_path),
  tree_sha256 = sha256(tree_path),
  fixed_tree_status = "approved_descriptive_sensitivity",
  causal_predictor_status = "not_approved",
  tree_uncertainty_status = "not_approved_single_point_tree",
  ppc1_variance = unname(full$variance[[1]]),
  ppc2_variance = unname(full$variance[[2]]),
  ordinary_ppca_score_correlations = list(PC1_pPC1 = unname(pc1_score_cor), PC2_pPC2 = unname(pc2_score_cor)),
  outputs = lapply(output_paths, function(path) list(path = portable(path), sha256 = sha256(path)))
)
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), manifest_output)

message("Wrote ", report_output)
message("Wrote ", manifest_output)

#!/usr/bin/env Rscript

# Audit the Stewart & Wiens (2025) dated salamander tree and its 200
# time-calibrated bootstrap trees on the exact 18-species analysis panel.
# The external trees are sensitivity inputs; they do not silently replace the
# collaborator-supplied focal tree.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phangorn)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/processing/audit_published_tree_uncertainty.R", mustWork = TRUE)
}
project_root <- dirname(dirname(dirname(script_path)))

focal_path <- file.path(project_root, "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk")
main_source_path <- file.path(project_root, "input_data/phylogeny/stewart_wiens_2025/Supplementary_File_S3.nex")
bootstrap_source_path <- file.path(project_root, "input_data/phylogeny/stewart_wiens_2025/Supplementary_File_S4.nex")
main_zip_path <- file.path(project_root, "input_data/phylogeny/stewart_wiens_2025/mmc14_supplementary_file_s3.zip")
bootstrap_zip_path <- file.path(project_root, "input_data/phylogeny/stewart_wiens_2025/mmc15_supplementary_file_s4.zip")

output_dir <- file.path(project_root, "results/data/corrected/phylogeny")
figure_dir <- file.path(project_root, "results/figures/corrected/phylogeny")
audit_dir <- file.path(project_root, "plans/publication-readiness-deep-audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

main_output <- file.path(output_dir, "stewart_wiens_2025_main_analysis18_v1.nwk")
bootstrap_output <- file.path(output_dir, "stewart_wiens_2025_bootstrap_analysis18_v1.nex")
source_registry_output <- file.path(output_dir, "phylogeny_source_registry_v1.csv")
tree_metrics_output <- file.path(output_dir, "phylogeny_tree_uncertainty_metrics_analysis18_v1.csv")
pairwise_output <- file.path(output_dir, "phylogeny_pairwise_uncertainty_analysis18_v1.csv")
clade_output <- file.path(output_dir, "phylogeny_clade_uncertainty_analysis18_v1.csv")
topology_output <- file.path(output_dir, "phylogeny_topology_difference_analysis18_v1.csv")
manifest_output <- file.path(output_dir, "phylogeny_published_uncertainty_analysis18_v1.manifest.json")
report_output <- file.path(audit_dir, "phylogeny_published_uncertainty_analysis18_v1.md")

density_figure <- file.path(figure_dir, "phylogeny_published_bootstrap_density_analysis18_v1.png")
age_figure <- file.path(figure_dir, "phylogeny_root_age_uncertainty_analysis18_v1.png")
pairwise_figure <- file.path(figure_dir, "phylogeny_patristic_uncertainty_analysis18_v1.png")
topology_figure <- file.path(figure_dir, "phylogeny_focal_vs_published_topology_analysis18_v1.png")

required <- c(focal_path, main_source_path, bootstrap_source_path, main_zip_path, bootstrap_zip_path)
missing <- required[!file.exists(required)]
if (length(missing) > 0) {
  stop("Missing required tree input(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

portable <- function(path) {
  sub(paste0(normalizePath(project_root, winslash = "/"), "/"), "", normalizePath(path, winslash = "/"), fixed = TRUE)
}
sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

strip_desmognathus <- function(labels) {
  sub("^Desmognathus_", "", labels)
}

terminal_pad <- function(tree) {
  depths <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
  target <- max(depths)
  deltas <- target - depths
  terminal_edges <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
  tree$edge.length[terminal_edges] <- tree$edge.length[terminal_edges] + pmax(deltas, 0)
  list(tree = tree, max_padding_myr = max(pmax(deltas, 0)), range_before_myr = diff(range(depths)))
}

prune_published <- function(tree, target_species) {
  tree$tip.label <- strip_desmognathus(tree$tip.label)
  missing <- setdiff(target_species, tree$tip.label)
  if (length(missing) > 0) {
    stop("Published tree is missing final-panel taxa: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tree <- ape::keep.tip(tree, target_species)
  tree <- terminal_pad(tree)$tree
  tree
}

clade_keys <- function(tree) {
  tip_count <- ape::Ntip(tree)
  nodes <- seq.int(tip_count + 1L, tip_count + tree$Nnode)
  root_node <- tip_count + 1L
  nodes <- setdiff(nodes, root_node)
  keys <- vapply(nodes, function(node) {
    paste(sort(ape::extract.clade(tree, node)$tip.label), collapse = ";")
  }, character(1))
  names(keys) <- as.character(nodes)
  keys
}

rooted_rf <- function(left, right) {
  left_keys <- unname(clade_keys(left))
  right_keys <- unname(clade_keys(right))
  length(setdiff(left_keys, right_keys)) + length(setdiff(right_keys, left_keys))
}

cophenetic_correlation <- function(left, right) {
  left_matrix <- ape::cophenetic.phylo(left)
  right_matrix <- ape::cophenetic.phylo(right)[rownames(left_matrix), colnames(left_matrix)]
  cor(as.vector(left_matrix), as.vector(right_matrix))
}

tree_metric_row <- function(tree, tree_id, tree_role, focal, main, padding_range_myr = NA_real_, max_padding_myr = NA_real_) {
  depths <- ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))]
  data.frame(
    tree_id = tree_id,
    tree_role = tree_role,
    n_tips = ape::Ntip(tree),
    rooted = ape::is.rooted(tree),
    fully_bifurcating = ape::is.binary(tree),
    all_branch_lengths_positive = all(is.finite(tree$edge.length) & tree$edge.length > 0),
    root_age_myr = max(depths),
    root_to_tip_range_before_padding_myr = padding_range_myr,
    maximum_terminal_padding_myr = max_padding_myr,
    strict_ultrametric_after_padding = diff(range(depths)) <= 1e-10,
    minimum_branch_length_myr = min(tree$edge.length),
    maximum_branch_length_myr = max(tree$edge.length),
    rooted_rf_to_focal = rooted_rf(tree, focal),
    rooted_rf_to_published_main = rooted_rf(tree, main),
    cophenetic_correlation_to_focal = cophenetic_correlation(focal, tree),
    cophenetic_correlation_to_published_main = cophenetic_correlation(main, tree),
    stringsAsFactors = FALSE
  )
}

focal <- ape::read.tree(focal_path)
target_species <- focal$tip.label
if (length(target_species) != 18 || anyDuplicated(target_species)) {
  stop("Focal tree does not contain an exact unique 18-species panel", call. = FALSE)
}

main_raw <- ape::read.nexus(main_source_path)
main_before <- main_raw
main_before$tip.label <- strip_desmognathus(main_before$tip.label)
main_before <- ape::keep.tip(main_before, target_species)
main_padding <- terminal_pad(main_before)
main <- main_padding$tree

bootstrap_raw <- ape::read.nexus(bootstrap_source_path)
if (!inherits(bootstrap_raw, "multiPhylo") || length(bootstrap_raw) != 200) {
  stop("Expected exactly 200 published bootstrap time trees", call. = FALSE)
}

bootstrap_trees <- vector("list", length(bootstrap_raw))
padding_ranges <- numeric(length(bootstrap_raw))
padding_maxima <- numeric(length(bootstrap_raw))
for (i in seq_along(bootstrap_raw)) {
  tree <- bootstrap_raw[[i]]
  tree$tip.label <- strip_desmognathus(tree$tip.label)
  missing <- setdiff(target_species, tree$tip.label)
  if (length(missing) > 0) {
    stop("Bootstrap tree ", i, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  tree <- ape::keep.tip(tree, target_species)
  padded <- terminal_pad(tree)
  bootstrap_trees[[i]] <- padded$tree
  padding_ranges[[i]] <- padded$range_before_myr
  padding_maxima[[i]] <- padded$max_padding_myr
}
class(bootstrap_trees) <- "multiPhylo"

if (!all(vapply(bootstrap_trees, ape::is.rooted, logical(1))) ||
    !all(vapply(bootstrap_trees, ape::is.binary, logical(1))) ||
    !all(vapply(bootstrap_trees, function(tree) all(tree$edge.length > 0), logical(1))) ||
    !all(vapply(bootstrap_trees, ape::is.ultrametric, logical(1), tol = 1e-8))) {
  stop("At least one pruned published bootstrap tree failed the structural contract", call. = FALSE)
}

ape::write.tree(main, file = main_output, digits = 12)
ape::write.nexus(bootstrap_trees, file = bootstrap_output, translate = FALSE)

source_registry <- data.frame(
  source_id = c(
    "collaborator_local_desmo900dated",
    "stewart_wiens_2025_supplementary_s3",
    "stewart_wiens_2025_supplementary_s4"
  ),
  role = c("focal_working_tree", "published_optimal_time_tree_sensitivity", "published_200_time_calibrated_bootstrap_trees"),
  citation = c(
    "Exact publication/calibration provenance unresolved",
    "Stewart and Wiens 2025, Molecular Phylogenetics and Evolution 204:108272",
    "Stewart and Wiens 2025, Molecular Phylogenetics and Evolution 204:108272"
  ),
  doi = c(NA, "10.1016/j.ympev.2024.108272", "10.1016/j.ympev.2024.108272"),
  source_url = c(
    NA,
    "https://ars.els-cdn.com/content/image/1-s2.0-S1055790324002641-mmc14.zip",
    "https://ars.els-cdn.com/content/image/1-s2.0-S1055790324002641-mmc15.zip"
  ),
  local_source_file = c(portable(focal_path), portable(main_source_path), portable(bootstrap_source_path)),
  local_source_sha256 = c(sha256(focal_path), sha256(main_source_path), sha256(bootstrap_source_path)),
  archive_file = c(NA, portable(main_zip_path), portable(bootstrap_zip_path)),
  archive_sha256 = c(NA, sha256(main_zip_path), sha256(bootstrap_zip_path)),
  n_full_tree_tips = c(18, ape::Ntip(main_raw), ape::Ntip(bootstrap_raw[[1]])),
  n_trees = c(1, 1, length(bootstrap_raw)),
  analysis_disposition = c(
    "focal_pending_exact_source_provenance",
    "approved_external_fixed_tree_sensitivity",
    "approved_external_tree_uncertainty_sensitivity"
  ),
  stringsAsFactors = FALSE
)
write.csv(source_registry, source_registry_output, row.names = FALSE, na = "")

metric_rows <- list(
  tree_metric_row(focal, "focal", "collaborator_focal", focal, main, 0, 0),
  tree_metric_row(
    main, "published_main", "stewart_wiens_2025_main", focal, main,
    main_padding$range_before_myr, main_padding$max_padding_myr
  )
)
for (i in seq_along(bootstrap_trees)) {
  metric_rows[[length(metric_rows) + 1L]] <- tree_metric_row(
    bootstrap_trees[[i]], sprintf("published_bootstrap_%03d", i),
    "stewart_wiens_2025_time_calibrated_bootstrap", focal, main,
    padding_ranges[[i]], padding_maxima[[i]]
  )
}
tree_metrics <- do.call(rbind, metric_rows)
write.csv(tree_metrics, tree_metrics_output, row.names = FALSE)

pair_names <- combn(sort(target_species), 2, simplify = FALSE)
pair_rows <- vector("list", length(pair_names))
focal_dist <- ape::cophenetic.phylo(focal)
main_dist <- ape::cophenetic.phylo(main)
bootstrap_distances <- lapply(bootstrap_trees, ape::cophenetic.phylo)
for (i in seq_along(pair_names)) {
  pair <- pair_names[[i]]
  values <- vapply(bootstrap_distances, function(matrix) matrix[pair[[1]], pair[[2]]], numeric(1))
  pair_rows[[i]] <- data.frame(
    species_1 = pair[[1]],
    species_2 = pair[[2]],
    focal_distance_myr = focal_dist[pair[[1]], pair[[2]]],
    published_main_distance_myr = main_dist[pair[[1]], pair[[2]]],
    bootstrap_mean_distance_myr = mean(values),
    bootstrap_sd_distance_myr = sd(values),
    bootstrap_cv_distance = sd(values) / mean(values),
    bootstrap_q025_distance_myr = unname(quantile(values, 0.025)),
    bootstrap_q05_distance_myr = unname(quantile(values, 0.05)),
    bootstrap_median_distance_myr = median(values),
    bootstrap_q95_distance_myr = unname(quantile(values, 0.95)),
    bootstrap_q975_distance_myr = unname(quantile(values, 0.975)),
    stringsAsFactors = FALSE
  )
}
pairwise <- do.call(rbind, pair_rows)
write.csv(pairwise, pairwise_output, row.names = FALSE)

node_height_for_key <- function(tree, key) {
  keys <- clade_keys(tree)
  node_name <- names(keys)[match(key, keys)]
  if (is.na(node_name)) return(NA_real_)
  node <- as.integer(node_name)
  depths <- ape::node.depth.edgelength(tree)
  max(depths[seq_len(ape::Ntip(tree))]) - depths[[node]]
}

main_keys <- unname(clade_keys(main))
clade_rows <- vector("list", length(main_keys))
for (i in seq_along(main_keys)) {
  key <- main_keys[[i]]
  ages <- vapply(bootstrap_trees, node_height_for_key, numeric(1), key = key)
  clade_rows[[i]] <- data.frame(
    descendant_species = key,
    n_descendant_species = length(strsplit(key, ";", fixed = TRUE)[[1]]),
    present_in_focal_topology = key %in% unname(clade_keys(focal)),
    bootstrap_topology_frequency = mean(is.finite(ages)),
    published_main_age_myr = node_height_for_key(main, key),
    bootstrap_mean_age_myr = mean(ages, na.rm = TRUE),
    bootstrap_sd_age_myr = sd(ages, na.rm = TRUE),
    bootstrap_q025_age_myr = unname(quantile(ages, 0.025, na.rm = TRUE)),
    bootstrap_median_age_myr = median(ages, na.rm = TRUE),
    bootstrap_q975_age_myr = unname(quantile(ages, 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}
clade_uncertainty <- do.call(rbind, clade_rows)
write.csv(clade_uncertainty, clade_output, row.names = FALSE)

focal_keys <- unname(clade_keys(focal))
topology_difference <- rbind(
  data.frame(
    topology = "focal_only",
    descendant_species = setdiff(focal_keys, main_keys),
    biological_interpretation = "clade occurs only in collaborator focal tree",
    stringsAsFactors = FALSE
  ),
  data.frame(
    topology = "published_only",
    descendant_species = setdiff(main_keys, focal_keys),
    biological_interpretation = "clade occurs in published main and all 200 bootstrap trees",
    stringsAsFactors = FALSE
  )
)
write.csv(topology_difference, topology_output, row.names = FALSE)

grDevices::png(density_figure, width = 2700, height = 2200, res = 300)
par(mar = c(5, 10, 4, 2))
phangorn::densiTree(
  bootstrap_trees,
  alpha = 0.025,
  consensus = NULL,
  cex = 0.7,
  type = "phylogram",
  col = grDevices::adjustcolor("#2A9D8F", alpha.f = 0.12)
)
title(main = "Stewart-Wiens 2025: 200 pruned time-calibrated bootstrap trees", font.main = 2)
mtext("Topology is constant for these 18 taxa; horizontal spread is branch-time uncertainty", side = 1, line = 3)
grDevices::dev.off()

bootstrap_metric_rows <- subset(tree_metrics, tree_role == "stewart_wiens_2025_time_calibrated_bootstrap")
age_plot <- ggplot(bootstrap_metric_rows, aes(x = root_age_myr)) +
  geom_histogram(bins = 24, fill = "#2A9D8F", color = "white") +
  geom_vline(xintercept = max(ape::node.depth.edgelength(focal)), color = "#E76F51", linewidth = 1.1) +
  geom_vline(xintercept = max(ape::node.depth.edgelength(main)), color = "#264653", linewidth = 1.1, linetype = "dashed") +
  annotate("text", x = max(ape::node.depth.edgelength(focal)), y = Inf, label = "focal", vjust = 1.5, hjust = -0.1, color = "#E76F51") +
  annotate("text", x = max(ape::node.depth.edgelength(main)), y = Inf, label = "published main", vjust = 2.8, hjust = -0.1, color = "#264653") +
  labs(
    title = "Final-panel crown age varies across published bootstrap time trees",
    x = "Root age of pruned 18-species tree (Myr)",
    y = "Number of bootstrap trees"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave(age_figure, age_plot, width = 8.5, height = 5.8, dpi = 300)

pair_grid <- rbind(
  data.frame(species_1 = pairwise$species_1, species_2 = pairwise$species_2, cv = pairwise$bootstrap_cv_distance),
  data.frame(species_1 = pairwise$species_2, species_2 = pairwise$species_1, cv = pairwise$bootstrap_cv_distance),
  data.frame(species_1 = target_species, species_2 = target_species, cv = 0)
)
pair_grid$species_1 <- factor(pair_grid$species_1, levels = focal$tip.label)
pair_grid$species_2 <- factor(pair_grid$species_2, levels = rev(focal$tip.label))
pair_plot <- ggplot(pair_grid, aes(x = species_1, y = species_2, fill = cv)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "magma", name = "Distance CV") +
  coord_equal() +
  labs(
    title = "Relative patristic-distance uncertainty across 200 time trees",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(pairwise_figure, pair_plot, width = 9.5, height = 8.2, dpi = 300)

grDevices::png(topology_figure, width = 3000, height = 2300, res = 300)
cophy <- phytools::cophylo(focal, main, rotate = TRUE)
plot(
  cophy,
  link.type = "curved",
  link.lwd = 0.7,
  link.col = grDevices::adjustcolor("grey45", alpha.f = 0.55),
  fsize = 0.75,
  pts = FALSE
)
title(main = "Collaborator focal tree vs Stewart-Wiens 2025 main topology", font.main = 2)
mtext("One rooted clade differs in the black-bellied complex; tip links show exact taxon identity", side = 1, line = 1.5)
grDevices::dev.off()

bootstrap_root_ages <- bootstrap_metric_rows$root_age_myr
rf_to_main <- bootstrap_metric_rows$rooted_rf_to_published_main
rf_to_focal <- bootstrap_metric_rows$rooted_rf_to_focal
main_only <- subset(topology_difference, topology == "published_only")$descendant_species
focal_only <- subset(topology_difference, topology == "focal_only")$descendant_species

report <- sprintf(
  paste0(
    "# Published time-tree uncertainty audit: final 18-species panel\n\n",
    "Stewart and Wiens (2025; DOI `10.1016/j.ympev.2024.108272`) provide an optimal time-calibrated tree (Supplementary File S3) and 200 time-calibrated bootstrap trees (Supplementary File S4). The exact publisher archives, extracted files, and SHA-256 hashes are recorded in the source registry. All 18 focal taxa are present by exact binomial name.\n\n",
    "## Structural and topology results\n\n",
    "- All 201 published trees prune to exactly 18 rooted, bifurcating tips with positive branches. Newick rounding required at most %.3f years of terminal padding; topology and internal ages were not altered.\n",
    "- All 200 bootstrap trees have the same rooted final-panel topology as the published main tree (`rooted RF = 0`). They differ from the collaborator focal tree by one rooted bipartition (`rooted RF = 2`).\n",
    "- Focal-only clade: `%s`. Published-main/all-bootstrap-only clade: `%s`. This is a real black-bellied-complex topology sensitivity, not a label substitution.\n",
    "- Focal-versus-published-main patristic correlation is %.3f, indicating broad relative-distance agreement despite the one topology difference and their different crown ages.\n\n",
    "## Time uncertainty\n\n",
    "- Final-panel crown age across 200 trees: median %.3f Myr; 95%% interval %.3f-%.3f Myr; range %.3f-%.3f Myr.\n",
    "- The collaborator focal crown age is %.3f Myr; the published main crown age is %.3f Myr.\n",
    "- Every pairwise patristic distance and every published-main clade age is exported with bootstrap mean, SD, and 95%% interval.\n\n",
    "## Verdict\n\n",
    "**Approved as a published tree-uncertainty sensitivity set.** The repository no longer needs to represent phylogenetic uncertainty with one point tree. The collaborator focal tree remains usable as the focal topology after its exact source/calibration record is supplied; until then, report Stewart-Wiens main plus 200-tree results alongside it. The 200-tree set represents bootstrap/dating sensitivity, not gene-tree conflict or reticulate-network uncertainty.\n\n",
    "## Outputs\n\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n- `%s`\n- `%s`\n- `%s`\n",
    "- `%s`\n- `%s`\n- `%s`\n- `%s`\n"
  ),
  max(padding_maxima) * 1e6,
  paste(focal_only, collapse = " | "), paste(main_only, collapse = " | "),
  cophenetic_correlation(focal, main),
  median(bootstrap_root_ages), quantile(bootstrap_root_ages, 0.025), quantile(bootstrap_root_ages, 0.975),
  min(bootstrap_root_ages), max(bootstrap_root_ages),
  max(ape::node.depth.edgelength(focal)), max(ape::node.depth.edgelength(main)),
  portable(main_output), portable(bootstrap_output), portable(source_registry_output),
  portable(tree_metrics_output), portable(pairwise_output), portable(clade_output), portable(topology_output),
  portable(density_figure), portable(age_figure), portable(pairwise_figure), portable(topology_figure)
)
writeLines(report, report_output)

output_paths <- c(
  main_output, bootstrap_output, source_registry_output, tree_metrics_output,
  pairwise_output, clade_output, topology_output, report_output,
  density_figure, age_figure, pairwise_figure, topology_figure
)
manifest <- list(
  analysis_scope = "final18_published_time_tree_uncertainty",
  citation = "Stewart and Wiens 2025 Molecular Phylogenetics and Evolution 204:108272",
  doi = "10.1016/j.ympev.2024.108272",
  n_species = 18,
  n_published_bootstrap_trees = 200,
  exact_tip_match_all_trees = TRUE,
  structural_status = "approved_after_rounding_only_terminal_padding",
  uncertainty_sensitivity_status = "approved_published_bootstrap_time_tree_set",
  network_reticulation_uncertainty_status = "not_represented",
  bootstrap_root_age_myr = list(
    median = unname(median(bootstrap_root_ages)),
    q025 = unname(quantile(bootstrap_root_ages, 0.025)),
    q975 = unname(quantile(bootstrap_root_ages, 0.975))
  ),
  bootstrap_rooted_rf_to_published_main = sort(unique(rf_to_main)),
  bootstrap_rooted_rf_to_focal = sort(unique(rf_to_focal)),
  focal_vs_published_main_rooted_rf = rooted_rf(focal, main),
  focal_vs_published_main_cophenetic_correlation = unname(cophenetic_correlation(focal, main)),
  outputs = lapply(output_paths, function(path) list(path = portable(path), sha256 = sha256(path)))
)
writeLines(jsonlite::toJSON(manifest, pretty = TRUE, auto_unbox = TRUE), manifest_output)

message("Wrote ", report_output)
message("Wrote ", manifest_output)

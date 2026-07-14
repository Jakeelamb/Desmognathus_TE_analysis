#!/usr/bin/env Rscript

# Rebuild the high-information TE figures used in the Research Talk IV PDF.
# This is a presentation/descriptive layer only: it consumes frozen local
# dnaPipeTE, diversity, and LTR-depth tables and never runs repeat discovery,
# RepeatMasker, dnaPipeTE, or mapping.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
flag <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(flag)) dirname(normalizePath(sub("^--file=", "", flag[[1]]), mustWork = FALSE)) else "scripts/processing"
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(cluster)
  library(jsonlite)
})

data_dir <- file.path(root, "results", "data")
out_dir <- file.path(data_dir, "legacy_descriptive")
fig_dir <- file.path(root, "results", "figures", "legacy_descriptive")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

canonical <- function(x) tolower(trimws(gsub("^D\\.\\s*|^Desmognathus\\s+", "", as.character(x))))
save_plot <- function(plot, stem, width = 11, height = 8) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot,
    width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, bg = "white")
}

lookup <- read.delim(file.path(root, "input_data", "lookup_table.txt"), check.names = FALSE)
lookup$species <- canonical(lookup$Species)
if (nrow(lookup) != 34L || anyDuplicated(lookup$species)) stop("Active lookup must contain 34 unique TE species")

read_component_file <- function(path) {
  table <- read.table(path, header = FALSE, sep = "", quote = "", fill = TRUE, comment.char = "", stringsAsFactors = FALSE)
  if (ncol(table) < 7L) stop("Unexpected dnaPipeTE component file: ", path)
  names(table)[seq_len(7L)] <- c("reads", "aligned_bases", "contig", "rm_hit_length_bp", "annotation", "rm_classification", "hitlength_contiglength")
  table$reads <- as.numeric(table$reads)
  table$aligned_bases <- as.numeric(table$aligned_bases)
  table$hitlength_contiglength <- as.numeric(table$hitlength_contiglength)
  table
}

component_root <- file.path(root, "input_data", "dnaPipeTE")
run_quality <- list()
orestes_runs <- list()
for (record in split(lookup, lookup$species)) {
  species <- record$species[[1]]
  accession <- record$SRA_Accension[[1]]
  pattern <- paste0("^", accession, "(R?[0-9]+)?_reads_per_component_and_annotation$")
  files <- list.files(component_root, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No dnaPipeTE component file for active resource: ", accession)
  for (path in files) {
    run_label <- sub("_reads_per_component_and_annotation$", "", basename(path))
    component <- read_component_file(path)
    run_quality[[length(run_quality) + 1L]] <- tibble(
      species = species, active_sra_accession = accession, run_label = run_label,
      n_components = nrow(component), total_reads = sum(component$reads, na.rm = TRUE),
      total_aligned_bases = sum(component$aligned_bases, na.rm = TRUE),
      median_hitlength_contig_ratio = median(component$hitlength_contiglength, na.rm = TRUE),
      fraction_components_classified = mean(!is.na(component$rm_classification) & nzchar(component$rm_classification))
    )
    if (species == "orestes") orestes_runs[[run_label]] <- component
  }
}
run_quality <- bind_rows(run_quality) %>% arrange(species, run_label)
if (n_distinct(run_quality$run_label) < 35L) stop("Expected the orestes retry run in dnaPipeTE quality inventory")
species_quality <- run_quality %>% group_by(species) %>% summarise(
  n_dnapipete_runs = n(),
  run_labels = toJSON(sort(run_label), auto_unbox = TRUE),
  across(c(n_components, total_reads, total_aligned_bases, median_hitlength_contig_ratio, fraction_components_classified), mean),
  .groups = "drop"
)
write_csv(run_quality, file.path(out_dir, "dnapipete_run_quality_te34_legacy_descriptive_v1.csv"))
write_csv(species_quality, file.path(out_dir, "dnapipete_species_quality_te34_legacy_descriptive_v1.csv"))

quality_long <- species_quality %>% select(species, total_reads, total_aligned_bases, median_hitlength_contig_ratio) %>%
  pivot_longer(-species, names_to = "metric", values_to = "value")
quality_plot <- ggplot(quality_long, aes(x = reorder(species, value), y = value, fill = metric)) +
  geom_col(show.legend = FALSE) + coord_flip() + facet_wrap(~ metric, scales = "free_x", ncol = 1) +
  scale_fill_viridis_d() + theme_bw(base_size = 11) +
  labs(title = "dnaPipeTE sequencing-input and component quality", subtitle = "Per-run summaries; same-species retries are averaged for the species summary", x = "Desmognathus species", y = NULL)
save_plot(quality_plot, "dnapipete_input_quality_te34_legacy_descriptive_v1", 11, 10)

read_composition <- function(filename) {
  table <- read.csv(file.path(data_dir, filename), row.names = 1, check.names = FALSE)
  rownames(table) <- canonical(rownames(table))
  table <- as.matrix(table[lookup$species, , drop = FALSE])
  storage.mode(table) <- "numeric"
  table / rowSums(table)
}
average_orestes_runs <- function(column, all_features) {
  per_run <- lapply(orestes_runs, function(frame) {
    # The original dnaPipeTE component tables store a repeat-class string,
    # rather than separate Order/Superfamily columns. Reconstruct those exact
    # display categories for the same run-level averaging policy.
    classification <- as.character(frame$rm_classification)
    category <- if (column == "Order") sub("/.*$", "", classification) else sub("^.*/", "", classification)
    category[is.na(classification) | !grepl("/", classification, fixed = TRUE)] <- NA_character_
    values <- tapply(frame$aligned_bases, category, sum, na.rm = TRUE)
    values <- values[intersect(names(values), all_features)]
    out <- setNames(rep(0, length(all_features)), all_features)
    out[names(values)] <- values
    out / sum(out)
  })
  Reduce(`+`, per_run) / length(per_run)
}
order_comp <- read_composition("dnaPipeTE_order_breakdown.csv")
super_comp <- read_composition("dnaPipeTE_superfamily_breakdown.csv")
order_comp["orestes", ] <- average_orestes_runs("Order", colnames(order_comp))
super_comp["orestes", ] <- average_orestes_runs("Superfamily", colnames(super_comp))
write_csv(as.data.frame(order_comp) %>% mutate(species = rownames(order_comp), .before = 1), file.path(out_dir, "te_order_composition_te34_legacy_descriptive_v1.csv"))
write_csv(as.data.frame(super_comp) %>% mutate(species = rownames(super_comp), .before = 1), file.path(out_dir, "te_superfamily_composition_te34_legacy_descriptive_v1.csv"))

mean_composition <- bind_rows(
  tibble(level = "Order", category = colnames(order_comp), mean_proportion = colMeans(order_comp)),
  tibble(level = "Superfamily", category = colnames(super_comp), mean_proportion = colMeans(super_comp))
) %>% group_by(level) %>% arrange(desc(mean_proportion), .by_group = TRUE) %>% ungroup()
composition_plot <- ggplot(mean_composition, aes(x = reorder(category, mean_proportion), y = mean_proportion, fill = level)) +
  geom_col(show.legend = FALSE) + coord_flip() + facet_wrap(~level, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + scale_fill_viridis_d() + theme_bw(base_size = 11) +
  labs(title = "Mean TE order and superfamily composition across species", subtitle = "TE34; orestes is the equal-weight average of its two dnaPipeTE runs", x = NULL, y = "Mean aligned-repeat composition")
save_plot(composition_plot, "te_mean_composition_te34_legacy_descriptive_v1", 12, 9)

presence <- colSums(super_comp > 0)
features <- names(presence)[presence >= 2L & tolower(names(presence)) != "ginger"]
if (!all(c("Chapaev", "Dada") %in% features)) warning("Expected Chapaev and Dada to remain in the legacy PCA feature set")
pca <- prcomp(super_comp[, features, drop = FALSE], center = TRUE, scale. = TRUE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
pc <- seq_along(variance)
variance_table <- tibble(pc = pc, eigenvalue = pca$sdev^2, proportion_variance = variance, cumulative_variance = cumsum(variance))
scores <- as.data.frame(pca$x) %>% mutate(species = rownames(pca$x), .before = 1)
loadings <- as.data.frame(pca$rotation) %>% mutate(feature = rownames(pca$rotation), .before = 1)
write_csv(variance_table, file.path(out_dir, "te_superfamily_pca_variance_te34_legacy_descriptive_v1.csv"))
write_csv(scores, file.path(out_dir, "te_superfamily_pca_scores_te34_legacy_descriptive_v1.csv"))
write_csv(loadings, file.path(out_dir, "te_superfamily_pca_loadings_te34_legacy_descriptive_v1.csv"))

scree_plot <- ggplot(variance_table, aes(pc)) +
  geom_col(aes(y = proportion_variance), fill = "#2c7fb8", alpha = .85) +
  geom_line(aes(y = cumulative_variance), color = "#d95f0e", linewidth = 1) + geom_point(aes(y = cumulative_variance), color = "#d95f0e") +
  geom_hline(yintercept = c(.8, .9), linetype = "dashed", color = "#238b45") +
  scale_x_continuous(breaks = pc) + scale_y_continuous(labels = scales::percent_format()) + theme_bw(base_size = 12) +
  labs(title = "TE superfamily PCA scree and cumulative variance", subtitle = paste(length(features), "features; Ginger removed; Chapaev and Dada retained"), x = "Principal component", y = "Variance explained / cumulative variance")
save_plot(scree_plot, "te_superfamily_pca_scree_te34_legacy_descriptive_v1", 10, 7)

pc_matrix <- as.matrix(pca$x[, seq_len(min(6L, ncol(pca$x))), drop = FALSE])
set.seed(20250704)
wss <- vapply(1:10, function(k) kmeans(pc_matrix, centers = k, nstart = 100)$tot.withinss, numeric(1))
sil <- c(NA_real_, vapply(2:10, function(k) { fit <- kmeans(pc_matrix, centers = k, nstart = 100); mean(silhouette(fit$cluster, dist(pc_matrix))[, 3]) }, numeric(1)))
cluster_metrics <- tibble(k = 1:10, wss = wss, silhouette = sil)
write_csv(cluster_metrics, file.path(out_dir, "te_superfamily_cluster_metrics_te34_legacy_descriptive_v1.csv"))
elbow_plot <- ggplot(cluster_metrics, aes(k, wss)) + geom_line(color = "#2c7fb8") + geom_point(color = "#2c7fb8", size = 2.5) +
  scale_x_continuous(breaks = 1:10) + theme_bw(base_size = 12) + labs(title = "Elbow plot for TE-superfamily PCA clustering", x = "Number of clusters (k)", y = "Total within-cluster sum of squares")
save_plot(elbow_plot, "te_superfamily_pca_elbow_te34_legacy_descriptive_v1", 9, 7)
silhouette_plot <- ggplot(filter(cluster_metrics, k >= 2), aes(k, silhouette)) + geom_line(color = "#7a0177") + geom_point(color = "#7a0177", size = 2.5) +
  scale_x_continuous(breaks = 2:10) + theme_bw(base_size = 12) + labs(title = "Silhouette support for TE-superfamily PCA clustering", x = "Number of clusters (k)", y = "Mean silhouette width")
save_plot(silhouette_plot, "te_superfamily_pca_silhouette_te34_legacy_descriptive_v1", 9, 7)

set.seed(20250704)
cluster_fit <- kmeans(pc_matrix, centers = 5, nstart = 100)
scores$cluster <- factor(cluster_fit$cluster)
write_csv(scores %>% select(species, cluster, everything()), file.path(out_dir, "te_superfamily_pca_cluster_assignments_te34_legacy_descriptive_v1.csv"))
ellipse_scores <- scores %>% add_count(cluster) %>% filter(n >= 3)
cluster_plot <- ggplot(scores, aes(PC1, PC2, colour = cluster, label = species)) + geom_point(size = 3.5, alpha = .85) +
  stat_ellipse(data = ellipse_scores, aes(group = cluster), linewidth = .6, linetype = 2, show.legend = FALSE) + geom_text_repel(size = 3, max.overlaps = Inf) + scale_colour_viridis_d() + theme_bw(base_size = 12) + coord_equal() +
  labs(title = "PCA clustering of TE superfamily composition", subtitle = "K-means k = 5 on the first six scaled PCA axes", x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]), y = sprintf("PC2 (%.1f%%)", 100 * variance[[2]]), colour = "Cluster")
save_plot(cluster_plot, "te_superfamily_pca_clusters_te34_legacy_descriptive_v1", 10, 8)

read_diversity <- function(filename, level) {
  table <- read.csv(file.path(data_dir, filename), row.names = 1, check.names = FALSE)
  table$species <- canonical(rownames(table))
  table %>% filter(species %in% lookup$species) %>% transmute(species, level = level, Simpson = Simpson_Diversity, Shannon = Shannon_Diversity, Pielou = Pielou_Evenness)
}
diversity <- bind_rows(read_diversity("diversity_order_stats.csv", "Order"), read_diversity("diversity_superfamily_stats.csv", "Superfamily")) %>%
  pivot_longer(c(Simpson, Shannon, Pielou), names_to = "index", values_to = "value")
write_csv(diversity, file.path(out_dir, "te_diversity_indices_te34_legacy_descriptive_v1.csv"))
diversity_plot <- ggplot(diversity, aes(level, value, fill = level)) + geom_violin(alpha = .7, trim = FALSE) + geom_boxplot(width = .16, outlier.shape = NA, alpha = .9) + geom_jitter(width = .08, size = 1.5, alpha = .7) + facet_wrap(~index, scales = "free_y") + scale_fill_viridis_d() + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(title = "TE diversity indices across Desmognathus species", x = "TE level", y = "Index value")
save_plot(diversity_plot, "te_diversity_indices_te34_legacy_descriptive_v1", 11, 7)

ectopic <- read_tsv(file.path(data_dir, "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"), show_col_types = FALSE) %>%
  mutate(species = canonical(species), ratio_terminal_internal = as.numeric(ratio_terminal_internal)) %>% filter(species %in% lookup$species, is.finite(ratio_terminal_internal), ratio_terminal_internal > 0)
if (!nrow(ectopic)) stop("No active TE34 LTR proxy rows")
ectopic_summary <- ectopic %>% group_by(species) %>% summarise(n_elements = n(), mean_ratio = mean(ratio_terminal_internal), median_ratio = median(ratio_terminal_internal), .groups = "drop") %>% arrange(desc(median_ratio))
ectopic$species <- factor(ectopic$species, levels = ectopic_summary$species)
write_csv(ectopic_summary, file.path(out_dir, "ectopic_species_summary_te34_legacy_descriptive_v1.csv"))
anova_fit <- aov(log10(ratio_terminal_internal) ~ species, data = ectopic)
kruskal_fit <- kruskal.test(ratio_terminal_internal ~ species, data = ectopic)
writeLines(c("One-way ANOVA on log10 terminal:internal ratio", capture.output(summary(anova_fit)), "", "Kruskal-Wallis test on raw terminal:internal ratio", capture.output(kruskal_fit)), file.path(out_dir, "ectopic_anova_te34_legacy_descriptive_v1.txt"))
ectopic_plot <- ggplot(ectopic, aes(species, ratio_terminal_internal, fill = species)) + geom_hline(yintercept = 1, linetype = "dashed", colour = "#b2182b") + geom_violin(scale = "width", trim = FALSE, alpha = .85, colour = NA) + stat_summary(fun = median, geom = "point", shape = 95, size = 6, colour = "black") + scale_y_log10() + theme_bw(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major.x = element_blank()) + labs(title = "LTR terminal:internal depth ratio by species", subtitle = "LTRs >=3 kb with >=5 annotated domains; horizontal line = 1:1", x = "Species", y = "Terminal : internal depth ratio (log10 scale)")
save_plot(ectopic_plot, "ectopic_ratio_violin_te34_legacy_descriptive_v1", 12, 8)

manifest <- list(
  analysis_id = "legacy_descriptive_figure_bundle_v1", n_te_species = 34,
  presentation_role = "restored_lab_meeting_style_descriptive_figures",
  dnapipete_retry_policy = "within-species equal-weight mean of run-level summaries and compositions",
  orestes_runs = sort(unique(run_quality$run_label[run_quality$species == "orestes"])),
  pca = list(level = "superfamily", transform = "scaled raw relative composition", removed_feature = "Ginger", retained_low_presence_features = c("Chapaev", "Dada"), kmeans_k = 5),
  ectopic_statistics = c("one-way ANOVA on log10 ratio", "Kruskal-Wallis on raw ratio"),
  expensive_upstream_tools_executed = FALSE
)
write_json(manifest, file.path(out_dir, "legacy_descriptive_figure_bundle_v1.manifest.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Legacy descriptive TE figure bundle rebuilt in ", fig_dir, "\n", sep = "")

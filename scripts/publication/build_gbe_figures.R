#!/usr/bin/env Rscript

# Build the first GBE-styled ggplot figure set from Publication/datasets.

## ---- figure-setup
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

ROOT <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(ROOT, "scripts/publication/gbe_theme.R"))
DATA <- file.path(ROOT, "Publication/datasets")
OUT <- file.path(ROOT, "Publication/figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_supp <- function(filename) {
  read_csv(file.path(DATA, filename), show_col_types = FALSE, progress = FALSE)
}

species_label <- function(x) {
  text <- as.character(x)
  ifelse(grepl("^D\\. ", text), text, paste("D.", text))
}

manifest_rows <- list()
record_figure <- function(id, basename, width_mm, height_mm, datasets, legend, alt_text, status) {
  paths <- save_gbe_figure(
    get(id, envir = .GlobalEnv),
    file.path(OUT, basename),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 300
  )
  manifest_rows[[length(manifest_rows) + 1L]] <<- tibble(
    figure_id = id,
    basename = basename,
    width_mm = width_mm,
    height_mm = height_mm,
    final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT,
    datasets = datasets,
    release_status = status,
    legend = legend,
    alt_text = alt_text,
    pdf = basename |> paste0(".pdf"),
    png = basename |> paste0(".png"),
    tiff = basename |> paste0(".tif")
  )
  invisible(paths)
}

## ---- figure-1
# Figure 1: focal time tree, analysis design, genomic resources, and dnaPipeTE QC.
te_panel <- read_supp("Supplementary_Data_S02_te_resource_panel34.csv")
path_panel <- read_supp("Supplementary_Data_S03_path_analysis_panel24.csv")
path_traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
tree_audit_edges <- read_supp("Supplementary_Data_S38_phylogeny_edges.csv")
tree_audit_nodes <- read_supp("Supplementary_Data_S39_phylogeny_nodes.csv")
tree_path <- file.path(DATA, "trees", "desmognathus_time_tree_path24_v1.nwk")
path_tree <- ladderize(read.tree(tree_path), right = FALSE)

tip_count <- Ntip(path_tree)
node_count <- Nnode(path_tree)
root_distance <- node.depth.edgelength(path_tree)
root_age <- max(root_distance[seq_len(tip_count)])
node_age <- root_age - root_distance
node_y <- node.height(path_tree)

tree_horizontal <- tibble(
  parent = path_tree$edge[, 1],
  child = path_tree$edge[, 2]
) %>%
  mutate(
    x = node_age[parent],
    xend = node_age[child],
    y = node_y[child],
    yend = node_y[child]
  )
tree_vertical <- tibble(parent = sort(unique(path_tree$edge[, 1]))) %>%
  rowwise() %>%
  mutate(
    x = node_age[parent],
    ymin = min(node_y[path_tree$edge[path_tree$edge[, 1] == parent, 2]]),
    ymax = max(node_y[path_tree$edge[path_tree$edge[, 1] == parent, 2]])
  ) %>%
  ungroup()
tree_tips <- tibble(
  node = seq_len(tip_count),
  tree_tip = path_tree$tip.label,
  x = 0,
  y = node_y[seq_len(tip_count)]
) %>%
  left_join(
    path_traits %>% select(species, tree_tip),
    by = "tree_tip",
    relationship = "one-to-one"
  )

audit_tip_labels <- tree_audit_nodes %>%
  filter(is_tip) %>%
  pull(tip_label)
if (
  tip_count != 24L ||
    node_count != 23L ||
    nrow(tree_horizontal) != 46L ||
    nrow(tree_vertical) != 23L ||
    nrow(tree_audit_edges) != 46L ||
    nrow(tree_audit_nodes) != 47L ||
    any(is.na(tree_tips$species)) ||
    !setequal(path_tree$tip.label, path_traits$tree_tip) ||
    !setequal(path_tree$tip.label, audit_tip_labels) ||
    !isTRUE(all.equal(
      tree_horizontal$x - tree_horizontal$xend,
      path_tree$edge.length,
      tolerance = 1e-12,
      check.attributes = FALSE
    )) ||
    abs(root_age - max(tree_audit_nodes$root_distance_myr)) > 1e-4
) {
  stop("Figure 1 requires the exact audited 24-tip Path24 time tree")
}
if (
  !is.rooted(path_tree) ||
    !is.binary(path_tree) ||
    !is.ultrametric(path_tree, tol = 1e-6) ||
    any(path_tree$edge.length < 0)
) {
  stop("Figure 1 Path24 tree must be rooted, binary, time-scaled, and have nonnegative edges")
}

phylogeny_a <- ggplot() +
  geom_segment(
    data = tree_horizontal,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = tree_vertical,
    aes(x = x, xend = x, y = ymin, yend = ymax),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = tree_tips,
    aes(x = -0.6, y = y, label = species),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    fontface = "italic"
  ) +
  scale_x_reverse(
    breaks = c(30, 20, 10, 0),
    limits = c(root_age, -14),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, tip_count + 0.5), expand = c(0, 0)) +
  labs(x = "Time before present (Ma; as supplied)", y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(3, 3, 3, 4, unit = "pt")
  )
assembly <- read_supp("Supplementary_Data_S05_assembly_quality.csv") %>%
  mutate(
    species_label = species_label(species),
    assembly_span_gb = total_sequence_length_bp / 1e9,
    contig_n50_kb = contig_n50_bp / 1e3,
    label = if_else(
      rank(-assembly_span_gb, ties.method = "first") <= 3L |
        rank(-contig_n50_kb, ties.method = "first") == 1L,
      species_label,
      NA_character_
    )
  )
dnapipete_qc <- read_supp("Supplementary_Data_S06_dnapipete_input_quality.csv") %>%
  mutate(species_label = species_label(species))

panel_counts <- tibble(
  panel = factor(
    c("Genomic TE resources", "Finalized phenotype/path panel", "Integrated TE × phenotype panel"),
    levels = rev(c("Genomic TE resources", "Finalized phenotype/path panel", "Integrated TE × phenotype panel"))
  ),
  species = c(
    nrow(te_panel),
    nrow(path_panel),
    nrow(inner_join(te_panel, path_panel, by = "species"))
  )
)
if (!identical(panel_counts$species, c(34L, 24L, 21L))) {
  stop("Figure 1 panel counts must remain TE34, Path24, and their 21-species overlap")
}

design_a <- ggplot(panel_counts, aes(species, panel)) +
  geom_col(width = 0.62, fill = GBE_COLORS[["blue"]]) +
  geom_text(aes(label = species), hjust = -0.35, size = GBE_LABEL_SIZE_MM) +
  scale_x_continuous(limits = c(0, 38), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Species", y = NULL) +
  theme_gbe()

design_b <- ggplot(assembly, aes(assembly_span_gb, contig_n50_kb)) +
  geom_point(shape = 21, size = 2.2, stroke = 0.35, fill = GBE_COLORS[["orange"]], colour = GBE_COLORS[["charcoal"]]) +
  geom_text_repel(data = filter(assembly, !is.na(label)), aes(label = label), size = GBE_LABEL_SIZE_MM, fontface = "italic", min.segment.length = 0, max.overlaps = Inf, seed = 20260713) +
  scale_y_log10(labels = label_number()) +
  labs(x = "Assembly span (Gb)", y = "Contig N50 (kb, log scale)") +
  theme_gbe()

design_c <- ggplot(dnapipete_qc, aes(repeat_aligned_fraction_of_sample, order_unresolved_fraction)) +
  geom_hline(yintercept = median(dnapipete_qc$order_unresolved_fraction, na.rm = TRUE), linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_point(aes(shape = n_dnapipete_runs > 1), size = 2.1, colour = GBE_COLORS[["green"]]) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("One run", "Retry-averaged"), name = NULL) +
  labs(x = "Repeat-aligned fraction of sampled bases", y = "Unresolved order mass") +
  theme_gbe()

figure_1_layout <- "
AAABB
AAACC
AAADD
"
Figure_1 <- phylogeny_a + design_a + design_b + design_c +
  plot_layout(design = figure_1_layout, guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")

record_figure(
  "Figure_1", "Figure_1_study_design_and_genomic_qc", 185, 142,
  "S02;S03;S04;S05;S06;S38;S39",
  "The focal 24-tip Path24 time tree is shown with branch time exactly as supplied by Alex Pyron; its publication or archive identifier, calibration method, and branch-length provenance remain pending. The declared panels contain 34 genomic TE species, 24 finalized phenotype/path species, and an exact 21-species overlap for integrated TE-phenotype summaries. Assembly contiguity and dnaPipeTE classification accounting are shown for TE34; the dashed line is the median unresolved order fraction.",
  "Panel A shows the supplied 24-tip time-calibrated Path24 phylogeny from approximately 30.6 million years before present to the present, with italic species labels. Panel B shows counts for the 34-species genomic panel, 24-species phenotype/path panel, and 21-species overlap. Panels C and D show assembly span versus contig N50 and repeat-aligned read fraction versus unresolved TE-order mass.",
  "main_candidate"
)

## ---- figure-2
# Figure 2: TE diversity and compositional ordination.
te_diversity <- read_supp("Supplementary_Data_S08_te_diversity.csv") %>%
  filter(te_level == "order") %>%
  select(
    species,
    shannon_entropy,
    gini_simpson
  ) %>%
  pivot_longer(
    cols = c(shannon_entropy, gini_simpson),
    names_to = "metric",
    values_to = "diversity_value"
  ) %>%
  mutate(
    metric = factor(
      recode(
        metric,
        shannon_entropy = "Shannon entropy",
        gini_simpson = "Gini-Simpson index"
      ),
      levels = c("Shannon entropy", "Gini-Simpson index")
    )
  )
diversity_counts <- te_diversity %>% count(metric)
if (
  n_distinct(te_diversity$species) != 34L ||
    nrow(diversity_counts) != 2L ||
    !all(diversity_counts$n == 34L)
) {
  stop("Figure 2 classified-only order diversity panel requires 34 species per metric")
}
te_scores <- read_supp("Supplementary_Data_S09_te_pca_scores.csv") %>%
  filter(te_level == "superfamily", composition_mode == "classified_conditional") %>%
  mutate(
    species_label = species_label(species),
    distance = sqrt(PC1^2 + PC2^2),
    label = if_else(rank(-distance, ties.method = "first") <= 10, species_label, NA_character_)
  )
te_variance <- read_supp("Supplementary_Data_S10_te_pca_variance.csv") %>%
  filter(PC <= 10)
pc1_variance <- te_variance %>% filter(PC == 1) %>% pull(variance_explained)
pc2_variance <- te_variance %>% filter(PC == 2) %>% pull(variance_explained)
if (length(pc1_variance) != 1L || length(pc2_variance) != 1L) {
  stop("Figure 2 requires one variance row for each of the first two PCs")
}

te_a <- ggplot(te_diversity, aes(x = "", y = diversity_value)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.35,
    fill = GBE_COLORS[["blue"]],
    alpha = 0.65
  ) +
  geom_point(
    position = position_jitter(width = 0.1, height = 0, seed = 20260713),
    shape = 21,
    size = 0.9,
    stroke = 0.2,
    fill = GBE_COLORS[["blue"]],
    colour = "white",
    alpha = 0.65
  ) +
  facet_wrap(vars(metric), ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Order-level diversity") +
  theme_gbe() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

te_b <- ggplot(te_scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = GBE_COLORS[["light_gray"]]) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = GBE_COLORS[["light_gray"]]) +
  geom_point(size = 1.8, colour = GBE_COLORS[["purple"]]) +
  geom_text_repel(data = filter(te_scores, !is.na(label)), aes(label = label), size = GBE_LABEL_SIZE_MM, fontface = "italic", max.overlaps = Inf, min.segment.length = 0, seed = 20260713) +
  labs(
    x = sprintf("Superfamily CLR-PC1 (%.1f%%)", 100 * pc1_variance),
    y = sprintf("Superfamily CLR-PC2 (%.1f%%)", 100 * pc2_variance)
  ) +
  theme_gbe()

te_c <- ggplot(te_variance, aes(PC, variance_explained)) +
  geom_col(width = 0.68, fill = GBE_COLORS[["green"]]) +
  scale_x_continuous(breaks = te_variance$PC) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Principal component", y = "Variance explained") +
  theme_gbe()

Figure_2 <- te_a + te_b + te_c +
  plot_layout(widths = c(1.05, 1.2, 0.85)) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")

record_figure(
  "Figure_2", "Figure_2_te_diversity_and_ordination", 185, 82,
  "S08;S09;S10;S11",
  "Panel A shows classified-only order-level Shannon entropy and the Gini-Simpson index (1 − Σᵢpᵢ²) across 34 species; classified categories are reclosed to one within each species. Boxes show medians and interquartile ranges, and points are species. Observed richness is not plotted as a diversity index. Unresolved aligned-base mass remains separate in S07 as dnaPipeTE classification-accounting QC and is not a diversity category. The retry-averaged 24-feature CLR-PCA summarizes species differences among superfamilies detected in all 34 species; only the ten most distant PCA scores are labeled.",
  "Panel A shows classified-only order-level Shannon entropy and Gini-Simpson index distributions for 34 species. Unresolved aligned-base mass is retained separately in S07 as dnaPipeTE quality-control accounting and is not included in diversity. Panel B shows species scores on the first two superfamily CLR principal components, and panel C shows variance explained by the first ten components.",
  "main_candidate_with_pca_supplementary"
)

## ---- figure-3
# Figure 3: RepeatMasker divergence landscape.
landscape <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
top_categories <- landscape %>%
  group_by(order) %>%
  summarise(hit_bp = sum(hit_bp), .groups = "drop") %>%
  slice_max(hit_bp, n = 6, with_ties = FALSE) %>%
  pull(order)
category_levels <- c(top_categories, "Other")
landscape_by_species <- landscape %>%
  mutate(repeat_category = if_else(order %in% top_categories, order, "Other")) %>%
  group_by(species, repeat_category, divergence_bin_start_pct) %>%
  summarise(percent_species_hit_bp = sum(percent_species_hit_bp), .groups = "drop") %>%
  right_join(
    expand_grid(
      species = sort(unique(landscape$species)),
      repeat_category = category_levels,
      divergence_bin_start_pct = 0:50
    ),
    by = c("species", "repeat_category", "divergence_bin_start_pct")
  ) %>%
  mutate(percent_species_hit_bp = replace_na(percent_species_hit_bp, 0))
if (
  n_distinct(landscape_by_species$species) != 34L ||
    nrow(landscape_by_species) != 34L * length(category_levels) * 51L
) {
  stop("Figure 3 landscape grid must zero-fill every TE34 category/bin combination")
}
landscape_summary <- landscape_by_species %>%
  group_by(repeat_category, divergence_bin_start_pct) %>%
  summarise(
    mean_percent = mean(percent_species_hit_bp),
    se_percent = sd(percent_species_hit_bp) / sqrt(n()),
    n_species = n(),
    .groups = "drop"
  )
if (!all(landscape_summary$n_species == 34L)) {
  stop("Figure 3 across-species means must give every TE34 species equal weight")
}
young_repeat <- landscape %>%
  group_by(species) %>%
  summarise(young_repeat_percent = sum(percent_species_hit_bp[divergence_bin_start_pct < 5]), .groups = "drop") %>%
  mutate(species_label = species_label(species), species_label = reorder(species_label, young_repeat_percent))

repeat_a <- ggplot(
  landscape_summary,
  aes(
    divergence_bin_start_pct,
    mean_percent,
    colour = repeat_category,
    linetype = repeat_category
  )
) +
  geom_line(linewidth = 0.55) +
  scale_color_gbe(
    name = "Repeat category",
    values = setNames(
      c(unname(GBE_COLORS[c("purple", "sky", "blue", "green", "yellow", "orange")]), GBE_COLORS[["gray"]]),
      category_levels
    )
  ) +
  scale_linetype_manual(
    name = "Repeat category",
    values = setNames(
      c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "42"),
      category_levels
    )
  ) +
  labs(x = "RepeatMasker divergence bin (%)", y = "Mean aligned repeat fraction (%)") +
  theme_gbe()

repeat_b <- ggplot(young_repeat, aes(young_repeat_percent, species_label)) +
  geom_segment(aes(x = 0, xend = young_repeat_percent, yend = species_label), linewidth = 0.3, colour = GBE_COLORS[["light_gray"]]) +
  geom_point(size = 1.5, colour = GBE_COLORS[["vermillion"]]) +
  labs(x = "Aligned repeat fraction below 5% divergence (%)", y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.y = element_text(face = "italic"))

Figure_3 <- repeat_a + repeat_b +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")

record_figure(
  "Figure_3", "Figure_3_repeat_divergence_landscape", 185, 125,
  "S12;S13",
  "RepeatMasker divergence profiles vary among reported repeat categories and species. Curves are equal-weight means across all 34 species after zero-filling absent category/bin combinations; `Unclassified`, when shown, is an annotation category rather than a biological TE order. The species panel sums repeat-aligned bases in bins below 5% divergence.",
  "The left panel shows mean repeat abundance across divergence bins for the six most abundant reported categories and all other categories combined. The right panel ranks species by the fraction of aligned repeats below five percent divergence.",
  "main_candidate"
)

## ---- figure-4
# Figure 4: terminal:internal LTR deletion-footprint proxy.
ltr <- read_supp("Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv") %>%
  filter(analysis_branch == "primary_iqr_filtered_coverage_ge_80pct", ratio_median > 0) %>%
  mutate(
    species_label = species_label(species),
    estimate = log2(ratio_median),
    ci_low = log2(median_bootstrap_ci_low),
    ci_high = log2(median_bootstrap_ci_high)
  )

species_order <- ltr %>%
  arrange(estimate, species) %>%
  pull(species)
species_levels <- species_label(species_order)
ltr <- ltr %>%
  mutate(species_label = factor(species_label, levels = species_levels))

ltr_elements <- read_supp("Supplementary_Data_S14_ltr_terminal_internal_elements.csv") %>%
  filter(primary_ltr_element) %>%
  mutate(
    species_label = factor(species_label(species), levels = species_levels),
    estimate = log2_ratio_terminal_internal_all_positions
  )
if (nrow(ltr_elements) != 380L || n_distinct(ltr_elements$species) != 30L ||
    any(!is.finite(ltr_elements$estimate))) {
  stop("Figure 4 must expose the exact 380-element, 30-species primary LTR branch")
}

ltr_element_checks <- ltr_elements %>%
  group_by(species) %>%
  summarise(
    element_n = n(),
    element_median = median(ratio_terminal_internal_all_positions),
    .groups = "drop"
  ) %>%
  right_join(
    ltr %>% select(species, n_elements, ratio_median),
    by = "species"
  ) %>%
  left_join(
    read_supp("Supplementary_Data_S16_ltr_resource_coverage.csv") %>%
      select(species, n_primary_iqr_filtered_elements),
    by = "species"
  )
if (any(is.na(ltr_element_checks$element_n)) ||
    !all(ltr_element_checks$element_n == ltr_element_checks$n_elements) ||
    !all(ltr_element_checks$n_elements == ltr_element_checks$n_primary_iqr_filtered_elements) ||
    !isTRUE(all.equal(
      ltr_element_checks$element_median,
      ltr_element_checks$ratio_median,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
  stop("Figure 4 element observations, medians, and displayed support must agree")
}

ltr_a <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_point(
    data = ltr_elements,
    aes(estimate, species_label),
    position = position_jitter(width = 0, height = 0.15, seed = 20260811),
    shape = 21,
    size = 1.05,
    stroke = 0.25,
    fill = "white",
    colour = GBE_COLORS[["gray"]],
    alpha = 0.72
  ) +
  geom_errorbar(
    data = ltr,
    aes(xmin = ci_low, xmax = ci_high, y = species_label),
    width = 0.32,
    linewidth = 0.55,
    colour = GBE_COLORS[["blue"]]
  ) +
  geom_point(
    data = ltr,
    aes(estimate, species_label),
    shape = 21,
    size = 2,
    stroke = 0.4,
    fill = GBE_COLORS[["blue"]],
    colour = GBE_COLORS[["charcoal"]]
  ) +
  annotate(
    "text",
    x = -0.04,
    y = Inf,
    label = "equal mean depth\n(1:1)",
    hjust = 1,
    vjust = 1.1,
    size = GBE_LABEL_SIZE_MM,
    colour = GBE_COLORS[["gray"]]
  ) +
  scale_x_continuous(
    breaks = c(-5, -4, -3, -2, -1, 0, 0.5),
    labels = c("1:32", "1:16", "1:8", "1:4", "1:2", "1:1", "1.4:1"),
    expand = expansion(mult = c(0.025, 0.04))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.45, 1.65))) +
  labs(x = "Terminal:internal mean-depth ratio (log2 scale)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_gbe(base_size = 7.5) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_text(face = "italic"),
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.35,
      colour = GBE_COLORS[["charcoal"]]
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 2, 4, 5, unit = "pt")
  )

ltr_n <- ggplot(ltr, aes(0, species_label)) +
  geom_text(aes(label = n_elements), size = GBE_LABEL_SIZE_MM, colour = GBE_COLORS[["charcoal"]]) +
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.45, 1.65))) +
  labs(title = "n", x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 7.5, base_family = resolve_gbe_font()) +
  theme(
    plot.title = element_text(
      size = 7.5,
      face = "bold",
      hjust = 0.5,
      colour = GBE_COLORS[["charcoal"]],
      margin = margin(b = 3, unit = "pt")
    ),
    plot.margin = margin(5, 1, 4, 0, unit = "pt")
  )

Figure_4 <- ltr_a + ltr_n +
  plot_layout(widths = c(1, 0.075))

record_figure(
  "Figure_4", "Figure_4_ltr_terminal_internal_proxy", 185, 165,
  "S14;S15;S16;S17",
  "The zero-aware terminal:internal LTR depth proxy is shown for the 380-element primary branch across 30 species. Light hollow points are individual element-level observations shown descriptively; blue points are species medians and blue intervals are 2,000-replicate element-bootstrap intervals. Elements are not biological replicates, and intervals are conditional on the recovered elements and sequencing resources. The right-hand n column gives the exact element support count for every species. Elements were retained inclusively within each species' two-sided 1.5-IQR fences and only when the left LTR, right LTR, and internal region each had at least 80% positive-depth coverage. The primary branch derives from 408 LTR/Gypsy candidates at least 3,000 bp long with at least five correctly parsed TEsorter domain annotations. Catahoula, kanawha, and valtos lacked tabout resources; lycos had no eligible LTR/Gypsy element with at least five domains. The dashed line marks equal terminal and internal mean depth (1:1). This mapping/deletion-footprint proxy is not an estimate of an ectopic-recombination, solo-LTR, or DNA-loss rate.",
  "Horizontal element-level observations and species median intervals show the log2 terminal-to-internal read-depth ratio for 30 Desmognathus species, ordered by median. A right-hand support count reports the number of retained elements per species, and a dashed vertical line marks equal terminal and internal depth.",
  "sensitivity_only_not_ectopic_rate"
)

## ---- figure-5
# Figure 5: reviewed upper-tail cell and nucleus morphology.
morphology <- read_supp("Supplementary_Data_S19_cell_nucleus_species_estimates.csv") %>%
  mutate(
    species_label = species_label(species),
    species_label = reorder(species_label, cell_area_um2),
    morphology_support = if_else(n_size_specimens == 1L, "One animal", "Two or more animals")
  )

morph_a <- ggplot(morphology, aes(cell_area_um2, species_label)) +
  geom_errorbar(aes(xmin = cell_area_um2_ci_low, xmax = cell_area_um2_ci_high), width = 0, linewidth = 0.4, colour = GBE_COLORS[["orange"]]) +
  geom_point(aes(shape = morphology_support), size = 1.8, colour = GBE_COLORS[["orange"]]) +
  scale_shape_manual(name = "Animal support", values = c("Two or more animals" = 16, "One animal" = 17)) +
  labs(x = expression("Cell area ("*mu*"m"^2*")"), y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.y = element_text(face = "italic"))

morph_b <- ggplot(morphology, aes(nucleus_area_um2, species_label)) +
  geom_errorbar(aes(xmin = nucleus_area_um2_ci_low, xmax = nucleus_area_um2_ci_high), width = 0, linewidth = 0.4, colour = GBE_COLORS[["blue"]]) +
  geom_point(aes(shape = morphology_support), size = 1.8, colour = GBE_COLORS[["blue"]]) +
  scale_shape_manual(name = "Animal support", values = c("Two or more animals" = 16, "One animal" = 17)) +
  guides(shape = "none") +
  labs(x = expression("Nucleus area ("*mu*"m"^2*")"), y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

Figure_5 <- morph_a + morph_b +
  plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")

record_figure(
  "Figure_5", "Figure_5_reviewed_cell_nucleus_morphology", 185, 120,
  "S18;S19;S20",
  "Reviewed upper-tail erythrocyte cell and corresponding nucleus areas differ among 24 finalized species. Points are species medians of all accepted reviewed pairs (up to 50 per species); triangles identify the seven species represented by one animal. Intervals are conditional paired-object bootstrap intervals and do not represent population-level among-individual uncertainty.",
  "Two aligned forest plots rank 24 Desmognathus species by cell area and show corresponding nucleus-area estimates with bootstrap intervals; point shape distinguishes one-animal from multiple-animal support.",
  "main_candidate_upper_tail_estimand"
)

## ---- figure-6
# Figure 6: relative nuclear IOD and pairwise Pagel-lambda PGLS sensitivity.
iod <- read_supp("Supplementary_Data_S22_relative_nuclear_iod_species.csv") %>%
  select(
    species,
    relative_iod_index,
    relative_iod_ci_low,
    relative_iod_ci_high,
    n_images,
    n_specimens
  ) %>%
  left_join(
    read_supp("Supplementary_Data_S24_nuclear_iod_quality_balance.csv") %>%
      select(species, passes_prespecified_balance_thresholds),
    by = "species"
  ) %>%
  mutate(
    iod_support = case_when(
      !passes_prespecified_balance_thresholds ~ "Quality-balance caution",
      n_images == 1L & n_specimens == 1L ~ "One image and animal",
      TRUE ~ "Balanced multiple-image support"
    )
  )
pair_data <- inner_join(iod, morphology %>% mutate(species = as.character(species)) %>% select(-species_label), by = "species")
pgls <- read_supp("Supplementary_Data_S26_pairwise_pgls_sensitivity.csv")

pgls_panel <- function(data, comparison_id, x, y, x_low, x_high, y_low, y_high, x_label, y_label, colour, shape_variable = NULL) {
  stat <- pgls %>% filter(.data$comparison == .env$comparison_id) %>% slice(1)
  x_values <- data[[x]]
  y_values <- data[[y]]
  x_log <- log10(x_values)
  y_log <- log10(y_values)
  x_sequence <- seq(min(x_values), max(x_values), length.out = 200)
  x_z <- (log10(x_sequence) - mean(x_log)) / sd(x_log)
  predicted_z <- stat$pgls_intercept_z + stat$pgls_standardized_beta * x_z
  model_variance <- stat$pgls_var_intercept + 2 * x_z * stat$pgls_cov_intercept_slope + x_z^2 * stat$pgls_var_slope
  critical <- qt(0.975, stat$pgls_degrees_freedom)
  model <- tibble(
    x = x_sequence,
    fit = 10^(mean(y_log) + sd(y_log) * predicted_z),
    low = 10^(mean(y_log) + sd(y_log) * (predicted_z - critical * sqrt(pmax(model_variance, 0)))),
    high = 10^(mean(y_log) + sd(y_log) * (predicted_z + critical * sqrt(pmax(model_variance, 0))))
  )
  lambda_text <- ifelse(stat$pagel_lambda < 0.001, "<0.001", sprintf("%.2f", stat$pagel_lambda))
  p_text <- ifelse(
    stat$pgls_p_value < 0.001,
    "<0.001",
    paste("=", sprintf("%.3f", stat$pgls_p_value))
  )
  annotation <- sprintf(
    "PGLS beta = %.2f (95%% CI %.2f-%.2f)\nP %s; lambda %s; n = %d",
    stat$pgls_standardized_beta, stat$pgls_ci_low, stat$pgls_ci_high,
    p_text, lambda_text, stat$n_species
  )
  panel <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_ribbon(data = model, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = colour, alpha = 0.12) +
    geom_line(data = model, aes(x = x, y = fit), inherit.aes = FALSE, colour = colour, linewidth = 0.65) +
    geom_errorbar(aes(ymin = .data[[y_low]], ymax = .data[[y_high]]), width = 0, linewidth = 0.3, colour = alpha(colour, 0.35)) +
    geom_errorbar(aes(xmin = .data[[x_low]], xmax = .data[[x_high]]), orientation = "y", width = 0, linewidth = 0.3, colour = alpha(colour, 0.35))
  if (is.null(shape_variable)) {
    panel <- panel + geom_point(shape = 21, size = 1.8, stroke = 0.3, fill = colour, colour = "white")
  } else {
    panel <- panel +
      geom_point(aes(shape = .data[[shape_variable]]), size = 2, stroke = 0.3, fill = colour, colour = "white") +
      scale_shape_manual(
        name = "IOD support",
        values = c(
          "Balanced multiple-image support" = 21,
          "One image and animal" = 24,
          "Quality-balance caution" = 23
        )
      )
  }
  panel +
    annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.02, vjust = 1.3, size = GBE_LABEL_SIZE_MM) +
    labs(x = x_label, y = y_label) +
    theme_gbe(base_size = 7.5)
}

iod_a <- pgls_panel(
  pair_data, "relative_iod_vs_nucleus_area", "relative_iod_index", "nucleus_area_um2",
  "relative_iod_ci_low", "relative_iod_ci_high", "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high",
  "Relative nuclear IOD index", expression("Nucleus area ("*mu*"m"^2*")"), GBE_COLORS[["blue"]], "iod_support"
)
iod_b <- pgls_panel(
  pair_data, "relative_iod_vs_cell_area", "relative_iod_index", "cell_area_um2",
  "relative_iod_ci_low", "relative_iod_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  "Relative nuclear IOD index", expression("Cell area ("*mu*"m"^2*")"), GBE_COLORS[["orange"]], "iod_support"
) + guides(shape = "none")
iod_c <- pgls_panel(
  morphology, "nucleus_area_vs_cell_area", "nucleus_area_um2", "cell_area_um2",
  "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  expression("Nucleus area ("*mu*"m"^2*")"), expression("Cell area ("*mu*"m"^2*")"), GBE_COLORS[["purple"]]
)

Figure_6 <- iod_a + iod_b + iod_c +
  plot_layout(widths = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")

record_figure(
  "Figure_6", "Figure_6_relative_iod_pairwise_pgls", 185, 92,
  "S19;S22;S24;S26",
  "Relative nuclear IOD, nucleus area, and cell area are positively associated in Pagel-lambda PGLS sensitivity fits. Circles identify balanced multiple-image IOD support, triangles identify balanced estimates from one image and animal, and diamonds identify four species that fail the declared quality-balance gate. Bars are conditional measurement-bootstrap intervals; lines and ribbons are PGLS fits and model confidence intervals. The IOD index is not an independently calibrated absolute genome size.",
  "Three scatterplots show relative nuclear IOD versus nucleus area, relative nuclear IOD versus cell area, and nucleus area versus cell area, each with measurement intervals and a phylogenetically corrected fitted line; point shape exposes single-image and quality-balance limitations in the IOD panels.",
  "sensitivity_only_relative_iod"
)

## ---- figure-7
# Figure 7: exploratory relative-IOD path comparison and calibration.
path_models <- read_supp("Supplementary_Data_S30_path_model_comparison.csv") %>%
  filter(analysis_type == "primary") %>%
  arrange(delta_CICc) %>%
  slice_head(n = 10) %>%
  mutate(
    model_label = gsub("_", " ", model),
    model_label = sub("^iod", "IOD", model_label),
    model_label = factor(model_label, levels = rev(model_label))
  )
simulation <- read_supp("Supplementary_Data_S36_path_simulation_calibration.csv") %>%
  mutate(
    generating_label = gsub("_", " ", generating_class),
    generating_label = sub("^iod", "IOD", generating_label),
    residual_label = recode(residual_regime, iid = "Independent residuals", brownian = "Brownian residuals")
  )

path_a <- ggplot(
  path_models,
  aes(
    delta_CICc,
    model_label,
    colour = globally_supported,
    shape = globally_supported
  )
) +
  geom_vline(xintercept = 2, linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_point(size = 1.7) +
  scale_color_manual(values = c(`TRUE` = GBE_COLORS[["blue"]], `FALSE` = GBE_COLORS[["vermillion"]]), name = "Global fit", labels = c("Rejected", "Passed")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 17), name = "Global fit", labels = c("Rejected", "Passed")) +
  labs(x = expression(Delta*"CICc from best Path24 model"), y = NULL) +
  theme_gbe(base_size = 7.5)

path_b <- ggplot(
  simulation,
  aes(
    generating_label,
    true_class_top_rate,
    colour = residual_label,
    shape = residual_label,
    group = residual_label
  )
) +
  geom_point(size = 1.7, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = unname(GBE_COLORS[c("blue", "orange")]), name = "Simulation regime") +
  scale_shape_manual(values = c(16, 17), name = "Simulation regime") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Generating equivalence class", y = "True-class top-model rate") +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

Figure_7 <- path_a + path_b +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") & theme(legend.position = "bottom")

record_figure(
  "Figure_7", "Figure_7_exploratory_phylogenetic_path_sensitivity", 185, 108,
  "S30-S36",
  "The finalized 24-species relative-IOD path-model ranking favors one Markov-equivalence class, while simulation recovery varies by generating class and residual regime. Delta CICc is relative to the best Path24 model; color and shape record the global d-separation fit gate. The IOD node is standardized log10 relative nuclear IOD, not absolute genome size, and these results do not orient causal arrows within an equivalence class.",
  "The left panel compares candidate relative-IOD equivalence classes in the finalized 24-species analysis by delta CICc and global-fit status. The right panel shows true-class recovery rates under independent and Brownian residual simulations.",
  "sensitivity_only_relative_iod_no_causal_claim"
)

## ---- path-dag-reference
# Unnumbered reference: one representative DAG for each scored equivalence class.
candidate_dags <- read_csv(
  file.path(
    ROOT, "analyses", "03_phylogenetic_path", "data", "candidate_dags.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
)
path_reference_models <- read_supp(
  "Supplementary_Data_S30_path_model_comparison.csv"
) %>%
  filter(analysis_type == "primary")

dag_representatives <- tribble(
  ~model, ~display_name, ~representative_dag, ~legacy_class, ~class_members,
  "independent", "Independent", "no_edges", "independent", 1L,
  "iod_nucleus_only", "IOD-nucleus only", "GS->NS", "genome_nucleus_only", 2L,
  "nucleus_cell_only", "Nucleus-cell only", "NS->CS", "nucleus_cell_only", 2L,
  "iod_cell_only", "IOD-cell only", "GS->CS", "genome_cell_only", 2L,
  "nucleus_bridge", "Nucleus bridge", "GS->NS;NS->CS", "nucleus_bridge", 3L,
  "nucleus_collider", "Nucleus collider", "GS->NS;CS->NS", "nucleus_collider", 1L,
  "cell_bridge", "Cell bridge", "GS->CS;CS->NS", "cell_bridge", 3L,
  "cell_collider", "Cell collider", "GS->CS;NS->CS", "cell_collider", 1L,
  "iod_bridge", "IOD bridge", "GS->CS;NS->GS", "genome_bridge", 3L,
  "iod_collider", "IOD collider", "NS->GS;CS->GS", "genome_collider", 1L
) %>%
  mutate(model_order = row_number())

if (
  nrow(candidate_dags) != 25L ||
    n_distinct(candidate_dags$equivalence_class) != 11L ||
    n_distinct(candidate_dags$equivalence_class[candidate_dags$testable_by_dsep]) != 10L ||
    nrow(path_reference_models) != 10L ||
    !setequal(path_reference_models$model, dag_representatives$model)
) {
  stop("Path-DAG reference requires the frozen 25-DAG, 11-class, 10-testable-class contract")
}

representative_audit <- dag_representatives %>%
  left_join(
    candidate_dags %>%
      select(
        representative_dag = dag_id,
        candidate_class = equivalence_class,
        candidate_members = class_member_count,
        testable_by_dsep
      ),
    by = "representative_dag",
    relationship = "one-to-one"
  )
if (
  any(is.na(representative_audit$candidate_class)) ||
    any(representative_audit$candidate_class != representative_audit$legacy_class) ||
    any(representative_audit$candidate_members != representative_audit$class_members) ||
    !all(representative_audit$testable_by_dsep)
) {
  stop("Path-DAG representative mapping drifted from candidate_dags.csv")
}

representative_audit <- representative_audit %>%
  mutate(
    panel_label = sprintf(
      "%s\n(%d equivalent %s)",
      display_name,
      class_members,
      if_else(class_members == 1L, "DAG", "DAGs")
    ),
    panel_label = factor(panel_label, levels = panel_label[order(model_order)])
  )

dag_nodes <- crossing(
  representative_audit %>% select(model, panel_label),
  tibble(
    node = c("relative_iod", "nucleus_size", "cell_size"),
    node_label = c("Relative\nIOD", "Nucleus\narea", "Cell\narea"),
    x = c(0.15, 0.85, 0.50),
    y = c(0.23, 0.23, 0.78),
    text_colour = c(
      GBE_COLORS[["charcoal"]],
      GBE_COLORS[["charcoal"]],
      "white"
    )
  )
)

segment_templates <- tribble(
  ~edge_key, ~x, ~y, ~xend, ~yend,
  "relative_iod->nucleus_size", 0.26, 0.23, 0.74, 0.23,
  "nucleus_size->relative_iod", 0.74, 0.23, 0.26, 0.23,
  "relative_iod->cell_size", 0.23, 0.34, 0.44, 0.67,
  "cell_size->relative_iod", 0.44, 0.67, 0.23, 0.34,
  "nucleus_size->cell_size", 0.77, 0.34, 0.56, 0.67,
  "cell_size->nucleus_size", 0.56, 0.67, 0.77, 0.34
)

dag_edges <- representative_audit %>%
  filter(representative_dag != "no_edges") %>%
  separate_rows(representative_dag, sep = ";") %>%
  separate(representative_dag, into = c("parent_legacy", "child_legacy"), sep = "->") %>%
  mutate(
    parent = recode(parent_legacy, GS = "relative_iod", NS = "nucleus_size", CS = "cell_size"),
    child = recode(child_legacy, GS = "relative_iod", NS = "nucleus_size", CS = "cell_size"),
    edge_key = paste(parent, child, sep = "->")
  ) %>%
  left_join(segment_templates, by = "edge_key", relationship = "many-to-one")
if (any(!complete.cases(dag_edges[, c("x", "y", "xend", "yend")]))) {
  stop("Path-DAG reference contains an unsupported directed edge")
}

independent_note <- representative_audit %>%
  filter(model == "independent") %>%
  mutate(x = 0.5, y = 0.49, label = "No directed paths")

Path_DAG_reference <- ggplot() +
  geom_segment(
    data = dag_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.65,
    colour = GBE_COLORS[["charcoal"]],
    arrow = grid::arrow(type = "closed", angle = 22, length = grid::unit(2.2, "mm"))
  ) +
  geom_point(
    data = dag_nodes,
    aes(x = x, y = y, fill = node),
    shape = 21,
    size = 15.5,
    stroke = 0.65,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = dag_nodes,
    aes(x = x, y = y, label = node_label, colour = text_colour),
    size = GBE_LABEL_SIZE_MM,
    fontface = "bold"
  ) +
  geom_text(
    data = independent_note,
    aes(x = x, y = y, label = label),
    size = GBE_LABEL_SIZE_MM,
    colour = GBE_COLORS[["gray"]]
  ) +
  facet_wrap(vars(panel_label), ncol = 2, drop = FALSE) +
  scale_fill_manual(
    values = c(
      relative_iod = GBE_COLORS[["blue"]],
      nucleus_size = GBE_COLORS[["orange"]],
      cell_size = GBE_COLORS[["purple"]]
    ),
    guide = "none"
  ) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    caption = paste0(
      "Representative DAGs of standardized log10 traits; arrows are not unique within classes.\n",
      "Saturated class (6 members) unscored: no d-separation claim."
    )
  ) +
  theme_void(base_size = 8, base_family = resolve_gbe_font()) +
  theme(
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.4,
      colour = GBE_COLORS[["charcoal"]]
    ),
    panel.spacing = grid::unit(5, "pt"),
    strip.background = element_rect(
      fill = "white",
      linewidth = 0.35,
      colour = GBE_COLORS[["charcoal"]]
    ),
    strip.text = element_text(
      size = 8,
      face = "bold",
      lineheight = 0.95,
      margin = margin(3, 2, 3, 2, unit = "pt")
    ),
    plot.caption = element_text(
      size = 7,
      hjust = 0,
      lineheight = 1.05,
      colour = GBE_COLORS[["gray"]],
      margin = margin(t = 5, unit = "pt")
    ),
    plot.margin = margin(4, 4, 4, 4, unit = "pt"),
    aspect.ratio = 0.47
  )

save_gbe_figure(
  Path_DAG_reference,
  file.path(OUT, "Path_DAG_equivalence_class_reference"),
  width_mm = 185,
  height_mm = 220,
  dpi = 300
)

## ---- figure-s1
# Supplementary Figure S1: exact integrated TE34 x Path24 overlap.
te_iod <- read_supp("Supplementary_Data_S40_te_relative_iod_overlap21.csv") %>%
  mutate(
    manual_review_inclusion = as.logical(manual_review_inclusion),
    inclusion_label = if_else(
      manual_review_inclusion,
      "Manual below-target inclusion",
      "Standard finalized support"
    )
  )
if (
  nrow(te_iod) != 21L ||
    n_distinct(te_iod$species) != 21L ||
    sum(te_iod$manual_review_inclusion) != 3L
) {
  stop("Supplementary Figure S1 requires the exact 21-species TE34 x Path24 overlap")
}

Figure_S1 <- ggplot(
  te_iod,
  aes(
    relative_iod_index,
    shannon_entropy,
    fill = inclusion_label,
    shape = inclusion_label
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.35,
    colour = GBE_COLORS[["gray"]]
  ) +
  geom_errorbar(
    aes(
      xmin = relative_iod_ci_low,
      xmax = relative_iod_ci_high
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.35,
    colour = GBE_COLORS[["blue"]],
    alpha = 0.55
  ) +
  geom_point(
    size = 2.4,
    stroke = 0.45,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text_repel(
    data = filter(te_iod, manual_review_inclusion),
    aes(label = scientific_name),
    size = GBE_LABEL_SIZE_MM,
    fontface = "italic",
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 20260729,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Standard finalized support" = GBE_COLORS[["sky"]],
      "Manual below-target inclusion" = GBE_COLORS[["orange"]]
    ),
    name = "Path24 support"
  ) +
  scale_shape_manual(
    values = c(
      "Standard finalized support" = 21,
      "Manual below-target inclusion" = 24
    ),
    name = "Path24 support"
  ) +
  labs(
    x = "Relative nuclear IOD index (Path24 median = 1)",
    y = "Superfamily Shannon entropy (classified composition)"
  ) +
  theme_gbe()

record_figure(
  "Figure_S1", "Figure_S1_te_relative_iod_overlap21", 120, 90,
  "S40",
  "Across the exact 21-species TE34-by-Path24 overlap, classified-superfamily Shannon entropy is shown against the relative nuclear-IOD index. Horizontal intervals are conditional 95% IOD intervals; the dashed line marks the Path24 median anchor of one. Triangles identify the three manually reviewed below-target Path24 inclusions. Relative IOD is an image-derived phenotype, not an independently validated absolute genome size.",
  "Scatterplot of 21 Desmognathus species comparing relative nuclear IOD with TE superfamily Shannon entropy. Three manually retained species use orange triangles; the remaining species use blue circles.",
  "supplementary_sensitivity_only_relative_iod"
)

## ---- figure-manifest
figure_manifest <- bind_rows(manifest_rows)
write_csv(figure_manifest, file.path(OUT, "FIGURE_MANIFEST.csv"))
legend_lines <- c(
  "# Figure legends and alt text",
  "",
  unlist(lapply(seq_len(nrow(figure_manifest)), function(i) {
    row <- figure_manifest[i, ]
    c(
      paste0("## ", row$figure_id),
      "",
      row$legend,
      "",
      paste0("**Alt text:** ", row$alt_text),
      ""
    )
  }))
)
write_lines(
  head(legend_lines, -1),
  file.path(OUT, "FIGURE_LEGENDS_AND_ALT_TEXT.md")
)

message("Wrote ", nrow(figure_manifest), " GBE-styled figures to ", OUT)

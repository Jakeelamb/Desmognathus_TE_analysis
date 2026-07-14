#!/usr/bin/env Rscript

# Build the first GBE-styled ggplot figure set from Publication/datasets.

suppressPackageStartupMessages({
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
    minimum_target_text_pt = GBE_BASE_SIZE - 1,
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

# Figure 1: analysis design, genomic resources, and dnaPipeTE QC.
availability <- read_supp("Supplementary_Data_S01_species_analysis_availability.csv")
te_panel <- read_supp("Supplementary_Data_S02_te_resource_panel34.csv")
cell_panel <- read_supp("Supplementary_Data_S03_cell_linked_panel21.csv")
path_panel <- read_supp("Supplementary_Data_S04_integrated_path_panel18.csv")
assembly <- read_supp("Supplementary_Data_S05_assembly_quality.csv") %>%
  mutate(
    species_label = species_label(species),
    assembly_span_gb = total_sequence_length_bp / 1e9,
    contig_n50_kb = contig_n50_bp / 1e3,
    label = if_else(
      assembly_span_gb >= quantile(assembly_span_gb, 0.9, na.rm = TRUE) |
        contig_n50_kb >= quantile(contig_n50_kb, 0.9, na.rm = TRUE),
      species_label,
      NA_character_
    )
  )
dnapipete_qc <- read_supp("Supplementary_Data_S06_dnapipete_input_quality.csv") %>%
  mutate(species_label = species_label(species))

panel_counts <- tibble(
  panel = factor(
    c("Genomic TE resources", "Linked microscopy", "Integrated phylogenetic"),
    levels = rev(c("Genomic TE resources", "Linked microscopy", "Integrated phylogenetic"))
  ),
  species = c(
    nrow(te_panel),
    nrow(cell_panel),
    nrow(path_panel)
  )
)

design_a <- ggplot(panel_counts, aes(species, panel)) +
  geom_col(width = 0.62, fill = GBE_COLORS[["blue"]]) +
  geom_text(aes(label = species), hjust = -0.35, size = 2.7) +
  scale_x_continuous(limits = c(0, 38), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Species", y = NULL) +
  theme_gbe()

design_b <- ggplot(assembly, aes(assembly_span_gb, contig_n50_kb)) +
  geom_point(shape = 21, size = 2.2, stroke = 0.35, fill = GBE_COLORS[["orange"]], colour = GBE_COLORS[["charcoal"]]) +
  geom_text_repel(data = filter(assembly, !is.na(label)), aes(label = label), size = 2.1, min.segment.length = 0, max.overlaps = Inf, seed = 20260713) +
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

Figure_1 <- design_a + design_b + design_c +
  plot_layout(widths = c(0.75, 1.15, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")

record_figure(
  "Figure_1", "Figure_1_study_design_and_genomic_qc", 185, 76,
  "S01;S05;S06",
  "Analysis-specific panels contain 34 genomic TE species, 21 linked-microscopy species, and 18 species in the integrated phylogenetic overlap. Assembly contiguity and dnaPipeTE classification accounting are shown for the genomic panel; the dashed line is the median unresolved order fraction.",
  "Three panels show sample sizes for the genomic, microscopy, and integrated datasets; assembly span versus contig N50; and repeat-aligned read fraction versus unresolved TE-order mass.",
  "main_candidate"
)

# Figure 2: TE diversity and compositional ordination.
te_diversity <- read_supp("Supplementary_Data_S08_te_diversity.csv") %>%
  mutate(
    level = recode(te_level, order = "Order", superfamily = "Superfamily"),
    mode = recode(
      composition_mode,
      classified_conditional = "Classified only",
      mass_aware_unresolved_bin = "Unresolved mass retained"
    )
  )
te_scores <- read_supp("Supplementary_Data_S09_te_pca_scores.csv") %>%
  filter(te_level == "superfamily", composition_mode == "classified_conditional") %>%
  mutate(
    species_label = species_label(species),
    distance = sqrt(PC1^2 + PC2^2),
    label = if_else(rank(-distance, ties.method = "first") <= 10, species_label, NA_character_)
  )
te_variance <- read_supp("Supplementary_Data_S10_te_pca_variance.csv") %>%
  filter(PC <= 10)

te_a <- ggplot(te_diversity, aes(mode, pielou_evenness, fill = level)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.35, alpha = 0.75) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.68), size = 0.8, alpha = 0.45) +
  scale_fill_gbe(name = "TE level", values = unname(GBE_COLORS[c("blue", "orange")])) +
  labs(x = NULL, y = "Pielou evenness") +
  theme_gbe() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

te_b <- ggplot(te_scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = GBE_COLORS[["light_gray"]]) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = GBE_COLORS[["light_gray"]]) +
  geom_point(size = 1.8, colour = GBE_COLORS[["purple"]]) +
  geom_text_repel(data = filter(te_scores, !is.na(label)), aes(label = label), size = 2.0, max.overlaps = Inf, min.segment.length = 0, seed = 20260713) +
  labs(x = "Superfamily CLR-PC1", y = "Superfamily CLR-PC2") +
  theme_gbe()

te_c <- ggplot(te_variance, aes(PC, VarExplained)) +
  geom_col(width = 0.68, fill = GBE_COLORS[["green"]]) +
  scale_x_continuous(breaks = te_variance$PC) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Principal component", y = "Variance explained") +
  theme_gbe()

Figure_2 <- te_a + te_b + te_c +
  plot_layout(widths = c(1.05, 1.2, 0.85), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")

record_figure(
  "Figure_2", "Figure_2_te_diversity_and_ordination", 185, 82,
  "S08;S09;S10;S11",
  "TE diversity depends on whether unresolved classification mass is retained, while CLR-PCA summarizes species differences among classified superfamilies. Boxes show medians and interquartile ranges; points are species; only the ten most distant PCA scores are labeled.",
  "Three panels compare Pielou evenness under two denominator choices, show species scores on the first two superfamily CLR principal components, and show variance explained by the first ten components.",
  "main_candidate_with_pca_supplementary"
)

# Figure 3: RepeatMasker divergence landscape.
landscape <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
top_orders <- landscape %>%
  group_by(order) %>%
  summarise(hit_bp = sum(hit_bp), .groups = "drop") %>%
  slice_max(hit_bp, n = 6, with_ties = FALSE) %>%
  pull(order)
landscape_summary <- landscape %>%
  mutate(order_group = if_else(order %in% top_orders, order, "Other")) %>%
  group_by(species, order_group, divergence_bin_start_pct) %>%
  summarise(percent_species_hit_bp = sum(percent_species_hit_bp), .groups = "drop") %>%
  group_by(order_group, divergence_bin_start_pct) %>%
  summarise(
    mean_percent = mean(percent_species_hit_bp),
    se_percent = sd(percent_species_hit_bp) / sqrt(n_distinct(species)),
    .groups = "drop"
  )
young_repeat <- landscape %>%
  group_by(species) %>%
  summarise(young_repeat_percent = sum(percent_species_hit_bp[divergence_bin_start_pct < 5]), .groups = "drop") %>%
  mutate(species_label = species_label(species), species_label = reorder(species_label, young_repeat_percent))

repeat_a <- ggplot(landscape_summary, aes(divergence_bin_start_pct, mean_percent, colour = order_group)) +
  geom_line(linewidth = 0.55) +
  scale_color_gbe(
    name = "TE order",
    values = setNames(
      c(unname(GBE_COLORS[c("blue", "orange", "green", "vermillion", "purple", "sky")]), GBE_COLORS[["gray"]]),
      c(top_orders, "Other")
    )
  ) +
  labs(x = "RepeatMasker divergence bin (%)", y = "Mean aligned repeat fraction (%)") +
  theme_gbe()

repeat_b <- ggplot(young_repeat, aes(young_repeat_percent, species_label)) +
  geom_segment(aes(x = 0, xend = young_repeat_percent, yend = species_label), linewidth = 0.3, colour = GBE_COLORS[["light_gray"]]) +
  geom_point(size = 1.5, colour = GBE_COLORS[["vermillion"]]) +
  labs(x = "Aligned repeat fraction below 5% divergence (%)", y = NULL) +
  theme_gbe(base_size = 7.5)

Figure_3 <- repeat_a + repeat_b +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")

record_figure(
  "Figure_3", "Figure_3_repeat_divergence_landscape", 185, 145,
  "S12;S13",
  "RepeatMasker divergence profiles vary among TE orders and species. Curves are across-species means of each species-normalized order/bin contribution; the species panel sums repeat-aligned bases in bins below 5% divergence.",
  "The left panel shows mean repeat abundance across divergence bins for the six most abundant TE orders and all other orders combined. The right panel ranks species by the fraction of aligned repeats below five percent divergence.",
  "main_candidate"
)

# Figure 4: terminal:internal LTR deletion-footprint proxy.
ltr <- read_supp("Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv") %>%
  filter(analysis_branch == "coverage_ge_80pct", ratio_median > 0) %>%
  mutate(
    species_label = species_label(species),
    estimate = log2(ratio_median),
    ci_low = log2(median_bootstrap_ci_low),
    ci_high = log2(median_bootstrap_ci_high),
    species_label = reorder(species_label, estimate)
  )
ltr_coverage <- read_supp("Supplementary_Data_S16_ltr_resource_coverage.csv") %>%
  filter(species %in% ltr$species) %>%
  mutate(species_label = factor(species_label(species), levels = levels(ltr$species_label)))

ltr_a <- ggplot(ltr, aes(estimate, species_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.4, colour = GBE_COLORS[["blue"]]) +
  geom_point(size = 1.5, colour = GBE_COLORS[["blue"]]) +
  labs(x = expression(log[2]~"terminal:internal depth ratio"), y = NULL) +
  theme_gbe(base_size = 7.5)

ltr_b <- ggplot(ltr_coverage, aes(n_usable_elements, species_label)) +
  geom_col(width = 0.62, fill = GBE_COLORS[["orange"]]) +
  labs(x = "Usable LTR elements", y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

Figure_4 <- ltr_a + ltr_b +
  plot_layout(widths = c(1.55, 0.8)) +
  plot_annotation(tag_levels = "A")

record_figure(
  "Figure_4", "Figure_4_ltr_terminal_internal_proxy", 185, 165,
  "S14;S15;S16;S17",
  "The zero-aware terminal:internal LTR depth proxy varies among species but has uneven element support. Points are species medians for elements with at least 80% reported positional coverage; intervals are 2,000-replicate bootstrap confidence intervals; the dashed line marks equal terminal and internal depth.",
  "A forest plot ranks species by median log2 terminal-to-internal LTR read depth and a matched bar plot shows the number of usable elements per species.",
  "sensitivity_only_not_ectopic_rate"
)

# Figure 5: reviewed upper-tail cell and nucleus morphology.
morphology <- read_supp("Supplementary_Data_S19_cell_nucleus_species_estimates.csv") %>%
  mutate(
    species_label = species_label(species),
    species_label = reorder(species_label, cell_area_um2)
  )

morph_a <- ggplot(morphology, aes(cell_area_um2, species_label)) +
  geom_errorbar(aes(xmin = cell_area_um2_ci_low, xmax = cell_area_um2_ci_high), width = 0, linewidth = 0.4, colour = GBE_COLORS[["orange"]]) +
  geom_point(size = 1.6, colour = GBE_COLORS[["orange"]]) +
  labs(x = expression("Cell area ("*mu*"m"^2*")"), y = NULL) +
  theme_gbe(base_size = 7.5)

morph_b <- ggplot(morphology, aes(nucleus_area_um2, species_label)) +
  geom_errorbar(aes(xmin = nucleus_area_um2_ci_low, xmax = nucleus_area_um2_ci_high), width = 0, linewidth = 0.4, colour = GBE_COLORS[["blue"]]) +
  geom_point(size = 1.6, colour = GBE_COLORS[["blue"]]) +
  labs(x = expression("Nucleus area ("*mu*"m"^2*")"), y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

Figure_5 <- morph_a + morph_b +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

record_figure(
  "Figure_5", "Figure_5_reviewed_cell_nucleus_morphology", 185, 132,
  "S18;S19;S20",
  "Reviewed upper-tail erythrocyte cell and corresponding nucleus areas differ among 21 species. Points are medians of the exact 50 accepted pairs per species; intervals are conditional paired-object bootstrap intervals and do not represent population-level among-individual uncertainty.",
  "Two aligned forest plots rank 21 Desmognathus species by cell area and show corresponding nucleus-area estimates with bootstrap intervals.",
  "main_candidate_upper_tail_estimand"
)

# Figure 6: relative nuclear IOD and pairwise Pagel-lambda PGLS sensitivity.
iod <- read_supp("Supplementary_Data_S22_relative_nuclear_iod_species.csv") %>%
  select(species, relative_iod_index, relative_iod_ci_low, relative_iod_ci_high)
pair_data <- inner_join(iod, morphology %>% mutate(species = as.character(species)) %>% select(-species_label), by = "species")
pgls <- read_supp("Supplementary_Data_S26_pairwise_pgls_sensitivity.csv")

pgls_panel <- function(data, comparison_id, x, y, x_low, x_high, y_low, y_high, x_label, y_label, colour) {
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
  annotation <- sprintf(
    "PGLS beta = %.2f (95%% CI %.2f-%.2f)\nP = %.3f; lambda %s; n = %d",
    stat$pgls_standardized_beta, stat$pgls_ci_low, stat$pgls_ci_high,
    stat$pgls_p_value, lambda_text, stat$n_species
  )
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_ribbon(data = model, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = colour, alpha = 0.12) +
    geom_line(data = model, aes(x = x, y = fit), inherit.aes = FALSE, colour = colour, linewidth = 0.65) +
    geom_errorbar(aes(ymin = .data[[y_low]], ymax = .data[[y_high]]), width = 0, linewidth = 0.3, colour = alpha(colour, 0.35)) +
    geom_errorbar(aes(xmin = .data[[x_low]], xmax = .data[[x_high]]), orientation = "y", width = 0, linewidth = 0.3, colour = alpha(colour, 0.35)) +
    geom_point(shape = 21, size = 1.8, stroke = 0.3, fill = colour, colour = "white") +
    annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.02, vjust = 1.3, size = 2.25) +
    labs(x = x_label, y = y_label) +
    theme_gbe(base_size = 7.5)
}

iod_a <- pgls_panel(
  pair_data, "genome_size_vs_nucleus_area", "relative_iod_index", "nucleus_area_um2",
  "relative_iod_ci_low", "relative_iod_ci_high", "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high",
  "Relative nuclear IOD index", expression("Nucleus area ("*mu*"m"^2*")"), GBE_COLORS[["blue"]]
)
iod_b <- pgls_panel(
  pair_data, "genome_size_vs_cell_area", "relative_iod_index", "cell_area_um2",
  "relative_iod_ci_low", "relative_iod_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  "Relative nuclear IOD index", expression("Cell area ("*mu*"m"^2*")"), GBE_COLORS[["orange"]]
)
iod_c <- pgls_panel(
  morphology, "nucleus_area_vs_cell_area", "nucleus_area_um2", "cell_area_um2",
  "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  expression("Nucleus area ("*mu*"m"^2*")"), expression("Cell area ("*mu*"m"^2*")"), GBE_COLORS[["purple"]]
)

Figure_6 <- iod_a + iod_b + iod_c + plot_layout(widths = c(1, 1, 1)) + plot_annotation(tag_levels = "A")

record_figure(
  "Figure_6", "Figure_6_relative_iod_pairwise_pgls", 185, 78,
  "S19;S22;S26;S27-S29",
  "Relative nuclear IOD, nucleus area, and cell area are positively associated in Pagel-lambda PGLS sensitivity fits. Points are species estimates; bars are conditional measurement-bootstrap intervals; lines and ribbons are PGLS fits and model confidence intervals. The IOD index is not an independently calibrated absolute genome size.",
  "Three scatterplots show relative nuclear IOD versus nucleus area, relative nuclear IOD versus cell area, and nucleus area versus cell area, each with measurement intervals and a phylogenetically corrected fitted line.",
  "sensitivity_only_relative_iod"
)

# Figure 7: exploratory path-model comparison and actual-tree calibration.
path_models <- read_supp("Supplementary_Data_S30_path_model_comparison.csv") %>%
  group_by(family) %>%
  arrange(delta_CICc, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    family_label = recode(
      family,
      integrated = "Integrated TE-IOD-morphology",
      iod_morphology = "Relative IOD-morphology",
      te_iod = "TE-relative IOD",
      terminal_internal_iod = "LTR proxy-relative IOD"
    ),
    model_key = paste(family_label, model_label, sep = "___"),
    model_key = factor(model_key, levels = rev(unique(model_key)))
  )
simulation <- read_supp("Supplementary_Data_S36_path_simulation_calibration.csv") %>%
  mutate(
    scenario_label = if_else(
      scenario == "independent_null",
      "Independent null",
      paste0("Observed chain x", effect_scale)
    ),
    family_label = recode(
      family,
      integrated = "Integrated",
      iod_morphology = "IOD-morphology",
      te_iod = "TE-IOD"
    )
  )

path_a <- ggplot(path_models, aes(delta_CICc, model_key, colour = global_fit_pass)) +
  geom_vline(xintercept = 2, linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_point(size = 1.7) +
  facet_grid(family_label ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_color_manual(values = c(`TRUE` = GBE_COLORS[["blue"]], `FALSE` = GBE_COLORS[["vermillion"]]), name = "Global fit", labels = c("Rejected", "Passed")) +
  labs(x = expression(Delta*"CICc from family-best model"), y = NULL) +
  theme_gbe(base_size = 7) +
  theme(strip.text = element_text(size = 7, face = "bold"))

path_b <- ggplot(simulation, aes(scenario_label, primary_calibration_rate, colour = family_label, group = family_label)) +
  geom_errorbar(aes(ymin = primary_rate_ci_low, ymax = primary_rate_ci_high), width = 0.15, linewidth = 0.4, position = position_dodge(width = 0.45)) +
  geom_point(size = 1.7, position = position_dodge(width = 0.45)) +
  scale_color_gbe(name = "Model family", values = unname(GBE_COLORS[c("blue", "orange", "green")])) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "Strict recovery or false-selection rate") +
  theme_gbe(base_size = 7.5) +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

Figure_7 <- path_a + path_b +
  plot_layout(widths = c(1.45, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")

record_figure(
  "Figure_7", "Figure_7_exploratory_phylogenetic_path_sensitivity", 185, 145,
  "S30-S36",
  "Phylogenetic path-model rankings are specification-dependent and actual-tree simulations limit causal discrimination. Delta CICc is relative within each family; color records the global d-separation fit gate. Simulation points and Wilson intervals show strict expected-model recovery under signal and false non-null selection under independent traits.",
  "The left panel compares candidate path models within four families by delta CICc and global-fit status. The right panel shows simulation recovery and false-selection rates across model families and effect sizes.",
  "sensitivity_only_no_causal_claim"
)

figure_manifest <- bind_rows(manifest_rows)
write_csv(figure_manifest, file.path(OUT, "FIGURE_MANIFEST.csv"))
write_lines(
  c(
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
  ),
  file.path(OUT, "FIGURE_LEGENDS_AND_ALT_TEXT.md")
)

message("Wrote ", nrow(figure_manifest), " GBE-styled figures to ", OUT)

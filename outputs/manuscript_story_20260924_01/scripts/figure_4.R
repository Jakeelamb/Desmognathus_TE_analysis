#!/usr/bin/env Rscript
# Display-only assembly of frozen LTR distributions and the existing PGLS fit.
source("scripts/common.R")

out <- file.path(RUN, "figures", "figure_4")
dir.create(file.path(out, "plotted_data"), recursive = TRUE, showWarnings = FALSE)
blue <- GBE_COLORS[["blue"]]
summaries <- read_supp("Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv") %>%
  filter(analysis_branch == "primary_iqr_filtered_coverage_ge_80pct")
rows <- tibble(species = phylogenetic_species_order(summaries$species),
               row = seq_len(nrow(summaries)))
summaries <- summaries %>% left_join(rows, by = "species", relationship = "one-to-one") %>%
  arrange(row) %>% mutate(species_label = species_label(species),
    estimate = log2(ratio_median), ci_low = log2(median_bootstrap_ci_low),
    ci_high = log2(median_bootstrap_ci_high))
elements <- read_supp("Supplementary_Data_S14_ltr_terminal_internal_elements.csv") %>%
  filter(primary_ltr_element) %>%
  left_join(rows, by = "species", relationship = "many-to-one") %>%
  mutate(estimate = log2_ratio_terminal_internal_all_positions)
genomes <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv") %>%
  transmute(species = tree_tip, genome_size_gbp = genome_size_pg_fuscus_anchored * .978,
    genome_size_pg = genome_size_gbp / .978)
summaries <- summaries %>% left_join(genomes, by = "species", relationship = "one-to-one") %>%
  mutate(support_label = paste0(if_else(is.na(genome_size_pg), "n/a",
    sprintf("%.1f", genome_size_pg)), " / ", n_elements))
checks <- elements %>% group_by(species) %>%
  summarise(observed_n = n(), observed_median = median(ratio_terminal_internal_all_positions), .groups = "drop") %>%
  left_join(summaries %>% select(species, n_elements, ratio_median), by = "species", relationship = "one-to-one")
stopifnot(nrow(elements) == 380L, nrow(summaries) == 30L,
  n_distinct(elements$species) == 30L, !anyNA(elements$row),
  all(checks$observed_n == checks$n_elements),
  max(abs(checks$observed_median - checks$ratio_median)) < 1e-12,
  all(elements$ratio_terminal_internal_all_positions > 0),
  all(elements$estimate >= log2(.025)), all(elements$estimate <= log2(1.6)))

bandwidth <- stats::bw.nrd0(elements$estimate)
distribution <- ggplot(elements, aes(estimate, row)) +
  geom_violin(aes(group = species), orientation = "y", position = "identity",
    bw = bandwidth, trim = TRUE, scale = "width", width = .8, quantiles = NULL,
    fill = "white", colour = "#B3B3B3", linewidth = .3) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = .4, colour = "#555555") +
  geom_point(position = position_jitter(width = 0, height = .15, seed = 20260828),
    shape = 21, size = .95, stroke = .25, fill = "white", colour = "#333333", alpha = .75) +
  geom_segment(data = summaries, aes(x = ci_low, xend = ci_high, y = row, yend = row),
    inherit.aes = FALSE, colour = blue, linewidth = .4) +
  geom_segment(data = summaries, aes(x = estimate, xend = estimate, y = row - .25, yend = row + .25),
    inherit.aes = FALSE, colour = blue, linewidth = .8) +
  scale_x_continuous(breaks = log2(c(.025, .05, .1, .25, .5, 1, 1.5)),
    labels = c("0.025", "0.05", "0.1", "0.25", "0.5", "1", "1.5"),
    limits = log2(c(.025, 1.6)), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.5, 30.5), breaks = summaries$row,
    labels = summaries$species_label, expand = c(0, 0),
    sec.axis = dup_axis(breaks = summaries$row, labels = summaries$support_label,
      name = "Genome size (pg/1C) / elements (n)")) +
  labs(x = "LTR terminal-to-internal depth ratio (log2 scale; 1 = equal depth)", y = NULL,
    subtitle = "LTR depth distributions across 30 species",
    caption = "Points: 380 elements. Blue: median and 95% bootstrap interval. n/a: genome size unavailable.") +
  theme_gbe() + theme(axis.text.y.left = element_text(face = "italic", size = 7.5),
    axis.text.y.right = element_text(size = 7.5), axis.title.y.right = element_text(size = 8),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#AAAAAA", linewidth = .3),
    plot.subtitle = element_text(size = 9, face = "bold", margin = margin(b = 5)),
    plot.caption = element_text(size = 7.5), plot.margin = margin(5, 5, 5, 5))

pgls_dir <- file.path(RUN, "inputs", "ltr_pgls")
pgls <- read_csv(file.path(pgls_dir, "species_data.csv"), show_col_types = FALSE)
prediction <- read_csv(file.path(pgls_dir, "fitted_relationship.csv"), show_col_types = FALSE)
model_summary <- read_csv(file.path(pgls_dir, "model_summary.csv"), show_col_types = FALSE)
primary <- filter(model_summary, model == "log_log_lambda")
joined <- pgls %>% left_join(summaries %>% select(species, n_elements, ratio_median,
  median_bootstrap_ci_low, median_bootstrap_ci_high, genome_size_gbp), by = "species", suffix = c("", "_source"))
for (field in c("n_elements", "ratio_median", "median_bootstrap_ci_low", "median_bootstrap_ci_high", "genome_size_gbp")) {
  stopifnot(!anyNA(joined[[paste0(field, "_source")]]),
    max(abs(joined[[field]] - joined[[paste0(field, "_source")]])) < 1e-12)
}
stopifnot(nrow(pgls) == 19L, sum(pgls$n_elements) == 239L,
  nrow(primary) == 1L, primary$n_species == 19L,
  abs(primary$slope - .024196124321278045) < 1e-14,
  abs(primary$ci_low - (-.5608651711893331)) < 1e-14,
  abs(primary$ci_high - .6092574198318892) < 1e-14,
  abs(primary$p_value - .93148859408955) < 1e-14)

# Fixed presentation-only label destinations preserve every observation.
layout <- tibble(
  species = c("campi", "orestes", "welteri", "perlapsus", "anicetus", "aeneus",
    "wrighti", "fuscus", "amphileucus", "auriculatus", "gvnigeusgwotli", "mavrokoilius",
    "apalachicolae", "tilleyi", "marmoratus", "intermedius", "monticola", "ocoee", "bairdi"),
  label_x = c(11.6, 11.6, 11.6, 11.6, 11.6, 9.4, 13.8, 15.7, 17.9, 20.9,
    13.3, 11.6, 16.2, 16.4, 19.3, 19.4, 19.4, 19.4, 19.4),
  label_y = c(.96, .84, .49, .41, .34, .70, 1.18, 1.07, 1.29, 1.09,
    .275, .58, .30, .40, .34, .83, .70, .59, .47),
  label_hjust = c(rep(1, 5), rep(.5, 6), 1, rep(.5, 3), rep(0, 4)))
pgls_display <- pgls %>% left_join(layout, by = "species", relationship = "one-to-one") %>%
  mutate(species_label = species_label(species),
    line_x = if_else(label_hjust == 1, label_x * 1.015, if_else(label_hjust == 0, label_x / 1.015, label_x)),
    line_y = if_else(label_hjust == .5, label_y * if_else(label_y > ratio_median, 1 / 1.035, 1.035), label_y),
    genome_size_pg = genome_size_gbp / .978,
    genome_size_ci_low_pg = genome_size_ci_low_gbp / .978,
    genome_size_ci_high_pg = genome_size_ci_high_gbp / .978,
    label_x_pg = label_x / .978, line_x_pg = line_x / .978)
prediction_display <- prediction %>% mutate(genome_size_pg = genome_size_gbp / .978)
stopifnot(!anyNA(pgls_display))
x_limits <- range(c(7.7, 23, pgls$genome_size_ci_low_gbp, pgls$genome_size_ci_high_gbp)) / .978
y_limits <- range(c(.23, 1.8, pgls$median_bootstrap_ci_low,
  pgls$median_bootstrap_ci_high, prediction$ratio_ci_low, prediction$ratio_ci_high))
annotation <- sprintf("PGLS: n = %d species; lambda < 0.001\nSlope = %.3f (95%% CI %.3f to %.3f); P = %.3f",
  primary$n_species, primary$slope, primary$ci_low, primary$ci_high, primary$p_value)
scatter <- ggplot(pgls_display, aes(genome_size_pg, ratio_median)) +
  geom_ribbon(data = prediction_display, aes(x = genome_size_pg, ymin = ratio_ci_low, ymax = ratio_ci_high),
    inherit.aes = FALSE, fill = blue, alpha = .12) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = .4, colour = "#555555") +
  geom_segment(aes(x = genome_size_ci_low_pg, xend = genome_size_ci_high_pg, yend = ratio_median),
    colour = "#BDBDBD", linewidth = .3) +
  geom_segment(aes(y = median_bootstrap_ci_low, yend = median_bootstrap_ci_high, xend = genome_size_pg),
    colour = "#BDBDBD", linewidth = .3) +
  geom_line(data = prediction_display, aes(x = genome_size_pg, y = ratio_fit),
    inherit.aes = FALSE, colour = blue, linewidth = .7) +
  geom_segment(data = pgls_display, aes(x = genome_size_pg, y = ratio_median, xend = line_x_pg, yend = line_y),
    inherit.aes = FALSE, colour = "#999999", linewidth = .2) +
  geom_point(shape = 21, fill = "white", colour = GBE_COLORS[["charcoal"]], size = 1.8, stroke = .45) +
  geom_text(data = pgls_display, aes(x = label_x_pg, y = label_y, label = species_label, hjust = label_hjust),
    inherit.aes = FALSE, family = GBE_FONT, fontface = "italic", size = GBE_LABEL_SIZE_MM) +
  annotate("text", x = 8.05 / .978, y = 1.69, label = annotation, hjust = 0, vjust = 1,
    family = GBE_FONT, size = GBE_LABEL_SIZE_MM, lineheight = 1.15) +
  scale_x_log10(breaks = c(9, 12, 15, 18, 22), limits = x_limits, expand = expansion(mult = .02)) +
  scale_y_log10(breaks = c(.25, .4, .6, .8, 1, 1.5), limits = y_limits,
    labels = label_number(), expand = expansion(mult = .02)) +
  labs(x = "Genome size (pg/1C; log10 scale)",
    y = "Species median LTR depth ratio\n(log10 scale)",
    subtitle = "No detectable association with genome size",
    caption = "Blue line and band: PGLS fit and 95% coefficient interval. Gray bars: measurement intervals.\nDotted line: equal terminal and internal depth. Genome sizes: assembly-anchored IOD estimates.") +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(plot.subtitle = element_text(size = 9, face = "bold", margin = margin(b = 5)),
    plot.caption = element_text(size = 7.5), plot.margin = margin(5, 5, 5, 5))

legend_a <- paste("The exact 380 primary LTR/Gypsy elements from 30 species are shown in the canonical phylogenetic row order.",
  "Gray horizontal outlines are log2-ratio densities with one pooled bandwidth, equal maximum widths, and tails trimmed to the observed range; violin widths do not represent sample size.",
  "Hollow points show all retained elements, jittered only vertically. Blue bars show released species medians of untransformed ratios; thin blue intervals are the released 95% intervals from 2,000 within-species element bootstraps.",
  "The right-hand labels give genome size (pg/1C) and element count; n/a denotes unavailable genome size without imputation.",
  "The depth ratio is length-weighted mean per-base depth across both terminal LTRs divided by mean depth across the entire intervening sequence, including zero-depth bases.",
  "The primary branch retains LTR/Gypsy elements at least 3,000 bp long with at least five correctly parsed TEsorter domains, inclusive within-species 1.5-IQR fences, and at least 80% positive-depth coverage in all three regions.",
  "The dotted ratio-one line denotes equal mapped depth; it is not an empirically calibrated no-deletion baseline. This proxy is not an ectopic-recombination or DNA-loss rate. Elements are not independent biological replicates.")
legend_b <- paste("Nineteen species with both genome-size estimates and primary LTR ratios contribute one median each (239 elements in total).",
  "The existing stored Pagel-lambda PGLS fit, coefficient interval, and predictions are reused without refitting: log10 median ratio is regressed on log10 genome size. Stored genome coordinates and intervals in Gbp are divided by 0.978 to display pg/1C consistently with the other manuscript figures. This constant unit conversion leaves the log-log slope, its interval, P value, and all fitted ratio values unchanged.",
  sprintf("The slope is %.6f (95%% t interval %.6f to %.6f; P = %.6f); lambda = %.3g at its lower boundary.",
    primary$slope, primary$ci_low, primary$ci_high, primary$p_value, primary$lambda),
  "The blue band is a coefficient interval for the fitted geometric mean, not an individual-species prediction interval. Gray bars show released genome-size intervals and element-bootstrap median intervals; these measurement uncertainties are not propagated through the regression.",
  "The interval permits associations in either direction; the nonsignificant test does not establish equivalence. Original source-tree branch lengths were retained in the stored fit; publication and calibration provenance remain unresolved.")
boundary <- paste("Genome sizes are IOD-derived estimates in pg/1C anchored to the D. fuscus assembly scale; they are not independently measured species C-values and their intervals exclude assembly-anchor uncertainty.",
  "This figure assembles existing results only; no filtering rules, data values, model estimates, or confidence intervals were changed. The excluded leave-one-species-out influence plot is not included.")
alt_a <- "Thirty italic species labels align with gray horizontal density outlines, 380 hollow points, blue medians and intervals, a dotted equal-depth line, and genome-size/count labels at right."
alt_b <- "Nineteen labelled species medians with gray measurement intervals lie around a nearly horizontal blue PGLS line and a broad confidence band on logarithmic axes."
collage <- (distribution + labs(tag = "a") + theme(plot.tag.position = "topleft")) /
  (scatter + labs(tag = "b") + theme(plot.tag.position = "topleft")) + plot_layout(heights = c(1.48, 1))
save_candidate(4, "LTR depth distributions and genome size", collage,
  panels = list(
    a = list(plot = distribution, title = "LTR depth distributions", width_mm = 180, height_mm = 180,
      legend = paste(legend_a, boundary), alt_text = alt_a),
    b = list(plot = scatter, title = "LTR ratio and genome size", width_mm = 180, height_mm = 120,
      legend = paste(legend_b, boundary), alt_text = alt_b)),
  width_mm = 180, height_mm = 280,
  datasets = "S04;S14;S15;stored LTR PGLS species_data/model_summary/fitted_relationship",
  legend = paste(legend_a, legend_b, boundary), alt_text = paste(alt_a, alt_b))
write_csv(elements, file.path(out, "plotted_data", "element_ratios.csv"))
write_csv(summaries, file.path(out, "plotted_data", "species_distributions_and_annotations.csv"))
write_csv(pgls_display, file.path(out, "plotted_data", "pgls_points_and_label_layout.csv"))
write_csv(prediction_display, file.path(out, "plotted_data", "pgls_predictions_display_pg.csv"))
file.copy(file.path(pgls_dir, c("model_summary.csv", "fitted_relationship.csv", "species_data.csv")),
  file.path(out, "plotted_data"), overwrite = TRUE)
write_csv(tibble(check = c("retained_elements", "distribution_species", "pgls_species", "pgls_elements",
  "median_replay_max_absolute_error", "pgls_coefficients_and_ci_exact", "pgls_input_summary_join_exact"),
  value = c(380, 30, 19, 239, max(abs(checks$observed_median - checks$ratio_median)), 1, 1),
  status = "PASS"), file.path(out, "validation.csv"))
write_lines(c("# Figure 4 source joins and validation", "",
  "All inputs are snapshotted under this isolated run. No model was refitted.",
  "S14 primary_ltr_element supplies all 380 raw ratios. S15 primary_iqr_filtered_coverage_ge_80pct supplies the 30 released medians, counts and bootstrap intervals.",
  "S04 tree_tip joins to LTR species. Display values are pg/1C; corresponding Gbp values remain in the exported table for exact comparison to the stored regression inputs. Missing values are displayed as n/a without imputation.",
  "Row order is the shared canonical source-tree tip-y order, matching the existing Figure 4 distribution.",
  "Stored inputs/ltr_pgls species_data.csv is joined by species against S15 and S04; counts, medians, both interval endpoints and genome sizes must match within 1e-12.",
  "Stored model_summary.csv and fitted_relationship.csv are copied byte-for-byte to plotted_data. They supply all PGLS inference and predictions.",
  "Plot-only genome coordinates, interval endpoints, and label positions in Gbp are divided by 0.978 to display pg/1C. pgls_predictions_display_pg.csv and pgls_points_and_label_layout.csv add explicit pg columns. Ratio predictions, ratio intervals, model coefficients, and original CSV copies remain unchanged. A constant predictor-unit conversion changes an intercept but leaves the log-log slope and its inference unchanged; plotting the transformed x coordinates preserves the stored fitted line exactly.",
  "Ratio density is estimated in log2 space for display, matching the prior figure. The PGLS axes are log10, matching the stored model. Both axes label the original ratios directly.",
  "The collage is 180 x 280 mm; isolated panels are 180 x 180 and 180 x 120 mm. All text is Arial and at least 7.5 pt except the shared panel-tag settings, which are larger."),
  file.path(out, "SOURCE_JOINS_AND_QA.md"))
message("Figure 4 complete: 380 elements / 30 species, stored 19-species fit preserved.")

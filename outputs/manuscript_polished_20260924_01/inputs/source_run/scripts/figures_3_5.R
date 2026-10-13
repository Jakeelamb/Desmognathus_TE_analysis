source("scripts/common.R")
suppressPackageStartupMessages(library(phylolm))

# This isolated candidate uses snapshot inputs and writes only inside this run.
path_traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
annotations <- read_supp("Supplementary_Data_S42_te_bonsai_tip_annotations.csv")
diversity <- path_traits %>%
  transmute(species = sub("^D\\. ", "", species),
            genome_size_pg_fuscus_anchored, genome_size_pg_ci_low, genome_size_pg_ci_high) %>%
  inner_join(read_supp("Supplementary_Data_S08_te_diversity.csv"), by = "species", relationship = "one-to-one") %>%
  inner_join(annotations %>% select(species, reproductive_strategy), by = "species", relationship = "one-to-one") %>%
  mutate(species_label = species_label(species),
         reproductive_strategy = factor(reproductive_strategy, levels = te_strategy_levels))
stopifnot(nrow(diversity) == 21L, all(diversity$n_superfamilies_included == 10L))

analysis3 <- file.path(RUN, "figures", "figure_3", "analysis")
dir.create(analysis3, recursive = TRUE, showWarnings = FALSE)
write_csv(diversity, file.path(analysis3, "plotted_species_data.csv"))

# Same full-sample formula, tree pruning, lambda bounds and data as existing Fig 9.
fit_diversity <- function(metric) {
  model_data <- as.data.frame(diversity, stringsAsFactors = FALSE)
  rownames(model_data) <- model_data$species
  fit <- phylolm(reformulate("genome_size_pg_fuscus_anchored", response = metric),
                 data = model_data, phy = prune_source_tree(model_data$species),
                 model = "lambda", lower.bound = 0, upper.bound = 1)
  coefficients <- summary(fit)$coefficients
  slope <- coefficients["genome_size_pg_fuscus_anchored", "Estimate"]
  se <- coefficients["genome_size_pg_fuscus_anchored", "StdErr"]
  t_value <- coefficients["genome_size_pg_fuscus_anchored", "t.value"]
  p_value <- coefficients["genome_size_pg_fuscus_anchored", "p.value"]
  critical <- qt(0.975, df = nrow(model_data) - 2L)
  x <- seq(min(model_data$genome_size_pg_fuscus_anchored),
           max(model_data$genome_size_pg_fuscus_anchored), length.out = 160)
  design <- cbind(`(Intercept)` = 1, genome_size_pg_fuscus_anchored = x)
  fitted <- as.numeric(design %*% coef(fit))
  fitted_se <- sqrt(pmax(rowSums((design %*% vcov(fit)) * design), 0))
  expected <- if (metric == "shannon_entropy") c(4.25, 0.053) else c(5.15, 0.035)
  stopifnot(round(t_value^2, 2) == expected[1], round(p_value, 3) == expected[2], slope < 0)
  list(summary = tibble(metric = metric, response = metric, predictor = "genome_size_pg_fuscus_anchored",
    n_species = nrow(model_data), model = "Pagel_lambda_PGLS", transformation = "untransformed",
    coefficient_scale = "diversity units per pg/1C", intercept = unname(coef(fit)[1]),
    slope = slope, standard_error = se, ci_low = slope - critical * se,
    ci_high = slope + critical * se, degrees_freedom = nrow(model_data) - 2L,
    t_value = t_value, f_statistic = t_value^2, p_value = p_value, lambda = fit$optpar,
    interval_method = "two-sided t coefficient interval", measurement_error_propagated = FALSE),
    predictions = tibble(response = metric, genome_size_pg_fuscus_anchored = x,
                         fit = fitted, low = fitted - critical * fitted_se,
                         high = fitted + critical * fitted_se))
}
diversity_fits <- lapply(c("shannon_entropy", "gini_simpson"), fit_diversity)
diversity_summary <- bind_rows(lapply(diversity_fits, `[[`, "summary"))
diversity_predictions <- bind_rows(lapply(diversity_fits, `[[`, "predictions"))
write_csv(diversity_summary, file.path(analysis3, "model_summary.csv"))
write_csv(diversity_predictions, file.path(analysis3, "fitted_relationships.csv"))
write.tree(prune_source_tree(diversity$species), file.path(analysis3, "pgls_21_species.tre"))

diversity_panel <- function(metric, y_label, label_all = FALSE) {
  statistic <- diversity_summary %>% filter(response == .env$metric)
  predictions <- diversity_predictions %>% filter(response == .env$metric)
  stopifnot(nrow(statistic) == 1L, nrow(predictions) == 160L)
  label_data <- diversity %>% mutate(display_label = if (label_all) species_label else
    ifelse(species %in% c("aeneus", "wrighti", "valtos", "bairdi", "amphileucus", "auriculatus"), species_label, ""))
  annotation <- sprintf("Slope = %.4f [%.4f, %.4f]\nP = %.3f; n = 21; lambda = %.2f",
                        statistic$slope, statistic$ci_low, statistic$ci_high,
                        statistic$p_value, statistic$lambda)
  ggplot(diversity, aes(x = genome_size_pg_fuscus_anchored, y = .data[[metric]])) +
    geom_ribbon(data = predictions, aes(y = NULL, ymin = low, ymax = high),
                fill = "#DDDDDD", alpha = 0.6) +
    geom_line(data = predictions, aes(y = fit), colour = GBE_COLORS[["charcoal"]], linewidth = 0.65) +
    geom_errorbar(aes(xmin = genome_size_pg_ci_low, xmax = genome_size_pg_ci_high,
                      colour = reproductive_strategy), orientation = "y", width = 0,
                  linewidth = 0.35, alpha = 0.55) +
    geom_point(aes(colour = reproductive_strategy, shape = reproductive_strategy),
               size = 2.1, fill = "white", stroke = 0.65) +
    geom_text_repel(data = label_data, aes(label = display_label), family = GBE_FONT, fontface = "italic",
                    size = 7 / ggplot2::.pt, seed = 20260924, max.overlaps = Inf,
                    box.padding = 0.35, point.padding = 0.22, min.segment.length = 0,
                    nudge_x = if (label_all) ifelse(label_data$genome_size_pg_fuscus_anchored < 12, 0,
                      ifelse(label_data$genome_size_pg_fuscus_anchored < median(diversity$genome_size_pg_fuscus_anchored), -3.4, 3.4)) else ifelse(label_data$species %in% c("amphileucus", "auriculatus"), 3.4, 0),
                    force = 3, force_pull = 0.15, xlim = c(8.4, 27.6),
                    segment.color = "#999999", segment.size = 0.2,
                    max.time = Inf, max.iter = 30000,
                    ylim = range(c(diversity[[metric]], predictions$low, predictions$high))) +
    annotate("text", x = Inf, y = Inf, label = annotation, hjust = 1.04, vjust = 1.15,
             family = GBE_FONT, size = 7.5 / ggplot2::.pt, lineheight = 1.1) +
    scale_colour_manual(name = NULL, values = te_strategy_colors,
                        breaks = te_strategy_levels, labels = te_strategy_labels) +
    scale_shape_manual(name = NULL, values = te_strategy_shapes,
                       breaks = te_strategy_levels, labels = te_strategy_labels) +
    scale_x_continuous(limits = c(8, 28), breaks = c(10, 15, 20, 25), expand = expansion(mult = 0.015)) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.23))) +
    labs(x = "Genome size (pg/1C)", y = y_label) + theme_gbe() +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(plot.margin = margin(5, 7, 5, 7), legend.position = "bottom")
}
p3a <- diversity_panel("shannon_entropy", "Shannon entropy")
p3b <- diversity_panel("gini_simpson", "Gini–Simpson index")
figure3 <- (p3a / p3b) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")
legend3 <- paste(
  "TE evenness tends to decrease with genome size across the 21 overlapping species.",
  "Shannon entropy (a) and Gini–Simpson index (b) summarize each species' ten most abundant positive classified dnaPipeTE superfamilies, reclosed to sum to one.",
  "Species points are coloured and shaped by reproductive strategy; a single full-sample Pagel-lambda PGLS is fitted per metric, without subgroup selection.",
  "Six reference species are labelled in the overview; all 21 are labelled in the isolated panels. Label selection changes no observations or fits.",
  "Lines and gray ribbons show the fitted mean and 95% coefficient confidence band. Annotations give unstandardized slopes per pg/1C, 95% t intervals, unadjusted P values, species n, and fitted lambda.",
  "Horizontal bars show released conditional measurement intervals. Genome-size estimates are IOD-derived and D. fuscus assembly-anchored, not independent species C-value measurements; anchor uncertainty is excluded.",
  "Predictor measurement error is displayed but not propagated into the fits. The supplied dated tree is preserved, with unresolved calibration/citation provenance.",
  "These companion diversity metrics are not independent replications and do not establish causal effects.")
save_candidate(3, "TE evenness and genome size", figure3,
  list(a = list(plot = diversity_panel("shannon_entropy", "Shannon entropy", TRUE), title = "Shannon entropy and genome size", width_mm = 180, height_mm = 115,
                legend = legend3, alt_text = "All 21 species have labelled Shannon points and genome-size intervals around a negative full-sample PGLS line."),
       b = list(plot = diversity_panel("gini_simpson", "Gini–Simpson index", TRUE), title = "Gini–Simpson index and genome size", width_mm = 180, height_mm = 115,
                legend = legend3, alt_text = "The same 21 species have labelled Gini–Simpson points and genome-size intervals around a negative full-sample PGLS line.")),
  180, 215, "S04;S08;S42;trees/desmo900dated_source46.tre", legend3,
  "Two labelled 21-species scatterplots show negative full-sample diversity–genome-size fits with intervals; colours and shapes identify reproductive strategy.")

# Figure 5 reads the released S26 model coefficients/covariances without refitting.
analysis5 <- file.path(RUN, "figures", "figure_5", "analysis")
dir.create(analysis5, recursive = TRUE, showWarnings = FALSE)
morphology <- read_supp("Supplementary_Data_S19_cell_nucleus_species_estimates.csv")
pair_data <- path_traits %>%
  select(species, genome_size_pg_fuscus_anchored, genome_size_pg_ci_low, genome_size_pg_ci_high,
         n_genome_size_nuclei, n_genome_size_images, n_genome_size_specimens) %>%
  inner_join(morphology, by = "species", relationship = "one-to-one") %>%
  left_join(read_supp("Supplementary_Data_S24_nuclear_iod_quality_balance.csv") %>%
              select(species, passes_prespecified_balance_thresholds), by = "species", relationship = "one-to-one") %>%
  mutate(genome_size_support = case_when(!passes_prespecified_balance_thresholds ~ "Quality-balance caution",
    n_genome_size_images == 1L & n_genome_size_specimens == 1L ~ "One image and animal",
    TRUE ~ "Balanced multiple-image support"))
stopifnot(nrow(pair_data) == 24L, !anyNA(pair_data), sum(pair_data$n_size_cells) == 1152L)
write_csv(pair_data, file.path(analysis5, "plotted_species_data.csv"))
pgls <- read_supp("Supplementary_Data_S26_pairwise_pgls_sensitivity.csv")
write_csv(pgls, file.path(analysis5, "model_summary.csv"))

cell_specs <- tibble(panel = letters[1:3],
  comparison = c("genome_size_vs_cell_area", "genome_size_vs_nucleus_area", "nucleus_area_vs_cell_area"),
  x = c("genome_size_pg_fuscus_anchored", "genome_size_pg_fuscus_anchored", "nucleus_area_um2"),
  y = c("cell_area_um2", "nucleus_area_um2", "cell_area_um2"),
  x_low = c("genome_size_pg_ci_low", "genome_size_pg_ci_low", "nucleus_area_um2_ci_low"),
  x_high = c("genome_size_pg_ci_high", "genome_size_pg_ci_high", "nucleus_area_um2_ci_high"),
  y_low = c("cell_area_um2_ci_low", "nucleus_area_um2_ci_low", "cell_area_um2_ci_low"),
  y_high = c("cell_area_um2_ci_high", "nucleus_area_um2_ci_high", "cell_area_um2_ci_high"),
  x_label = c("Genome size (pg/1C)", "Genome size (pg/1C)", "Nucleus area (µm²)"),
  y_label = c("Cell area (µm²)", "Nucleus area (µm²)", "Cell area (µm²)"))
predictions5 <- lapply(seq_len(nrow(cell_specs)), function(i) {
  spec <- cell_specs[i, ]; stat <- pgls %>% filter(comparison == spec$comparison)
  x_log <- log10(pair_data[[spec$x]]); y_log <- log10(pair_data[[spec$y]])
  stopifnot(nrow(stat) == 1L, stat$n_species == 24L,
    abs(mean(x_log) - stat$predictor_log10_mean) < 1e-7,
    abs(sd(x_log) - stat$predictor_log10_sd) < 1e-7,
    abs(mean(y_log) - stat$outcome_log10_mean) < 1e-7,
    abs(sd(y_log) - stat$outcome_log10_sd) < 1e-7)
  x_sequence <- seq(min(pair_data[[spec$x]]), max(pair_data[[spec$x]]), length.out = 200)
  x_z <- (log10(x_sequence) - mean(x_log)) / sd(x_log)
  predicted_z <- stat$pgls_intercept_z + stat$pgls_standardized_beta * x_z
  variance <- stat$pgls_var_intercept + 2 * x_z * stat$pgls_cov_intercept_slope + x_z^2 * stat$pgls_var_slope
  critical <- qt(0.975, stat$pgls_degrees_freedom)
  tibble(comparison = spec$comparison, x = x_sequence,
    fit = 10^(mean(y_log) + sd(y_log) * predicted_z),
    low = 10^(mean(y_log) + sd(y_log) * (predicted_z - critical * sqrt(pmax(variance, 0)))),
    high = 10^(mean(y_log) + sd(y_log) * (predicted_z + critical * sqrt(pmax(variance, 0)))))
}) %>% bind_rows()
write_csv(predictions5, file.path(analysis5, "fitted_relationships.csv"))

cell_panel <- function(i, label_all = FALSE) {
  spec <- cell_specs[i, ]; stat <- pgls %>% filter(comparison == spec$comparison)
  model <- predictions5 %>% filter(comparison == spec$comparison)
  point_labels <- pair_data %>% mutate(display_label = if (label_all) species else
    ifelse(species %in% c("D. aeneus", "D. valtos", "D. wrighti", "D. campi", "D. bairdi"), species, ""))
  p_label <- if (stat$pgls_p_value < 0.001) formatC(stat$pgls_p_value, format = "e", digits = 2) else sprintf("%.3f", stat$pgls_p_value)
  annotation <- sprintf("β = %.3f [%.3f, %.3f]\nP = %s; n = 24",
                        stat$pgls_standardized_beta, stat$pgls_ci_low, stat$pgls_ci_high, p_label)
  label_y_range <- range(c(pair_data[[spec$y_low]], pair_data[[spec$y_high]], model$low, model$high))
  plot <- ggplot(pair_data, aes(x = .data[[spec$x]], y = .data[[spec$y]])) +
    geom_ribbon(data = model, aes(x = x, y = NULL, ymin = low, ymax = high), inherit.aes = FALSE,
                fill = "#DDDDDD", alpha = 0.6) +
    geom_line(data = model, aes(x = x, y = fit), inherit.aes = FALSE, linewidth = 0.65, colour = "#222222") +
    geom_errorbar(aes(ymin = .data[[spec$y_low]], ymax = .data[[spec$y_high]]), width = 0,
                  colour = "#AAAAAA", linewidth = 0.3) +
    geom_errorbar(aes(xmin = .data[[spec$x_low]], xmax = .data[[spec$x_high]]), orientation = "y", width = 0,
                  colour = "#AAAAAA", linewidth = 0.3) +
    geom_point(shape = 21, fill = "white", colour = "#222222", stroke = 0.6, size = 2) +
    geom_text_repel(data = point_labels, aes(label = display_label), family = GBE_FONT, fontface = "italic",
                    size = 7 / ggplot2::.pt, seed = 20260924 + i, max.overlaps = Inf,
                    box.padding = 0.35, point.padding = 0.2, min.segment.length = 0,
                    nudge_x = if (label_all) ifelse(point_labels[[spec$x]] < 12, 0,
                      ifelse(point_labels[[spec$x]] < median(pair_data[[spec$x]]), -3.4, 3.4)) else 0,
                    force = 3, force_pull = 0.15,
                    segment.color = "#999999", segment.size = 0.2,
                    max.time = Inf, max.iter = 30000, ylim = label_y_range) +
    annotate("text", x = -Inf, y = Inf, label = annotation, hjust = -0.03, vjust = 1.1,
             family = GBE_FONT, size = 7.5 / ggplot2::.pt, lineheight = 1.1) +
    labs(x = spec$x_label, y = spec$y_label) + theme_gbe() +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
    theme(plot.margin = margin(5, 7, 5, 7), legend.position = "none")
  if (i %in% c(1, 2)) plot <- plot + scale_x_continuous(limits = c(8, 28), breaks = c(10, 15, 20, 25), expand = expansion(mult = 0.015))
  if (i == 3) plot <- plot + scale_x_continuous(expand = expansion(mult = c(0.12, 0.12)))
  # Identical area scale in the two cell-area panels, retaining the model bands and intervals.
  if (i %in% c(1, 3)) plot <- plot + coord_cartesian(ylim = c(185, 395), expand = FALSE)
  plot
}
p5a <- cell_panel(1, TRUE); p5b <- cell_panel(2); p5c <- cell_panel(3)
figure5 <- (p5a / (p5b + p5c)) + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "a")
legend5 <- paste(
  "Genome size is positively associated with selected erythrocyte area across the sampled species.",
  "All 24 species appear in each panel: cell area versus genome size (a), nucleus area versus genome size (b), and cell area versus nucleus area (c).",
  "The larger first panel labels all species; panels b and c label five reference species for orientation, while all species are labelled in the isolated exports.",
  "The exact released S26 Pagel-lambda PGLS coefficients and covariance matrices are reused without refitting. Fits use log10 traits followed by sample standardization, then are back-transformed for display.",
  "Annotations give standardized slopes, 95% t intervals and unadjusted P values; fitted lambda is approximately zero for all three fits.",
  "Black lines and gray ribbons show fitted means and 95% model confidence bands. Horizontal and vertical bars show released conditional measurement intervals.",
  "Cellular summaries describe selected large erythrocytes and their nuclei, based on 1,152 accepted pairs. Seven species have one-animal morphology support.",
  "Genome sizes are IOD-derived D. fuscus assembly-anchored estimates, not independent species C-values; three estimates have one-image/one-animal support and four fail the recorded IOD quality-balance gate (S24).",
  "Uniform circles retain every observation; support flags remain in the plotted-data export. Measurement intervals exclude anchor uncertainty and do not represent population-level sampling uncertainty. Associations do not establish causal direction.")
save_candidate(5, "Cellular dimensions and genome size", figure5,
  list(a = list(plot = p5a, title = "Cell area and genome size", width_mm = 180, height_mm = 125,
                legend = legend5, alt_text = "All 24 labelled species show cell-area and genome-size intervals, with a positive fitted relationship and confidence band."),
       b = list(plot = cell_panel(2, TRUE), title = "Nucleus area and genome size", width_mm = 180, height_mm = 125,
                legend = legend5, alt_text = "All 24 labelled species show nucleus-area and genome-size intervals around the weaker positive fitted relationship."),
       c = list(plot = cell_panel(3, TRUE), title = "Cell area and nucleus area", width_mm = 180, height_mm = 125,
                legend = legend5, alt_text = "All 24 labelled species show paired cell-area and nucleus-area intervals around a positive fitted relationship.")),
  180, 220, "S04;S19;S24;S26", legend5,
  "A large 24-species genome–cell panel sits above genome–nucleus and nucleus–cell panels, showing all observations, intervals, and stored PGLS fits.")
write_lines(c("Checks passed: full-sample diversity model F/P match the former Figure 9 after rounding; no subgroup model selection.",
              "Figure 3: 21 species in each metric, both raw-data and fitted-model numerical exports.",
              "Figure 5: 24 species per panel and 1,152 accepted paired objects; S26 coefficients reused, no model refit.",
              "Figure 5 predictor/outcome log10 means and SDs match stored S26 to 1e-7.",
              "S26 t intervals are used; S28 normal-approximation intervals are not substituted."),
            file.path(RUN, "qa", "figures_3_5_numerical_checks.txt"))
message("Completed isolated figures 3 and 5 with all panels and numerical exports.")

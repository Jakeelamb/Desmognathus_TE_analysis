source("scripts/common.R")

# Plot the frozen path-analysis release; no models are refitted here.
out <- file.path(RUN, "figures", "figure_6")
analysis_dir <- file.path(out, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
input_files <- c(
  file.path(DATA, c("Supplementary_Data_S30_path_model_comparison.csv",
    "Supplementary_Data_S32_path_model_stability.csv",
    "Supplementary_Data_S33_path_tree_sensitivity.csv",
    "Supplementary_Data_S36_path_simulation_calibration.csv",
    "Supplementary_Data_S37_path_dsep_tests.csv")),
  file.path(RUN, "inputs", "path", "candidate_dags.csv"),
  file.path(DATA, "Supplementary_Data_S28_pairwise_phylogenetic_regressions.csv"))
hash_inputs <- function() system2("sha256sum", shQuote(input_files), stdout = TRUE)
input_hashes <- hash_inputs()
primary <- read_csv(input_files[1], show_col_types = FALSE)
stability <- read_csv(input_files[2], show_col_types = FALSE)
trees <- read_csv(input_files[3], show_col_types = FALSE)
calibration <- read_csv(input_files[4], show_col_types = FALSE)
dsep <- read_csv(input_files[5], show_col_types = FALSE)
dags <- read_csv(input_files[6], show_col_types = FALSE)
regressions <- read_csv(input_files[7], show_col_types = FALSE)

stopifnot(nrow(primary) == 10L, all(primary$n_species == 24L),
  primary$model[which.min(primary$CICc)] == "cell_bridge",
  abs(sum(primary$w) - 1) < 1e-12,
  sum(primary$globally_supported) == 1L,
  all(primary$globally_supported == (primary$p >= 0.05)))
winner <- primary %>% filter(model == "cell_bridge")
stopifnot(abs(winner$w - 0.9167423772846516) < 1e-12,
  abs(winner$p - 0.975509706887934) < 1e-12)
bridge_claim <- dsep %>% filter(candidate_model == "cell_bridge")
stopifnot(nrow(bridge_claim) == 1L, bridge_claim$independence_claim == "n ~ c + g",
  abs(bridge_claim$component_p - winner$p) < 1e-12)
winner_dags <- dags %>% filter(equivalence_class == "cell_bridge")
stopifnot(nrow(winner_dags) == 3L,
  setequal(winner_dags$dag_id, c("GS->CS;CS->NS", "NS->CS;CS->GS", "CS->GS;CS->NS")),
  all(winner_dags$class_member_count == 3L))

blue <- GBE_COLORS[["blue"]]
gray <- GBE_COLORS[["gray"]]
black <- GBE_COLORS[["charcoal"]]
text_size <- 7.5 / ggplot2::.pt
panel_theme <- theme_gbe() + theme(
  plot.title = element_text(size = 9, face = "bold", margin = margin(b = 4)),
  plot.subtitle = element_text(size = 7.5, lineheight = 1.1, margin = margin(b = 5)),
  plot.caption = element_text(size = 7, lineheight = 1.1, hjust = 0, margin = margin(t = 6)),
  plot.margin = margin(6, 8, 6, 8),
  axis.text.y = element_text(size = 7.5),
  legend.text = element_text(size = 7.5))

# Each row is one source-enumerated member of the same equivalence class.
dag_order <- c("GS->CS;CS->NS", "NS->CS;CS->GS", "CS->GS;CS->NS")
node_positions <- tibble(node = c("GS", "CS", "NS"), x = c(1.6, 5, 8.4),
  label = c("Estimated\ngenome size", "Cell area", "Nucleus area"))
dag_nodes <- bind_rows(lapply(seq_along(dag_order), function(i) {
  node_positions %>% mutate(dag_id = dag_order[i], y = 4 - i)
}))
dag_edges <- bind_rows(lapply(seq_along(dag_order), function(i) {
  e <- strsplit(strsplit(dag_order[i], ";", fixed = TRUE)[[1]], "->", fixed = TRUE)
  bind_rows(lapply(e, function(pair) {
    x_start <- node_positions$x[match(pair[1], node_positions$node)]
    x_end <- node_positions$x[match(pair[2], node_positions$node)]
    direction <- sign(x_end - x_start)
    tibble(dag_id = dag_order[i], from = pair[1], to = pair[2],
      x = x_start + direction * 1.18, xend = x_end - direction * 1.18,
      y = 4 - i, yend = 4 - i)
  }))
}))
p6a <- ggplot() +
  geom_rect(data = dag_nodes, aes(xmin = x - 1.10, xmax = x + 1.10,
    ymin = y - 0.40, ymax = y + 0.40, fill = node == "CS"),
    colour = "#9A9A9A", linewidth = 0.35) +
  geom_segment(data = dag_edges, aes(x = x, xend = xend, y = y, yend = yend),
    arrow = grid::arrow(length = grid::unit(1.9, "mm"), type = "closed"),
    linewidth = 0.6, colour = blue) +
  geom_text(data = dag_nodes, aes(x = x, y = y, label = label),
    family = GBE_FONT, size = 8 / ggplot2::.pt, lineheight = 1.05) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#E5F1F8"), guide = "none") +
  scale_x_continuous(limits = c(0.1, 9.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, 3.5), expand = c(0, 0)) +
  labs(title = "Alternative directions within the best-supported path class",
    subtitle = "All three predict genome–nucleus independence conditional on cell area.",
    caption = "Equally compatible orientations: the analysis supports a conditional structure, while arrow direction remains unresolved.") +
  panel_theme + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title = element_blank())

model_names <- c(cell_bridge = "Cell bridge", cell_collider = "Cell collider",
  genome_bridge = "Genome bridge", genome_cell_only = "Genome–cell only",
  genome_collider = "Genome collider", nucleus_bridge = "Nucleus bridge",
  nucleus_cell_only = "Nucleus–cell only", nucleus_collider = "Nucleus collider",
  genome_nucleus_only = "Genome–nucleus only", independent = "Independence")
primary_plot <- primary %>% arrange(rank) %>%
  mutate(y = 11 - rank, display_label = unname(model_names[model]),
    highlight = ifelse(model == "cell_bridge", "winner", "other"),
    fit_status = ifelse(globally_supported, "Pass (P ≥ 0.05)", "Rejected (P < 0.05)"),
    weight_label = ifelse(w >= 0.001, sprintf("%.3f", w), "<0.001"))
p6b <- ggplot(primary_plot, aes(x = delta_CICc, y = y)) +
  geom_vline(xintercept = 2, colour = "#A0A0A0", linetype = "dashed", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = delta_CICc, yend = y, colour = highlight), linewidth = 0.45) +
  geom_point(aes(shape = fit_status, colour = highlight), size = 2.25, stroke = 0.65) +
  geom_text(aes(x = 34.8, label = weight_label), hjust = 1, size = text_size, family = GBE_FONT) +
  annotate("text", x = 34.8, y = 10.8, label = "Weight", hjust = 1,
    family = GBE_FONT, size = text_size, fontface = "bold") +
  scale_colour_manual(values = c(winner = blue, other = gray), guide = "none") +
  scale_shape_manual(name = "Global fit", values = c(`Pass (P ≥ 0.05)` = 16, `Rejected (P < 0.05)` = 4),
    breaks = c("Pass (P ≥ 0.05)", "Rejected (P < 0.05)")) +
  scale_y_continuous(breaks = primary_plot$y, labels = primary_plot$display_label,
    limits = c(0.6, 11.1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(-0.4, 35.2), expand = c(0, 0)) +
  labs(title = "Primary comparison of ten path classes", subtitle =
    sprintf("24 species; cell bridge: weight = %.3f, global-fit P = %.4f", winner$w, winner$p),
    x = "ΔCICc from the best-supported class", y = NULL,
    caption = "Dashed line: ΔCICc = 2. The nucleus-bridge class is 13.77 CICc units behind the cell-bridge class.") +
  panel_theme + theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = "bottom", legend.key.width = grid::unit(3, "mm"))

stability_order <- c("measurement_bootstrap_primary", "leave_one_species_out_24",
  "evolutionary_model_choice_primary", "published_bootstrap_trees_18")
stability_plot <- stability %>% filter(model == "cell_bridge") %>%
  mutate(order = match(analysis, stability_order)) %>% arrange(order) %>%
  mutate(y = 5 - order,
    label = c("Measurement resampling\n24 species", "Leave one species out\n23 species per fit",
      "Evolution models\n24 species", "Published trees\n18 species"),
    count_label = paste0(top_count, "/", n_fits), top_percent = 100 * top_rate)
stopifnot(nrow(stability_plot) == 4L,
  identical(stability_plot$top_count, c(246, 24, 4, 200)),
  identical(stability_plot$n_fits, c(250, 24, 4, 200)))
tree_nucleus <- stability %>% filter(analysis == "published_bootstrap_trees_18", model == "nucleus_bridge")
stopifnot(tree_nucleus$admissible_competitive_count == 200L,
  tree_nucleus$admissible_competitive_rate == 1)
tree_subset <- trees %>% filter(replicate_id %in% c("focal_tree_18", "published_main_18"),
  model %in% c("cell_bridge", "nucleus_bridge"))
stopifnot(nrow(tree_subset) == 4L, all(tree_subset$n_species == 18L),
  all(tree_subset$delta_CICc <= 2), all(tree_subset$globally_supported))
p6stability <- ggplot(stability_plot, aes(x = top_percent, y = y)) +
  geom_segment(aes(x = 0, xend = top_percent, yend = y), colour = blue, linewidth = 5.5, lineend = "butt") +
  geom_text(aes(x = 50, label = count_label), family = GBE_FONT, colour = "white", size = text_size) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.6, 4.4), breaks = stability_plot$y, labels = stability_plot$label,
    expand = c(0, 0)) +
  labs(title = "Top-rank stability", subtitle = "How often the cell-bridge class ranks first",
    x = "Fits ranked first (%)", y = NULL,
    caption = "All 200 tree fits also support nucleus bridge\nwithin ΔCICc ≤ 2; 18-species support is shared.") +
  panel_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7, lineheight = 1.05))

calibration_order <- c("cell_bridge", "nucleus_bridge", "genome_bridge", "independent")
calibration_plot <- calibration %>%
  mutate(class_y = 5 - match(generating_class, calibration_order),
    regime_label = ifelse(residual_regime == "brownian", "Brownian", "Independent"),
    y = class_y + ifelse(residual_regime == "brownian", 0.15, -0.15),
    recovery_percent = 100 * true_class_top_rate,
    display_label = paste0(round(recovery_percent), "%"))
stopifnot(nrow(calibration_plot) == 8L, all(calibration_plot$n_simulations == 100L),
  all(calibration_plot$recovery_percent >= 0 & calibration_plot$recovery_percent <= 100),
  setequal(calibration_plot$generating_class, calibration_order),
  setequal(calibration_plot$residual_regime, c("brownian", "iid")))
p6calibration <- ggplot(calibration_plot, aes(x = recovery_percent, y = y, colour = regime_label)) +
  geom_segment(aes(x = 0, xend = recovery_percent, yend = y), linewidth = 0.4, alpha = 0.65) +
  geom_point(aes(shape = regime_label), size = 1.7, stroke = 0.6) +
  geom_text(aes(x = recovery_percent + 4, label = display_label), colour = black,
    hjust = 0, family = GBE_FONT, size = 7 / ggplot2::.pt) +
  scale_colour_manual(name = "Residuals", values = c(Brownian = blue, Independent = GBE_COLORS[["vermillion"]])) +
  scale_shape_manual(name = "Residuals", values = c(Brownian = 16, Independent = 15)) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.6, 4.4), breaks = 4:1,
    labels = model_names[calibration_order], expand = c(0, 0)) +
  labs(title = "Recovery of the generating class", subtitle = "100 simulations per class and residual regime",
    x = "Generating class ranked first (%)", y = NULL,
    caption = "Calibration measures recovery of the path class;\nit does not resolve directions within that class.") +
  panel_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7), legend.text = element_text(size = 7),
    legend.title = element_text(size = 7), legend.key.width = grid::unit(3, "mm"))

contrast_order <- c("genome_to_nucleus", "genome_and_cell_to_nucleus",
  "genome_to_cell", "genome_and_nucleus_to_cell")
contrast_plot <- regressions %>% filter(predictor == "genome_size", relation_id %in% contrast_order) %>%
  mutate(order = match(relation_id, contrast_order)) %>% arrange(order) %>%
  mutate(panel = factor(ifelse(outcome == "nucleus_size", "Outcome: nucleus area", "Outcome: cell area"),
      levels = c("Outcome: nucleus area", "Outcome: cell area")),
    estimate_type = ifelse(grepl("and", relation_id), "Adjusted", "Marginal"),
    row_label = c("Marginal", "Adjusted for\ncell area", "Marginal", "Adjusted for\nnucleus area"),
    row_label = factor(row_label, levels = c("Adjusted for\ncell area", "Adjusted for\nnucleus area", "Marginal")))
stopifnot(nrow(contrast_plot) == 4L, all(contrast_plot$n_species == 24L),
  all(abs(contrast_plot$ci_low - (contrast_plot$standardized_coefficient - 1.96 * contrast_plot$standard_error)) < 1e-12),
  all(abs(contrast_plot$ci_high - (contrast_plot$standardized_coefficient + 1.96 * contrast_plot$standard_error)) < 1e-12),
  abs(contrast_plot$standardized_coefficient[2] - 0.0079231983550237) < 1e-12,
  abs(contrast_plot$standardized_coefficient[4] - 0.5748580416809171) < 1e-12)
p6contrast <- ggplot(contrast_plot, aes(x = standardized_coefficient, y = row_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#A0A0A0") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = estimate_type),
    orientation = "y", width = 0.12, linewidth = 0.5) +
  geom_point(aes(colour = estimate_type, shape = estimate_type), size = 2.3, stroke = 0.6) +
  scale_colour_manual(values = c(Marginal = gray, Adjusted = blue), guide = "none") +
  scale_shape_manual(values = c(Marginal = 16, Adjusted = 17), guide = "none") +
  scale_x_continuous(limits = c(-0.6, 1.1), breaks = c(-0.5, 0, 0.5, 1), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  labs(title = "Conditional associations distinguish the bridge classes",
    subtitle = "Estimated genome size is the predictor in both comparisons; 24 species.",
    x = "Standardized coefficient (95% interval)", y = NULL,
    caption = "Genome–nucleus point estimate is near zero after cell-area adjustment.\nGenome–cell association remains positive after nucleus-area adjustment.") +
  panel_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5, lineheight = 1.05), panel.spacing.x = grid::unit(5, "mm"))

figure6 <- (p6a / p6b / p6contrast / p6stability) +
  plot_layout(heights = c(0.85, 1.4, 1.0, 0.9)) + plot_annotation(tag_levels = "a")
legend_a <- paste(
  "The cell-bridge class contains three distinct directed acyclic graphs, reproduced from the frozen candidate-DAG enumeration.",
  "Estimated genome size -> cell area -> nucleus area, its reverse chain, and cell area -> both traits all imply genome–nucleus conditional independence given cell area.",
  "These orientations are statistically indistinguishable in this three-trait comparison; their common conditional structure is the supported result.")
legend_b <- paste(
  "Ten testable path-equivalence classes are ranked by CICc using 24 species, standardized log10 traits, the supplied dated tree and Pagel-lambda regressions.",
  "The cell-bridge class has CICc = 13.3829, Fisher C = 0.04959, one independence claim, global-fit P = 0.97551, and relative weight = 0.91674.",
  "Its tested claim is the residual genome–nucleus association after conditioning on cell area (P = 0.97551).",
  "All nine other tested classes fail the P >= 0.05 global-fit criterion; the nucleus-bridge class has delta CICc = 13.77171 and global-fit P = 0.0009971.",
  "Weights compare only the supplied ten classes; six saturated DAGs impose no testable d-separation claim and are not ranked.",
  "All primary component regressions have fitted lambda approximately 1e-7. No direction-specific coefficients are assigned to the winning equivalence class.")
legend_stability <- paste(
  "Cell bridge ranks first in 246/250 measurement resamples, 24/24 leave-one-species-out fits, 4/4 evolutionary-model specifications and 200/200 published-tree fits.",
  "The first and third checks retain 24 species, species omission retains 23 per fit, and the published-tree analysis contains only 18 species.",
  "All 200 published-tree fits also retain nucleus bridge as globally admissible and within delta CICc <= 2, so top ranking does not indicate unique support.",
  "On the focal tree pruned to the same 18 species, cell-bridge weight is 0.43357 and nucleus-bridge weight is 0.32240 (delta CICc = 0.59252); the published main tree gives essentially the same values.",
  "These reduced-species results primarily expose dependence of discrimination on taxon coverage, rather than establishing strong topology robustness.")
legend_calibration <- paste(
  "All eight released calibration settings are shown: four generating classes under Brownian and independent residual regimes, with 100 simulations each.",
  "True-class top-recovery rates are cell bridge 78%/80%, nucleus bridge 49%/41%, genome bridge 29%/33%, and independence 40%/36% (Brownian/independent).",
  "These are recovery frequencies under the specified simulation settings, not posterior probabilities for the observed causal directions.")
legend_contrast <- paste(
  "Released S28 Pagel-lambda coefficients compare marginal and conditional associations, with standardized log10 traits and estimated genome size as predictor in all four regressions.",
  "The genome–nucleus coefficient is 0.41902 (95% normal interval 0.03960 to 0.79844) marginally, and 0.00792 (-0.49195 to 0.50780) after adjusting for cell area.",
  "The genome–cell coefficient is 0.71666 (0.42523 to 1.00810) marginally, and 0.57486 (0.27994 to 0.86977) after adjusting for nucleus area.",
  "Intervals are the stored S28 normal approximations (coefficient +/- 1.96 standard errors), not the t intervals used for the Figure 5 pairwise displays.",
  "The adjusted genome–nucleus P is 0.97551; the adjusted genome–cell P is 0.0009971. These estimates explain the contrasting independence tests but are not estimates of indirect causal effects.",
  "The broad adjusted genome–nucleus interval does not establish an exactly zero effect.",
  "No models were refitted and these coefficients are not labelled on the three equivalent DAG orientations.")
legend <- paste(
  "Phylogenetic path comparisons favor the class connecting estimated genome size and nucleus area through cell area in the 24-species dataset.",
  "The winning cell-bridge equivalence class supports a mechanistic question about how the three traits covary, while retaining all three compatible causal orientations.",
  "(a)", legend_a, "(b)", legend_b, "(c)", legend_contrast, "(d)", legend_stability,
  "An accompanying simulation-calibration figure preserves all eight released recovery settings.",
  "All results are taken from the frozen path-analysis release; this figure introduces no model refitting or data filtering.",
  "Genome size is an IOD-derived D. fuscus assembly-anchored estimate, and cellular traits describe selected large erythrocytes and their nuclei. Measurement uncertainty and taxon coverage are evaluated by the displayed sensitivity analyses.")
save_candidate(6, "Path analysis of genome size, cell area and nucleus area", figure6,
  list(a = list(plot = p6a, title = "Alternative causal orientations in the cell-bridge class", width_mm = 180, height_mm = 90,
    legend = legend_a, alt_text = "Three equally compatible DAGs show genome to cell to nucleus, the reverse chain, and cell area pointing to both genome and nucleus."),
    b = list(plot = p6b, title = "Primary comparison of ten path classes", width_mm = 180, height_mm = 135,
      legend = legend_b, alt_text = "Ten models are ranked by delta CICc; only cell bridge passes global fit, with weight 0.917. Nucleus bridge is 13.77 units behind."),
    c = list(plot = p6contrast, title = "Marginal and conditional genome–cellular associations", width_mm = 180, height_mm = 100,
      legend = legend_contrast, alt_text = "Genome–nucleus coefficient decreases from 0.419 to 0.008 after adjusting for cell area, whereas genome–cell decreases from 0.717 to 0.575 after adjusting for nucleus area."),
    d = list(plot = p6stability, title = "Path-class ranking across sensitivity analyses", width_mm = 180, height_mm = 100,
      legend = legend_stability, alt_text = "Cell bridge ranks first in 246 of 250 resamples, all 24 omissions, all four evolutionary models, and all 200 trees. The 18-species tree fits also support nucleus bridge.")),
  width_mm = 180, height_mm = 270, datasets = "S28;S30;S32;S33;S36;S37;path/candidate_dags.csv",
  legend = legend,
  alt_text = "Four panels combine equivalent cell-bridge DAGs, all ten primary path-class rankings, marginal-versus-conditional genome associations and four sensitivity checks.")

save_gbe_figure(p6calibration, file.path(out, "validation", "Simulation_calibration"), 180, 110)
write_lines(c("# Simulation calibration of path-class recovery", "", legend_calibration),
  file.path(out, "validation", "SIMULATION_CALIBRATION_LEGEND.md"))

write_csv(winner_dags, file.path(analysis_dir, "cell_bridge_equivalent_dags.csv"))
write_csv(dag_nodes, file.path(analysis_dir, "dag_nodes.csv"))
write_csv(dag_edges, file.path(analysis_dir, "dag_edges.csv"))
write_csv(primary_plot, file.path(analysis_dir, "primary_class_comparison.csv"))
write_csv(stability_plot, file.path(analysis_dir, "cell_bridge_stability.csv"))
write_csv(tree_subset, file.path(analysis_dir, "reduced_taxon_comparison.csv"))
write_csv(calibration_plot, file.path(analysis_dir, "simulation_recovery.csv"))
write_csv(contrast_plot, file.path(analysis_dir, "marginal_conditional_coefficients.csv"))
write_csv(bridge_claim, file.path(analysis_dir, "cell_bridge_independence_test.csv"))
stopifnot(identical(input_hashes, hash_inputs()))
write_lines(input_hashes, file.path(analysis_dir, "input_checksums.sha256"))
write_lines(c(
  "PASS: only frozen inputs were read; all seven SHA256 digests are unchanged after rendering.",
  "PASS: 10 primary classes; 24 species; relative weights sum to 1; only cell bridge passes global fit.",
  "PASS: three DAG members exactly match the source cell-bridge equivalence class; no borrowed edge coefficients.",
  "PASS: cell-bridge d-separation P equals primary global-fit P with one independence claim.",
  "PASS: stability counts are 246/250, 24/24, 4/4, 200/200; all 200 tree fits also support nucleus bridge.",
  "PASS: focal and published-main 18-species fits both retain cell and nucleus bridges within delta CICc <= 2.",
  "PASS: all eight simulation settings retained, with 100 simulations each and recovery rates in [0,100].",
  "PASS: four S28 coefficient contrasts retained with unchanged 1.96-SE normal intervals and 24-species sample size.",
  "No analysis was refitted and no source data or original figures were edited."),
  file.path(RUN, "qa", "figure_6_numerical_checks.txt"))
message("Completed isolated Figure 6: one collage, four panels, one calibration support plot and nine numerical exports.")

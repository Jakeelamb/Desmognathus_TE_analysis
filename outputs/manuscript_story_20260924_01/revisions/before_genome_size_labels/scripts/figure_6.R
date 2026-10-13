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

# Each arrow uses the exact S28 regression whose predictor and outcome match
# that orientation. Every child here has one parent; no coefficients are averaged.
dag_order <- c("GS->CS;CS->NS", "NS->CS;CS->GS", "CS->GS;CS->NS")
node_positions <- tibble(node = c("GS", "CS", "NS"), x = c(1.45, 5, 8.55),
  label = c("Estimated\ngenome size", "Cell area", "Nucleus area"))
dag_rows <- tibble(dag_id = dag_order, y = c(3.6, 2.0, 0.55),
  header = c("Genome-first hypothesis", "Equivalent reverse-chain hypothesis", "Equivalent cell-first hypothesis"))
edge_lookup <- tibble(from = c("GS", "CS", "NS", "CS"),
  to = c("CS", "NS", "CS", "GS"),
  relation_id = c("genome_to_cell", "cell_to_nucleus", "nucleus_to_cell", "cell_to_genome"),
  expected_predictor = c("genome_size", "cell_size", "nucleus_size", "cell_size"),
  expected_outcome = c("cell_size", "nucleus_size", "cell_size", "genome_size")) %>%
  left_join(regressions, by = "relation_id")
stopifnot(nrow(edge_lookup) == 4L, all(edge_lookup$predictor == edge_lookup$expected_predictor),
  all(edge_lookup$outcome == edge_lookup$expected_outcome), all(edge_lookup$n_species == 24L),
  all(edge_lookup$evolutionary_model == "lambda"))
dag_nodes <- bind_rows(lapply(seq_along(dag_order), function(i) {
  node_positions %>% mutate(dag_id = dag_order[i], y = dag_rows$y[i])
}))
dag_edges <- bind_rows(lapply(seq_along(dag_order), function(i) {
  e <- strsplit(strsplit(dag_order[i], ";", fixed = TRUE)[[1]], "->", fixed = TRUE)
  bind_rows(lapply(e, function(pair) {
    x_start <- node_positions$x[match(pair[1], node_positions$node)]
    x_end <- node_positions$x[match(pair[2], node_positions$node)]
    direction <- sign(x_end - x_start)
    tibble(dag_id = dag_order[i], from = pair[1], to = pair[2],
      x = x_start + direction * 1.18, xend = x_end - direction * 1.18,
      y = dag_rows$y[i], yend = dag_rows$y[i])
  }))
})) %>% left_join(edge_lookup, by = c("from", "to")) %>%
  mutate(coefficient_label = sprintf("β = %.2f", standardized_coefficient),
    arrow_width_mm = 1.4 * abs(standardized_coefficient))
stopifnot(nrow(dag_edges) == 6L, !anyNA(dag_edges$standardized_coefficient),
  all(abs(dag_edges$ci_low - (dag_edges$standardized_coefficient - 1.96 * dag_edges$standard_error)) < 1e-12),
  all(abs(dag_edges$ci_high - (dag_edges$standardized_coefficient + 1.96 * dag_edges$standard_error)) < 1e-12))
p6a <- ggplot() +
  annotate("rect", xmin = 0.12, xmax = 9.88, ymin = 2.94, ymax = 4.37,
    fill = "#F1F7FB", colour = NA) +
  geom_rect(data = dag_nodes, aes(xmin = x - 1.10, xmax = x + 1.10,
    ymin = y - 0.32, ymax = y + 0.32, fill = node == "CS"),
    colour = "#9A9A9A", linewidth = 0.35) +
  geom_segment(data = dag_edges, aes(x = x, xend = xend, y = y, yend = yend, linewidth = arrow_width_mm),
    arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed"), colour = blue) +
  scale_linewidth_identity() +
  geom_text(data = dag_edges, aes(x = (x + xend) / 2, y = y + 0.32, label = coefficient_label),
    family = GBE_FONT, size = 8.5 / ggplot2::.pt, colour = blue, fontface = "bold") +
  geom_text(data = dag_rows, aes(x = 0.3, y = y + 0.65, label = header), hjust = 0,
    family = GBE_FONT, size = 7.5 / ggplot2::.pt, colour = gray) +
  geom_text(data = dag_nodes, aes(x = x, y = y, label = label),
    family = GBE_FONT, size = 8 / ggplot2::.pt, lineheight = 1.05) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#E5F1F8"), guide = "none") +
  scale_x_continuous(limits = c(0.1, 9.9), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.05, 4.5), expand = c(0, 0)) +
  labs(title = "Cell area connects genome-size and nuclear-size variation",
    subtitle = sprintf("Cell-bridge class: shared weight w = %.3f  |  24 species", winner$w),
    caption = "Numbers and arrow widths show standardized path coefficients. All three directions fit equally well.") +
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

path_intervals <- bind_rows(
  edge_lookup %>% filter(relation_id %in% c("genome_to_cell", "cell_to_nucleus")) %>%
    mutate(role = "Displayed path"),
  contrast_plot %>% filter(estimate_type == "Adjusted") %>% mutate(role = "Conditional check")) %>%
  mutate(order = match(relation_id, c("genome_to_cell", "cell_to_nucleus",
    "genome_and_cell_to_nucleus", "genome_and_nucleus_to_cell"))) %>% arrange(order) %>%
  mutate(y = 5 - order,
    row_label = c("Genome → cell", "Cell → nucleus", "Genome → nucleus\n(adjusted for cell area)",
      "Genome → cell\n(adjusted for nucleus area)"),
    interval_label = sprintf("%.2f [%.2f, %.2f]", standardized_coefficient, ci_low, ci_high))
stopifnot(nrow(path_intervals) == 4L, all(path_intervals$n_species == 24L),
  identical(path_intervals$relation_id, c("genome_to_cell", "cell_to_nucleus",
    "genome_and_cell_to_nucleus", "genome_and_nucleus_to_cell")))
p6intervals <- ggplot(path_intervals, aes(x = standardized_coefficient, y = y)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#A0A0A0", linewidth = 0.35) +
  geom_hline(yintercept = 2.5, colour = "#E0E0E0", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = role),
    orientation = "y", width = 0.13, linewidth = 0.6) +
  geom_point(aes(colour = role, shape = role), size = 2.5) +
  geom_text(aes(x = 1.18, label = interval_label), hjust = 0,
    family = GBE_FONT, size = 7.5 / ggplot2::.pt, colour = black) +
  annotate("text", x = 1.18, y = 4.65, label = "β [95% interval]", hjust = 0,
    family = GBE_FONT, size = 7.5 / ggplot2::.pt, fontface = "bold") +
  scale_colour_manual(values = c(`Displayed path` = blue, `Conditional check` = gray), guide = "none") +
  scale_shape_manual(values = c(`Displayed path` = 16, `Conditional check` = 17), guide = "none") +
  scale_x_continuous(limits = c(-0.6, 2.05), breaks = c(-0.5, 0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, 4.9), breaks = path_intervals$y,
    labels = path_intervals$row_label, expand = c(0, 0)) +
  labs(title = "Path strengths and the remaining genome–nucleus association",
    subtitle = "Upper rows: paths in the genome-first hypothesis. Lower rows: conditional comparisons.",
    x = "Standardized coefficient, β", y = NULL,
    caption = "Intervals: estimate ± 1.96 SE. The adjusted genome–nucleus estimate is near zero, with substantial uncertainty.") +
  panel_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5, lineheight = 1.05))

figure6 <- (p6a / p6b / p6intervals / p6stability) +
  plot_layout(heights = c(1.2, 1.25, 1.0, 0.8)) + plot_annotation(tag_levels = "a")
legend_a <- paste(
  "The cell-bridge class contains three distinct directed acyclic graphs, reproduced from the frozen candidate-DAG enumeration.",
  "Estimated genome size -> cell area -> nucleus area, its reverse chain, and cell area -> both traits all imply genome–nucleus conditional independence given cell area.",
  "Arrow labels are direction-specific standardized phylogenetic regression coefficients from S28; every child in these diagrams has one parent, so its displayed coefficient comes from the matching two-variable regression.",
  "Arrow width is proportional to the absolute coefficient. Traits are log10-transformed and standardized to unit sample standard deviation across 24 species.",
  "The genome-first orientation has genome -> cell beta = 0.71666 and cell -> nucleus beta = 0.57930. Reverse and fork coefficients are retrieved separately, rather than copied or averaged; they agree after rounding to two decimals.",
  "These orientations are statistically indistinguishable in this three-trait comparison and share one class weight of 0.91674; their common conditional structure is the supported result. The highlighted orientation is a hypothesis, not a uniquely identified causal direction.")
legend_b <- paste(
  "Ten testable path-equivalence classes are ranked by CICc using 24 species, standardized log10 traits, the supplied dated tree and Pagel-lambda regressions.",
  "The cell-bridge class has CICc = 13.3829, Fisher C = 0.04959, one independence claim, global-fit P = 0.97551, and relative weight = 0.91674.",
  "Its tested claim is the residual genome–nucleus association after conditioning on cell area (P = 0.97551).",
  "All nine other tested classes fail the P >= 0.05 global-fit criterion; the nucleus-bridge class has delta CICc = 13.77171 and global-fit P = 0.0009971.",
  "Weights compare only the supplied ten classes; six saturated DAGs impose no testable d-separation claim and are not ranked.",
  "All primary component regressions have fitted lambda approximately 1e-7. The class weight measures relative model support and is distinct from the direction-specific path coefficients in panel a.")
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
  "No models were refitted; these conditional comparisons are distinct from the single-parent paths shown in the equivalent DAG orientations.")
legend_intervals <- paste(
  "Upper rows show the two component regressions in the genome-first hypothesis: genome -> cell beta = 0.71666 (95% normal interval 0.42523 to 1.00810), and cell -> nucleus beta = 0.57930 (0.23869 to 0.91992).",
  "Lower rows show separate conditional regressions: genome -> nucleus adjusted for cell area, beta = 0.00792 (-0.49195 to 0.50780); genome -> cell adjusted for nucleus area, beta = 0.57486 (0.27994 to 0.86977).",
  "Circles identify displayed single-parent paths and triangles identify conditional comparisons. All four fits use 24 species, standardized log10 traits and Pagel-lambda phylogenetic regression.",
  "Intervals are normal approximations from S28 (beta +/- 1.96 SE), not bootstrap intervals or estimates of indirect effects. The broad conditional genome–nucleus interval does not establish an exactly zero effect.")
legend <- paste(
  "Phylogenetic path comparisons favor the class connecting estimated genome size and nucleus area through cell area in the 24-species dataset.",
  "The winning cell-bridge equivalence class supports a mechanistic question about how the three traits covary, while retaining all three compatible causal orientations.",
  "(a)", legend_a, "(b)", legend_b, "(c)", legend_intervals, "(d)", legend_stability,
  "Accompanying figures retain all eight simulation-recovery settings and the original marginal-versus-conditional coefficient comparison.",
  "All results are taken from the frozen path-analysis release; this figure introduces no model refitting or data filtering.",
  "Genome size is an IOD-derived D. fuscus assembly-anchored estimate, and cellular traits describe selected large erythrocytes and their nuclei. Measurement uncertainty and taxon coverage are evaluated by the displayed sensitivity analyses.")
save_candidate(6, "Path analysis of genome size, cell area and nucleus area", figure6,
  list(a = list(plot = p6a, title = "Weighted causal hypotheses linking genome, cell and nucleus", width_mm = 180, height_mm = 100,
    legend = legend_a, alt_text = "Three equally compatible causal hypotheses have shared class weight 0.917. Direction-specific genome–cell arrows have coefficient 0.72 and cell–nucleus arrows 0.58."),
    b = list(plot = p6b, title = "Primary comparison of ten path classes", width_mm = 180, height_mm = 135,
      legend = legend_b, alt_text = "Ten models are ranked by delta CICc; only cell bridge passes global fit, with weight 0.917. Nucleus bridge is 13.77 units behind."),
    c = list(plot = p6intervals, title = "Path coefficients and uncertainty intervals", width_mm = 180, height_mm = 100,
      legend = legend_intervals, alt_text = "Genome-to-cell and cell-to-nucleus path coefficients are 0.72 and 0.58. Conditional genome-to-nucleus is 0.008 with a wide interval; conditional genome-to-cell remains positive."),
    d = list(plot = p6stability, title = "Path-class ranking across sensitivity analyses", width_mm = 180, height_mm = 100,
      legend = legend_stability, alt_text = "Cell bridge ranks first in 246 of 250 resamples, all 24 omissions, all four evolutionary models, and all 200 trees. The 18-species tree fits also support nucleus bridge.")),
  width_mm = 180, height_mm = 270, datasets = "S28;S30;S32;S33;S36;S37;path/candidate_dags.csv",
  legend = legend,
  alt_text = "Four panels show weighted alternative causal hypotheses, model-class support, coefficient uncertainty, and sensitivity checks for the genome–cell–nucleus relationship.")

save_gbe_figure(p6contrast, file.path(out, "validation", "Marginal_conditional_comparison"), 180, 100)
write_lines(c("# Marginal and conditional associations", "", legend_contrast),
  file.path(out, "validation", "MARGINAL_CONDITIONAL_LEGEND.md"))

save_gbe_figure(p6calibration, file.path(out, "validation", "Simulation_calibration"), 180, 110)
write_lines(c("# Simulation calibration of path-class recovery", "", legend_calibration),
  file.path(out, "validation", "SIMULATION_CALIBRATION_LEGEND.md"))

write_csv(winner_dags, file.path(analysis_dir, "cell_bridge_equivalent_dags.csv"))
write_csv(dag_nodes, file.path(analysis_dir, "dag_nodes.csv"))
write_csv(dag_edges, file.path(analysis_dir, "dag_edges.csv"))
write_csv(path_intervals, file.path(analysis_dir, "path_coefficient_intervals.csv"))
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
  "PASS: three DAG members exactly match the source cell-bridge equivalence class; each of six arrows uses its matching direction-specific single-parent S28 regression.",
  "PASS: all six arrow widths are proportional to absolute beta; intervals match beta +/- 1.96 SE. Shared class weight is distinct from coefficients.",
  "PASS: cell-bridge d-separation P equals primary global-fit P with one independence claim.",
  "PASS: stability counts are 246/250, 24/24, 4/4, 200/200; all 200 tree fits also support nucleus bridge.",
  "PASS: focal and published-main 18-species fits both retain cell and nucleus bridges within delta CICc <= 2.",
  "PASS: all eight simulation settings retained, with 100 simulations each and recovery rates in [0,100].",
  "PASS: four S28 coefficient contrasts retained with unchanged 1.96-SE normal intervals and 24-species sample size.",
  "No analysis was refitted and no source data or original figures were edited."),
  file.path(RUN, "qa", "figure_6_numerical_checks.txt"))
message("Completed isolated Figure 6: one collage, four panels, one calibration support plot and nine numerical exports.")

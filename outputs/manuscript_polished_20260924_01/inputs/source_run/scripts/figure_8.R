source("scripts/common.R")

# Rendering only: the separate body-size analysis writes all fitted estimates.
out <- file.path(RUN, "figures", "figure_8")
analysis_dir <- file.path(out, "analysis")
display_dir <- file.path(out, "display")
dir.create(display_dir, recursive = TRUE, showWarnings = FALSE)

blue <- GBE_COLORS[["blue"]]
gray <- "#777777"
black <- GBE_COLORS[["charcoal"]]
red <- GBE_COLORS[["vermillion"]]
panel_style <- theme_gbe() + theme(
  plot.title = element_text(size = 10, face = "bold", margin = margin(b = 4)),
  plot.subtitle = element_text(size = 8, lineheight = 1.12, margin = margin(b = 6)),
  plot.caption = element_text(size = 7.5, lineheight = 1.12, hjust = 0, margin = margin(t = 7)),
  plot.margin = margin(6, 9, 6, 9),
  axis.text.y = element_text(size = 8, lineheight = 1.08),
  strip.text = element_text(size = 8, face = "bold"),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8))

coefficient_label <- function(x) ifelse(abs(x) < 0.05, sprintf("%.3f", x), sprintf("%.2f", x))

# The SVL node is drawn above all three traits. Direction is a stated hypothesis;
# coefficient magnitude and interval status have separate visual encodings.
nodes <- tibble(node = c("b", "g", "c", "n"),
  x = c(5, 1.4, 5, 8.6), y = c(5.05, 1.65, 1.65, 1.65),
  label = c("Body size (SVL)", "Genome size", "Cell area", "Nucleus area"))
edge_geometry <- tribble(
  ~from, ~to, ~x, ~y, ~xend, ~yend, ~label_x, ~label_y,
  "b", "g", 3.90, 4.57, 1.90, 2.13, 1.95, 3.60,
  "b", "c", 5.00, 4.57, 5.00, 2.13, 6.14, 3.50,
  "b", "n", 6.10, 4.57, 8.10, 2.13, 8.05, 3.60,
  "g", "c", 2.70, 1.65, 3.70, 1.65, 3.20, 0.74,
  "c", "n", 6.30, 1.65, 7.30, 1.65, 6.80, 0.74,
  "n", "c", 7.30, 1.65, 6.30, 1.65, 6.80, 0.74,
  "g", "n", 1.40, 1.17, 8.60, 1.17, 5.00, -0.49)

make_network <- function(edges, title, subtitle) {
  stopifnot(nrow(edges) %in% c(3L, 5L), !anyDuplicated(paste(edges$from, edges$to)))
  edges <- edges %>% left_join(edge_geometry, by = c("from", "to")) %>%
    mutate(interval_includes_zero = ci_low <= 0 & ci_high >= 0,
      edge_colour = ifelse(interval_includes_zero, gray, ifelse(estimate < 0, red, black)),
      edge_style = ifelse(interval_includes_zero, "dashed", "solid"),
      arrow_width_mm = pmax(0.25, 1.7 * abs(estimate)),
      beta_label = coefficient_label(estimate),
      interval_label = sprintf("[%.2f, %.2f]", ci_low, ci_high),
      long_edge = from == "g" & to == "n")
  stopifnot(!anyNA(edges$x), !anyNA(edges$estimate))
  short <- edges %>% filter(!long_edge)
  long <- edges %>% filter(long_edge)
  p <- ggplot() +
    geom_segment(data = short, aes(x = x, y = y, xend = xend, yend = yend,
      colour = edge_colour, linetype = edge_style, linewidth = arrow_width_mm),
      arrow = grid::arrow(length = grid::unit(2.8, "mm"), type = "closed"), lineend = "butt")
  if (nrow(long)) {
    # Route the non-adjacent genome–nucleus edge below the middle node, with the
    # arrowhead only at its nucleus endpoint; no crossing or hidden node occurs.
    p <- p +
      geom_segment(data = long, aes(x = x, y = y, xend = x, yend = -0.08,
        colour = edge_colour, linetype = edge_style, linewidth = arrow_width_mm)) +
      geom_segment(data = long, aes(x = x, y = -0.08, xend = xend, yend = -0.08,
        colour = edge_colour, linetype = edge_style, linewidth = arrow_width_mm)) +
      geom_segment(data = long, aes(x = xend, y = -0.08, xend = xend, yend = yend,
        colour = edge_colour, linetype = edge_style, linewidth = arrow_width_mm),
        arrow = grid::arrow(length = grid::unit(2.8, "mm"), type = "closed"))
  }
  p <- p + scale_colour_identity() + scale_linewidth_identity() + scale_linetype_identity() +
    geom_rect(data = nodes, aes(xmin = x - 1.30, xmax = x + 1.30,
      ymin = y - 0.48, ymax = y + 0.48), fill = "white", colour = black, linewidth = 0.5) +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
      family = GBE_FONT, size = 9.5 / ggplot2::.pt) +
    geom_label(data = edges, aes(x = label_x, y = label_y + 0.15, label = beta_label,
      colour = edge_colour), family = GBE_FONT, size = 9.5 / ggplot2::.pt,
      fontface = "bold", fill = "white", linewidth = 0, label.padding = grid::unit(0.4, "mm")) +
    geom_label(data = edges, aes(x = label_x, y = label_y - 0.23, label = interval_label),
      family = GBE_FONT, size = 7.5 / ggplot2::.pt, colour = gray,
      fill = "white", linewidth = 0, label.padding = grid::unit(0.3, "mm")) +
    scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ifelse(nrow(long), -0.88, 0.10), 5.75), expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle,
      caption = "Full-panel standardized coefficients [95% t intervals]; widths scale with |coefficient|.\nGray dashed arrows: interval includes zero. Arrow directions are specified hypotheses.") +
    panel_style + theme(axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
      axis.title = element_blank())
  list(plot = p, displayed_edges = edges)
}

model_labels <- c(body_only = "SVL alone", cell_bridge = "Genome → cell\n→ nucleus",
  nucleus_bridge = "Genome → nucleus\n→ cell", genome_bridge = "Cell ← genome\n→ nucleus")
cohort_labels <- c(full24 = "All species (n = 24)", without_small_direct = "Without aeneus + wrighti (n = 22)")

make_comparison <- function(comparison) {
  model_order <- comparison %>% filter(cohort == "full24") %>% arrange(CICc) %>% pull(model)
  comparison <- comparison %>% mutate(y = 5 - match(model, model_order),
    cohort_label = factor(cohort_labels[cohort], levels = cohort_labels),
    fit_status = ifelse(p >= 0.05, "P ≥ 0.05", "P < 0.05"),
    highlight = ifelse(rank == 1, "First ranked", "Other candidate"),
    weight_label = ifelse(w < 0.001, "<0.001", sprintf("%.3f", w)))
  xmax <- max(3, comparison$delta_CICc)
  weight_x <- xmax * 1.35
  p <- ggplot(comparison, aes(x = delta_CICc, y = y)) +
    geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.35, colour = "#9B9B9B") +
    geom_segment(aes(x = 0, xend = delta_CICc, yend = y, colour = highlight), linewidth = 0.45) +
    geom_point(aes(colour = highlight, shape = fit_status), size = 2.5, stroke = 0.7) +
    geom_text(aes(x = weight_x, label = weight_label), hjust = 1,
      family = GBE_FONT, size = 7.5 / ggplot2::.pt) +
    annotate("text", x = weight_x, y = 4.8, label = "Weight", hjust = 1,
      family = GBE_FONT, size = 7.5 / ggplot2::.pt, fontface = "bold") +
    scale_colour_manual(values = c(`First ranked` = blue, `Other candidate` = gray), guide = "none") +
    scale_shape_manual(name = "Global-fit test", values = c(`P ≥ 0.05` = 16, `P < 0.05` = 4)) +
    scale_x_continuous(limits = c(-0.02 * xmax, xmax * 1.40),
      breaks = pretty(c(0, xmax), n = 4), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 5), breaks = 4:1,
      labels = unname(model_labels[model_order]), expand = c(0, 0)) +
    facet_wrap(~cohort_label, nrow = 1) +
    labs(title = "Support among four SVL-adjusted hypotheses",
      subtitle = "Every candidate includes SVL → genome size, cell area and nucleus area.",
      x = "ΔCICc within each species panel", y = NULL,
      caption = "Dashed line: ΔCICc = 2. Weights compare only these four candidates within each panel.\nExclusion is a sensitivity analysis; body size and development mode differ together.") +
    panel_style + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      panel.spacing.x = grid::unit(7, "mm"), legend.position = "bottom")
  list(plot = p, displayed_comparison = comparison)
}

contrast_labels <- c(genome_cell = "Genome → cell", cell_nucleus = "Cell → nucleus",
  direct_genome_nucleus = "Genome → nucleus\n(conditioned on cell)")
make_contrasts <- function(contrasts) {
  contrasts <- contrasts %>% mutate(base_y = 4 - match(relation, names(contrast_labels)),
    y = base_y + ifelse(adjustment == "Without SVL", 0.14, -0.14),
    adjustment = factor(adjustment, levels = c("Without SVL", "With SVL")),
    cohort_label = factor(cohort_labels[cohort], levels = cohort_labels),
    interval_includes_zero = ci_low <= 0 & ci_high >= 0)
  p <- ggplot(contrasts, aes(x = estimate, y = y, colour = adjustment)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#9B9B9B") +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.11, linewidth = 0.5) +
    geom_point(aes(shape = adjustment), size = 2.5, stroke = 0.7) +
    scale_colour_manual(values = c(`Without SVL` = gray, `With SVL` = blue), name = NULL) +
    scale_shape_manual(values = c(`Without SVL` = 16, `With SVL` = 17), name = NULL) +
    scale_y_continuous(limits = c(0.5, 3.5), breaks = 3:1,
      labels = unname(contrast_labels), expand = c(0, 0)) +
    scale_x_continuous(breaks = scales::breaks_pretty(5), expand = expansion(mult = 0.04)) +
    facet_wrap(~cohort_label, nrow = 1) +
    labs(title = "Do cellular associations persist after SVL adjustment?",
      subtitle = "Same species and response variables in each within-panel comparison.",
      x = "Full-panel standardized coefficient (95% t interval)", y = NULL,
      caption = "The genome–nucleus comparison also conditions on cell area. All traits are log10 transformed.\nThe 22-species coefficients retain the means and SDs of the full 24-species panel.") +
    panel_style + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      panel.spacing.x = grid::unit(7, "mm"), legend.position = "bottom")
  list(plot = p, displayed_contrasts = contrasts)
}

input_files <- file.path(analysis_dir, c("path_model_comparison.csv", "path_edge_estimates.csv",
  "regression_coefficients.csv", "candidate_edges.csv"))
stopifnot(all(file.exists(input_files)))
hash_inputs <- function() system2("sha256sum", shQuote(input_files), stdout = TRUE)
input_hashes <- hash_inputs()
source_comparison <- read_csv(input_files[1], show_col_types = FALSE)
source_edges <- read_csv(input_files[2], show_col_types = FALSE)
source_regressions <- read_csv(input_files[3], show_col_types = FALSE)
candidate_edges <- read_csv(input_files[4], show_col_types = FALSE)
cohort_key <- c(full24 = "full24", subset22 = "without_small_direct")

comparison <- source_comparison %>% filter(candidate_set == "body_adjusted_four") %>%
  mutate(cohort = unname(cohort_key[panel]))
stopifnot(nrow(comparison) == 8L, !anyNA(comparison$cohort),
  setequal(comparison$model, names(model_labels)),
  all(comparison$n_species == ifelse(comparison$panel == "full24", 24L, 22L)),
  all(abs(comparison %>% group_by(panel) %>% summarise(total = sum(w)) %>% pull(total) - 1) < 1e-12))
winner <- comparison %>% filter(panel == "full24", rank == 1)
stopifnot(nrow(winner) == 1L)
network_edges <- source_edges %>%
  filter(candidate_set == "body_adjusted_four", panel == "full24", model == winner$model) %>%
  rename(from = parent, to = child, estimate = standardized_coefficient)
expected_edges <- candidate_edges %>% filter(candidate_set == "body_adjusted_four", model == winner$model)
stopifnot(nrow(network_edges) == nrow(expected_edges),
  setequal(paste(network_edges$from, network_edges$to), paste(expected_edges$parent, expected_edges$child)),
  all(network_edges$n_species == 24L),
  all(abs(network_edges$ci_low - (network_edges$estimate - qt(.975, network_edges$residual_df) * network_edges$standard_error)) < 1e-12),
  all(abs(network_edges$ci_high - (network_edges$estimate + qt(.975, network_edges$residual_df) * network_edges$standard_error)) < 1e-12))
network_names <- c(body_only = "SVL-only hypothesis", cell_bridge = "SVL-adjusted cell-bridge hypothesis",
  nucleus_bridge = "SVL-adjusted nucleus-bridge hypothesis", genome_bridge = "SVL-adjusted genome-bridge hypothesis")
network <- make_network(network_edges, network_names[[winner$model]],
  sprintf("Highest-ranked of four candidates; n = 24; weight = %.3f; global-fit P = %.3f", winner$w, winner$p))
comparison_plot <- make_comparison(comparison)

contrast_specification <- tribble(
  ~relation_id, ~term, ~relation, ~adjustment,
  "genome_to_cell", "g", "genome_cell", "Without SVL",
  "genome_to_cell_body_adjusted", "g", "genome_cell", "With SVL",
  "cell_to_nucleus", "c", "cell_nucleus", "Without SVL",
  "cell_to_nucleus_body_adjusted", "c", "cell_nucleus", "With SVL",
  "genome_and_cell_to_nucleus", "g", "direct_genome_nucleus", "Without SVL",
  "genome_and_cell_to_nucleus_body_adjusted", "g", "direct_genome_nucleus", "With SVL")
contrasts <- source_regressions %>% inner_join(contrast_specification, by = c("relation_id", "term")) %>%
  rename(estimate = standardized_coefficient) %>% mutate(cohort = unname(cohort_key[panel]))
stopifnot(nrow(contrasts) == 12L, !anyNA(contrasts$cohort),
  !anyDuplicated(paste(contrasts$panel, contrasts$relation_id, contrasts$term)),
  all(abs(contrasts$ci_low - (contrasts$estimate - qt(.975, contrasts$residual_df) * contrasts$standard_error)) < 1e-12),
  all(abs(contrasts$ci_high - (contrasts$estimate + qt(.975, contrasts$residual_df) * contrasts$standard_error)) < 1e-12),
  all(contrasts$n_species == ifelse(contrasts$panel == "full24", 24L, 22L)))
contrast_plot <- make_contrasts(contrasts)

methods <- paste(
  "All regressions use Pagel-lambda phylogenetic regression on log10 traits standardized using the full 24-species means and sample SDs.",
  "The 22-species panel retains the same scaling and excludes only D. aeneus and D. wrighti, which are both the smallest species and the only sampled direct developers.",
  "This exclusion is a joint influence/developmental sensitivity analysis, not an independently replicated small-bodied group or an outlier filter.",
  "The focal supplied dated phylogeny is pruned to the respective species panel.",
  "Intervals are coefficient +/- qt(0.975, n - p) times SE, where p counts all regression coefficients including the intercept; they are conditional approximations, not bootstrap intervals integrating lambda uncertainty.",
  "SVL is the collaborator-provided species maximum snout-vent length, not the length of the imaged individual.",
  "Genome size is an IOD-derived D. fuscus assembly-anchored estimate; cell and nucleus areas are species medians among retained upper-tail paired erythrocyte objects.")
legend_a <- paste(
  sprintf("The displayed network is the highest-ranked full24 candidate, %s, with CICc %.5f, relative weight %.5f and global-fit P %.5f.",
    winner$model, winner$CICc, winner$w, winner$p),
  "All its structural paths are retained, including paths with broad or zero-crossing intervals. Labels give standardized coefficients and 95% t intervals.",
  "Arrow width equals max(0.25 mm, 1.7 times absolute coefficient). Gray dashed arrows indicate intervals including zero; other positive paths are black and negative paths vermillion.",
  "SVL is specified as a common predictor. The figure does not identify the direction of causation, demonstrate that SVL is exogenous, or model-average effects.",
  "The ranking concerns a limited set of four hypotheses and their associated independence claims; equivalent reorientations cannot be separated by these observations.", methods)
legend_b <- paste(
  "Four candidates are compared separately in the full24 and subset22 panels. All include body size -> genome size, cell area and nucleus area.",
  "SVL alone adds no genome/cell/nucleus links; cell bridge adds genome -> cell -> nucleus; nucleus bridge adds genome -> nucleus -> cell; genome bridge adds genome -> cell and genome -> nucleus.",
  "The plotted delta CICc is relative to the best candidate in the same panel and candidate set. Relative weights sum to one separately within each panel; they are not causal probabilities.",
  "Blue identifies the first-ranked candidate; a circle denotes a global-fit P >= 0.05 and a cross P < 0.05. The dashed threshold is delta CICc = 2.",
  "This is a finite, non-exhaustive candidate comparison. Saturated four-trait graphs have no testable d-separation claim and are not ranked.",
  "Absolute CICc values and weights should not be compared across species panels or to the earlier ten-class three-trait set.", methods)
legend_c <- paste(
  "Within each species panel, paired coefficients compare the same response and focal predictor with and without SVL adjustment.",
  "The first two contrasts are cell ~ genome versus cell ~ genome + SVL, and nucleus ~ cell versus nucleus ~ cell + SVL.",
  "The third is nucleus ~ genome + cell versus nucleus ~ genome + cell + SVL, displaying the direct genome-nucleus association after conditioning on cell area.",
  "The fully connected last regression is a separate coefficient diagnostic; it is not assigned the reduced cell-bridge candidate's model weight.",
  "Gray circles omit SVL and blue triangles include it; horizontal bars are 95% t intervals and the vertical dashed line is zero.", methods)
legend <- paste("Body size is incorporated into exploratory phylogenetic path comparisons and targeted species-exclusion sensitivity.",
  "(a)", legend_a, "(b)", legend_b, "(c)", legend_c)

figure8 <- network$plot / comparison_plot$plot / contrast_plot$plot +
  plot_layout(heights = c(1.16, 1, 1)) + plot_annotation(tag_levels = "a")
save_candidate(8, "Body size in the genome–cell–nucleus path analysis", figure8,
  list(a = list(plot = network$plot, title = "Weighted four-trait hypothesis with body size",
      width_mm = 180, height_mm = 110, legend = legend_a,
      alt_text = "A four-trait path diagram includes maximum snout-vent length, genome size, cell area and nucleus area. All displayed path estimates have their 95% intervals; dashed arrows have intervals crossing zero."),
    b = list(plot = comparison_plot$plot, title = "Four body-size hypotheses across two species panels",
      width_mm = 180, height_mm = 110, legend = legend_b,
      alt_text = "Four body-size-adjusted candidate structures are ranked separately among all 24 species and among 22 species excluding aeneus and wrighti. Each candidate has a delta CICc, a relative weight and a global-fit marker."),
    c = list(plot = contrast_plot$plot, title = "Cellular path coefficients before and after SVL adjustment",
      width_mm = 180, height_mm = 110, legend = legend_c,
      alt_text = "Paired phylogenetic coefficients compare genome-to-cell, cell-to-nucleus, and genome-to-nucleus conditioned on cell area, before and after SVL adjustment, in both species panels.")),
  width_mm = 180, height_mm = 285,
  datasets = "body_size/traits24.csv;body_size/organismal_traits.csv;body_size/path24_time_tree.nwk;figure_8/analysis",
  legend = legend,
  alt_text = "Three panels integrate body size into a weighted four-trait network, a comparison of four candidate hypotheses, and matched coefficient contrasts across full24 and targeted subset22 species panels.")

write_csv(network$displayed_edges, file.path(display_dir, "network_edges.csv"))
write_csv(nodes, file.path(display_dir, "network_nodes.csv"))
write_csv(comparison_plot$displayed_comparison, file.path(display_dir, "candidate_comparison.csv"))
write_csv(contrast_plot$displayed_contrasts, file.path(display_dir, "coefficient_contrasts.csv"))
write_csv(contrast_specification, file.path(display_dir, "contrast_specification.csv"))
stopifnot(identical(input_hashes, hash_inputs()))
write_lines(input_hashes, file.path(display_dir, "input_checksums.sha256"))
write_lines(c(
  "PASS: all four input SHA256 digests unchanged after rendering; no model refitted by this script.",
  "PASS: eight candidate rows identify exactly four candidates in each of two species panels; weights sum to one separately.",
  "PASS: the displayed network is selected from the first-ranked full24 candidate, retaining every matching structural edge.",
  "PASS: all displayed network and contrast intervals match coefficient +/- qt(.975,residual_df)*SE.",
  "PASS: 12 contrast rows contain paired SVL adjustments for three focal associations in the full24 and subset22 panels.",
  "PASS: all source-matched sample sizes are 24 or 22 as specified; original figures, inputs and notes are untouched."),
  file.path(display_dir, "numerical_checks.txt"))
message("Completed Figure 8: one collage, three isolated panels, source-matched network and coefficient displays.")

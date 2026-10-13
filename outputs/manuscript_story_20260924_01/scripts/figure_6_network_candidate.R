source("scripts/common.R")

# Separate coefficient display under a specified genome-first orientation.
# The fully connected DAG has no testable d-separation claim. Do not attach
# the reduced cell-bridge class's CICc weight to this different diagram.
out <- file.path(RUN, "figures", "figure_6", "alternatives", "direct_indirect")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source_path <- file.path(DATA, "Supplementary_Data_S28_pairwise_phylogenetic_regressions.csv")
before <- system2("sha256sum", shQuote(source_path), stdout = TRUE)
source_data <- read_csv(source_path, show_col_types = FALSE)
edges <- tibble(
  from = c("GS", "CS", "GS"), to = c("CS", "NS", "NS"),
  relation_id = c("genome_to_cell", "genome_and_cell_to_nucleus", "genome_and_cell_to_nucleus"),
  predictor = c("genome_size", "cell_size", "genome_size"),
  expected_outcome = c("cell_size", "nucleus_size", "nucleus_size"),
  x = c(3.0, 7.7, 1.6), y = c(4.25, 3.72, 3.72),
  xend = c(6.3, 7.7, 6.3), yend = c(4.25, 1.23, 0.7),
  label_x = c(4.65, 8.25, 3.6), label_y = c(4.64, 2.5, 2.0)) %>%
  left_join(source_data, by = c("relation_id", "predictor")) %>%
  mutate(coefficient_label = ifelse(abs(standardized_coefficient) < 0.05,
      sprintf("%.3f", standardized_coefficient), sprintf("%.2f", standardized_coefficient)),
    interval_includes_zero = ci_low <= 0 & ci_high >= 0,
    edge_colour = ifelse(interval_includes_zero, "#6F6F6F", "#222222"),
    edge_style = ifelse(interval_includes_zero, "dashed", "solid"),
    arrow_width_mm = pmax(0.25, 1.8 * abs(standardized_coefficient)),
    row_label = c("Genome → cell", "Cell → nucleus\n(adjusted for genome size)",
      "Genome → nucleus\n(adjusted for cell area)"),
    interval_label = sprintf("%s [%.2f, %.2f]", coefficient_label, ci_low, ci_high),
    row = 3:1)
stopifnot(nrow(edges) == 3L, !anyNA(edges$standardized_coefficient),
  all(edges$outcome == edges$expected_outcome), all(edges$n_species == 24L),
  all(edges$evolutionary_model == "lambda"),
  all(abs(edges$ci_low - (edges$standardized_coefficient - 1.96 * edges$standard_error)) < 1e-12),
  all(abs(edges$ci_high - (edges$standardized_coefficient + 1.96 * edges$standard_error)) < 1e-12),
  identical(edges$interval_includes_zero, c(FALSE, FALSE, TRUE)))
nodes <- tibble(node = c("GS", "CS", "NS"), x = c(1.6, 7.7, 7.7), y = c(4.25, 4.25, 0.7),
  label = c("Genome size", "Cell area", "Nucleus area"))
style <- theme_gbe() + theme(
  plot.title = element_text(size = 10, face = "bold", hjust = 0, margin = margin(b = 5)),
  plot.subtitle = element_text(size = 8, hjust = 0, margin = margin(b = 7)),
  plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.1, margin = margin(t = 7)),
  plot.margin = margin(8, 10, 8, 10))
network <- ggplot() +
  geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend,
    colour = edge_colour, linetype = edge_style, linewidth = arrow_width_mm),
    arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"), lineend = "butt") +
  scale_colour_identity() + scale_linewidth_identity() + scale_linetype_identity() +
  geom_rect(data = nodes, aes(xmin = x - 1.3, xmax = x + 1.3, ymin = y - 0.48, ymax = y + 0.48),
    fill = "white", colour = "#222222", linewidth = 0.55) +
  geom_text(data = nodes, aes(x = x, y = y, label = label),
    family = GBE_FONT, size = 10 / ggplot2::.pt, lineheight = 1.1) +
  geom_label(data = edges, aes(x = label_x, y = label_y, label = coefficient_label,
    colour = edge_colour), family = GBE_FONT, size = 10 / ggplot2::.pt,
    fontface = "bold", fill = "white", linewidth = 0, label.padding = grid::unit(1, "mm")) +
  scale_x_continuous(limits = c(0, 9.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 5.1), expand = c(0, 0)) +
  labs(title = "Direct and indirect paths among genome, cell and nucleus",
    subtitle = "Specified genome-first hypothesis · standardized coefficients · 24 species",
    caption = "Arrow width scales with coefficient magnitude, with a minimum visible width.\nDashed gray path: its 95% interval includes zero. Directions are specified by the hypothesis.") +
  style + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), axis.title = element_blank())

intervals <- ggplot(edges, aes(x = standardized_coefficient, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#B0B0B0", linewidth = 0.35) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, colour = edge_colour),
    orientation = "y", width = 0.12, linewidth = 0.6) +
  geom_point(aes(colour = edge_colour, shape = interval_includes_zero), size = 2.5, stroke = 0.7) +
  geom_text(aes(x = 1.23, label = interval_label), hjust = 0,
    family = GBE_FONT, size = 8 / ggplot2::.pt) +
  annotate("text", x = 1.23, y = 3.62, label = "β [95% interval]", hjust = 0,
    family = GBE_FONT, size = 8 / ggplot2::.pt, fontface = "bold") +
  scale_colour_identity() + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(-0.6, 2.3), breaks = c(-0.5, 0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 3.9), breaks = edges$row, labels = edges$row_label,
    expand = c(0, 0)) +
  labs(title = "Uncertainty in the three displayed paths", x = "Standardized coefficient, β", y = NULL,
    caption = "Nucleus area is regressed on genome size and cell area together.\nIntervals are normal approximations (estimate ± 1.96 SE).") +
  style + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8, lineheight = 1.05))

legend <- paste(
  "Reference-style candidate depicting direct and indirect paths under a specified genome-first causal hypothesis: genome -> cell, cell -> nucleus and genome -> nucleus.",
  "Component equations are cell ~ genome and nucleus ~ genome + cell, fitted with Pagel-lambda phylogenetic regressions on log10-transformed, sample-standardized traits from 24 species.",
  "Coefficients and normal 95% intervals come unchanged from frozen S28: genome -> cell beta = 0.716664 (0.425232 to 1.008097), cell -> nucleus conditional on genome beta = 0.573623 (0.073748 to 1.073498), and genome -> nucleus conditional on cell beta = 0.007923 (-0.491952 to 0.507798).",
  "The cell-to-nucleus coefficient differs from the previous two-edge chain because this model additionally adjusts nucleus area for genome size.",
  "Arrow widths equal max(0.25 mm, 1.8 times absolute beta); gray dashed arrows and open interval points identify coefficients whose normal interval includes zero. All three point estimates are positive, so no negative-effect color is needed.",
  "The weak direct path is retained; its broad interval does not establish an exactly zero effect. Arrow directions are specified, not identified by this figure.",
  "This fully connected three-variable DAG is saturated: it has no d-separation claim and was not ranked by the current CICc workflow. The cell-bridge class weight of 0.917 is not assigned to this diagram.",
  "No model-averaging, refitting, new filtering, indirect-effect significance test or new measurement was performed.",
  "Genome size is an IOD-derived D. fuscus assembly-anchored estimate; cell and nucleus traits are areas of selected erythrocytes.")
alt <- "A triangular genome-first hypothesis retains all three paths: genome to cell 0.72, cell to nucleus 0.57, and genome to nucleus 0.008. The weak direct path is dashed; intervals appear below."
collage <- network / intervals + plot_layout(heights = c(1.35, 1)) + plot_annotation(tag_levels = "a")
save_gbe_figure(collage, file.path(out, "Direct_indirect_candidate"), 180, 180, dpi = 600)
save_gbe_figure(network + labs(tag = "a") + theme(plot.tag.position = "topleft"),
  file.path(out, "Direct_indirect_network"), 180, 110, dpi = 600)
save_gbe_figure(intervals + labs(tag = "b") + theme(plot.tag.position = "topleft"),
  file.path(out, "Direct_indirect_intervals"), 180, 95, dpi = 600)
write_csv(edges, file.path(out, "displayed_path_estimates.csv"))
write_csv(nodes, file.path(out, "node_positions.csv"))
write_lines(c("# Direct and indirect path candidate", "", legend, "", paste0("Alt text: ", alt)),
  file.path(out, "LEGEND_AND_ALT_TEXT.md"))
write_lines(before, file.path(out, "input_checksum.sha256"))
stopifnot(identical(before, system2("sha256sum", shQuote(source_path), stdout = TRUE)))
message("PASS: three source-matched coefficients, unchanged normal intervals, no class weight assigned; candidate collage and two panels rendered.")

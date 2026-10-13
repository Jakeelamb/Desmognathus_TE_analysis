source("scripts/common.R")

# Display only: read all fitted values from the frozen source run.
out <- file.path(RUN, "figures", "figure_6")
dir.create(file.path(out, "display"), recursive = TRUE, showWarnings = FALSE)
coef13 <- read_analysis(13, "structural_coefficients.csv")
coef8 <- read_analysis(8, "regression_coefficients.csv")
scores <- read_analysis(13, "model_comparison.csv")
edges <- read_analysis(13, "candidate_edges.csv")
equiv <- read_analysis(13, "equivalent_orientations.csv")
classes <- read_analysis(13, "equivalence_classes.csv")
pt <- function(x) x / ggplot2::.pt
ink <- "#202020"; blue <- "#0072B2"; orange <- "#D55E00"; grey <- "#777777"
diagram_theme <- function() theme_void(base_family = resolve_gbe_font(), base_size = 8.5) +
  theme(plot.title = element_text(size = 9, face = "bold", margin = margin(b = 4)),
        plot.subtitle = element_text(size = 7.5, lineheight = 1.08, margin = margin(b = 5)),
        plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.08, margin = margin(t = 5)),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.margin = margin(5, 6, 5, 6))

# A: One representative of H5, with its own structural equations.
network <- coef13 %>% filter(panel == "full24", model == "cell_body_branch", term != "(Intercept)") %>%
  mutate(parent = term, child = response, beta_label = sprintf("%.2f", estimate),
         arrow_width_mm = 1.6 * abs(estimate))
stopifnot(nrow(network) == 3, all(network$estimate > 0),
  setequal(paste(network$parent, network$child), c("g c", "c n", "c b")))
nodes_a <- tibble(node = c("g", "c", "n", "b"),
  x = c(1.4, 5, 8.6, 5), y = c(5, 5, 5, 1.25),
  label = c("Genome\nsize", "Cell\narea", "Nuclear\narea", "Body\nlength"),
  half_w = c(1.18, 1.02, 1.18, 1.20), half_h = .63)
network <- network %>% left_join(nodes_a %>% select(parent = node, x, y, half_w, half_h), by = "parent") %>%
  left_join(nodes_a %>% select(child = node, xend = x, yend = y, end_w = half_w, end_h = half_h), by = "child") %>%
  mutate(dx = xend - x, dy = yend - y,
    xs = x + ifelse(dx != 0, sign(dx) * (half_w + .08), 0),
    xe = xend - ifelse(dx != 0, sign(dx) * (end_w + .09), 0),
    ys = y + ifelse(dy != 0, sign(dy) * (half_h + .05), 0),
    ye = yend - ifelse(dy != 0, sign(dy) * (end_h + .06), 0),
    label_x = ifelse(dy != 0, x - .63, (x + xend) / 2),
    label_y = ifelse(dy != 0, (y + yend) / 2, y + .63))
a <- ggplot() +
  geom_segment(data = network, aes(xs, ys, xend = xe, yend = ye, linewidth = arrow_width_mm),
    colour = ink, arrow = grid::arrow(length = grid::unit(2.1, "mm"), type = "closed")) +
  scale_linewidth_identity() +
  geom_rect(data = nodes_a, aes(xmin = x - half_w, xmax = x + half_w, ymin = y - half_h, ymax = y + half_h),
    fill = "#F3F5F6", colour = "#484848", linewidth = .4) +
  geom_text(data = nodes_a, aes(x, y, label = label), family = resolve_gbe_font(), size = pt(8), lineheight = 1) +
  geom_text(data = network, aes(label_x, label_y, label = beta_label), family = resolve_gbe_font(), size = pt(8.5), fontface = "bold") +
  coord_cartesian(xlim = c(0, 10), ylim = c(.2, 6.25), expand = FALSE, clip = "off") +
  labs(title = "Representative hypothesis", subtitle = "H5 class; 24 species",
    caption = "Arrows assume genome-first ordering.\nWidth and labels show standardized coefficients.") + diagram_theme()

# B: Separate conditional checks, explicitly not the equations of A.
conditional <- coef8 %>% filter(panel == "full24") %>%
  filter((relation_id == "genome_to_cell_body_adjusted" & term == "g") |
         (relation_id == "genome_and_cell_to_nucleus_body_adjusted" & term %in% c("c", "g"))) %>%
  mutate(y = case_when(response == "c" ~ 3, term == "c" ~ 2, TRUE ~ 1),
    row_label = case_when(response == "c" ~ "Genome → cell\n+ body length",
                         term == "c" ~ "Cell → nucleus\n+ genome, body length",
                         TRUE ~ "Genome → nucleus\n+ cell, body length"),
    includes_zero = ci_low <= 0 & ci_high >= 0,
    beta_label = sprintf("%.2f", standardized_coefficient))
stopifnot(nrow(conditional) == 3, all(conditional$n_species == 24),
  isTRUE(all.equal(conditional$standardized_coefficient[conditional$y == 3], .5453794365491734)))
b <- ggplot(conditional, aes(standardized_coefficient, y)) +
  geom_vline(xintercept = 0, colour = "#AAAAAA", linewidth = .35, linetype = "dashed") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = y, colour = includes_zero), linewidth = .7) +
  geom_point(aes(colour = includes_zero), size = 2.1) +
  geom_text(aes(label = beta_label), nudge_y = .25, family = resolve_gbe_font(), size = pt(7.5)) +
  scale_colour_manual(values = c("FALSE" = blue, "TRUE" = grey), guide = "none") +
  scale_y_continuous(breaks = conditional$y, labels = conditional$row_label, limits = c(.55, 3.6), expand = expansion(mult = 0)) +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1), limits = c(-.62, 1.30)) +
  labs(title = "Conditional associations", subtitle = "Separate checks; 24 species",
    x = "Standardized coefficient (95% CI)", y = NULL,
    caption = "+ indicates included covariates.\nIntervals condition on fitted phylogenetic λ.") +
  theme_gbe() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5, lineheight = 1), axis.title.x = element_text(size = 7.5),
    plot.caption = element_text(size = 7.5, colour = "#333333"), panel.grid = element_blank())

# C: Complete eight-class comparison, with undirected skeletons.
model_order <- scores %>% filter(panel == "full24", test_order == "primary") %>% arrange(rank) %>% pull(model)
names_short <- c(cell_body_branch = "H5  Cell–body link", genome_body_branch = "H7  Genome–body link",
  body_upstream_cell = "H1  Cell bridge + body", body_disconnected = "H8  Body disconnected",
  body_upstream_genome = "H3  Genome fork + body", body_upstream_nucleus = "H2  Nucleus bridge + body",
  nucleus_cell_body_chain = "H6  Nucleus–cell–body", body_only = "H4  Body only")
comparisons <- scores %>% filter(test_order == "primary") %>%
  mutate(y = 9 - match(model, model_order), x_start = ifelse(panel == "full24", 0, .84),
    x_end = x_start + w, weight_label = ifelse(w < .001, "<0.001", sprintf("%.3f", w)),
    label_x = x_end + .035, group_colour = ifelse(panel == "full24", blue, orange))
stopifnot(nrow(comparisons) == 16, max(abs(tapply(comparisons$w, comparisons$panel, sum) - 1)) < 1e-10)
key_nodes <- tibble(node = c("g", "c", "n", "b"),
  node_x = c(-.59, -.405, -.22, -.405), node_y = c(0, .31, 0, -.31))
skeleton_nodes <- tidyr::crossing(model = model_order, node = key_nodes$node) %>%
  left_join(key_nodes, by = "node") %>% mutate(row_y = 9 - match(model, model_order),
    x = node_x, y = row_y + node_y, label = toupper(node))
skeleton_edges <- edges %>%
  left_join(key_nodes %>% select(parent = node, x = node_x, y_offset = node_y), by = "parent") %>%
  left_join(key_nodes %>% select(child = node, xend = node_x, yend_offset = node_y), by = "child") %>%
  mutate(row_y = 9 - match(model, model_order), y = row_y + y_offset, yend = row_y + yend_offset)
model_labels <- tibble(model = model_order, y = 9 - seq_along(model_order), label = unname(names_short[model_order]))
ticks <- tidyr::crossing(panel = c("full24", "subset22"), tick = c(0, .25, .5)) %>%
  mutate(x = tick + ifelse(panel == "full24", 0, .84), label = ifelse(tick == 0, "0", sprintf("%.2f", tick)))
c <- ggplot() +
  annotate("rect", xmin = -1.33, xmax = 1.62, ymin = seq(1.5, 7.5, 2), ymax = seq(2.5, 8.5, 2), fill = "#F5F6F7", colour = NA) +
  geom_segment(data = ticks, aes(x = x, xend = x, y = .7, yend = 8.4), colour = "#E1E4E6", linewidth = .25) +
  geom_segment(data = comparisons, aes(x = x_start, xend = x_end, y = y, yend = y, colour = group_colour), linewidth = 1.9) +
  geom_point(data = comparisons, aes(x = x_end, y = y, colour = group_colour), size = 1.1) +
  scale_colour_identity() +
  geom_text(data = comparisons, aes(label_x, y, label = weight_label), hjust = 0, family = resolve_gbe_font(), size = pt(7.5)) +
  geom_text(data = model_labels, aes(-1.30, y, label = label), hjust = 0, family = resolve_gbe_font(), size = pt(7.5)) +
  geom_segment(data = skeleton_edges, aes(x = x, xend = xend, y = y, yend = yend), linewidth = .35, colour = "#666666") +
  geom_point(data = skeleton_nodes, aes(x, y), shape = 21, fill = "white", colour = "white", size = 3.3, stroke = 0) +
  geom_text(data = skeleton_nodes, aes(x, y, label = label), family = resolve_gbe_font(), size = pt(7.5)) +
  geom_text(data = ticks, aes(x, .35, label = label), family = resolve_gbe_font(), size = pt(7.5)) +
  annotate("text", x = -.405, y = 9.02, label = "Structure", family = resolve_gbe_font(), size = pt(7.5), fontface = "bold") +
  annotate("text", x = .29, y = 9.02, label = "All 24 species", family = resolve_gbe_font(), size = pt(8), fontface = "bold", colour = blue) +
  annotate("text", x = 1.16, y = 9.02, label = "22 species*", family = resolve_gbe_font(), size = pt(8), fontface = "bold", colour = orange) +
  annotate("text", x = .70, y = -.20, label = "Relative model weight", family = resolve_gbe_font(), size = pt(8)) +
  coord_cartesian(xlim = c(-1.33, 1.62), ylim = c(-.38, 9.2), expand = FALSE, clip = "off") +
  labs(title = "Support across eight candidate classes",
    caption = "Weights compare candidates within each species set; they do not identify causal direction.\n*Excludes D. aeneus and D. wrighti. G = genome, C = cell, N = nucleus, B = body.") + diagram_theme()

# D: Four exactly enumerated orientations of H5, not a simplified 3-node graph.
orientations <- equiv %>% filter(model == "cell_body_branch")
roots <- orientations %>% group_by(orientation) %>%
  summarise(root = setdiff(unique(c(parent, child)), unique(child)), .groups = "drop") %>%
  mutate(order = match(root, c("g", "c", "n", "b")),
    label = c(g = "Genome first", c = "Cell first", n = "Nucleus first", b = "Body first")[root])
stopifnot(nrow(roots) == 4, setequal(roots$root, c("g", "c", "n", "b")),
  classes$equivalent_orientations[classes$model == "cell_body_branch"] == 4)
tiny_nodes <- tibble(node = c("g", "c", "n", "b"), nx = c(.3, 1.17, 2.04, 1.17), ny = c(1.95, 1.95, 1.95, .60))
equiv_nodes <- tidyr::crossing(orientation = roots$orientation, node = tiny_nodes$node) %>%
  left_join(roots, by = "orientation") %>% left_join(tiny_nodes, by = "node") %>%
  mutate(x = nx + (order - 1) * 2.55, y = ny, label_node = toupper(node))
equiv_edges <- orientations %>% left_join(roots, by = "orientation") %>%
  left_join(tiny_nodes %>% select(parent = node, x = nx, y = ny), by = "parent") %>%
  left_join(tiny_nodes %>% select(child = node, xend = nx, yend = ny), by = "child") %>%
  mutate(dx = xend - x, dy = yend - y,
    xs = x + (order - 1) * 2.55 + ifelse(dx != 0, sign(dx) * .19, 0),
    xe = xend + (order - 1) * 2.55 - ifelse(dx != 0, sign(dx) * .19, 0),
    ys = y + ifelse(dy != 0, sign(dy) * .19, 0),
    ye = yend - ifelse(dy != 0, sign(dy) * .19, 0))
d <- ggplot() +
  geom_segment(data = equiv_edges, aes(xs, ys, xend = xe, yend = ye), linewidth = .5, colour = ink,
    arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")) +
  geom_point(data = equiv_nodes, aes(x, y), shape = 21, size = 5, fill = "#F3F5F6", colour = "#666666", stroke = .4) +
  geom_text(data = equiv_nodes, aes(x, y, label = label_node), family = resolve_gbe_font(), size = pt(8)) +
  geom_text(data = roots, aes(x = 1.17 + (order - 1) * 2.55, y = 2.80, label = label), family = resolve_gbe_font(), size = pt(8), fontface = "bold") +
  coord_cartesian(xlim = c(0, 10), ylim = c(.15, 3.1), expand = FALSE, clip = "off") +
  labs(title = "Equivalent directions within H5",
    subtitle = "All four orientations imply the same conditional independences.",
    caption = "These arrows are schematic; the data do not select one orientation.") + diagram_theme()

write_csv(network, file.path(out, "display", "representative_network_coefficients.csv"))
write_csv(conditional, file.path(out, "display", "conditional_coefficients.csv"))
write_csv(comparisons, file.path(out, "display", "expanded_candidate_weights.csv"))
write_csv(orientations, file.path(out, "display", "h5_equivalent_orientations.csv"))
legend_a <- paste("One genome-first representative of H5, the cell-linked body-size class, for all 24 species.",
  "Equations are standardized log10 cell area ~ genome size; nuclear area ~ cell area; maximum body length ~ cell area.",
  "Arrow width is proportional to the absolute standardized regression coefficient; labels are rounded to two decimals.",
  "These are the Figure 13 structural coefficients, not the body-adjusted coefficients in panel b.",
  "Directions are model assumptions. H5 and H7 have similar relative support, and H5 contains four observationally equivalent orientations.")
legend_b <- paste("Separate conditional-regression checks on all 24 species, using the Figure 8 exports:",
  "cell area ~ genome size + body length; nuclear area ~ cell area + genome size + body length.",
  "Points give standardized coefficients and segments the 95% t intervals, with n minus regression coefficient count degrees of freedom.",
  "Plus signs in row labels identify included covariates. Blue intervals exclude zero; grey intervals include zero.",
  "These checks are distinct from the representative equations in panel a and do not identify causal effects.")
legend_c <- paste("Relative CICc weights for the expanded eight-class candidate set, using primary independence-test order b,g,c,n.",
  "Full24 and the subset excluding D. aeneus and D. wrighti (22 species) are fitted and normalized separately.",
  "Rows are ordered by full24 support. Miniature undirected graphs show class skeletons; they do not assert that every orientation belongs to that class.",
  "No class contains an unshielded collider. H1–H4 correspond to the former body-upstream representatives.",
  "H5 links genome, nucleus and body through cell; H6 links genome–nucleus–cell–body; H7 links body–genome–cell–nucleus; H8 disconnects body.",
  "Weights are support relative to this candidate set, not probabilities of causal truth. H5 and H7 are separated by only 0.624 CICc in full24.",
  "Reversing the independence-test order changes relative support, and the three leading subset22 classes are closely ranked.")
legend_d <- paste("The four enumerated acyclic orientations of H5 share the same skeleton and have no unshielded colliders.",
  "They therefore encode the same graphical conditional-independence claims. Arrow widths in this panel have no effect-size meaning.",
  "Graphical equivalence does not require separately fitted reversed phylogenetic regressions to return identical numerical scores when covariance parameters are re-estimated.")
legend <- paste("Figure 6. Cellular associations and alternative causal hypotheses.",
  "(a) One genome-first orientation of H5 for 24 species. Arrow widths and labels encode absolute standardized coefficients from",
  "cell area ~ genome size, nuclear area ~ cell area, and body length ~ cell area. Directions are hypotheses.",
  "(b) Separate conditional checks: cell area ~ genome size + body length and nuclear area ~ cell area + genome size + body length.",
  "Points and segments show coefficients and 95% t intervals; grey intervals include zero. These are different equations from panel a.",
  "(c) CICc weights across eight candidate classes for 24 species and after excluding D. aeneus and D. wrighti (22 species).",
  "Both excluded species are small direct developers; the comparison changes body-size range, development and taxon composition together.",
  "Weights sum to one within each cohort and candidate set; they are not probabilities of causal truth.",
  "H5 and H7 differ by only 0.624 CICc in full24. Skeletons show adjacency, not arrow direction.",
  "The primary independence-test order is b,g,c,n; reversing it changes support, with the leading subset22 classes closely ranked.",
  "(d) All four acyclic H5 orientations share a skeleton and no unshielded colliders, hence the same graphical conditional independences.",
  "These schematic arrows do not encode effect size. Re-estimated residual covariance can yield numerical differences between reversed regressions.",
  "G = genome size; C = erythrocyte area; N = nuclear area; B = species maximum snout–vent length.",
  "All log10 traits use fixed full24 means and sample SDs. Regressions use Pagel-lambda covariance;",
  "intervals condition on fitted lambda and omit tree and measurement uncertainty.",
  "The expanded eight-class comparison has not undergone measurement or tree resampling.",
  "Genome-size calibration and selection of large eligible erythrocytes limit biological interpretation.")
alt <- "Four panels show a weighted genome-first hypothesis, three body-adjusted regression coefficients, support for eight candidate classes in 24 and 22 species, and four equivalent causal orientations. Genome–cell and cell–nucleus coefficients are positive after body adjustment; the additional genome–nucleus interval spans zero. H5 and H7 lead in 24 species, whereas H8, H7 and H5 are close in 22. The four H5 orientations cannot be distinguished by graphical conditional independences."
combined <- ((a | b) / c / d) + plot_layout(heights = c(1, 1.64, .66), widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(family = resolve_gbe_font(), size = 11, face = "bold")))
save_candidate(6, "Cellular associations and alternative causal hypotheses", combined,
  panels = list(
    a = list(plot = a, title = "Representative genome-first hypothesis", width_mm = 90, height_mm = 79, legend = legend_a,
      alt_text = "Genome size points to cell area with coefficient 0.72, cell area points to nuclear area with 0.58, and cell area points to maximum body length with 0.36. These directions are assumed."),
    b = list(plot = b, title = "Conditional associations", width_mm = 105, height_mm = 79, legend = legend_b,
      alt_text = "Body-adjusted genome–cell coefficient is 0.55; the cell–nucleus coefficient adjusted for genome and body is 0.65; the additional genome–nucleus coefficient is 0.05 with a wide interval spanning zero."),
    c = list(plot = c, title = "Support across eight candidate classes", width_mm = 180, height_mm = 118, legend = legend_c,
      alt_text = "Eight candidate skeletons and their model weights are displayed. Full24 weights for H5 and H7 are 0.524 and 0.384. Subset22 weights for H8, H7 and H5 are 0.387, 0.295 and 0.283."),
    d = list(plot = d, title = "Equivalent directions within H5", width_mm = 180, height_mm = 51, legend = legend_d,
      alt_text = "Four H5 star diagrams have genome, cell, nucleus or body at the root. All retain a cell-centered skeleton without a collider and share the same graphical conditional independences.")),
  width_mm = 180, height_mm = 235,
  datasets = "Frozen Figure 8 regression_coefficients.csv; Figure 13 structural_coefficients.csv, model_comparison.csv, candidate_edges.csv, equivalence_classes.csv, equivalent_orientations.csv",
  legend = legend, alt_text = alt)
writeLines(c("Display-only assertions passed:",
  "Three H5 edges and three explicitly conditioned regression checks retain all 24 species.",
  "All 16 primary-order candidate weights displayed; weights sum to one within each cohort.",
  "All four enumerated H5 orientations displayed; four distinct roots verified.",
  "No fitting or resampling functions called; all fitted quantities read from frozen inputs."),
  file.path(out, "display", "numerical_checks.txt"))

#!/usr/bin/env Rscript
# Isolated manuscript candidates; source data and the existing gallery are read-only.
source("scripts/common.R")

annotations <- read_supp("Supplementary_Data_S42_te_bonsai_tip_annotations.csv") %>%
  mutate(reproductive_strategy = factor(reproductive_strategy, levels = te_strategy_levels))
sizes <- read_supp("Supplementary_Data_S40_te_genome_size_overlap21.csv")
stopifnot(nrow(annotations) == 34L, nrow(sizes) == 21L,
          identical(as.integer(table(annotations$reproductive_strategy)), c(7L, 24L, 3L)))
species_order <- phylogenetic_species_order(annotations$species)

# Figure 1: topology-only tree with directly aligned trait and composition tracks.
tree <- prune_source_tree(annotations$species)
n_tips <- Ntip(tree)
node_y <- node.height(tree)
edges <- reorder.phylo(tree, "cladewise")$edge
depth <- numeric(n_tips + Nnode(tree))
for (i in seq_len(nrow(edges))) depth[edges[i, 2]] <- depth[edges[i, 1]] + 1
tip_x <- max(depth) + 0.5
horizontal <- tibble(parent = tree$edge[, 1], child = tree$edge[, 2]) %>%
  mutate(x = depth[parent] / tip_x,
         xend = if_else(child <= n_tips, 1, depth[child] / tip_x),
         y = node_y[child], yend = node_y[child])
vertical <- tibble(parent = sort(unique(tree$edge[, 1]))) %>%
  rowwise() %>% mutate(x = depth[parent] / tip_x,
    ymin = min(node_y[tree$edge[tree$edge[, 1] == parent, 2]]),
    ymax = max(node_y[tree$edge[tree$edge[, 1] == parent, 2]])) %>% ungroup()
tracks <- tibble(species = tree$tip.label, y = node_y[seq_len(n_tips)]) %>%
  left_join(annotations, by = "species", relationship = "one-to-one") %>%
  left_join(sizes %>% select(species, genome_size_pg_fuscus_anchored,
    genome_size_pg_ci_low, genome_size_pg_ci_high, n_nuclei, n_images, n_specimens,
    sample_size_support, quality_match_status), by = "species", relationship = "one-to-one") %>%
  mutate(genome_size_available = !is.na(genome_size_pg_fuscus_anchored),
         genome_label = if_else(genome_size_available,
           sprintf("%.1f", genome_size_pg_fuscus_anchored), "n/a"))
stopifnot(identical(tracks$species[order(tracks$y)], species_order),
          sum(tracks$genome_size_available) == 21L,
          sum(!tracks$genome_size_available) == 13L)
composition <- read_csv(file.path(RUN, "inputs", "te_order_composition.csv"),
                        show_col_types = FALSE) %>%
  pivot_longer(-species, names_to = "source_order", values_to = "fraction") %>%
  mutate(order = if_else(source_order %in% c("LTR", "LINE", "TIR"), source_order, "Other")) %>%
  group_by(species, order) %>% summarise(fraction = sum(fraction), .groups = "drop") %>%
  mutate(order = factor(order, levels = c("LTR", "LINE", "TIR", "Other"))) %>%
  arrange(species, order) %>% group_by(species) %>%
  mutate(xmax = cumsum(fraction), xmin = lag(xmax, default = 0)) %>% ungroup() %>%
  left_join(tracks %>% select(species, y), by = "species", relationship = "many-to-one")
stopifnot(nrow(composition) == 136L, all(is.finite(composition$y)),
          all(abs((composition %>% group_by(species) %>% summarise(total = sum(fraction)))$total - 1) < 1e-10))

strategy_guide <- guide_legend(nrow = 1, override.aes = list(size = 2))
f1a <- ggplot() +
  geom_segment(data = horizontal, aes(x, y, xend = xend, yend = yend), linewidth = .25) +
  geom_segment(data = vertical, aes(x, ymin, xend = x, yend = ymax), linewidth = .25) +
  geom_point(data = tracks, aes(1, y, fill = reproductive_strategy, shape = reproductive_strategy),
             colour = GBE_COLORS[["charcoal"]], size = 1.85, stroke = .35) +
  geom_text(data = tracks, aes(1.12, y, label = display_species), hjust = 0,
            family = GBE_FONT, fontface = "italic", size = 7.5 / ggplot2::.pt) +
  scale_fill_manual(values = te_strategy_colors, breaks = te_strategy_levels,
    labels = te_strategy_labels, name = "Reproduction", guide = strategy_guide) +
  scale_shape_manual(values = te_strategy_shapes, breaks = te_strategy_levels,
    labels = te_strategy_labels, name = "Reproduction", guide = strategy_guide) +
  scale_x_continuous(limits = c(-.02, 2.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.5, 34.5), expand = c(0, 0)) +
  labs(title = "Phylogeny and\nreproductive strategy", x = NULL, y = NULL) +
  theme_gbe() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text = element_blank(), plot.title = element_text(size = 8.5, face = "bold", hjust = .5),
    legend.text = element_text(size = 7.5), legend.title = element_text(size = 7.5))

f1b <- ggplot(tracks, aes(genome_size_pg_fuscus_anchored, y)) +
  geom_segment(data = filter(tracks, genome_size_available),
    aes(x = genome_size_pg_ci_low, xend = genome_size_pg_ci_high, yend = y),
    colour = GBE_COLORS[["gray"]], linewidth = .4) +
  geom_point(data = filter(tracks, genome_size_available), shape = 21,
    fill = "white", size = 1.8, stroke = .45) +
  geom_text(data = filter(tracks, genome_size_available), aes(y = y + .34, label = genome_label),
    size = 7.5 / ggplot2::.pt, family = GBE_FONT, vjust = .5) +
  geom_text(data = filter(tracks, !genome_size_available), aes(x = 17.5, label = "n/a"),
    size = 7.5 / ggplot2::.pt, colour = GBE_COLORS[["gray"]], family = GBE_FONT) +
  scale_x_continuous(limits = c(8, 27), breaks = c(10, 18, 26), expand = expansion(mult = .04)) +
  scale_y_continuous(limits = c(.5, 34.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Estimated\ngenome size", x = "pg/1C", y = NULL) +
  theme_gbe() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(size = 8.5, face = "bold", hjust = .5))

f1c <- ggplot(composition) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - .34, ymax = y + .34, fill = order)) +
  scale_fill_manual(values = GBE_TE_ORDER_COLORS, breaks = c("LTR", "LINE", "TIR", "Other"),
    name = "TE order", guide = guide_legend(nrow = 1)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
    labels = label_percent(accuracy = 1), oob = squish, expand = c(0, 0)) +
  scale_y_continuous(limits = c(.5, 34.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Broad TE\ncomposition", x = "Classified aligned bases", y = NULL) +
  theme_gbe() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.title = element_text(size = 8.5, face = "bold", hjust = .5),
    legend.text = element_text(size = 7.5), legend.title = element_text(size = 7.5),
    plot.margin = margin(3, 13, 3, 4, unit = "pt"))

f1legend <- paste(
  "A shared LTR-rich repeat composition spans 34 sampled Desmognathus species.",
  "(a) The supplied species tree is pruned and shown as topology only; branch lengths are omitted because citation and calibration provenance remain unresolved.",
  "Tip symbols identify reproductive strategy (7 aquatic-egg species, 24 terrestrial-egg species with aquatic larvae, and 3 direct developers).",
  "(b) Estimated genome size is available for 21 species; n/a marks the 13 unavailable estimates.",
  "Points and one-decimal labels show D. fuscus assembly-anchored, nuclear-IOD-derived pg/1C estimates; bars show conditional 95% measurement-bootstrap intervals.",
  "These are not independently measured species C-values; intervals exclude anchor uncertainty and are conditional on observed images and animals. Species-level support remains in the accompanying tables.",
  "(c) LTR, LINE, TIR and Other are shares of classified dnaPipeTE aligned bases after the released technical-run aggregation.",
  "Other combines the remaining classified orders; unresolved mass is outside this denominator. Fractions are not whole-genome occupancy.",
  "All panels use one fixed phylogenetic row order.")
f1alt <- "Thirty-four species align across a topology-only tree with reproductive-strategy tip symbols, 21 estimated genome sizes with intervals and 13 explicit missing values, and stacked LTR, LINE, TIR and other classified repeat proportions."
f1 <- (f1a + f1b + f1c) + plot_layout(widths = c(1.32, .80, 1.05), guides = "collect") +
  plot_annotation(tag_levels = "a", caption = "Estimated genome size: IOD-derived, D. fuscus anchored.  n/a: no estimate.  Composition: classified dnaPipeTE aligned bases.") &
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.caption = element_text(size = 7, colour = GBE_COLORS[["gray"]], hjust = 0))

standalone_rows <- function(p) p +
  scale_y_continuous(limits = c(.5, 34.5), breaks = tracks$y,
    labels = tracks$display_species, expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 7.5, face = "italic"),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())
panels1 <- list(
  a = list(plot = f1a + guides(fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)),
    title = "Species phylogeny and reproductive strategy", width_mm = 105, height_mm = 240,
    legend = "Topology-only supplied-tree prune for all 34 TE species, with reproductive-strategy tip symbols. No branch-length or ancestral-state inference is shown.",
    alt_text = "Topology-only phylogeny with all 34 species labels and strategy symbols."),
  b = list(plot = standalone_rows(f1b), title = "Estimated genome size and missingness", width_mm = 110, height_mm = 230,
    legend = "The same 34 species retain their phylogenetic row order. Twenty-one IOD-derived, D. fuscus-anchored genome-size estimates have conditional measurement intervals; 13 are labelled n/a. Intervals omit anchor uncertainty and population-level variation.",
    alt_text = "Labelled genome-size forest plot with all 34 species and explicit missing values."),
  c = list(plot = standalone_rows(f1c), title = "Shared broad TE composition", width_mm = 125, height_mm = 230,
    legend = "Thirty-four species' LTR, LINE, TIR and remaining classified-order proportions from dnaPipeTE. Every row sums to 100% of classified aligned bases; this denominator excludes unresolved mass and does not estimate genome occupancy.",
    alt_text = "Labelled stacked bars show LTR-rich classified repeat composition across all 34 species."))
save_candidate(1, "Shared repeat architecture across the radiation", f1, panels1, 180, 245,
  "S40; S42; te_order_composition.csv; trees/desmo900dated_source46.tre", f1legend, f1alt)
dir.create("figures/figure_1/analysis", showWarnings = FALSE)
write_csv(tracks %>% arrange(y), "figures/figure_1/analysis/species_tracks.csv")
write_csv(composition %>% arrange(y, order), "figures/figure_1/analysis/broad_composition.csv")
write_csv(horizontal, "figures/figure_1/analysis/tree_horizontal_segments.csv")
write_csv(vertical, "figures/figure_1/analysis/tree_vertical_segments.csv")

# Figure 2: unchanged all-hit denominator, zero-filled complete profile grid.
landscape <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
orders <- c("LTR", "LINE", "TIR")
profiles <- expand_grid(species = species_order, order = orders, divergence_bin_start_pct = 0:49) %>%
  left_join(landscape %>% filter(order %in% orders, divergence_bin_start_pct < 50) %>%
    select(species, order, divergence_bin_start_pct, percent_species_hit_bp),
    by = c("species", "order", "divergence_bin_start_pct"), relationship = "one-to-one") %>%
  mutate(percent_species_hit_bp = replace_na(percent_species_hit_bp, 0),
         order = factor(order, levels = orders)) %>%
  left_join(annotations %>% select(species, reproductive_strategy), by = "species", relationship = "many-to-one")
means <- profiles %>% group_by(reproductive_strategy, order, divergence_bin_start_pct) %>%
  summarise(mean_percent = mean(percent_species_hit_bp), n_species = n(), .groups = "drop")
low <- profiles %>% filter(divergence_bin_start_pct < 5) %>%
  group_by(species, order, reproductive_strategy) %>%
  summarise(low_divergence_percent = sum(percent_species_hit_bp), .groups = "drop") %>%
  mutate(display_species = factor(species_label(species), levels = species_label(species_order)))
stopifnot(nrow(profiles) == 5100L, nrow(means) == 450L, nrow(low) == 102L,
  identical(sort(unique(means$n_species)), c(3L, 7L, 24L)),
  all(profiles$percent_species_hit_bp >= 0 & profiles$percent_species_hit_bp <= 3),
  all(low$low_divergence_percent >= 0 & low$low_divergence_percent <= 6),
  !anyNA(profiles$reproductive_strategy))
strategy_n_labels <- c("Aquatic eggs (n = 7)", "Terrestrial eggs + aquatic larvae (n = 24)", "Direct development (n = 3)")
line_types <- setNames(c("solid", "longdash", "dotdash"), te_strategy_levels)
f2a <- ggplot(profiles, aes(divergence_bin_start_pct, percent_species_hit_bp, group = species)) +
  geom_line(colour = "#B9B9B9", linewidth = .20, alpha = .45) +
  geom_line(data = means, aes(y = mean_percent, group = reproductive_strategy,
    colour = reproductive_strategy, linetype = reproductive_strategy), linewidth = .70) +
  facet_grid(cols = vars(order)) +
  scale_colour_manual(values = te_strategy_colors, breaks = te_strategy_levels, labels = strategy_n_labels, name = NULL) +
  scale_linetype_manual(values = line_types, breaks = te_strategy_levels, labels = strategy_n_labels, name = NULL) +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3, expand = expansion(mult = c(0, .025))) +
  labs(x = "RepeatMasker sequence divergence (%)", y = "All aligned hit bp per bin (%)") +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(panel.spacing = unit(16, "pt"), legend.position = "none")
f2b <- ggplot(low, aes(low_divergence_percent, display_species)) +
  geom_segment(aes(x = 0, xend = low_divergence_percent, yend = display_species),
    colour = "#DEDEDE", linewidth = .32) +
  geom_point(aes(fill = reproductive_strategy, colour = reproductive_strategy,
    shape = reproductive_strategy), size = 1.85, stroke = .4) +
  facet_grid(cols = vars(order)) +
  scale_colour_manual(values = te_strategy_colors, breaks = te_strategy_levels, labels = strategy_n_labels, name = NULL) +
  scale_fill_manual(values = te_strategy_colors, breaks = te_strategy_levels, labels = strategy_n_labels, name = NULL) +
  scale_shape_manual(values = te_strategy_shapes, breaks = te_strategy_levels, labels = strategy_n_labels, name = NULL) +
  scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6), expand = expansion(mult = c(.025, .035))) +
  labs(x = "All aligned hit bp below 5% divergence (%)", y = NULL) +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(axis.text.y = element_text(size = 7.5, face = "italic"),
    panel.spacing = unit(16, "pt"), legend.position = "none")
# One shared key records both the strategy line and point encodings.
key <- ggplot(tibble(reproductive_strategy = factor(te_strategy_levels, levels = te_strategy_levels),
                     y = c(3, 2, 1), label = strategy_n_labels)) +
  geom_segment(aes(x = 0, xend = .10, y = y, yend = y, colour = reproductive_strategy,
                   linetype = reproductive_strategy), linewidth = .65) +
  geom_point(aes(x = .05, y = y, colour = reproductive_strategy,
                 fill = reproductive_strategy, shape = reproductive_strategy), size = 1.85, stroke = .4) +
  geom_text(aes(x = .12, y = y, label = label), family = GBE_FONT, size = 7.5 / ggplot2::.pt, hjust = 0) +
  scale_colour_manual(values = te_strategy_colors) + scale_fill_manual(values = te_strategy_colors) +
  scale_shape_manual(values = te_strategy_shapes) + scale_linetype_manual(values = line_types) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.5, 3.5), expand = c(0, 0)) +
  theme_void(base_family = GBE_FONT) + theme(legend.position = "none", plot.margin = margin(1, 0, 1, 0))
key_caption <- ggplot() + annotate("text", x = 0, y = .5,
  label = "Gray curves: individual species\nColoured curves: equal-weight group means",
  family = GBE_FONT, size = 7.5 / ggplot2::.pt, colour = GBE_COLORS[["gray"]], hjust = 0, lineheight = 1.25) +
  xlim(0, 1) + ylim(0, 1) + theme_void()
key_row <- key + key_caption + plot_layout(widths = c(1.22, 1))
key_points <- key
key_points$layers <- key_points$layers[-1]
point_caption <- ggplot() + annotate("text", x = 0, y = .5,
  label = "Points: individual species\nRows: shared phylogenetic order",
  family = GBE_FONT, size = 7.5 / ggplot2::.pt, colour = GBE_COLORS[["gray"]], hjust = 0, lineheight = 1.25) +
  xlim(0, 1) + ylim(0, 1) + theme_void()
key_row_points <- key_points + point_caption + plot_layout(widths = c(1.22, 1))
f2 <- ((f2a + labs(tag = "a")) / key_row / (f2b + labs(tag = "b"))) +
  plot_layout(heights = c(1, .25, 2.0))
f2legend <- paste(
  "Representation of low-divergence repeats differs across sampled reproductive strategies, with the clearest separation among LTRs.",
  "(a) LTR, LINE and TIR profiles show all 34 species as thin gray curves and equal-weight means for 7 aquatic-egg species, 24 terrestrial-egg species with aquatic larvae, and 3 direct developers as coloured curves.",
  "Each 1-percentage-point bin is expressed as its share of all RepeatMasker-aligned hit bp for that species; missing order/bin combinations are zero-filled across bins 0–49.",
  "All orders share 0–50% divergence and 0–3% aligned-hit-bp axes. No within-order normalization is applied.",
  "(b) For every species, the same all-hit denominator is used for the sum of bins 0–4, shown on matched 0–6% axes in the shared phylogenetic order.",
  "RepeatMasker alignments come from dnaPipeTE repeat-contig assemblies. Strategy comparisons are descriptive and not phylogenetically corrected.",
  "Sequence divergence is not insertion age or a measured transposition rate; aligned-hit fractions are not whole-genome occupancy. No biological confidence intervals are inferred from the aligned hits.")
f2alt <- "Three equal-scale divergence panels show individual species and reproductive-strategy means. Below, a 34-species dot plot compares LTR, LINE and TIR fractions below 5% divergence. The three direct developers have the lowest low-divergence LTR fractions; aquatic-egg species are intermediate; terrestrial-egg species with aquatic larvae are higher."
panels2 <- list(
  a = list(plot = ((f2a + labs(tag = "a")) / key_row) + plot_layout(heights = c(1, .25)), title = "Repeat-divergence profiles", width_mm = 180, height_mm = 105,
    legend = f2legend, alt_text = "Shared-axis LTR, LINE and TIR sequence-divergence profiles: gray species curves and coloured equal-weight reproductive-strategy means."),
  b = list(plot = ((f2b + labs(tag = "b")) / key_row_points) + plot_layout(heights = c(1, .13)), title = "Low-divergence repeat representation", width_mm = 180, height_mm = 190,
    legend = "All 34 species retain the phylogenetic row order. Points sum bins 0–4 of the released RepeatMasker landscape for LTR, LINE and TIR and retain each species' full aligned-hit-bp denominator. The three axes are identically scaled 0–6%. Reproductive-strategy groupings are descriptive, with no new test or biological interval.",
    alt_text = "All-species dot plot with three matched columns for LTR, LINE and TIR fractions below 5% divergence."))
save_candidate(2, "Contrasting repeat-divergence profiles", f2, panels2, 180, 270,
  "S12; S42; trees/desmo900dated_source46.tre", f2legend, f2alt)
dir.create("figures/figure_2/analysis", showWarnings = FALSE)
write_csv(profiles, "figures/figure_2/analysis/species_divergence_profiles.csv")
write_csv(means, "figures/figure_2/analysis/strategy_mean_profiles.csv")
write_csv(low, "figures/figure_2/analysis/low_divergence_fractions.csv")
write_lines(c("All 34 species retained; strategy counts 7/24/3.",
  "Profiles: 5,100 species-order-bin rows; 450 strategy-order-bin means; 102 low-divergence rows.",
  "Primary denominator remains all species RepeatMasker aligned-hit bp. Missing bins zero-filled, not reclosed.",
  "Figure 1 has 21 available genome estimates and 13 explicit missing values; every composition row closes to one.",
  "No inferential analyses were performed."), "figures/figure_2/analysis/QA.txt")
message("Completed isolated figures 1 and 2")

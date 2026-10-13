source('scripts/common.R')
analysis <- file.path(RUN, 'figures/figure_8/analysis')
out <- file.path(RUN, 'figures/figure_7/analysis')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
traits <- read_csv(file.path(analysis, 'traits_joined.csv'), show_col_types = FALSE)
fits <- read_csv(file.path(analysis, 'regression_coefficients.csv'), show_col_types = FALSE) %>%
  rename(beta = standardized_coefficient, se = standard_error, lambda = fitted_phylogenetic_parameter)
curves <- read_csv(file.path(analysis, 'pairwise_predictions.csv'), show_col_types = FALSE) %>%
  rename(svl_mm = max_svl_mm, fit = fitted, low = ci_low, high = ci_high)
morph <- read_supp('Supplementary_Data_S19_cell_nucleus_species_estimates.csv') %>%
  select(species, cell_area_um2_ci_low, cell_area_um2_ci_high,
         nucleus_area_um2_ci_low, nucleus_area_um2_ci_high)
plot_data <- traits %>% mutate(species = species_label(tree_tip),
    genome_size_pg_fuscus_anchored_ci_low = genome_size_pg_ci_low,
    genome_size_pg_fuscus_anchored_ci_high = genome_size_pg_ci_high) %>%
  left_join(morph, by = 'species', relationship = 'one-to-one') %>%
  mutate(development = if_else(tree_tip %in% c('aeneus', 'wrighti'), 'Direct development', 'Aquatic-larval development'),
         development = factor(development, levels = c('Aquatic-larval development', 'Direct development')),
         species_order = factor(species, levels = rev(phylogenetic_species_levels(tree_tip))))
stopifnot(nrow(plot_data) == 24L, !anyDuplicated(plot_data$tree_tip),
          all(is.finite(plot_data$max_svl_mm)), all(plot_data$max_svl_mm > 0),
          all(is.finite(plot_data$cell_area_um2_ci_low)))
write_csv(plot_data, file.path(out, 'plotted_species_data.csv'))
write_csv(fits %>% filter(relation_id %in% c('body_to_genome','body_to_cell','body_to_nucleus')),
          file.path(out, 'body_size_regressions.csv'))
write_csv(curves, file.path(out, 'fitted_relationships.csv'))
colours <- c('Aquatic-larval development' = GBE_COLORS[['gray']], 'Direct development' = GBE_COLORS[['purple']])
shapes <- c('Aquatic-larval development' = 21, 'Direct development' = 24)
line_labels <- c(full24 = 'All 24 species', subset22 = '22 species (without aeneus + wrighti)')
style <- theme_gbe() + theme(plot.title = element_text(size = 9, face = 'bold', margin = margin(b = 5)),
  plot.margin = margin(6, 8, 5, 6), axis.text = element_text(size = 7.5),
  legend.text = element_text(size = 7.5), legend.position = 'bottom')
body_panel <- function() {
  ggplot(plot_data, aes(x = max_svl_mm, y = species_order)) +
    geom_segment(aes(x = 0, xend = max_svl_mm, yend = species_order), colour = '#D9D9D9', linewidth = .45) +
    geom_point(aes(colour = development, shape = development), fill = 'white', size = 2.1, stroke = .65) +
    geom_text(aes(label = format(max_svl_mm), x = max_svl_mm + 4), hjust = 0, family = GBE_FONT, size = 7 / ggplot2::.pt) +
    scale_colour_manual(values = colours, drop = FALSE, name = NULL) +
    scale_shape_manual(values = shapes, drop = FALSE, name = NULL) +
    scale_x_continuous(limits = c(0, 136), breaks = c(0, 40, 80, 120), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = .5)) +
    labs(title = 'Body size across the 24-species panel', x = 'Maximum SVL (mm)', y = NULL) + style +
    theme(axis.text.y = element_text(face = 'italic', size = 7), axis.line.y = element_blank(),
          axis.ticks.y = element_blank(), legend.position = 'none') +
    guides(colour = 'none', shape = 'none')
}
format_p <- function(x) ifelse(x < .001, formatC(x, format = 'e', digits = 1), sprintf('%.3f', x))
scatter <- function(id, column, y_title, title, label_all = FALSE) {
  estimates <- fits %>% filter(relation_id == id, predictor == 'body_size')
  prediction <- curves %>% filter(relation_id == id) %>%
    mutate(panel = factor(panel, levels = c('full24','subset22')))
  stopifnot(nrow(estimates) == 2L, setequal(estimates$panel, c('full24','subset22')),
    all(prediction$fit > 0), all(prediction$low > 0))
  estimates <- estimates[match(c('full24','subset22'), estimates$panel), ]
  annotation <- paste(sprintf('%s species: β = %.2f; P = %s', estimates$n_species,
                              estimates$beta, format_p(estimates$p_value)), collapse = '\n')
  labels <- plot_data %>% mutate(display_label = if (label_all) species else
    if_else(tree_tip %in% c('wrighti','aeneus','amphileucus','valtos','ochrophaeus'), species, ''))
  repel_layer <- if (label_all) function(...) geom_label_repel(..., fill = 'white',
    label.size = 0, label.padding = grid::unit(.6, 'mm')) else geom_text_repel
  p <- ggplot(plot_data, aes(x = max_svl_mm, y = .data[[column]])) +
    geom_ribbon(data = prediction %>% filter(panel == 'full24'),
      aes(x = svl_mm, y = NULL, ymin = low, ymax = high), inherit.aes = FALSE,
      fill = '#E1E1E1', alpha = .7) +
    geom_errorbar(aes(ymin = .data[[paste0(column, '_ci_low')]], ymax = .data[[paste0(column, '_ci_high')]], colour = development),
      width = 0, linewidth = .3, alpha = .55) +
    geom_line(data = prediction, aes(x = svl_mm, y = fit, linetype = panel),
      colour = GBE_COLORS[['charcoal']], linewidth = .6, inherit.aes = FALSE) +
    geom_point(aes(colour = development, shape = development), fill = 'white', size = 2.1, stroke = .7) +
    repel_layer(data = labels, aes(label = display_label), family = GBE_FONT, fontface = 'italic',
      size = 7 / ggplot2::.pt, seed = 20260924, max.overlaps = Inf, min.segment.length = 0,
      box.padding = .4, point.padding = .3, force = 3, force_pull = .15,
      max.time = Inf, max.iter = 30000, segment.size = .2, segment.color = '#AAAAAA') +
    annotate('text', x = Inf, y = Inf, label = annotation, hjust = 1.03, vjust = 1.15,
      family = GBE_FONT, size = 7 / ggplot2::.pt, lineheight = 1.15) +
    scale_colour_manual(values = colours, drop = FALSE, name = NULL) +
    scale_shape_manual(values = shapes, drop = FALSE, name = NULL) +
    scale_linetype_manual(values = c(full24 = 'solid', subset22 = '22'), labels = line_labels, name = NULL) +
    scale_x_log10(breaks = c(30, 50, 80, 120), limits = c(25, 140), expand = expansion(mult = .01)) +
    scale_y_log10(expand = expansion(mult = c(.15, .34))) +
    labs(title = title, x = 'Maximum SVL (mm; log scale)', y = y_title) + style +
    guides(colour = guide_legend(order = 1, ncol = 2), shape = guide_legend(order = 1, ncol = 2),
           linetype = guide_legend(order = 2, ncol = 1))
  p
}
panels <- list(a = body_panel(),
  b = scatter('body_to_genome','genome_size_pg_fuscus_anchored','Genome size (pg/1C; log scale)','Body size and genome size'),
  c = scatter('body_to_cell','cell_area_um2','Cell area (µm²; log scale)','Body size and cell area'),
  d = scatter('body_to_nucleus','nucleus_area_um2','Nucleus area (µm²; log scale)','Body size and nucleus area'))
legend <- paste(
  'Maximum snout–vent length (SVL) and cellular traits for all 24 species with genome-size and morphology data.',
  '(a) Species maximum SVL in the established phylogenetic display order; values are collaborator-provided species maxima, not measurements of the microscopy animals.',
  '(b–d) Genome size, cell area and nucleus area versus SVL. Both axes are logarithmic but tick values retain the measurement units.',
  'The solid line is the full 24-species Pagel-lambda PGLS; the dashed line refits the relationship after excluding only D. aeneus and D. wrighti (22 species). Each line is restricted to its fitted sample SVL range.',
  'Gray shading is the full-panel 95% confidence band for the fitted log-scale mean, back-transformed to original units; it is not a prediction interval.',
  'Annotations give coefficients standardized using full24 log10 means and SDs for both panels and unadjusted two-sided t P values. Full coefficient intervals and lambda estimates are exported.',
  'Vertical intervals show released species measurement intervals; their uncertainty and SVL uncertainty are not propagated into PGLS.',
  'The two direct developers are purple triangles; other species are gray circles. Five reference species are labeled in the collage and all species in the isolated scatterplots.',
  'The two excluded species are simultaneously the smallest taxa and the only direct developers. This sensitivity cannot identify body size or development as the cause of any coefficient change, and differing significance levels are not a test of slope differences.',
  'Genome size is IOD-derived and D. fuscus assembly-anchored; cell and nucleus areas summarize selected erythrocytes. These are exploratory interspecific associations, without multiplicity adjustment.')
collage <- (panels$a | panels$b) / (panels$c | panels$d) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom', legend.box = 'vertical')
panel_list <- list(
  a = list(plot = panels$a, title = 'Maximum SVL in phylogenetic species order', width_mm = 160, height_mm = 160,
    legend = legend, alt_text = 'Twenty-four species have recorded maximum SVL; wrighti is 29 mm, aeneus 31 mm and amphileucus 122 mm.'),
  b = list(plot = scatter('body_to_genome','genome_size_pg_fuscus_anchored','Genome size (pg/1C; log scale)','Body size and genome size',TRUE),
    title = 'Body size and genome size', width_mm = 180, height_mm = 135, legend = legend,
    alt_text = 'All 24 species labeled by maximum SVL and genome size; full-panel and 22-species sensitivity PGLS lines.'),
  c = list(plot = scatter('body_to_cell','cell_area_um2','Cell area (µm²; log scale)','Body size and cell area',TRUE),
    title = 'Body size and cell area', width_mm = 180, height_mm = 135, legend = legend,
    alt_text = 'All 24 species labeled by maximum SVL and cell area, with measurement intervals and two PGLS fits.'),
  d = list(plot = scatter('body_to_nucleus','nucleus_area_um2','Nucleus area (µm²; log scale)','Body size and nucleus area',TRUE),
    title = 'Body size and nucleus area', width_mm = 180, height_mm = 135, legend = legend,
    alt_text = 'All 24 species labeled by maximum SVL and nucleus area, with measurement intervals and two PGLS fits.'))
save_candidate(7, 'Body size, genome size and cellular dimensions', collage, panel_list, 180, 230,
  'S04;S19;inputs/body_size;figure_8/analysis', legend,
  'One phylogenetically ordered body-size panel and three trait relationships compare all 24 species with a sensitivity excluding the two smallest direct developers.')
message('Figure 7 rendered: four-panel collage and four isolated panels.')

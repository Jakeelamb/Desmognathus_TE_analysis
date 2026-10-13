source('scripts/common.R')

# Display-only refinement. Every observation, fit, and interval is read from the
# frozen prior run; this script contains no model-fitting or resampling call.
strategy <- read_supp('Supplementary_Data_S42_te_bonsai_tip_annotations.csv') %>%
  select(tree_tip = species, reproductive_strategy)
display_levels <- c(te_strategy_levels, 'unassigned')
display_colours <- c(te_strategy_colors[te_strategy_levels], unassigned = '#777777')
display_labels <- c(te_strategy_labels[te_strategy_levels], unassigned = 'Egg placement unassigned')
display_shapes <- c(te_strategy_shapes[te_strategy_levels], unassigned = 23)
strategy_scale <- function() scale_colour_manual(values = display_colours,
  breaks = display_levels, labels = display_labels, name = NULL, drop = FALSE)
add_strategy <- function(dat, id = 'tree_tip') {
  dat %>% select(-any_of('reproductive_strategy')) %>%
    left_join(strategy, by = setNames('tree_tip', id), relationship = 'many-to-one') %>%
    mutate(reproductive_strategy = factor(coalesce(reproductive_strategy, 'unassigned'), levels = display_levels))
}
style <- theme_gbe() + theme(legend.position = 'bottom', legend.key.width = grid::unit(3, 'mm'),
  legend.spacing.x = grid::unit(2, 'mm'), plot.margin = margin(6, 6, 5, 6))
shape_scale <- function() scale_shape_manual(values = display_shapes, breaks = display_levels,
  labels = display_labels, name = NULL, drop = FALSE)
strategy_guide <- guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2)),
  shape = guide_legend(nrow = 2, byrow = TRUE))
label_layer <- function(dat, all = FALSE, seed = 20260924) {
  # Identical observations; overview labels select fixed reference species only.
  geom_label_repel(data = dat, aes(label = display_label), family = GBE_FONT,
    fontface = 'italic', size = GBE_LABEL_SIZE_MM, seed = seed, max.overlaps = Inf,
    box.padding = .35, point.padding = .22, min.segment.length = 0,
    force = if (all) 5 else 2, force_pull = .10, max.time = Inf, max.iter = 75000,
    segment.size = .18, segment.color = '#999999', show.legend = FALSE,
    fill = 'white', linetype = 0, label.size = NA, label.padding = grid::unit(.35, 'mm'), label.r = grid::unit(0, 'mm'))
}
copy_analysis <- function(number, dat, filename) {
  folder <- file.path(RUN, 'figures', paste0('figure_', number), 'analysis')
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  write_csv(dat, file.path(folder, filename))
}

# FIGURE 5: two genome relationships, conditional allometry, paired allocation.
pair_data <- read_analysis(5, 'plotted_species_data.csv') %>%
  mutate(tree_tip = sub('^D\\. ', '', species)) %>% add_strategy()
pair_models <- read_analysis(5, 'model_summary.csv')
pair_predictions <- read_analysis(5, 'fitted_relationships.csv')
allometry <- read_analysis(10, 'released_species_data.csv') %>%
  rename(tree_tip = species) %>% add_strategy() %>% mutate(species = morphology_label)
allometry_models <- read_analysis(10, 'allometry_model_summary.csv')
allometry_predictions <- read_analysis(10, 'fitted_relationship.csv')
allocation <- read_analysis(10, 'equal_animal_species_data.csv') %>% add_strategy(id = 'species')
animals <- read_analysis(10, 'paired_animal_summaries.csv') %>% add_strategy(id = 'species')
display_order <- phylogenetic_species_order(allocation$species)
allocation <- allocation %>% mutate(y = match(species, display_order))
animals <- animals %>% select(-y, -animal_y) %>%
  left_join(allocation %>% select(species, y), by = 'species', relationship = 'many-to-one') %>%
  group_by(species) %>% mutate(animal_y = y + (row_number() - (n() + 1) / 2) * .13) %>% ungroup()
stopifnot(nrow(pair_data) == 24, nrow(allometry) == 24, nrow(allocation) == 24,
  nrow(animals) == 42, sum(allocation$n_paired_cells) == 1152,
  sum(pair_data$reproductive_strategy == 'unassigned') == 3,
  identical(sort(pair_data$tree_tip), sort(allometry$tree_tip)),
  all(pair_data$cell_area_um2[match(allometry$tree_tip, pair_data$tree_tip)] == allometry$cell_area_um2))
copy_analysis(5, pair_data, 'genome_cellular_species.csv')
copy_analysis(5, pair_models, 'frozen_genome_cellular_models.csv')
copy_analysis(5, pair_predictions, 'frozen_genome_cellular_predictions.csv')
copy_analysis(5, allometry_models, 'frozen_allometry_models.csv')
copy_analysis(5, allometry_predictions, 'frozen_allometry_predictions.csv')
copy_analysis(5, allocation, 'paired_nuclear_fraction_species.csv')
copy_analysis(5, animals, 'paired_nuclear_fraction_animals.csv')

genome_panel <- function(response, all = FALSE) {
  comparison_id <- if (response == 'cell_area_um2') 'genome_size_vs_cell_area' else 'genome_size_vs_nucleus_area'
  m <- pair_models %>% filter(comparison == comparison_id)
  pred <- pair_predictions %>% filter(comparison == comparison_id)
  dat <- pair_data %>% mutate(display_label = if (all) species else
    if_else(tree_tip %in% c('aeneus', 'wrighti', 'valtos'), species, ''))
  title <- if (response == 'cell_area_um2') 'Genome size and cell area' else 'Genome size and nucleus area'
  ylab <- if (response == 'cell_area_um2') 'Cell area (µm²)' else 'Nucleus area (µm²)'
  ggplot(dat, aes(genome_size_pg_fuscus_anchored, .data[[response]])) +
    geom_ribbon(data = pred, aes(x = x, y = NULL, ymin = low, ymax = high), inherit.aes = FALSE, fill = '#E4E4E4') +
    geom_line(data = pred, aes(x = x, y = fit), inherit.aes = FALSE, linewidth = .55, colour = '#222222') +
    geom_errorbar(aes(ymin = .data[[paste0(response, '_ci_low')]], ymax = .data[[paste0(response, '_ci_high')]]),
      width = 0, linewidth = .25, colour = '#ADADAD') +
    geom_errorbar(aes(xmin = genome_size_pg_ci_low, xmax = genome_size_pg_ci_high), orientation = 'y',
      width = 0, linewidth = .25, colour = '#ADADAD') +
    geom_point(aes(colour = reproductive_strategy, shape = reproductive_strategy), fill = 'white', size = 1.95, stroke = .75) +
    label_layer(dat, all) + strategy_scale() + shape_scale() + strategy_guide +
    scale_x_continuous(breaks = c(10, 15, 20, 25), expand = expansion(mult = c(.06, .07))) +
    scale_y_continuous(expand = expansion(mult = c(.12, .12))) +
    labs(title = title, subtitle = sprintf('β = %.2f [%.2f, %.2f]; P = %s; n = 24',
      m$pgls_standardized_beta, m$pgls_ci_low, m$pgls_ci_high, format_p(m$pgls_p_value)),
      x = 'Genome size (pg/1C)', y = ylab) + style
}

allometry_panel <- function(all = FALSE) {
  m <- allometry_models %>% filter(analysis_id == 'released24')
  dat <- allometry %>% mutate(display_label = if (all) species else
    if_else(tree_tip %in% c('aeneus', 'wrighti', 'ochrophaeus', 'brimleyorum'), species, ''))
  ggplot(dat, aes(cell_area_um2, nucleus_area_um2)) +
    geom_ribbon(data = allometry_predictions, aes(y = NULL, ymin = low, ymax = high), fill = '#E4E4E4') +
    geom_line(data = allometry_predictions, aes(y = fitted), linewidth = .55) +
    geom_line(data = allometry_predictions, aes(y = isometry), linetype = 'dashed', colour = '#777777', linewidth = .45) +
    geom_errorbar(aes(xmin = cell_area_um2_ci_low, xmax = cell_area_um2_ci_high), orientation = 'y', width = 0, colour = '#ADADAD', linewidth = .25) +
    geom_errorbar(aes(ymin = nucleus_area_um2_ci_low, ymax = nucleus_area_um2_ci_high), width = 0, colour = '#ADADAD', linewidth = .25) +
    geom_point(aes(colour = reproductive_strategy, shape = reproductive_strategy), fill = 'white', size = 1.95, stroke = .75) +
    label_layer(dat, all, 20260925) + strategy_scale() + shape_scale() + strategy_guide +
    scale_x_log10(breaks = c(200, 250, 300, 350), expand = expansion(mult = c(.12, .15))) +
    scale_y_log10(breaks = c(30, 40, 50, 60), expand = expansion(mult = c(.13, .15))) +
    labs(title = 'Nuclear–cell area scaling', subtitle = sprintf('Slope = %.2f [%.2f, %.2f]; n = 24\nP(slope = 1) = %s',
      m$slope, m$ci_low, m$ci_high, format_p(m$p_vs_isometry)),
      x = 'Cell area (µm²; log scale)', y = 'Nucleus area (µm²; log scale)',
      caption = 'Dashed line: proportional scaling (slope = 1).') + style
}

allocation_panel <- function(isolated = FALSE) {
  p <- ggplot(allocation, aes(nuclear_fraction, y)) +
    geom_errorbar(aes(xmin = fraction_low, xmax = fraction_high, colour = reproductive_strategy), orientation = 'y', width = .12, linewidth = .5) +
    geom_point(data = animals, aes(nuclear_fraction, animal_y, colour = reproductive_strategy), shape = 1, size = 1.45, stroke = .55) +
    geom_point(aes(colour = reproductive_strategy), shape = 18, size = 2) +
    strategy_scale() + guides(colour = 'none') +
    scale_y_continuous(breaks = allocation$y, labels = allocation$morphology_label,
      limits = c(.4, 24.6), expand = c(0, 0)) +
    scale_x_continuous(labels = label_percent(accuracy = 1), breaks = c(.12, .16, .20, .24),
      expand = expansion(mult = c(.03, .06))) +
    labs(title = 'Paired nuclear fractions', subtitle = '24 species; 42 animals; 1,152 paired cells',
      x = 'Nucleus / cell area', y = NULL,
      caption = 'Circles: animals; diamonds: equal-animal means.\nBars: 95% paired resampling intervals.') + style +
    theme(axis.text.y = element_text(face = 'italic', size = 7.5), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = 'none')
  if (isolated) p <- p +
    geom_text(aes(x = .255, label = support), hjust = .5, size = GBE_LABEL_SIZE_MM, family = GBE_FONT, colour = '#555555') +
    annotate('text', x = .255, y = 25.15, label = 'Animals / cells', size = GBE_LABEL_SIZE_MM, family = GBE_FONT) +
    scale_y_continuous(breaks = allocation$y, labels = allocation$morphology_label, limits = c(.4, 25.8), expand = c(0, 0)) +
    scale_x_continuous(labels = label_percent(accuracy = 1), breaks = c(.12, .16, .20, .24), limits = c(.11, .28), expand = c(0, 0)) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = 21, linetype = 0))) +
    theme(legend.position = 'bottom')
  p
}

legend5 <- paste(
  'Genome size and cellular architecture in all 24 sampled species. (a,b) Cell area and nucleus area versus genome size, with log10-standardized Pagel-lambda PGLS fits. β denotes a standardized coefficient; brackets give 95% t intervals; P values are unadjusted.',
  '(c) Conditional PGLS of log10 nucleus area on log10 cell area. Its unstandardized slope is compared with slope one; the dashed proportional reference passes through the log-scale species centroid. This regression does not identify an errors-in-variables structural exponent.',
  'Solid lines and gray bands in a–c show back-transformed fitted log-scale means and 95% coefficient bands; point bars show released measurement intervals. Measurement errors are displayed but not propagated into these regressions.',
  '(d) Fractions are first computed within each actual paired cell and nucleus, then summarized as animal medians (open circles) and the species arithmetic mean of animal medians (diamonds). Bars are the central 95% range of 1,000 animal-then-paired-cell resamples. All 42 animals and 1,152 pairs remain. Species follow the fixed phylogenetic display order. Seven morphology species have only one animal; their intervals omit among-animal uncertainty.',
  'All cellular traits describe selected large erythrocytes and two-dimensional nuclear and cell areas. Genome size is IOD-derived and D. fuscus assembly-anchored; calibration and staining-run uncertainty are not represented by the plotted intervals.',
  'Shared colors and shapes identify recorded reproductive strategy in a–c. Gray diamonds mark the three species with no frozen egg-placement assignment (D. brimleyorum, D. folkertsi, D. ochrophaeus); their aquatic-larval development is known and they remain in every fit. The overview labels a few fixed reference species; isolated scatterplots label all 24. Label choice changes no data or fits. The association results do not determine causal direction.')
fig5 <- ((genome_panel('cell_area_um2') | genome_panel('nucleus_area_um2')) /
  (allometry_panel() | allocation_panel())) + plot_layout(heights = c(.88, 1.25), guides = 'collect') +
  plot_annotation(tag_levels = 'a') & theme(legend.position = 'bottom')
save_candidate(5, 'Genome size and cellular architecture', fig5,
  list(a = list(plot = genome_panel('cell_area_um2', TRUE), title = 'Genome size and cell area', width_mm = 180, height_mm = 135, legend = legend5, alt_text = 'All 24 species labelled; cell area rises with genome size under the stored PGLS fit.'),
    b = list(plot = genome_panel('nucleus_area_um2', TRUE), title = 'Genome size and nucleus area', width_mm = 180, height_mm = 135, legend = legend5, alt_text = 'All 24 species labelled; nucleus area has a weaker positive genome-size association.'),
    c = list(plot = allometry_panel(TRUE), title = 'Conditional nuclear–cell area scaling', width_mm = 180, height_mm = 140, legend = legend5, alt_text = 'All 24 species labelled on logarithmic area axes with fitted slope 0.45 and a slope-one reference.'),
    d = list(plot = allocation_panel(TRUE), title = 'Paired nuclear fractions across species and animals', width_mm = 180, height_mm = 165, legend = legend5, alt_text = 'All 24 species in phylogenetic order show animal medians, species means, uncertainty intervals and animal/cell counts.')),
  180, 240, 'Frozen prior Figures5/10 analysis; S04;S18;S19;S24;S26;S42', legend5,
  'Four panels integrate genome–cell and genome–nucleus relationships, conditional nuclear–cell scaling and paired nuclear fractions, retaining all 24 species.')

# FIGURE 7: body-size associations plus weighting/cohort sensitivity.
body_data <- read_analysis(7, 'plotted_species_data.csv') %>% add_strategy()
body_models <- read_analysis(7, 'body_size_regressions.csv')
body_predictions <- read_analysis(7, 'fitted_relationships.csv')
forest <- read_analysis(9, 'displayed_coefficients.csv') %>%
  mutate(relationship = recode(relation_id, genome_to_cell_body_adjusted = 'Cell ~ genome + SVL',
    cell_to_nucleus_body_adjusted = 'Nucleus ~ cell + SVL',
    genome_and_cell_to_nucleus_body_adjusted = 'Nucleus ~ genome + cell + SVL'),
    relationship = factor(relationship, levels = c('Cell ~ genome + SVL', 'Nucleus ~ cell + SVL', 'Nucleus ~ genome + cell + SVL')),
    summary = if_else(aggregation == 'Released summaries', 'Released', 'Equal animal'),
    row = paste0(summary, ' · ', n_species),
    row = factor(row, levels = rev(c('Released · 24', 'Equal animal · 24', 'Released · 22', 'Equal animal · 22'))))
stopifnot(nrow(body_data) == 24, nrow(forest) == 12,
  setequal(forest$n_species, c(22L, 24L)), all(body_predictions$fit > 0))
copy_analysis(7, body_data, 'body_size_species.csv')
copy_analysis(7, body_models, 'frozen_body_size_models.csv')
copy_analysis(7, body_predictions, 'frozen_body_size_predictions.csv')
copy_analysis(7, forest, 'frozen_conditional_sensitivity_coefficients.csv')
body_panel <- function(id, response, ylab, title, all = FALSE) {
  m <- body_models %>% filter(relation_id == id, predictor == 'body_size') %>% arrange(desc(n_species))
  pred <- body_predictions %>% filter(relation_id == id)
  dat <- body_data %>% mutate(display_label = if (all) species else
    if_else(tree_tip %in% c('aeneus', 'wrighti', 'amphileucus'), species, ''))
  ggplot(dat, aes(max_svl_mm, .data[[response]])) +
    geom_ribbon(data = pred %>% filter(panel == 'full24'), aes(x = svl_mm, y = NULL, ymin = low, ymax = high), inherit.aes = FALSE, fill = '#E4E4E4') +
    geom_line(data = pred, aes(x = svl_mm, y = fit, linetype = panel), inherit.aes = FALSE, colour = '#222222', linewidth = .55) +
    geom_errorbar(aes(ymin = .data[[paste0(response, '_ci_low')]], ymax = .data[[paste0(response, '_ci_high')]]), width = 0, linewidth = .25, colour = '#ADADAD') +
    geom_point(aes(colour = reproductive_strategy, shape = reproductive_strategy), fill = 'white', size = 1.95, stroke = .75) +
    label_layer(dat, all, 20260927) + strategy_scale() + shape_scale() + strategy_guide +
    scale_linetype_manual(values = c(full24 = 'solid', subset22 = 'dashed'), name = NULL,
      labels = c(full24 = 'All 24 species', subset22 = '22 species; aeneus + wrighti excluded')) +
    guides(linetype = guide_legend(ncol = 1, order = 2), colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
      shape = guide_legend(nrow = 2, byrow = TRUE, order = 1)) +
    scale_x_log10(breaks = c(30, 50, 80, 120), expand = expansion(mult = c(.07, .09))) +
    scale_y_log10(expand = expansion(mult = c(.13, .13))) +
    labs(title = title, subtitle = paste(sprintf('n = %s: β = %.2f; P = %s', m$n_species, m$beta, format_p(m$p_value)), collapse = '\n'),
      x = 'Maximum SVL (mm; log scale)', y = ylab) + style
}
forest_panel <- function() ggplot(forest, aes(estimate, row)) +
  geom_vline(xintercept = 0, colour = '#AAAAAA', linewidth = .35, linetype = 'dashed') +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = 'y', width = .16, linewidth = .45) +
  geom_point(aes(shape = summary), fill = 'white', size = 1.9, stroke = .6) +
  facet_wrap(~relationship, ncol = 1) + scale_shape_manual(values = c(Released = 21, 'Equal animal' = 18), guide = 'none') +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1), limits = c(-.75, 1.3), expand = expansion(mult = .015)) +
  labs(title = 'Sampling and species sensitivity', subtitle = 'Focal predictor: genome, cell, genome\nWhiskers: 95% t intervals',
    x = 'Standardized coefficient', y = NULL) + style +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5), strip.text = element_text(size = 7.5, hjust = 0),
    panel.spacing.y = grid::unit(3, 'mm'), legend.position = 'none')
legend7 <- paste(
  'Body size and sampling sensitivity of the 24-species cellular panel. (a–c) Genome size, cell area and nucleus area versus species maximum snout–vent length (SVL). SVL values are species maxima supplied by the collaborator, not body lengths of the microscopy animals.',
  'Both axes are logarithmic with original-unit tick labels. Solid lines are the full 24-species Pagel-lambda PGLS fits; dashed lines are fits excluding D. aeneus and D. wrighti (22 species), restricted to the observed body-size range of each fitted panel. The gray band is the full 24-species 95% confidence band for the log-scale fitted mean, back-transformed; it is not a prediction interval. Vertical bars are released trait measurement intervals.',
  'β values and the conditional coefficients in d use identical full 24-species log10 means and SDs. Species, not cells, are the model observations. Displayed P values are unadjusted. Neither SVL uncertainty nor predictor measurement error is propagated.',
  '(d) Three conditional regressions show the coefficient on genome size in cell ~ genome + SVL, cell area in nucleus ~ cell + SVL, and genome size in nucleus ~ genome + cell + SVL. Four rows compare released summaries and equal-animal summaries in 24 and 22 species. Whiskers are 95% t coefficient intervals conditional on the fitted phylogenetic parameter.',
  'Released morphology pools selected paired objects within species before taking the median. The equal-animal sensitivity takes the median of within-animal medians for morphology; IOD first averages image medians within animal, then averages animals. This median-of-medians morphology definition differs from the arithmetic mean-of-animal-medians nuclear-fraction summaries in Figure5d.',
  'The 22-species sensitivity removes exactly the two smallest species and only direct developers; their body size, developmental mode and lineage identity are confounded, so altered significance does not identify a causal mechanism or test differences between slopes.',
  'Shared strategy colors retain all 24 points, with three unassigned egg-placement records shown gray. Isolated scatterplots label all species. Cell and nuclear areas refer to selected large erythrocytes; genome size is IOD-derived and D. fuscus assembly-anchored. Seven morphology species and four IOD species have one animal. These sensitivity analyses do not resolve staining calibration, causal direction, or population representativeness.')
fig7 <- ((body_panel('body_to_genome', 'genome_size_pg_fuscus_anchored', 'Genome size (pg/1C; log scale)', 'Body size and genome size') |
  body_panel('body_to_cell', 'cell_area_um2', 'Cell area (µm²; log scale)', 'Body size and cell area')) /
  (body_panel('body_to_nucleus', 'nucleus_area_um2', 'Nucleus area (µm²; log scale)', 'Body size and nucleus area') | forest_panel())) +
  plot_layout(guides = 'collect', heights = c(1, 1.15)) + plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom', legend.box = 'vertical')
save_candidate(7, 'Body size and sampling sensitivity', fig7,
  list(a = list(plot = body_panel('body_to_genome', 'genome_size_pg_fuscus_anchored', 'Genome size (pg/1C; log scale)', 'Body size and genome size', TRUE), title = 'Body size and genome size', width_mm = 180, height_mm = 145, legend = legend7, alt_text = 'All 24 species labeled by species maximum SVL and genome size, with full 24-species and subset 22-species fitted lines.'),
    b = list(plot = body_panel('body_to_cell', 'cell_area_um2', 'Cell area (µm²; log scale)', 'Body size and cell area', TRUE), title = 'Body size and cell area', width_mm = 180, height_mm = 145, legend = legend7, alt_text = 'All 24 species labeled by maximum SVL and selected cell area, with two fitted lines.'),
    c = list(plot = body_panel('body_to_nucleus', 'nucleus_area_um2', 'Nucleus area (µm²; log scale)', 'Body size and nucleus area', TRUE), title = 'Body size and nucleus area', width_mm = 180, height_mm = 145, legend = legend7, alt_text = 'All 24 species labeled by maximum SVL and selected nucleus area, with full 24-species and subset 22-species fits.'),
    d = list(plot = forest_panel(), title = 'Conditional coefficient sensitivity', width_mm = 180, height_mm = 145, legend = legend7, alt_text = 'Twelve coefficients compare three conditional relationships under released or equal-animal summaries and 24 or 22 species.')),
  180, 235, 'Frozen prior Figures7/9 analysis; S04;S19;S42; collaborator maximum SVL', legend7,
  'Three body-size scatterplots retain every species and compare 24-species with 22-species fits; a coefficient forest adds animal-weighting and species-panel sensitivity.')

write_lines(c('PASS: 24 species retained in each scatterplot and in nuclear-fraction panel.',
  'PASS: 42 animals and 1152 paired cells represented in frozen nuclear-fraction summaries.',
  'PASS: three unassigned reproductive-strategy records retained explicitly.',
  'PASS: all model coefficients, predictions and intervals read from frozen prior exports; no models refitted.',
  'PASS:12 conditional coefficients; exactly2 aggregations x2 species panels x3 relationships.',
  'Arial target; minimum plot text7.5pt; collage widths180mm, heights240mm and235mm.'),
  file.path(RUN, 'qa', 'figures_5_7_checks.txt'))
message('Figures5 and7 polished without refitting any model.')

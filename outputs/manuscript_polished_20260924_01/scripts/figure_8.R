# Render the frozen exact21 integration; no models are fitted here.
source('scripts/common.R')
d <- read_analysis(11, 'integrated_species21.csv')
pred <- read_analysis(11, 'genome_predictions.csv') %>% filter(endpoint == 'within_ltr')
coef <- read_analysis(11, 'regression_coefficients.csv') %>% filter(endpoint == 'within_ltr')
matched <- read_analysis(11, 'cellular_coefficients_matched_comparison.csv')
stopifnot(nrow(d) == 21, !anyDuplicated(d$species), nrow(matched) == 9)
d <- d %>% mutate(strategy = factor(strategy, levels = te_strategy_levels),
                  label = species_label(species))
out <- file.path(RUN, 'figures/figure_8')
dir.create(file.path(out, 'plotted_data'), recursive = TRUE, showWarnings = FALSE)
write_csv(d, file.path(out, 'plotted_data/species21.csv'))
write_csv(pred, file.path(out, 'plotted_data/predictions.csv'))
write_csv(matched, file.path(out, 'plotted_data/matched_coefficients.csv'))
selected <- c('aeneus', 'wrighti', 'fuscus', 'auriculatus', 'amphileucus')
scatter <- function(all_labels = FALSE) {
  labels <- filter(d, all_labels | species %in% selected)
  p <- ggplot(d, aes(within_ltr_fraction, genome_size_pg_fuscus_anchored)) +
    geom_ribbon(data = filter(pred, fit_id == 'genome_body_adjusted'),
      aes(x = te_fraction, ymin = ci_low, ymax = ci_high), inherit.aes = FALSE,
      fill = '#E2E2E2') +
    geom_line(data = pred, aes(x = te_fraction, y = fitted, linetype = fit_id),
      inherit.aes = FALSE, colour = '#333333', linewidth = .55) +
    geom_errorbar(aes(ymin = genome_size_pg_ci_low, ymax = genome_size_pg_ci_high),
      width = 0, colour = '#999999', linewidth = .3) +
    geom_point(aes(colour = strategy, shape = strategy), fill = 'white', size = 2.2, stroke = .7) +
    geom_label_repel(data = labels, aes(label = label),
      family = GBE_FONT, fontface = 'italic', size = 7.5 / ggplot2::.pt,
      fill = 'white', linetype = 0, label.size = NA, label.padding = grid::unit(.07, 'lines'),
      label.r = grid::unit(0, 'mm'),
      seed = 240908, box.padding = .35, point.padding = .25, min.segment.length = 0,
      nudge_y = ifelse(labels$species %in% c('auriculatus', 'amphileucus', 'wrighti'), 2.4, 0),
      max.overlaps = Inf, max.iter = 30000, max.time = Inf, segment.size = .2,
      segment.color = '#888888') +
    strategy_colours() + strategy_shapes() +
    scale_linetype_manual(name = NULL, values = c(genome_body_adjusted = 'solid', genome_unadjusted = 'dashed'),
      breaks = c('genome_body_adjusted', 'genome_unadjusted'), labels = c('SVL-adjusted', 'Unadjusted')) +
    scale_x_continuous(labels = label_percent(accuracy = 1), breaks = seq(.02, .12, .02),
      expand = expansion(mult = c(.06, .06))) +
    scale_y_continuous(breaks = seq(5, 30, 5), expand = expansion(mult = c(.03, .10))) +
    labs(title = 'Repeat divergence and genome size',
      x = 'Low-divergence LTR fraction', y = 'Genome size (pg/1C)') +
    theme_gbe() + theme(legend.position = 'bottom', legend.box = 'vertical') +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, order = 1),
      shape = guide_legend(nrow = 2, byrow = TRUE, order = 1), linetype = guide_legend(order = 2))
  p
}
te <- coef %>% filter(term == 't', fit_id %in% c('genome_body_adjusted', 'cell_te_adjusted', 'nucleus_te_adjusted')) %>%
  mutate(row = match(fit_id, c('nucleus_te_adjusted','cell_te_adjusted','genome_body_adjusted')))
stopifnot(nrow(te) == 3, max(abs(te$estimate[c(3,2,1)])) < 2)
write_csv(te, file.path(out, 'plotted_data/te_coefficients.csv'))
pb <- ggplot(te, aes(estimate, row)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = '#999999', linewidth = .3) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = row), linewidth = .6) +
  geom_point(shape = 21, fill = 'white', size = 2.5, stroke = .7) +
  scale_y_continuous(breaks = 1:3, labels = c('Nucleus area\n| Cell, genome, SVL',
    'Cell area\n| Genome, SVL', 'Genome size\n| SVL'), limits = c(.5, 3.5)) +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1), limits = c(-.6, 1.1)) +
  labs(title = 'TE association with each trait', x = 'Standardized TE coefficient', y = NULL) +
  theme_gbe() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7.5), plot.title = element_text(size = 9))
models <- c('24 species: SVL', '21 species: SVL', '21 species: SVL + TE')
link_order <- c('Genome to cell', 'Cell to nucleus', 'Genome to nucleus')
matched <- matched %>% mutate(model = factor(model, levels = models), link = factor(link, levels = link_order),
  row = 4 - as.integer(model), panel = factor(as.character(link), levels = link_order,
    labels = c('Genome–cell\n| SVL', 'Cell–nucleus\n| Genome, SVL', 'Genome–nucleus\n| Cell, SVL')))
pc <- ggplot(matched, aes(estimate, row, shape = model, colour = model)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = '#999999', linewidth = .3) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = row), linewidth = .6) +
  geom_point(size = 2.3, stroke = .7) +
  facet_wrap(~panel, nrow = 1) +
  scale_shape_manual(name = NULL, values = c(16, 1, 18), labels = c('24 species + SVL', '21 species + SVL', '21 species + SVL + TE')) +
  scale_colour_manual(name = NULL, values = c('#888888', '#222222', '#0072B2'), labels = c('24 species + SVL', '21 species + SVL', '21 species + SVL + TE')) +
  scale_x_continuous(breaks = c(-.5, 0, .5, 1), limits = c(-.65, 1.4)) +
  scale_y_continuous(limits = c(.5, 3.5), breaks = NULL) +
  labs(title = 'Cellular coefficients change little with TE adjustment', x = 'Standardized coefficient [95% interval]', y = NULL) +
  theme_gbe() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.spacing.x = grid::unit(4, 'mm'), legend.position = 'bottom') +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
legend <- paste(
  'The measured LTR divergence signature does not account for the genome–cell association in the fitted models.',
  '(a) All 21 species with complete repeat, cellular and maximum-SVL data are shown. Five species are labelled in the collage; all 21 are named in the isolated panel.',
  'The LTR score is aligned-hit bp in divergence bins 0–4 divided by all LTR aligned-hit bp; it is not genomic TE load, an insertion rate or calibrated age.',
  'Solid and dashed lines are the frozen SVL-adjusted and marginal Pagel-lambda PGLS predictions. The adjusted line holds SVL at the 21-species geometric mean (62.13 mm); its gray band is a conditional 95% coefficient interval.',
  'Vertical bars show genome-size measurement intervals, which are not propagated into these fits.',
  '(b) TE-score coefficients for genome, cell and nuclear outcomes. The vertical bar means adjusted for the listed covariates. The genome coefficient is 0.311 [−0.381, 1.003], P = 0.358.',
  '(c) Existing 24-species, matched 21-species and TE-adjusted 21-species coefficients use the same original cellular scaling. Both nuclear equations include genome and cell area; SVL appears in every model and TE is added only to the indicated series.',
  'The matched genome–cell coefficient changes from 0.549 to 0.527 with TE adjustment. Wider intervals do not establish absent effects or mechanistic independence.',
  'Intervals are t intervals conditional on fitted lambda. Cellular traits are standardized log10 values using the full24 means and SDs; the logit TE fraction uses full34 scaling.',
  'D. brimleyorum, D. folkertsi and D. ochrophaeus lack the TE endpoint and are not imputed. Genome sizes are nuclear-IOD-derived and D. fuscus anchored; unresolved staining calibration, selected large-cell sampling and tree provenance remain limitations.',
  'Shapes and colors identify reproductive strategy in (a), and species-set/model adjustment in (c). This exploratory comparison does not identify a TE-to-genome-to-cell causal sequence.')
pa <- scatter()
collage <- (pa / pb / pc) + plot_layout(heights = c(1.25, .68, .75)) +
  plot_annotation(tag_levels = 'a')
panels <- list(
  a = list(plot = scatter(TRUE), title = 'LTR divergence and genome size', width_mm = 180, height_mm = 130,
    legend = legend, alt_text = 'All 21 species are labelled around marginal and body-adjusted genome-size regressions on within-LTR low-divergence fraction.'),
  b = list(plot = pb, title = 'Conditional associations of the TE score', width_mm = 180, height_mm = 85,
    legend = legend, alt_text = 'Three TE coefficients have broad intervals crossing zero for genome, cell and nuclear outcomes.'),
  c = list(plot = pc, title = 'Matched-species cellular coefficients', width_mm = 180, height_mm = 95,
    legend = legend, alt_text = 'Three sets of coefficients show little change in the genome–cell association when the measured TE score is included.'))
save_candidate(8, 'Repeat divergence and cellular integration', collage, panels, 180, 235,
  'Frozen original Figure 11: integrated_species21; genome_predictions; regression_coefficients; cellular_coefficients_matched_comparison',
  legend, 'A 21-species scatterplot and conditional-coefficient panels show uncertain TE associations and a persistent genome–cell relationship after TE adjustment.')
write_lines(c('# Figure 8 refinement', '',
  'Consolidates original Figure 11; numerical inputs and predictions are frozen, with no refitting.',
  'Uses exact21 species, shared strategy encodings, selected collage labels plus fully labelled isolated scatter.',
  'Separates the measured TE association from the matched-sample cellular coefficient comparison.',
  'Source rows and predictions are exported in plotted_data. All observations and intervals are retained.'), file.path(out, 'REFINEMENT_NOTES.md'))

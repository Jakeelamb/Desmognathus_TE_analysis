#!/usr/bin/env Rscript
# Presentation-only refinement. All observations, inference and predictions are frozen.
source('scripts/common.R')
selected <- commandArgs(trailingOnly = TRUE)
if (!length(selected)) selected <- as.character(1:4)
blue <- GBE_COLORS[['blue']]
grey <- '#B8B8B8'
label_size <- 7.5 / ggplot2::.pt
anno <- read_supp('Supplementary_Data_S42_te_bonsai_tip_annotations.csv')
strategy_fill <- function() scale_fill_manual(values = te_strategy_colors, breaks = te_strategy_levels,
  labels = te_strategy_labels, name = NULL, drop = FALSE)
strategy_shape <- function() scale_shape_manual(values = te_strategy_shapes, breaks = te_strategy_levels,
  labels = te_strategy_labels, name = NULL, drop = FALSE)
key_theme <- theme(legend.position = 'bottom', legend.box = 'vertical',
  legend.key.height = unit(9, 'pt'), legend.key.width = unit(16, 'pt'))
row_theme <- theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
panel <- function(p, title, width, height, legend, alt) list(plot = p, title = title,
  width_mm = width, height_mm = height, legend = legend, alt_text = alt)
write_notes <- function(n, text) write_lines(c(paste0('# Figure ', n, ' refinement'), '', text,
  '', 'Display only: no model refitting, additional filtering, species exclusions, or inferred data.',
  'Arial; minimum planned text size 7.5 pt. All main collages 180 mm wide and no taller than 245 mm.',
  'Original source figures and notes are retained in the earlier galleries and frozen input snapshot.'),
  file.path(RUN, 'figures', paste0('figure_', n), 'REFINEMENT_NOTES.md'))
copy_exports <- function(n, source_n, folder = 'analysis') {
  destination <- file.path(RUN, 'figures', paste0('figure_', n), paste0('source_figure_', source_n))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  f <- list.files(file.path(SOURCE, 'figures', paste0('figure_', source_n), folder),
    full.names = TRUE, pattern = '\\.(csv|nwk|tre|json)$')
  stopifnot(all(file.copy(f, destination, overwrite = TRUE)))
}
genome_boundary <- paste('Genome size is based on nuclear integrated optical density (IOD), anchored to the D. fuscus assembly scale.',
  'Measurement intervals exclude anchor uncertainty; run-specific calibration remains unresolved. Sizes are species summaries, not independently validated C-values.')

if ('1' %in% selected) {
  tracks <- read_analysis(1, 'species_tracks.csv')
  horizontal <- read_analysis(1, 'tree_horizontal_segments.csv')
  vertical <- read_analysis(1, 'tree_vertical_segments.csv')
  composition <- read_analysis(1, 'broad_composition.csv') %>%
    mutate(order = factor(order, levels = c('LTR', 'LINE', 'TIR', 'Other')))
  stopifnot(nrow(tracks) == 34, sum(tracks$genome_size_available) == 21,
    identical(tracks$species[order(tracks$y)], phylogenetic_species_order(tracks$species)))
  a <- ggplot() +
    geom_segment(data = horizontal, aes(x, y, xend = xend, yend = yend), linewidth = .28) +
    geom_segment(data = vertical, aes(x, ymin, xend = x, yend = ymax), linewidth = .28) +
    geom_point(data = tracks, aes(1, y, fill = reproductive_strategy, shape = reproductive_strategy),
      colour = '#222222', size = 1.9, stroke = .3) +
    geom_text(data = tracks, aes(1.13, y, label = display_species), hjust = 0,
      family = GBE_FONT, fontface = 'italic', size = label_size) +
    strategy_fill() + strategy_shape() +
    scale_x_continuous(limits = c(-.02, 2.95), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 34.5), expand = c(0, 0)) +
    labs(title = 'Species phylogeny', x = NULL, y = NULL) + theme_gbe() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()) + key_theme +
    guides(fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
  # Numeric annotations sit in a dedicated right column instead of competing with intervals.
  b <- ggplot(tracks, aes(genome_size_pg_fuscus_anchored, y)) +
    geom_segment(data = filter(tracks, genome_size_available),
      aes(x = genome_size_pg_ci_low, xend = genome_size_pg_ci_high, yend = y), linewidth = .45, colour = grey) +
    geom_point(data = filter(tracks, genome_size_available), shape = 21, fill = 'white', size = 1.65, stroke = .4) +
    geom_text(aes(x = 29.3, label = genome_label), family = GBE_FONT, size = label_size, hjust = 1) +
    scale_x_continuous(limits = c(8, 30), breaks = c(10, 18, 26), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 34.5), breaks = NULL, expand = c(0, 0)) +
    labs(title = 'Genome size', x = 'pg/1C', y = NULL) + theme_gbe() + row_theme
  c <- ggplot(composition) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = y - .34, ymax = y + .34, fill = order)) +
    scale_fill_manual(values = GBE_TE_ORDER_COLORS, breaks = c('LTR', 'LINE', 'TIR', 'Other'), name = NULL) +
    scale_x_continuous(breaks = c(0, .5, 1), labels = label_percent(accuracy = 1), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 1)) +
    scale_y_continuous(limits = c(.5, 34.5), breaks = NULL, expand = c(0, 0)) +
    labs(title = 'Repeat composition', x = 'Classified aligned bases', y = NULL) + theme_gbe() + row_theme + key_theme +
    guides(fill = guide_legend(nrow = 1))
  # Some frozen cumulative sums exceed 1 by floating-point rounding. Coordinate
  # clipping preserves those complete rectangles; scale limits would censor them.
  stopifnot(nrow(ggplot_build(c)$data[[1]]) == nrow(composition),
    !anyNA(ggplot_build(c)$data[[1]]$xmax))
  stand_rows <- function(p) p + scale_y_continuous(limits = c(.5, 34.5), breaks = tracks$y,
    labels = tracks$display_species, expand = c(0, 0)) + theme(axis.text.y = element_text(face = 'italic', size = 7.5))
  legend <- paste('(a) Topology-only source species tree for 34 Desmognathus species; branch lengths are omitted.',
    'Reproductive strategies are aquatic eggs (7 species), terrestrial eggs with aquatic larvae (24), and direct development (3).',
    '(b) Twenty-one species have genome sizes, numeric labels and conditional 95% measurement-bootstrap intervals; n/a marks 13 unavailable estimates.',
    '(c) dnaPipeTE classified aligned bases are grouped into LTR, LINE, TIR and remaining classified orders (Other). Each row sums to 100% of this denominator.',
    'Unresolved sequence is outside the denominator; these proportions are not whole-genome occupancy. All three panels share the same phylogenetic row order.',
    genome_boundary, 'Source-tree publication and calibration provenance require resolution before submission.')
  collage <- (a + b + c) + plot_layout(widths = c(1.43, .88, 1.03), guides = 'collect') +
    plot_annotation(tag_levels = 'a') & key_theme
  save_candidate(1, 'Repeat architecture across Desmognathus', collage,
    list(a = panel(a, 'Species phylogeny', 115, 220, legend, 'Topology-only tree with 34 named species and reproductive strategies.'),
      b = panel(stand_rows(b), 'Genome sizes', 120, 205, legend, 'Twenty-one species genome sizes, numeric labels, intervals, and 13 missing values.'),
      c = panel(stand_rows(c) + theme(plot.margin = margin(5, 16, 5, 6)), 'Repeat composition', 140, 210, legend, 'Classified LTR, LINE, TIR and other repeat shares for all 34 species.')),
    180, 225, 'Frozen Figure 1 species_tracks, tree segments and broad_composition; S40; S42', legend,
    'A fixed 34-species row order aligns phylogeny, genome sizes, and LTR-rich repeat composition.')
  copy_exports(1, 1)
  write_notes(1, 'Retains prior Figure 1a–c, with a dedicated numeric genome-size column, consistent circular genome-size points, simpler titles, and compact aligned rows. Software and denominator are explicit in the caption. No heatmap is used.')
}

if ('2' %in% selected) {
  profiles <- read_analysis(2, 'species_divergence_profiles.csv') %>% mutate(order = factor(order, levels = c('LTR', 'LINE', 'TIR')))
  means <- read_analysis(2, 'strategy_mean_profiles.csv') %>% mutate(order = factor(order, levels = c('LTR', 'LINE', 'TIR')))
  scores <- read_analysis(12, 'repeat_score_species34.csv') %>%
    mutate(display_species = factor(display_species, levels = phylogenetic_species_levels(species)),
      strategy = factor(strategy, levels = te_strategy_levels))
  preds <- read_analysis(12, 'strategy_predictions.csv') %>% filter(endpoint == 'within_ltr')
  contrasts <- read_analysis(12, 'plotted_strategy_contrasts.csv') %>%
    mutate(contrast = factor(contrast, levels = c('Direct development', 'Terrestrial eggs + aquatic larvae')),
      model = factor(model, levels = c('Unadjusted', 'SVL adjusted')))
  omnibus <- read_analysis(12, 'strategy_omnibus_tests.csv') %>% filter(endpoint == 'within_ltr', svl_adjusted)
  stopifnot(n_distinct(profiles$species) == 34, nrow(scores) == 34, nrow(contrasts) == 4,
    all(abs(scores$within_ltr_fraction * scores$ltr_fraction_all_hits - scores$low_ltr_fraction_all_hits) < 1e-14))
  line_types <- setNames(c('solid', 'longdash', 'dotdash'), te_strategy_levels)
  a <- ggplot(profiles, aes(divergence_bin_start_pct, percent_species_hit_bp, group = species)) +
    geom_line(colour = '#B9B9B9', linewidth = .20, alpha = .48) +
    geom_line(data = means, aes(y = mean_percent, group = reproductive_strategy, colour = reproductive_strategy,
      linetype = reproductive_strategy), linewidth = .67) +
    facet_grid(cols = vars(order)) + strategy_colours() +
    scale_linetype_manual(values = line_types, breaks = te_strategy_levels, labels = te_strategy_labels, name = NULL) +
    scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 3), breaks = 0:3, expand = expansion(mult = c(0, .025))) +
    labs(title = 'Repeat-divergence profiles', x = 'Sequence divergence (%)', y = 'All hit bp per bin (%)') +
    theme_gbe() + theme(panel.spacing = unit(14, 'pt'), legend.position = 'none')
  b <- ggplot(scores, aes(within_ltr_fraction * 100, display_species)) +
    geom_segment(aes(x = 0, xend = within_ltr_fraction * 100, yend = display_species), colour = '#DDDDDD', linewidth = .3) +
    geom_point(aes(fill = strategy, shape = strategy), size = 1.8, stroke = .35) + strategy_fill() + strategy_shape() +
    scale_x_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15), expand = expansion(mult = c(.01, .025))) +
    labs(title = 'Within-LTR low divergence', x = 'LTR hit bp below 5% divergence (%)', y = NULL) +
    theme_gbe() + row_theme + theme(axis.text.y = element_text(face = 'italic', size = 7.5), legend.position = 'none')
  c <- ggplot(scores, aes(max_svl_mm, within_ltr_fraction * 100)) +
    geom_ribbon(data = preds, aes(x = max_svl_mm, y = NULL, ymin = ci_low * 100, ymax = ci_high * 100, fill = strategy),
      inherit.aes = FALSE, alpha = .13) +
    geom_line(data = preds, aes(x = max_svl_mm, y = fitted * 100, colour = strategy), inherit.aes = FALSE, linewidth = .6) +
    geom_point(aes(fill = strategy, shape = strategy), size = 1.9, stroke = .35) + strategy_fill() + strategy_shape() + strategy_colours() +
    scale_x_log10(breaks = c(30, 50, 80, 120)) +
    scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    labs(title = 'Body size and strategy', x = 'Maximum SVL (mm; log scale)', y = 'Within-LTR fraction (%)') +
    theme_gbe() + theme(legend.position = 'none')
  # Vertical contrasts keep both group labels readable at a half-page width.
  d <- ggplot(contrasts, aes(contrast, estimate, colour = model, shape = model)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = '#999999', linewidth = .35) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .12, linewidth = .45, position = position_dodge(.3)) +
    geom_point(size = 2, position = position_dodge(.3)) +
    scale_colour_manual(values = c('Unadjusted' = '#999999', 'SVL adjusted' = '#222222'), name = NULL) +
    scale_shape_manual(values = c('Unadjusted' = 1, 'SVL adjusted' = 16), name = NULL) +
    scale_x_discrete(labels = c('Direct development' = 'Direct\ndevelopment', 'Terrestrial eggs + aquatic larvae' = 'Terrestrial eggs +\naquatic larvae')) +
    labs(title = 'Strategy contrasts', x = NULL, y = 'Logit difference vs aquatic eggs') + theme_gbe() + key_theme +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  # Shared strategy key uses points and line styles, and carries group counts.
  keydata <- tibble(strategy = te_strategy_levels, y = c(3, 2, 1),
    label = c('Aquatic eggs (n = 7)', 'Terrestrial eggs + aquatic larvae (n = 24)', 'Direct development (n = 3)'))
  key <- ggplot(keydata) + geom_segment(aes(x = 0, xend = .07, y = y, yend = y, colour = strategy, linetype = strategy), linewidth = .6) +
    geom_point(aes(x = .035, y = y, fill = strategy, shape = strategy), size = 1.8, stroke = .35) +
    geom_text(aes(.085, y, label = label), hjust = 0, family = GBE_FONT, size = label_size) +
    strategy_fill() + strategy_shape() + strategy_colours() + scale_linetype_manual(values = line_types) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(.5, 3.5), expand = c(0, 0)) +
    theme_void(base_family = GBE_FONT) + theme(legend.position = 'none', plot.margin = margin(0, 5, 0, 5))
  legend <- paste('(a) All 34 species are shown as gray divergence profiles; colored curves are equal-weight reproductive-group means.',
    'Each 1-percentage-point bin contains its percentage of all RepeatMasker-aligned hit bp; bins 0–49 are displayed with missing order/bin combinations zero-filled.',
    '(b) All 34 named species retain the shared phylogenetic order. The within-LTR fraction is LTR hit bp in bins 0–4 divided by all LTR hit bp across all divergence bins.',
    'It separates relative profile shape from total LTR representation; it is not genome occupancy, insertion age, or transposition rate.',
    '(c) Stored Pagel-lambda PGLS predictions from logit(within-LTR fraction) ~ strategy + log10(maximum SVL), restricted to each group\'s observed SVL range; bands are conditional 95% intervals.',
    sprintf('(d) Stored raw-logit contrasts and 95%% t intervals compare each group with aquatic eggs, before and after SVL adjustment. The adjusted omnibus likelihood-ratio P is %.6f (2 df).', omnibus$p_value),
    'The adjusted terrestrial-versus-aquatic contrast is 0.439 (95% CI 0.061–0.818); the direct-development contrast is −0.448 (−1.037–0.142).',
    'Strategies were expert-coded after PCA inspection; group sizes are not counts of independent evolutionary transitions. Contrasts remain exploratory and conditional on the source tree and fitted lambda.',
    'The within-LTR fraction covaries with dnaPipeTE unresolved fraction (unadjusted Spearman rho = −0.452, P = 0.0072), but not clearly with aligned sample fraction (rho = −0.161, P = 0.362). These cross-pipeline quality associations do not establish the origin of the group pattern.')
  collage <- ((a + labs(tag = 'a')) / key / ((b + labs(tag = 'b')) | ((c + labs(tag = 'c')) / (d + labs(tag = 'd'))))) +
    plot_layout(heights = c(.92, .30, 2.26))
  save_candidate(2, 'LTR divergence, body size and reproductive strategy', collage,
    list(a = panel(a / key + plot_layout(heights = c(1, .26)), 'Repeat-divergence profiles', 180, 100, legend, 'Gray curves show all 34 species; colored means compare three strategies across LTR, LINE and TIR.'),
      b = panel(b / key + plot_layout(heights = c(1, .15)), 'Within-LTR fraction by species', 135, 225, legend, 'All 34 species are labelled in a dot plot of low-divergence hit bp within LTRs.'),
      c = panel(c / key + plot_layout(heights = c(1, .22)), 'Body size and reproductive strategy', 160, 110, legend, 'Thirty-four species and frozen group PGLS curves on maximum body length.'),
      d = panel(d, 'Phylogenetic strategy contrasts', 150, 100, legend, 'Unadjusted and body-adjusted contrasts and intervals relative to aquatic eggs.')),
    180, 245, 'Frozen Figure 2 profile exports; Figure 12 scores, predictions, contrasts and omnibus tests', legend,
    'All-species repeat profiles and named within-LTR fractions accompany body-size predictions and phylogenetic strategy contrasts.')
  copy_exports(2, 2); copy_exports(2, 12)
  write_notes(2, 'Integrates prior Figure 2a with Figure 12a–c. All 34 species remain visible in the divergence overlay and named within-LTR track. The 102-point three-order low-divergence track and full three-denominator decomposition remain available in the earlier gallery. Primary within-LTR and original all-hit denominators are explicitly distinguished; the QC association is retained in the caption.')
}

if ('3' %in% selected) {
  diversity <- read_analysis(3, 'plotted_species_data.csv')
  predictions <- read_analysis(3, 'fitted_relationships.csv')
  stats <- read_analysis(3, 'model_summary.csv')
  stopifnot(nrow(diversity) == 21, nrow(stats) == 2, all(diversity$n_superfamilies_included == 10))
  make_diversity <- function(metric, title, label_all = FALSE) {
    pr <- filter(predictions, response == metric)
    p <- ggplot(diversity, aes(genome_size_pg_fuscus_anchored, .data[[metric]])) +
      geom_ribbon(data = pr, aes(y = NULL, ymin = low, ymax = high), fill = '#DDDDDD', alpha = .65) +
      geom_line(data = pr, aes(y = fit), colour = '#222222', linewidth = .6) +
      geom_segment(aes(x = genome_size_pg_ci_low, xend = genome_size_pg_ci_high, yend = .data[[metric]]), colour = grey, linewidth = .35) +
      geom_point(aes(fill = reproductive_strategy, shape = reproductive_strategy), colour = '#222222', size = 2, stroke = .4) +
      strategy_fill() + strategy_shape() +
      scale_x_continuous(limits = c(8, 27), breaks = c(10, 15, 20, 25), expand = expansion(mult = .015)) +
      scale_y_continuous(expand = expansion(mult = if (label_all) c(.20, .20) else c(.07, .07))) +
      labs(title = title, x = 'Genome size (pg/1C)', y = title) + theme_gbe() + key_theme +
      guides(fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    if (!label_all) {
      label_data <- filter(diversity, species %in% c('aeneus', 'wrighti', 'valtos'))
      p <- p + geom_text_repel(data = label_data, aes(label = species_label), family = GBE_FONT, fontface = 'italic',
        size = label_size, seed = 20260924, max.overlaps = Inf, box.padding = .2, point.padding = .2,
        nudge_x = ifelse(label_data$species == 'aeneus', 1.7, 0),
        nudge_y = diff(range(diversity[[metric]])) * .075,
        min.segment.length = 0, segment.color = '#999999', segment.size = .2, max.time = Inf, max.iter = 30000)
    }
    if (label_all) {
      left <- diversity$genome_size_pg_fuscus_anchored < median(diversity$genome_size_pg_fuscus_anchored)
      p <- p + geom_label_repel(aes(label = species_label), family = GBE_FONT, fontface = 'italic',
        fill = 'white', linewidth = 0, label.padding = unit(.06, 'lines'),
        size = label_size, seed = 20260924, max.overlaps = Inf, box.padding = .38, point.padding = .25,
        nudge_x = ifelse(left, 12.0, 20.5) - diversity$genome_size_pg_fuscus_anchored,
        hjust = ifelse(left, 1, 0), direction = 'y', force = 10, force_pull = .01,
        min.segment.length = 0, segment.color = '#999999', segment.size = .2,
        max.time = Inf, max.iter = 80000)
    }
    p
  }
  a <- make_diversity('shannon_entropy', 'Shannon entropy')
  b <- make_diversity('gini_simpson', 'Gini–Simpson index')
  legend <- paste('Twenty-one species with both genome-size and repeat-diversity data are shown in each panel.',
    '(a) Shannon entropy and (b) Gini–Simpson index summarize each species\' ten most abundant positive classified dnaPipeTE superfamilies, reclosed to sum to one.',
    'These related measures represent diversity/evenness of the selected superfamily distribution, not total TE burden or independent replications.',
    'Points identify reproductive strategy; horizontal bars show genome-size measurement intervals. All observations contribute to the same full-sample fit per panel.',
    'Lines and bands reuse stored Pagel-lambda PGLS fitted means and 95% coefficient intervals without refitting.',
    'The Shannon slope is −0.007558 per pg/1C (95% CI −0.015234–0.000118; P = 0.0533), and the Gini–Simpson slope is −0.003623 (−0.006964–−0.000283; P = 0.0350). Both fitted lambdas are 1.',
    'P values are unadjusted. Predictor measurement error and tree uncertainty are not propagated. D. aeneus, D. wrighti and D. valtos are labelled in the collage to identify the genome-size extremes and the high-diversity direct developer; label selection does not change the 21-species sample. All 21 species are individually named in isolated panels.', genome_boundary)
  collage <- (a | b) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'a') & key_theme
  save_candidate(3, 'Repeat diversity and genome size', collage,
    list(a = panel(make_diversity('shannon_entropy', 'Shannon entropy', TRUE), 'Shannon entropy and genome size', 180, 125, legend, 'Twenty-one named species with genome-size intervals and a negative fitted Shannon relationship.'),
      b = panel(make_diversity('gini_simpson', 'Gini–Simpson index', TRUE), 'Gini–Simpson index and genome size', 180, 125, legend, 'The same 21 named species and intervals around a negative fitted Gini–Simpson relationship.')),
    180, 110, 'Frozen Figure 3 plotted_species_data, model_summary and fitted_relationships', legend,
    'Two compact 21-species plots show repeat diversity against genome size with unchanged full-sample PGLS fits and intervals.')
  copy_exports(3, 3)
  write_notes(3, 'Retains both prior Figure 3 metrics and all 21 species; replaces tall repeated panels with a compact side-by-side manuscript collage. The collage labels aeneus, wrighti and valtos; the expanded isolated panels name all species using readable two-column leader labels. Label selection does not select observations or alter fits. All fitted predictions and model summaries are copied unchanged. Full inference is in the caption, freeing plotting space.')
}

if ('4' %in% selected) {
  elements <- read_analysis(4, 'element_ratios.csv', 'plotted_data')
  summaries <- read_analysis(4, 'species_distributions_and_annotations.csv', 'plotted_data')
  pgls <- read_analysis(4, 'pgls_points_and_label_layout.csv', 'plotted_data') %>%
    left_join(anno %>% select(species, reproductive_strategy), by = 'species', relationship = 'one-to-one')
  prediction <- read_analysis(4, 'pgls_predictions_display_pg.csv', 'plotted_data')
  ms <- read_analysis(4, 'model_summary.csv', 'plotted_data') %>% filter(model == 'log_log_lambda')
  stopifnot(nrow(elements) == 380, nrow(summaries) == 30, nrow(pgls) == 19, sum(pgls$n_elements) == 239,
    abs(ms$slope - .024196124321278045) < 1e-14,
    identical(summaries$species[order(summaries$row)], phylogenetic_species_order(summaries$species)))
  a <- ggplot(elements, aes(estimate, row)) +
    geom_violin(aes(group = species), orientation = 'y', position = 'identity', bw = stats::bw.nrd0(elements$estimate),
      trim = TRUE, scale = 'width', width = .82, quantiles = NULL, fill = 'white', colour = '#B3B3B3', linewidth = .3) +
    geom_vline(xintercept = 0, linetype = 'dotted', linewidth = .4, colour = '#555555') +
    geom_point(position = position_jitter(width = 0, height = .15, seed = 20260828), shape = 21,
      size = .92, stroke = .23, fill = 'white', colour = '#333333', alpha = .75) +
    geom_segment(data = summaries, aes(x = ci_low, xend = ci_high, y = row, yend = row), inherit.aes = FALSE, colour = blue, linewidth = .4) +
    geom_segment(data = summaries, aes(x = estimate, xend = estimate, y = row - .25, yend = row + .25), inherit.aes = FALSE, colour = blue, linewidth = .8) +
    scale_x_continuous(breaks = log2(c(.025, .05, .1, .25, .5, 1, 1.5)), labels = c('0.025', '0.05', '0.1', '0.25', '0.5', '1', '1.5'),
      limits = log2(c(.025, 1.6)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(.5, 30.5), breaks = summaries$row, labels = summaries$species_label, expand = c(0, 0),
      sec.axis = dup_axis(breaks = summaries$row, labels = summaries$support_label, name = 'Genome size (pg/1C) / elements (n)')) +
    labs(title = 'LTR terminal-to-internal depth distributions', x = 'Terminal-to-internal depth ratio (log2 scale)', y = NULL) +
    theme_gbe() + row_theme + theme(axis.text.y.left = element_text(face = 'italic', size = 7.5),
      axis.text.y.right = element_text(size = 7.5), axis.title.y.right = element_text(size = 8),
      panel.border = element_rect(fill = NA, colour = '#AAAAAA', linewidth = .3))
  make_scatter <- function(label_all = FALSE) {
    b <- ggplot(pgls, aes(genome_size_pg, ratio_median)) +
      geom_ribbon(data = prediction, aes(x = genome_size_pg, ymin = ratio_ci_low, ymax = ratio_ci_high), inherit.aes = FALSE, fill = blue, alpha = .12) +
      geom_hline(yintercept = 1, linetype = 'dotted', linewidth = .4, colour = '#555555') +
      geom_segment(aes(x = genome_size_ci_low_pg, xend = genome_size_ci_high_pg, yend = ratio_median), colour = grey, linewidth = .3) +
      geom_segment(aes(y = median_bootstrap_ci_low, yend = median_bootstrap_ci_high, xend = genome_size_pg), colour = grey, linewidth = .3) +
      geom_line(data = prediction, aes(x = genome_size_pg, y = ratio_fit), inherit.aes = FALSE, colour = blue, linewidth = .7) +
      geom_point(aes(fill = reproductive_strategy, shape = reproductive_strategy), colour = '#222222', size = 2, stroke = .4) +
      strategy_fill() + strategy_shape() +
      scale_x_log10(breaks = c(9, 12, 15, 18, 22), limits = range(c(7.7, 23, pgls$genome_size_ci_low_gbp, pgls$genome_size_ci_high_gbp)) / .978,
        expand = expansion(mult = .02)) +
      scale_y_log10(breaks = c(.25, .5, 1, 1.5), labels = label_number(),
        limits = range(c(.23, 1.55, pgls$median_bootstrap_ci_low, pgls$median_bootstrap_ci_high, prediction$ratio_ci_low, prediction$ratio_ci_high)),
        expand = expansion(mult = .02)) +
      labs(title = 'Species median ratio and genome size', x = 'Genome size (pg/1C; log10 scale)', y = 'Median depth ratio\n(log10 scale)') +
      theme_gbe() + key_theme + guides(fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    if (label_all) b <- b + geom_segment(aes(xend = line_x_pg, yend = line_y), colour = '#999999', linewidth = .2) +
      geom_text(aes(x = label_x_pg, y = label_y, label = species_label, hjust = label_hjust),
        family = GBE_FONT, fontface = 'italic', size = label_size)
    b
  }
  b <- make_scatter()
  legend <- paste('(a) All 380 retained LTR/Gypsy elements from 30 species are shown in the fixed phylogenetic row order; no extra outliers were removed for this refinement.',
    'Gray density outlines use a pooled bandwidth in log2-ratio space, equal maximum widths and observed-range trimming. Hollow points show every retained element, jittered vertically only.',
    'Blue marks show species medians of raw ratios and released 95% intervals from 2,000 element bootstraps. Right labels give genome size and element count; n/a marks unavailable genome size.',
    'The ratio is length-weighted mean per-base depth across both terminal LTRs divided by mean depth across the entire intervening sequence, including zero-depth bases.',
    'The frozen primary filter requires elements at least 3,000 bp long, at least five parsed TEsorter domains, inclusive within-species 1.5-IQR fences, and at least 80% positive-depth coverage in each of the three regions.',
    'Ratio 1 indicates equal mapped depth; it is not a calibrated no-deletion baseline. This statistic is a mapping-depth proxy, not an ectopic-recombination or DNA-loss rate. Elements are not independent biological replicates.',
    '(b) Nineteen overlapping species (239 elements) contribute one median each. Stored log10–log10 Pagel-lambda PGLS predictions are reused unchanged, with slope 0.024196 (95% t interval −0.560865–0.609257; P = 0.931489), lambda approximately zero.',
    'The band is a 95% coefficient interval, not a species prediction interval. Gray bars show released measurement intervals, which are not propagated into the regression. The broad slope interval allows associations in either direction.',
    'All 19 species are named in the isolated regression panel. Stored Gbp coordinates are divided by 0.978 to display pg/1C; this unit conversion leaves the log–log slope and inference unchanged.', genome_boundary)
  collage <- (a / b) + plot_layout(heights = c(1.88, 1), guides = 'collect') + plot_annotation(tag_levels = 'a') & key_theme
  save_candidate(4, 'LTR depth ratios and genome size', collage,
    list(a = panel(a, 'LTR depth distributions', 180, 190, legend, 'Thirty named species, all 380 element ratios, blue medians and intervals, and genome-size/count annotations.'),
      b = panel(make_scatter(TRUE), 'LTR ratio and genome size', 180, 125, legend, 'Nineteen named species and measurement intervals around a nearly flat stored PGLS line with a broad band.')),
    180, 245, 'Frozen Figure 4 element ratios, species summaries, stored PGLS data, summary and predictions', legend,
    'Horizontal violins show all retained elements above a 19-species PGLS comparison with genome size.')
  copy_exports(4, 4, 'plotted_data')
  write_notes(4, 'Retains the user-requested reference-style horizontal violins and all existing extremes. The collage uses uncluttered species points in the single 19-species regression; the expanded isolated regression labels all species. All medians, intervals, counts and stored predictions are unchanged. The previously excluded PGLSb influence plot is not restored. No stronger deletion interpretation is introduced.')
}
message('Selected polished repeat figures rendered from frozen exports: ', paste(selected, collapse = ', '))

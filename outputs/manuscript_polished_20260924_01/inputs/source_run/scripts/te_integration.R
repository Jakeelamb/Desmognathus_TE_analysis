#!/usr/bin/env Rscript
# Fixed exploratory endpoint and finite models; isolated files only.
source('scripts/common.R')
suppressPackageStartupMessages({library(phylolm); library(jsonlite)})
options(warn = 1)
set.seed(20260925)
ROOT <- normalizePath(file.path(RUN, '../..'), mustWork = TRUE)
SNAP <- file.path(RUN, 'inputs/te_integration')
OUT11 <- file.path(RUN, 'figures/figure_11/analysis')
OUT12 <- file.path(RUN, 'figures/figure_12/analysis')
for (p in c(SNAP, OUT11, OUT12)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
sha <- function(p) sub(' .*', '', system2('sha256sum', shQuote(p), stdout = TRUE))
sources <- c(
  landscape = file.path(ROOT, 'analyses/01_transposable_elements/data/repeatmasker_divergence_landscape.csv'),
  frozen_landscape = file.path(DATA, 'Supplementary_Data_S12_repeatmasker_divergence_landscape.csv'),
  frozen_traits24 = file.path(DATA, 'Supplementary_Data_S04_path_analysis_traits24.csv'),
  traits24 = file.path(RUN, 'figures/figure_8/analysis/traits_joined.csv'),
  scaling24 = file.path(RUN, 'figures/figure_8/analysis/scaling.csv'),
  body_regressions24 = file.path(RUN, 'figures/figure_8/analysis/regression_coefficients.csv'),
  organismal_traits = file.path(ROOT, 'analyses/03_phylogenetic_path/data/organismal_traits.csv'),
  annotations34 = file.path(DATA, 'Supplementary_Data_S42_te_bonsai_tip_annotations.csv'),
  quality34 = file.path(DATA, 'Supplementary_Data_S06_dnapipete_input_quality.csv'),
  mass34 = file.path(DATA, 'Supplementary_Data_S07_dnapipete_mass_accounting.csv'),
  overlap21 = file.path(DATA, 'Supplementary_Data_S40_te_genome_size_overlap21.csv'),
  tree46 = file.path(DATA, 'trees/desmo900dated_source46.tre'),
  tree24 = file.path(RUN, 'inputs/body_size/path24_time_tree.nwk'))
dest <- setNames(file.path(SNAP, paste0(names(sources), ifelse(grepl('\\.(nwk|tre)$', sources), '.nwk', '.csv'))), names(sources))
stopifnot(all(file.exists(sources)))
for (id in names(sources)) {
  if (!file.exists(dest[[id]])) stopifnot(file.copy(sources[[id]], dest[[id]], overwrite = FALSE))
  stopifnot(identical(sha(sources[[id]]), sha(dest[[id]])))
}
manifest <- tibble(source_id = names(sources), source_path = unname(sources), snapshot_path = unname(dest), sha256 = vapply(dest, sha, character(1)))
write_csv(manifest, file.path(SNAP, 'MANIFEST.csv'))
# This design is written before any model fit, and must not change on a refit.
design <- c('# Fixed exploratory TE integration design', '',
 'Primary endpoint: RepeatMasker LTR aligned hit bp in bins 0–4 divided by all LTR aligned hit bp across every available divergence bin. No threshold search. Each species is one observation. Fraction endpoints must be strictly between zero and one; otherwise stop rather than add an arbitrary pseudocount.',
 'Context and denominator sensitivity: low-divergence LTR aligned hit bp / all RepeatMasker aligned hit bp. Decomposition: this fraction = (all LTR bp / all hit bp) * (low-divergence LTR bp / all LTR bp). These are hit-weighted pipeline summaries, not genomic TE occupancy, insertion rate, absolute age or direct activity measurements.',
 'Figure 11 uses exactly the 21 species common to the released TE and cellular panels, with no missing-species imputation. Scores are logit transformed and standardized across all 34 TE species. Genome, cell, nucleus and maximum SVL retain the existing full24 log10 mean and SD. Pagel-lambda PGLS fits g~t; g~t+b; c~g+b; c~g+t+b; n~c+g+b; n~c+g+t+b. The last equation is a descriptive conditional diagnostic, not a selected causal graph. Repeat this finite set for the all-hit denominator sensitivity; do not select the preferred denominator by its P value.',
 'Figure 12 uses all 34 TE species and the released, expert post-hoc three-level reproductive strategies. Fit t~1, t~log10(SVL), t~strategy, and t~strategy+log10(SVL), estimating lambda and ML residual variance independently per model. Report an approximate chi-square likelihood-ratio omnibus test (2 df) for strategy, unadjusted and adjusted for SVL. Report contrasts with aquatic eggs as reference using t confidence intervals. No stepwise or threshold search; no pairwise Tukey screening.',
 'Two predefined descriptive quality checks: Spearman association of the primary score with order-unresolved fraction and repeat-aligned fraction of the configured dnaPipeTE sample. This is a cross-pipeline quality flag, not an absolute-genome-load correction; retain both correlations regardless of sign or strength.',
 'All tests are exploratory. Strategy assignments were made after initial PCA review; group sizes, phylogenetic clustering and small direct-developer sample limit inference. Conditional intervals omit tree, measurement and lambda-estimation uncertainty. Existing figures, data and notes are unchanged.')
design_path <- file.path(SNAP, 'DESIGN_BEFORE_FITTING.md')
if (file.exists(design_path)) stopifnot(identical(readLines(design_path), design)) else writeLines(design, design_path)
read_snapshot <- function(id) read_csv(dest[[id]], show_col_types = FALSE, progress = FALSE)
land <- read_snapshot('landscape')
stopifnot(identical(sha(dest[['landscape']]), sha(dest[['frozen_landscape']])))
ann <- read_snapshot('annotations34')
org <- read_snapshot('organismal_traits')
quality <- read_snapshot('quality34')
traits24 <- read_snapshot('traits24')
frozen24 <- read_snapshot('frozen_traits24')
stopifnot(!anyDuplicated(ann$species), !anyDuplicated(org$species), !anyDuplicated(quality$species),
  nrow(ann) == 34L, nrow(traits24) == 24L, nrow(frozen24) == 24L,
  !anyDuplicated(land[, c('species','order','divergence_bin_start_pct')]))
raw_trait_columns <- c('genome_size_pg_fuscus_anchored','cell_area_um2','nucleus_area_um2')
stopifnot(setequal(traits24$tree_tip, frozen24$tree_tip),
  isTRUE(all.equal(as.matrix(traits24[, raw_trait_columns]),
    as.matrix(frozen24[match(traits24$tree_tip, frozen24$tree_tip), raw_trait_columns]),
    check.attributes = FALSE, tolerance = 0)))
scaling24 <- read_snapshot('scaling24')
for (v in scaling24$variable) {
  s <- filter(scaling24, variable == v)
  stopifnot(max(abs(traits24[[v]] - (log10(traits24[[s$raw_column]]) - s$full24_mean) / s$full24_sd)) < 1e-12)
}
denoms <- land %>% group_by(species) %>% summarise(all_hit_bp = first(species_total_hit_bp),
  n_denominators = n_distinct(species_total_hit_bp), summed_hit_bp = sum(hit_bp), .groups = 'drop')
stopifnot(all(denoms$n_denominators == 1L), all(denoms$all_hit_bp == denoms$summed_hit_bp),
  all(land$hit_bp >= 0), all(land$divergence_bin_start_pct >= 0),
  all(land$divergence_bin_start_pct %% 1 == 0), setequal(denoms$species, ann$species))
score34 <- land %>% filter(order == 'LTR') %>% group_by(species) %>%
  summarise(ltr_hit_bp = sum(hit_bp), low_ltr_hit_bp = sum(hit_bp[divergence_bin_start_pct < 5]), .groups = 'drop') %>%
  left_join(denoms, by = 'species', relationship = 'one-to-one') %>%
  mutate(within_ltr_fraction = low_ltr_hit_bp / ltr_hit_bp,
    ltr_fraction_all_hits = ltr_hit_bp / all_hit_bp,
    low_ltr_fraction_all_hits = low_ltr_hit_bp / all_hit_bp) %>%
  left_join(ann, by = 'species', relationship = 'one-to-one') %>%
  left_join(org %>% select(species, max_svl_mm, max_svl_source_id, max_svl_source_taxon), by = 'species', relationship = 'one-to-one') %>%
  left_join(quality %>% select(species, order_unresolved_fraction, repeat_aligned_fraction_of_sample,
      configured_quantification_sample_bp, run_aggregation), by = 'species', relationship = 'one-to-one') %>%
  mutate(strategy = factor(reproductive_strategy, levels = te_strategy_levels),
    log10_svl = log10(max_svl_mm),
    display_species = factor(display_species, levels = rev(phylogenetic_species_levels(species))))
stopifnot(nrow(score34) == 34L, all(score34$within_ltr_fraction > 0 & score34$within_ltr_fraction < 1),
  all(score34$low_ltr_fraction_all_hits > 0 & score34$low_ltr_fraction_all_hits < 1),
  !anyNA(score34$max_svl_mm), all(score34$max_svl_mm > 0), !anyNA(score34$strategy),
  identical(as.integer(table(score34$strategy)), c(7L, 24L, 3L)),
  max(abs(score34$low_ltr_fraction_all_hits - score34$ltr_fraction_all_hits * score34$within_ltr_fraction)) < 1e-15)
te_scaling <- bind_rows(lapply(c(within_ltr = 'within_ltr_fraction', all_hit = 'low_ltr_fraction_all_hits'), function(column) {
  x <- qlogis(score34[[column]])
  tibble(raw_column = column, transformation = 'logit', mean34 = mean(x), sd34 = sd(x), n_scaling = 34L)
})) %>% mutate(endpoint = c('within_ltr','all_hit'))
for (e in te_scaling$endpoint) {
  s <- filter(te_scaling, endpoint == e)
  score34[[paste0('t_', e)]] <- (qlogis(score34[[s$raw_column]]) - s$mean34) / s$sd34
  stopifnot(abs(mean(score34[[paste0('t_',e)]])) < 1e-12, abs(sd(score34[[paste0('t_',e)]]) - 1) < 1e-12)
}
write_csv(te_scaling, file.path(OUT11, 'te_scaling.csv'))
write_csv(score34, file.path(OUT12, 'repeat_score_species34.csv'))
joined <- traits24 %>% rename(path_species = species, species = tree_tip) %>%
  inner_join(score34 %>% select(species, within_ltr_fraction, low_ltr_fraction_all_hits,
      ltr_fraction_all_hits, low_ltr_hit_bp, ltr_hit_bp, all_hit_bp, t_within_ltr, t_all_hit,
      strategy, order_unresolved_fraction, repeat_aligned_fraction_of_sample), by = 'species', relationship = 'one-to-one')
overlap <- read_snapshot('overlap21')
stopifnot(nrow(joined) == 21L, setequal(joined$species, overlap$species),
  setequal(setdiff(traits24$tree_tip, joined$species), c('brimleyorum','folkertsi','ochrophaeus')),
  max(abs(joined$genome_size_pg_fuscus_anchored - overlap$genome_size_pg_fuscus_anchored[match(joined$species, overlap$species)])) < 1e-8,
  !anyNA(joined[, c('g','c','n','b','t_within_ltr','t_all_hit')]))
write_csv(joined, file.path(OUT11, 'integrated_species21.csv'))
coverage <- traits24 %>% transmute(species = tree_tip, genome_size_pg_fuscus_anchored, cell_area_um2,
  nucleus_area_um2, max_svl_mm, te_included = tree_tip %in% joined$species,
  status = if_else(te_included, 'Complete TE + cellular + SVL', 'No TE data; not imputed'))
write_csv(coverage, file.path(OUT11, 'species_coverage24.csv'))
tree21 <- keep.tip(read.tree(dest[['tree24']]), joined$species)
tree34 <- keep.tip(read.tree(dest[['tree46']]), score34$species)
stopifnot(Ntip(tree21) == 21L, Ntip(tree34) == 34L, all(tree21$edge.length >= 0), all(tree34$edge.length >= 0))
write.tree(tree21, file.path(OUT11, 'analysis_tree21.nwk'))
write.tree(tree34, file.path(OUT12, 'analysis_tree34.nwk'))
warning_rows <- list(); independent_gls_errors <- numeric()
fit_one <- function(formula, data, tree, fit_id, endpoint, analysis) {
  d <- as.data.frame(data); rownames(d) <- d$species; d <- d[tree$tip.label, , drop = FALSE]
  stopifnot(identical(rownames(d), tree$tip.label))
  fit <- withCallingHandlers(phylolm(formula, data = d, phy = tree, model = 'lambda', lower.bound = 1e-7, upper.bound = 1),
    warning = function(w) warning_rows[[length(warning_rows) + 1L]] <<- tibble(analysis = analysis,
      endpoint = endpoint, fit_id = fit_id, warning = conditionMessage(w)))
  sm <- as.data.frame(summary(fit)$coefficients); df <- fit$n - fit$d; crit <- qt(.975, df)
  rows <- tibble(analysis = analysis, endpoint = endpoint, fit_id = fit_id,
    equation = paste(deparse(formula), collapse = ' '), term = rownames(sm),
    estimate = sm$Estimate, standard_error = sm$StdErr, ci_low = sm$Estimate - crit * sm$StdErr,
    ci_high = sm$Estimate + crit * sm$StdErr, t_value = sm$t.value, p_value = sm$p.value,
    lambda = as.numeric(fit$optpar), lambda_boundary = lambda <= 1.01e-7 | lambda >= .999999,
    n_species = fit$n, residual_df = df, n_coefficients = fit$d, log_likelihood = fit$logLik,
    interval_method = 't, n minus regression coefficients; conditional on fitted lambda')
  X <- model.matrix(formula, d); y <- model.response(model.frame(formula, d))
  V <- vcv(tree)[tree$tip.label, tree$tip.label]; diag0 <- diag(V); V <- V * fit$optpar; diag(V) <- diag0
  gls_beta <- solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, y)))
  err <- max(abs(as.numeric(gls_beta) - as.numeric(coef(fit))))
  independent_gls_errors <<- c(independent_gls_errors, err)
  stopifnot(err < 1e-6, max(abs(rows$p_value - 2 * pt(-abs(rows$t_value), df))) < 1e-12,
    max(abs(as.numeric(X %*% coef(fit)) - as.numeric(fit$fitted.values))) < 1e-10,
    all(rows$ci_low <= rows$estimate), all(rows$ci_high >= rows$estimate))
  list(fit = fit, rows = rows, data = d, formula = formula)
}
spec11 <- list(genome_unadjusted = g ~ t, genome_body_adjusted = g ~ t + b,
  cell_baseline = c ~ g + b, cell_te_adjusted = c ~ g + t + b,
  nucleus_baseline = n ~ c + g + b, nucleus_te_adjusted = n ~ c + g + t + b)
fits11 <- list(); rows11 <- list(); covariance11 <- list()
for (e in te_scaling$endpoint) for (id in names(spec11)) {
  d <- joined; d$t <- d[[paste0('t_', e)]]
  result <- fit_one(spec11[[id]], d, tree21, id, e, 'integration21')
  fits11[[paste(e, id, sep = '/')]] <- result
  rows11[[length(rows11) + 1L]] <- result$rows
  vc <- result$fit$vcov; grid <- expand.grid(term1 = rownames(vc), term2 = colnames(vc), stringsAsFactors = FALSE)
  covariance11[[length(covariance11) + 1L]] <- tibble(endpoint = e, fit_id = id, term1 = grid$term1, term2 = grid$term2,
    covariance = mapply(function(a,b) vc[a,b], grid$term1, grid$term2))
}
reg11 <- bind_rows(rows11)
write_csv(reg11, file.path(OUT11, 'regression_coefficients.csv'))
write_csv(bind_rows(covariance11), file.path(OUT11, 'coefficient_covariance.csv'))
scaling24 <- read_snapshot('scaling24')
gs <- filter(scaling24, variable == 'g'); bs <- filter(scaling24, variable == 'b')
pred11 <- list()
for (e in te_scaling$endpoint) for (id in c('genome_unadjusted','genome_body_adjusted')) {
  result <- fits11[[paste(e, id, sep = '/')]]; s <- filter(te_scaling, endpoint == e)
  fraction <- plogis(seq(min(qlogis(joined[[s$raw_column]])), max(qlogis(joined[[s$raw_column]])), length.out = 121))
  d <- tibble(t = (qlogis(fraction) - s$mean34) / s$sd34, b = mean(joined$b))
  X <- model.matrix(delete.response(terms(result$formula)), d)
  eta <- as.numeric(X %*% coef(result$fit)); se <- sqrt(rowSums((X %*% result$fit$vcov) * X))
  crit <- qt(.975, result$fit$n - result$fit$d)
  raw <- function(z) 10^(z * gs$full24_sd + gs$full24_mean)
  pred11[[length(pred11) + 1L]] <- tibble(endpoint = e, fit_id = id, te_fraction = fraction,
    fitted = raw(eta), ci_low = raw(eta - crit*se), ci_high = raw(eta + crit*se),
    svl_held_mm = 10^(mean(joined$b) * bs$full24_sd + bs$full24_mean), n_species = 21L,
    interval_type = 'Fixed-effect conditional median and 95% confidence interval; not individual prediction interval')
  stopifnot(max(abs(log10(raw(eta)) - (eta * gs$full24_sd + gs$full24_mean))) < 1e-12)
}
pred11 <- bind_rows(pred11); write_csv(pred11, file.path(OUT11, 'genome_predictions.csv'))
baseline24 <- read_snapshot('body_regressions24') %>% filter(panel == 'full24',
  relation_id %in% c('genome_to_cell_body_adjusted', 'genome_and_cell_to_nucleus_body_adjusted')) %>%
  mutate(link = case_when(relation_id == 'genome_to_cell_body_adjusted' & term == 'g' ~ 'Genome to cell',
    relation_id == 'genome_and_cell_to_nucleus_body_adjusted' & term == 'c' ~ 'Cell to nucleus',
    relation_id == 'genome_and_cell_to_nucleus_body_adjusted' & term == 'g' ~ 'Genome to nucleus')) %>%
  filter(!is.na(link)) %>% transmute(link, estimate = standardized_coefficient, ci_low, ci_high,
    model = '24 species: SVL', n_species)
core21 <- reg11 %>% filter(endpoint == 'within_ltr') %>% mutate(link = case_when(
  fit_id %in% c('cell_baseline','cell_te_adjusted') & term == 'g' ~ 'Genome to cell',
  fit_id %in% c('nucleus_baseline','nucleus_te_adjusted') & term == 'c' ~ 'Cell to nucleus',
  fit_id %in% c('nucleus_baseline','nucleus_te_adjusted') & term == 'g' ~ 'Genome to nucleus')) %>%
  filter(!is.na(link)) %>% transmute(link, estimate, ci_low, ci_high,
    model = if_else(grepl('te_adjusted',fit_id), '21 species: SVL + TE', '21 species: SVL'), n_species)
core <- bind_rows(baseline24, core21)
stopifnot(nrow(core) == 9L)
write_csv(core, file.path(OUT11, 'cellular_coefficients_matched_comparison.csv'))
spec12 <- list(null = t ~ 1, body = t ~ log10_svl, strategy = t ~ strategy, strategy_body = t ~ strategy + log10_svl)
fits12 <- list(); rows12 <- list(); omnibus <- list(); pred12 <- list()
for (e in te_scaling$endpoint) {
  d <- score34; d$t <- qlogis(d[[filter(te_scaling, endpoint == e)$raw_column]])
  for (id in names(spec12)) {
    result <- fit_one(spec12[[id]], d, tree34, id, e, 'strategy34')
    fits12[[paste(e,id,sep = '/')]] <- result; rows12[[length(rows12)+1L]] <- result$rows
  }
  for (adjusted in c(FALSE,TRUE)) {
    alternative <- if (adjusted) 'strategy_body' else 'strategy'; null <- if (adjusted) 'body' else 'null'
    full <- fits12[[paste(e,alternative,sep = '/')]]$fit; reduced <- fits12[[paste(e,null,sep = '/')]]$fit
    lr <- 2*(full$logLik - reduced$logLik); df <- full$d - reduced$d
    stopifnot(lr >= -1e-8, df == 2L)
    omnibus[[length(omnibus)+1L]] <- tibble(endpoint = e, comparison = paste(alternative, 'vs', null),
      svl_adjusted = adjusted, alternative_logLik = full$logLik, null_logLik = reduced$logLik,
      likelihood_ratio = lr, df = df, p_value = pchisq(max(0,lr), df = df, lower.tail = FALSE),
      n_species = 34L, test = 'Exploratory approximate chi-square ML likelihood-ratio test; lambda separately fitted')
  }
  result <- fits12[[paste(e,'strategy_body',sep = '/')]]
  for (strategy_id in te_strategy_levels) {
    ds <- filter(score34, strategy == strategy_id)
    new <- tibble(log10_svl = seq(min(ds$log10_svl), max(ds$log10_svl), length.out = 81),
      strategy = factor(strategy_id, levels = te_strategy_levels))
    X <- model.matrix(~ strategy + log10_svl, new)
    eta <- as.numeric(X %*% coef(result$fit)); se <- sqrt(rowSums((X %*% result$fit$vcov)*X)); crit <- qt(.975, result$fit$n-result$fit$d)
    pred12[[length(pred12)+1L]] <- tibble(endpoint = e, strategy = strategy_id, max_svl_mm = 10^new$log10_svl,
      fitted = plogis(eta), ci_low = plogis(eta-crit*se), ci_high = plogis(eta+crit*se),
      n_group = nrow(ds), n_model = 34L, interval_type = 'Inverse-logit conditional median and 95% confidence interval; group observed SVL range')
  }
}
reg12 <- bind_rows(rows12); omnibus <- bind_rows(omnibus); pred12 <- bind_rows(pred12)
write_csv(reg12, file.path(OUT12, 'regression_coefficients.csv'))
write_csv(omnibus, file.path(OUT12, 'strategy_omnibus_tests.csv'))
write_csv(pred12, file.path(OUT12, 'strategy_predictions.csv'))
qc <- bind_rows(lapply(c('order_unresolved_fraction','repeat_aligned_fraction_of_sample'), function(column) {
  test <- cor.test(score34$within_ltr_fraction, score34[[column]], method = 'spearman', exact = FALSE)
  tibble(quality_variable = column, n_species = 34L, spearman_rho = as.numeric(test$estimate), p_value = test$p.value,
    purpose = 'Descriptive cross-pipeline quality association; not phylogenetic, not a genome-load correction')
}))
write_csv(qc, file.path(OUT12, 'predefined_quality_associations.csv'))
warnings <- if (length(warning_rows)) bind_rows(warning_rows) else tibble(analysis = character(),endpoint = character(),fit_id = character(),warning = character())
write_csv(warnings, file.path(OUT11, 'warnings_all_fits.csv'))
write_csv(filter(warnings, analysis == 'strategy34'), file.path(OUT12, 'warnings.csv'))
stopifnot(identical(manifest$sha256, vapply(dest, sha, character(1))))
verification <- list(passed = TRUE, n_te = 34L, n_integration = 21L,
  no_imputation = TRUE, excluded_missing_te = c('brimleyorum','folkertsi','ochrophaeus'),
  source_snapshots_unchanged = TRUE, canonical_landscape_matches_frozen = TRUE,
  all_bins_in_ltr_denominator = TRUE, no_pseudocount_needed = TRUE,
  maximum_denominator_identity_error = max(abs(score34$low_ltr_fraction_all_hits - score34$ltr_fraction_all_hits * score34$within_ltr_fraction)),
  maximum_independent_GLS_coefficient_error = max(independent_gls_errors),
  all_coefficients_and_predictions_verified = TRUE,
  group_n = as.list(table(score34$strategy)), configured_sample_bp = unique(score34$configured_quantification_sample_bp),
  phylolm_version = as.character(packageVersion('phylolm')), n_warning_occurrences = nrow(warnings))
write_json(verification, file.path(OUT11,'verification.json'), pretty = TRUE, auto_unbox = TRUE)
write_json(verification, file.path(OUT12,'verification.json'), pretty = TRUE, auto_unbox = TRUE)
methods <- c('# Focused TE integration and divergence decomposition', '', design[-c(1,2)], '',
 '## Implementation', '',
 'Exact species identifiers link the frozen Path24/body-size table to the released TE34 tables. The TE strategy assignments are taken from S42, including its accession-specific fuscus assignment; they are not inferred anew from morphology labels. Fitting uses the same dated source-tree branch lengths as existing analyses; publication/calibration provenance remains unresolved. Figure 11 uses the frozen Path24 tree pruned to 21; Figure 12 uses the full source tree pruned to 34. These are species phylogenies, not Bonsai composition trees.',
 'Genome, cell, nucleus and SVL use existing full24 standardized log10 variables. Figure 11 TE scores use the logit fraction standardized over all 34 TE species, so the TE coefficient is an outcome full24 SD per one full34 SD of logit TE fraction. Figure 12 strategy coefficients are raw logit-fraction differences (not standardized effect sizes), and the SVL slope is per log10 mm. Aquatic eggs is the reference group. Both zero and one are absent, so no correction is used. Displayed percentages always invert the transformation.',
 'phylolm estimates Pagel lambda by ML on [1e-7,1]. Coefficient confidence intervals use t_(n-p) * SE, consistent with package coefficient P values. Every displayed P value is unadjusted exploratory evidence. Confidence intervals condition on estimated lambda and exclude trait, tree, calibration and within-species uncertainty. Raw prediction curves represent the conditional median after inverse log/logit transformation, not arithmetic means or individual prediction intervals.',
 'No model is selected among the fixed regression equations. Cell and nuclear comparisons use a common21 panel, with the previous full24 result shown separately to distinguish membership from additional adjustment. Conditional coefficient arrows in textual labels denote formula direction, not identified causation. The most complex nuclear diagnostic uses four predictors plus intercept (n = 21, residual df = 16); it is not evidence for a five-node causal network.',
 'Strategy tests compare ML fits with two extra group parameters; the chi-square reference is approximate in this small phylogenetic sample. The seven aquatic-egg,24 terrestrial-egg and three direct-developer species are not seven/24/three independent evolutionary transitions. Strategy was expert-coded post hoc after PCA review, and strategy/body size are strongly structured together. Fitted group curves are limited to each group\'s observed SVL range. Neither a strategy contrast nor low divergence establishes TE activation or insertion timing.',
 'RepeatMasker fractions are computed from aligned hit bp in that pipeline. Overlaps/hit weighting are retained from the source. dnaPipeTE unresolved mass and repeat-aligned sample fraction come from a separate pipeline and are only descriptive QC covariates. S06 configured_quantification_sample_bp is 1,500,000,000 (1.5 Gb); it is not a whole-genome size or occupancy denominator. No genome-wide TE-load inference is made.',
 'Verification checks snapshot SHA256, exact21 membership, no imputation, all 34 SVL/strategy coverage, no duplicate species or bins, positive fractions, exact hit-denominator sums and multiplicative identity. Every regression coefficient is independently reconstructed with GLS using the lambda-transformed ape::vcv tree covariance; t intervals, fitted values and raw-scale predictions are checked. Boundary lambda warnings remain in the build log and CSV.', '',
 'Reproduce from this isolated run with:', '',
 '    PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH Rscript --no-environ scripts/te_integration.R')
writeLines(methods, file.path(OUT11, 'README.md')); writeLines(methods, file.path(OUT12, 'README.md'))

# Render only the two new figures. All published and earlier candidates stay intact.
style <- theme_gbe() + theme(plot.title = element_text(size = 9, face = 'bold', margin = margin(b = 5)),
  plot.subtitle = element_text(size = 7.5, margin = margin(b = 4)),
  axis.text = element_text(size = 7.5), plot.margin = margin(6,8,5,6))
fmtp <- function(x) ifelse(x < .001, formatC(x, format = 'e', digits = 1), sprintf('%.3f',x))
strategy_scales <- list(scale_fill_manual(values = te_strategy_colors, breaks = te_strategy_levels,
  labels = te_strategy_labels, name = NULL), scale_colour_manual(values = te_strategy_colors,
  breaks = te_strategy_levels, labels = te_strategy_labels, name = NULL),
  scale_shape_manual(values = te_strategy_shapes, breaks = te_strategy_levels, labels = te_strategy_labels, name = NULL))
scatter11 <- function(all_labels = FALSE) {
  curves <- filter(pred11, endpoint == 'within_ltr')
  est <- filter(reg11, endpoint == 'within_ltr', fit_id == 'genome_body_adjusted', term == 't')
  lab <- joined %>% mutate(label = if (all_labels) species_label(species) else
    if_else(species %in% c('aeneus','wrighti','valtos','orestes'), species_label(species), ''))
  ggplot(joined, aes(within_ltr_fraction*100, genome_size_pg_fuscus_anchored)) +
    geom_ribbon(data = filter(curves,fit_id == 'genome_body_adjusted'),
      aes(x = te_fraction*100, y = NULL, ymin = ci_low, ymax = ci_high), inherit.aes = FALSE, fill = '#E3E3E3') +
    geom_line(data = curves, aes(x = te_fraction*100, y = fitted, linetype = fit_id), inherit.aes = FALSE, linewidth = .65) +
    geom_errorbar(aes(ymin = genome_size_pg_ci_low, ymax = genome_size_pg_ci_high), linewidth = .28, width = 0, colour = '#9E9E9E') +
    geom_point(aes(fill = strategy, shape = strategy), colour = '#222222', size = 2.2, stroke = .45) +
    geom_text_repel(data = lab, aes(label = label), family = GBE_FONT, fontface = 'italic', size = 7/ggplot2::.pt,
      seed = 20260925, max.overlaps = Inf, max.time = Inf, max.iter = 20000, box.padding = .35,
      point.padding = .25, force = 2, min.segment.length = 0, segment.size = .2) +
    strategy_scales[c(1,3)] + scale_linetype_manual(values = c(genome_unadjusted = '22', genome_body_adjusted = 'solid'),
      labels = c(genome_unadjusted = 'Unadjusted', genome_body_adjusted = 'SVL adjusted'), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(.06,.15))) +
    scale_y_log10(expand = expansion(mult = c(.10,.18))) +
    labs(title = 'LTR divergence signature and genome size',
      subtitle = sprintf('21 species; SVL-adjusted β = %.2f; P = %s', est$estimate,fmtp(est$p_value)),
      x = 'LTR hit bp below 5% divergence (% of all LTR hit bp)', y = 'Genome size (pg/1C; log scale)') + style +
    guides(fill = guide_legend(ncol = 1, order = 1), shape = guide_legend(ncol = 1, order = 1), linetype = guide_legend(nrow = 1, order = 2))
}
effect11 <- reg11 %>% filter(endpoint == 'within_ltr', term == 't',
  fit_id %in% c('genome_body_adjusted','cell_te_adjusted','nucleus_te_adjusted')) %>%
  mutate(label = factor(case_when(fit_id == 'genome_body_adjusted' ~ 'Genome size | SVL',
    fit_id == 'cell_te_adjusted' ~ 'Cell area | genome + SVL', TRUE ~ 'Nucleus area | cell + genome + SVL'),
    levels = rev(c('Genome size | SVL','Cell area | genome + SVL','Nucleus area | cell + genome + SVL'))))
p11b <- ggplot(effect11, aes(estimate,label)) + geom_vline(xintercept = 0, linetype = '22', colour = '#999999', linewidth = .4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = 'y', width = .15, linewidth = .55) +
  geom_point(size = 2.2, shape = 21, fill = GBE_COLORS[['purple']], stroke = .5) +
  labs(title = 'Conditional association with the TE signature', subtitle = 'Same 21 species; fixed regression equations',
    x = 'Standardized TE coefficient (95% CI)', y = NULL) + style +
  scale_y_discrete(expand = expansion(add = .6))
core$model <- factor(core$model,levels = c('24 species: SVL','21 species: SVL','21 species: SVL + TE'))
core$link <- factor(core$link, levels = rev(c('Genome to cell','Cell to nucleus','Genome to nucleus')))
p11c <- ggplot(core, aes(estimate,link, colour = model,shape = model)) +
  geom_vline(xintercept = 0, linetype = '22', colour = '#999999', linewidth = .4) +
  geom_errorbar(aes(xmin = ci_low,xmax = ci_high), orientation = 'y', width = .14,
    linewidth = .5, position = position_dodge(width = .55)) +
  geom_point(size = 2,position = position_dodge(width = .55)) +
  scale_colour_manual(values = c('#999999',GBE_COLORS[['blue']],GBE_COLORS[['purple']]), name = NULL) +
  scale_shape_manual(values = c(1,16,17), name = NULL) +
  labs(title = 'Do cellular associations depend on TE adjustment?', subtitle = 'Separate taxon membership from adding the TE score',
    x = 'Standardized cellular coefficient (95% CI)', y = NULL) + style +
  guides(colour = guide_legend(ncol = 1),shape = guide_legend(ncol = 1))
legend11 <- paste('Focused exploratory integration of one prespecified repeat-divergence signature with genome size, cell area, nucleus area and maximum SVL.',
 'The primary TE score is RepeatMasker LTR hit bp below 5% divergence divided by all LTR hit bp, across all available divergence bins; it is not TE load, insertion rate or age.',
 '(a) All 21 complete species with genome-size measurement intervals, a marginal PGLS curve (dashed) and an SVL-adjusted curve (solid) at the geometric mean SVL of these 21 species. Gray shading is a conditional 95% confidence band. Adjusted curves need not follow the marginal point cloud.',
 '(b) TE coefficients from g~t+b, c~g+t+b, and n~c+g+t+b, where b is maximum SVL. The vertical bar denotes conditioning, not a causal direction. All use the same 21 species.',
 '(c) Genome-to-cell and the joint cell/genome-to-nucleus coefficients compare existing full24 body-adjusted fits, matched21 body-adjusted fits, and matched21 fits additionally including the TE score. Joint nuclear models retain both cell and genome predictors.',
 'Traits use the existing full24 log10 standardization; the logit TE fraction uses full34 standardization. Pagel-lambda ML regressions use t_(n-p) confidence intervals conditional on fitted lambda. Missing TE data for D. brimleyorum, D. folkertsi and D. ochrophaeus are not imputed.',
 'The all-hit denominator sensitivity and every equation, coefficient,lambda, warning and species inclusion are exported. Genome-size calibration, individual-animal uncertainty, SVL uncertainty and tree provenance remain as defined for prior figures. Associations do not establish an ordered TE-to-genome-to-cell causal chain.')
p11a <- scatter11()
collage11 <- p11a / p11b / p11c + plot_layout(heights = c(1.5,.72,1.05)) + plot_annotation(tag_levels = 'a')
save_candidate(11,'Repeat divergence and the genome–cell relationship',collage11,
  list(a = list(plot = scatter11(TRUE),title = 'LTR divergence signature and genome size',width_mm = 180,height_mm = 140,legend = legend11,alt_text = 'Twenty-one species with observed genome size plotted against their within-LTR low-divergence fraction, with marginal and SVL-adjusted PGLS curves.'),
    b = list(plot = p11b,title = 'Conditional TE associations',width_mm = 180,height_mm = 85,legend = legend11,alt_text = 'Conditional TE score coefficients for genome, cell and nucleus outcomes with 95% intervals.'),
    c = list(plot = p11c,title = 'Cellular associations before and after TE adjustment',width_mm = 180,height_mm = 115,legend = legend11,alt_text = 'Three cellular coefficients compared across full24, matched21, and TE-adjusted21 analyses.')),
  180,285,'inputs/te_integration;S04;S12;S42;figure8analysis',legend11,
  'A 21-species genome-size scatterplot and two coefficient panels link the fixed LTR-divergence signature with existing cellular relationships.')

decomposition <- score34 %>% select(species,display_species,strategy,ltr_fraction_all_hits,within_ltr_fraction,low_ltr_fraction_all_hits) %>%
  pivot_longer(c(ltr_fraction_all_hits,within_ltr_fraction,low_ltr_fraction_all_hits),names_to = 'component',values_to = 'fraction') %>%
  mutate(component = factor(component,levels = c('ltr_fraction_all_hits','within_ltr_fraction','low_ltr_fraction_all_hits'),
    labels = c('LTR / all hits','Low-divergence / LTR','Low-divergence LTR / all hits')))
write_csv(decomposition, file.path(OUT12,'plotted_decomposition.csv'))
p12a <- ggplot(decomposition,aes(fraction*100,display_species)) +
  geom_segment(aes(x = 0,xend = fraction*100,yend = display_species),colour = '#DDDDDD',linewidth = .4) +
  geom_point(aes(fill = strategy,shape = strategy),size = 1.8,stroke = .35) +
  facet_grid(cols = vars(component),scales = 'free_x') + strategy_scales[c(1,3)] +
  scale_x_continuous(expand = expansion(mult = c(.01,.08))) +
  labs(title = 'Separate LTR representation from its divergence profile',
    subtitle = 'Left fraction × middle fraction = right fraction; each panel has its own percentage axis',
    x = 'Aligned hit bp (%)',y = NULL) + style + theme_gbe_framed_facets() +
  theme(axis.text.y = element_text(face = 'italic',size = 7.2), strip.text = element_text(size = 7.5),
    panel.spacing = unit(9,'pt')) + guides(fill = guide_legend(ncol = 1),shape = guide_legend(ncol = 1))
omni_primary <- filter(omnibus,endpoint == 'within_ltr',svl_adjusted)
p12b <- ggplot(score34,aes(max_svl_mm,within_ltr_fraction*100)) +
  geom_ribbon(data = filter(pred12,endpoint == 'within_ltr'),
    aes(x = max_svl_mm,y = NULL,ymin = ci_low*100,ymax = ci_high*100,fill = strategy),
    inherit.aes = FALSE,alpha = .12) +
  geom_line(data = filter(pred12,endpoint == 'within_ltr'),
    aes(x = max_svl_mm,y = fitted*100,colour = strategy),inherit.aes = FALSE,linewidth = .6) +
  geom_point(aes(fill = strategy,shape = strategy),size = 2,stroke = .4) + strategy_scales +
  scale_x_log10(breaks = c(30,50,80,120)) +
  labs(title = 'Within-LTR signature, body size and strategy',
    subtitle = sprintf('34 species; strategy + SVL; omnibus P = %s',fmtp(omni_primary$p_value)),
    x = 'Maximum SVL (mm; log scale)',y = 'Low-divergence LTR hit bp (% of LTR hit bp)') + style +
  guides(fill = 'none',shape = 'none',colour = 'none')
contrast <- reg12 %>% filter(endpoint == 'within_ltr',fit_id %in% c('strategy','strategy_body'),grepl('^strategy',term)) %>%
  mutate(contrast = factor(if_else(term == 'strategydirect_development','Direct development','Terrestrial eggs + aquatic larvae'),
    levels = c('Direct development','Terrestrial eggs + aquatic larvae')),
    model = factor(if_else(fit_id == 'strategy','Unadjusted','SVL adjusted'),levels = c('Unadjusted','SVL adjusted')))
stopifnot(nrow(contrast) == 4L)
write_csv(contrast,file.path(OUT12,'plotted_strategy_contrasts.csv'))
p12c <- ggplot(contrast,aes(estimate,contrast,colour = model,shape = model)) +
  geom_vline(xintercept = 0,colour = '#999999',linetype = '22',linewidth = .4) +
  geom_errorbar(aes(xmin = ci_low,xmax = ci_high),orientation = 'y',width = .12,linewidth = .5,
    position = position_dodge(width = .4)) + geom_point(size = 2.2,position = position_dodge(width = .4)) +
  scale_colour_manual(values = c('#999999','#222222'),name = NULL) + scale_shape_manual(values = c(1,16),name = NULL) +
  labs(title = 'Phylogenetic strategy contrasts',subtitle = 'Reference: aquatic eggs (n = 7); other groups: n = 24 and 3',
    x = 'Difference in logit within-LTR fraction (95% CI)',y = NULL) + style +
  guides(colour = guide_legend(nrow = 1),shape = guide_legend(nrow = 1))
legend12 <- paste('Prespecified decomposition and phylogenetic analysis of the low-divergence LTR signature for all 34 TE species.',
 '(a) Each row shows total LTR hit bp / all RepeatMasker hit bp, below 5% LTR hit bp / all LTR hit bp, and below 5% LTR hit bp / all hit bp. Fractions satisfy the exact multiplicative identity; percentage values must be divided by 100 before multiplying. Axes differ among components. All available divergence bins enter the LTR denominator.',
 'The primary endpoint is the middle within-LTR fraction. It isolates profile shape from overall LTR representation in this pipeline, without estimating genomic occupancy or recent activity. The previous all-hit endpoint remains in the right panel and exported denominator sensitivities.',
 '(b) Primary endpoint against maximum SVL, with Pagel-lambda PGLS group curves from logit(fraction)~strategy+log10(SVL), limited to each group\'s observed SVL range; shading shows conditional 95% confidence intervals.',
 '(c) Raw-logit contrasts relative to aquatic eggs, with and without SVL. Groups are 7 aquatic-egg species, 24 terrestrial-egg species with aquatic larvae, and 3 direct developers. Intervals use t_(n-p) and condition on fitted lambda. The omnibus test compares strategy+SVL versus SVL by approximate ML likelihood ratio with 2 df; the unadjusted test is also exported.',
 'The three reproductive strategies were expert-coded after PCA inspection, and a group count is not a count of independent evolutionary transitions. No slope interaction, species filtering, divergence cutoff search or TE family screening was performed. All associations are exploratory; neither group separation nor low divergence establishes insertion timing or causation.',
 'Colors and species order match earlier TE figures; all 34 have maximum SVL. Trait uncertainty and uncertain source-tree calibration are not propagated. RepeatMasker aligned-hit denominators and dnaPipeTE quality denominators are kept distinct.')
collage12 <- p12a / p12b / p12c + plot_layout(heights = c(2.6,1,.8)) + plot_annotation(tag_levels = 'a')
save_candidate(12,'What drives the low-divergence LTR pattern?',collage12,
  list(a = list(plot = p12a,title = 'LTR fraction decomposition',width_mm = 180,height_mm = 230,legend = legend12,alt_text = 'Three aligned phylogenetic species tracks decompose the all-hit low-divergence LTR share into total LTR representation and within-LTR low-divergence fraction.'),
    b = list(plot = p12b,title = 'Within-LTR divergence signature and SVL',width_mm = 180,height_mm = 115,legend = legend12,alt_text = 'Thirty-four species and strategy-specific PGLS curves relate within-LTR low divergence to body size.'),
    c = list(plot = p12c,title = 'Strategy contrasts with phylogeny and SVL',width_mm = 180,height_mm = 90,legend = legend12,alt_text = 'Two reproductive-strategy contrasts versus aquatic eggs with and without SVL adjustment.')),
  180,340,'inputs/te_integration;S06;S12;S42;organismal_traits',legend12,
  'Thirty-four species show the exact LTR-profile decomposition alongside phylogenetic body-size and reproductive-strategy analyses.')
for (out in c(OUT11,OUT12)) capture.output(sessionInfo(),file = file.path(out,'session_info.txt'))
print(reg11 %>% filter(endpoint == 'within_ltr',term == 't') %>% select(fit_id,estimate,ci_low,ci_high,p_value,lambda))
print(omnibus); print(qc)
message('PASS: Figures 11 and 12, exact21 integration, full34 decomposition, all numerical assertions and independent GLS checks.')

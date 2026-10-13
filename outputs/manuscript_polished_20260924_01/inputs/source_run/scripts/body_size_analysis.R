#!/usr/bin/env Rscript
# Isolated, finite candidate comparison: organismal size and cellular traits.
# Run from manuscript_story_20260924_01; this script writes only new body-size files.
suppressPackageStartupMessages({
  library(ape); library(phylolm); library(phylopath)
  library(dplyr); library(readr); library(tibble); library(jsonlite)
})
options(warn = 1)
set.seed(20260924)
RUN <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(basename(RUN) == "manuscript_story_20260924_01")
ROOT <- normalizePath(file.path(RUN, "../.."), mustWork = TRUE)
SNAP <- file.path(RUN, "inputs", "body_size")
OUT <- file.path(RUN, "figures", "figure_8", "analysis")
dir.create(SNAP, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sha <- function(path) sub(" .*", "", system2("sha256sum", shQuote(path), stdout = TRUE))
sources <- c(
  traits24 = file.path(RUN, "inputs/datasets/Supplementary_Data_S04_path_analysis_traits24.csv"),
  canonical_traits24 = file.path(ROOT, "analyses/03_phylogenetic_path/data/path24_traits.csv"),
  organismal_traits = file.path(ROOT, "analyses/03_phylogenetic_path/data/organismal_traits.csv"),
  path24_time_tree = file.path(ROOT, "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk"),
  frozen_path24_tree = file.path(RUN, "inputs/datasets/trees/desmognathus_time_tree_path24_v1.nwk"),
  frozen_three_trait_ranking = file.path(RUN, "inputs/datasets/Supplementary_Data_S30_path_model_comparison.csv")
)
stopifnot(all(file.exists(sources)))
dest <- file.path(SNAP, paste0(names(sources), ifelse(grepl("\\.nwk$", sources), ".nwk", ".csv")))
names(dest) <- names(sources)
for (id in names(sources)) {
  if (!file.exists(dest[[id]])) stopifnot(file.copy(sources[[id]], dest[[id]], overwrite = FALSE))
  stopifnot(identical(sha(sources[[id]]), sha(dest[[id]])))
}
manifest <- tibble(source_id = names(sources), source_path = unname(sources),
  snapshot_path = unname(dest), sha256 = vapply(dest, sha, character(1)))
write_csv(manifest, file.path(SNAP, "MANIFEST.csv"))
read_snapshot <- function(id) read_csv(dest[[id]], show_col_types = FALSE, progress = FALSE)
traits <- read_snapshot("traits24")
canonical <- read_snapshot("canonical_traits24")
organismal <- read_snapshot("organismal_traits")
raw_columns <- c(b = "max_svl_mm", g = "genome_size_pg_fuscus_anchored", c = "cell_area_um2", n = "nucleus_area_um2")
labels <- c(b = "body_size", g = "genome_size", c = "cell_size", n = "nucleus_size")
trait_columns <- c("tree_tip", unname(raw_columns[-1]))
stopifnot(isTRUE(all.equal(as.data.frame(traits[, trait_columns]),
  as.data.frame(canonical[, trait_columns]), tolerance = 0, check.attributes = FALSE)))
stopifnot(nrow(traits) == 24L, !anyDuplicated(traits$tree_tip), !anyDuplicated(organismal$species))
joined <- traits %>% left_join(organismal %>% rename(tree_tip = species), by = "tree_tip")
stopifnot(nrow(joined) == 24L, !anyNA(joined[, unname(raw_columns)]),
  all(is.finite(as.matrix(joined[, unname(raw_columns)]))), all(as.matrix(joined[, unname(raw_columns)]) > 0))
tree <- read.tree(dest[["path24_time_tree"]])
frozen_tree <- read.tree(dest[["frozen_path24_tree"]])
stopifnot(setequal(tree$tip.label, joined$tree_tip), Ntip(tree) == 24L,
  isTRUE(all.equal.phylo(tree, frozen_tree, use.edge.length = TRUE)),
  all(is.finite(tree$edge.length)), all(tree$edge.length >= 0))
excluded <- c("aeneus", "wrighti")
stopifnot(setequal(joined$tree_tip[joined$development_mode == "direct_development"], excluded))
scaling <- bind_rows(lapply(names(raw_columns), function(v) {
  x <- log10(joined[[raw_columns[[v]]]])
  tibble(variable = v, trait = labels[[v]], raw_column = raw_columns[[v]],
    transformation = "log10", full24_mean = mean(x), full24_sd = sd(x), n_scaling = 24L)
}))
stopifnot(all(scaling$full24_sd > 0))
for (v in scaling$variable) {
  row <- scaling[scaling$variable == v, ]
  joined[[v]] <- (log10(joined[[row$raw_column]]) - row$full24_mean) / row$full24_sd
  stopifnot(abs(mean(joined[[v]])) < 1e-12, abs(sd(joined[[v]]) - 1) < 1e-12)
}
joined$included_full24 <- TRUE
joined$included_subset22 <- !joined$tree_tip %in% excluded
write_csv(joined, file.path(OUT, "traits_joined.csv"))
write_csv(scaling, file.path(OUT, "scaling.csv"))
panels <- list(full24 = joined, subset22 = joined[joined$included_subset22, ])
stopifnot(nrow(panels$subset22) == 22L,
  setequal(setdiff(panels$full24$tree_tip, panels$subset22$tree_tip), excluded))
warning_rows <- list()
capture_fit <- function(expr, panel, fit_id, source = "regression") {
  withCallingHandlers(expr, warning = function(w) {
    warning_rows[[length(warning_rows) + 1L]] <<- tibble(panel = panel,
      fit_id = fit_id, warning_source = source, warning = conditionMessage(w))
    # Deliberately do not muffle warnings: the build log retains every occurrence.
  })
}
prepared_panel <- function(panel) {
  raw <- panels[[panel]]
  phy <- keep.tip(tree, raw$tree_tip)
  d <- as.data.frame(raw[, names(raw_columns)])
  rownames(d) <- raw$tree_tip
  d <- d[phy$tip.label, , drop = FALSE]
  stopifnot(identical(rownames(d), phy$tip.label),
    identical(as.numeric(d$b), as.numeric(joined$b[match(phy$tip.label, joined$tree_tip)])))
  list(data = d, tree = phy)
}
equation_text <- function(f) paste(deparse(f), collapse = " ")
fit_regression <- function(formula, panel, relation_id, purpose = "regression") {
  prep <- prepared_panel(panel)
  fit <- capture_fit(phylolm(formula, data = prep$data, phy = prep$tree,
    model = "lambda", lower.bound = 1e-7, upper.bound = 1), panel, relation_id, purpose)
  sm <- as.data.frame(summary(fit)$coefficients)
  terms <- rownames(sm)
  df <- fit$n - fit$d
  crit <- qt(.975, df)
  response <- all.vars(formula)[[1]]
  rows <- tibble(panel = panel, relation_id = relation_id, equation = equation_text(formula),
    outcome = labels[[response]], response = response, term = terms,
    predictor = ifelse(terms == "(Intercept)", "(Intercept)", labels[terms]),
    standardized_coefficient = sm$Estimate, standard_error = sm$StdErr,
    ci_low = sm$Estimate - crit * sm$StdErr, ci_high = sm$Estimate + crit * sm$StdErr,
    t_value = sm$t.value, p_value = sm$p.value, fitted_phylogenetic_parameter = as.numeric(fit$optpar),
    lambda_boundary = as.numeric(fit$optpar) <= 1.01e-7 | as.numeric(fit$optpar) >= .999999,
    n_species = fit$n, n_regression_coefficients = fit$d, residual_df = df,
    r_squared = fit$r.squared, interval_method = "t: n minus regression coefficient count",
    standardization = "full24 log10 mean and sample SD")
  stopifnot(all(rows$ci_low <= rows$standardized_coefficient),
    all(rows$ci_high >= rows$standardized_coefficient),
    max(abs(rows$p_value - 2 * pt(-abs(rows$t_value), df))) < 1e-12,
    all.equal(as.numeric(fit$fitted.values), as.numeric(model.matrix(formula, prep$data) %*% coef(fit)), tolerance = 1e-10) == TRUE)
  list(fit = fit, rows = rows, data = prep$data)
}
specifications <- list(
  body_to_genome = g ~ b, body_to_cell = c ~ b, body_to_nucleus = n ~ b,
  genome_to_cell = c ~ g, cell_to_nucleus = n ~ c, genome_to_nucleus = n ~ g,
  genome_to_cell_body_adjusted = c ~ g + b,
  cell_to_nucleus_body_adjusted = n ~ c + b,
  genome_to_nucleus_body_adjusted = n ~ g + b,
  genome_and_cell_to_nucleus = n ~ g + c,
  genome_and_cell_to_nucleus_body_adjusted = n ~ g + c + b
)
regression_rows <- covariance_rows <- prediction_rows <- list()
for (panel in names(panels)) for (id in names(specifications)) {
  result <- fit_regression(specifications[[id]], panel, id)
  regression_rows[[length(regression_rows) + 1L]] <- result$rows
  fit <- result$fit
  vc <- fit$vcov
  idx <- expand.grid(term1 = rownames(vc), term2 = colnames(vc), stringsAsFactors = FALSE)
  covariance_rows[[length(covariance_rows) + 1L]] <- tibble(panel = panel, relation_id = id,
    term1 = idx$term1, term2 = idx$term2,
    covariance = mapply(function(a, b) vc[a, b], idx$term1, idx$term2))
  if (id %in% c("body_to_genome", "body_to_cell", "body_to_nucleus")) {
    response <- all.vars(specifications[[id]])[[1]]
    s_b <- scaling[scaling$variable == "b", ]
    s_y <- scaling[scaling$variable == response, ]
    body <- 10^seq(log10(min(panels[[panel]]$max_svl_mm)),
      log10(max(panels[[panel]]$max_svl_mm)), length.out = 121)
    b <- (log10(body) - s_b$full24_mean) / s_b$full24_sd
    X <- model.matrix(~b)
    eta <- as.numeric(X %*% coef(fit))
    se <- sqrt(rowSums((X %*% vc) * X))
    crit <- qt(.975, fit$n - fit$d)
    to_raw <- function(z) 10^(z * s_y$full24_sd + s_y$full24_mean)
    prediction_rows[[length(prediction_rows) + 1L]] <- tibble(panel = panel, relation_id = id,
      outcome = labels[[response]], response = response, max_svl_mm = body,
      standardized_body_size = b, standardized_fit = eta, standardized_se = se,
      fitted = to_raw(eta), ci_low = to_raw(eta - crit * se), ci_high = to_raw(eta + crit * se),
      n_species = fit$n, interval_type = "95% confidence interval for fixed-effect conditional median on raw scale")
    stopifnot(max(abs(log10(to_raw(eta)) - (eta * s_y$full24_sd + s_y$full24_mean))) < 1e-12)
  }
}
regressions <- bind_rows(regression_rows)
write_csv(regressions, file.path(OUT, "regression_coefficients.csv"))
write_csv(bind_rows(covariance_rows), file.path(OUT, "coefficient_covariance.csv"))
write_csv(bind_rows(prediction_rows), file.path(OUT, "pairwise_predictions.csv"))

four <- define_model_set(
  body_only = c(g ~ b, c ~ b, n ~ b),
  cell_bridge = c(g ~ b, c ~ b + g, n ~ b + c),
  nucleus_bridge = c(g ~ b, n ~ b + g, c ~ b + n),
  genome_bridge = c(g ~ b, c ~ b + g, n ~ b + g)
)
three <- define_model_set(
  independent = c(), genome_nucleus_only = c(n ~ g), nucleus_cell_only = c(c ~ n),
  genome_cell_only = c(c ~ g), nucleus_bridge = c(n ~ g, c ~ n),
  nucleus_collider = c(n ~ g + c), cell_bridge = c(c ~ g, n ~ c),
  cell_collider = c(c ~ g + n), genome_bridge = c(g ~ n, c ~ g),
  genome_collider = c(g ~ n + c)
)
sets <- list(body_adjusted_four = four, original_three = three)
candidate_rows <- list()
for (set_id in names(sets)) for (model in names(sets[[set_id]])) {
  dag <- sets[[set_id]][[model]]
  edges <- which(dag != 0, arr.ind = TRUE)
  if (nrow(edges)) candidate_rows[[length(candidate_rows) + 1L]] <- tibble(candidate_set = set_id,
    model = model, parent = rownames(dag)[edges[, 1]], child = colnames(dag)[edges[, 2]],
    parent_label = labels[rownames(dag)[edges[, 1]]], child_label = labels[colnames(dag)[edges[, 2]]])
}
write_csv(bind_rows(candidate_rows), file.path(OUT, "candidate_edges.csv"))
ranking_rows <- dsep_rows <- edge_rows <- list()
for (panel in names(panels)) for (set_id in names(sets)) {
  prep <- prepared_panel(panel)
  if (set_id == "original_three") prep$data <- prep$data[, c("g", "n", "c"), drop = FALSE]
  path <- capture_fit(phylo_path(sets[[set_id]], data = prep$data, tree = prep$tree,
    model = "lambda", lower.bound = 1e-7, upper.bound = 1), panel, set_id, "phylopath_runtime")
  if (length(path$warnings)) for (warning_text in unlist(path$warnings, use.names = FALSE)) {
    warning_rows[[length(warning_rows) + 1L]] <- tibble(panel = panel,
      fit_id = set_id, warning_source = "phylopath_stored", warning = warning_text)
    message("Stored phylopath warning [", panel, "/", set_id, "]: ", warning_text)
  }
  ranking <- as_tibble(as.data.frame(summary(path))) %>%
    mutate(candidate_set = set_id, panel = panel, rank = row_number(), n_species = nrow(prep$data),
      globally_supported = p >= .05, delta_competitive = delta_CICc <= 2,
      common_predictor_assumption = ifelse(set_id == "body_adjusted_four", "b -> g; b -> c; b -> n", "none"))
  stopifnot(abs(sum(ranking$w) - 1) < 1e-12, all(ranking$k > 0),
    all(ranking$q + 1 < ranking$n_species), all(diff(ranking$CICc) >= -1e-12),
    max(abs(ranking$delta_CICc - (ranking$CICc - min(ranking$CICc)))) < 1e-12)
  ranking_rows[[length(ranking_rows) + 1L]] <- ranking
  for (model in names(path$d_sep)) {
    ds <- path$d_sep[[model]]
    dsep_rows[[length(dsep_rows) + 1L]] <- tibble(candidate_set = set_id, panel = panel,
      model = model, independence_claim = ds$d_sep, component_p = ds$p,
      fitted_phylogenetic_parameter = ds$phylo_par, n_species = nrow(prep$data))
    # Fit component equations directly. phylopath::est_DAG would re-standardize
    # the subset and violate the fixed full24 scaling contract.
    dag <- sets[[set_id]][[model]]
    for (child in colnames(dag)) {
      parents <- rownames(dag)[dag[, child] != 0]
      if (!length(parents)) next
      formula <- reformulate(parents, response = child)
      component <- fit_regression(formula, panel, paste(set_id, model, child, sep = "/"), "path_edge")
      edge_rows[[length(edge_rows) + 1L]] <- component$rows %>%
        filter(term != "(Intercept)") %>% mutate(candidate_set = set_id, model = model,
          parent = term, child = child)
    }
    r <- ranking[ranking$model == model, ]
    stopifnot(abs(-2 * sum(log(ds$p)) - r$C) < 1e-10,
      abs(pchisq(r$C, df = 2 * r$k, lower.tail = FALSE) - r$p) < 1e-10)
  }
}
rankings <- bind_rows(ranking_rows)
write_csv(rankings, file.path(OUT, "path_model_comparison.csv"))
write_csv(bind_rows(dsep_rows), file.path(OUT, "path_dsep_tests.csv"))
write_csv(bind_rows(edge_rows), file.path(OUT, "path_edge_estimates.csv"))
warnings <- if (length(warning_rows)) bind_rows(warning_rows) else tibble(
  panel = character(), fit_id = character(), warning_source = character(), warning = character())
write_csv(warnings, file.path(OUT, "warnings.csv"))
frozen <- read_snapshot("frozen_three_trait_ranking") %>%
  filter(analysis_type == "primary", evolutionary_model == "lambda") %>%
  select(model, frozen_CICc = CICc, frozen_weight = w)
reproduction <- rankings %>% filter(panel == "full24", candidate_set == "original_three") %>%
  left_join(frozen, by = "model") %>%
  mutate(CICc_difference = CICc - frozen_CICc, weight_difference = w - frozen_weight)
stopifnot(nrow(reproduction) == 10L, !anyNA(reproduction$frozen_CICc),
  max(abs(reproduction$CICc_difference)) < 1e-5,
  max(abs(reproduction$weight_difference)) < 1e-6)
write_csv(reproduction, file.path(OUT, "full24_reproduction_check.csv"))
stopifnot(identical(manifest$sha256, vapply(dest, sha, character(1))))
verification <- list(passed = TRUE, full24_n = 24L, subset22_n = 22L,
  excluded_exactly = excluded, log10_scaling_reference = "full24 for both panels",
  all_species_positive_complete = TRUE, frozen_and_canonical_traits_identical = TRUE,
  frozen_and_canonical_trees_identical = TRUE, all_inputs_unchanged = TRUE,
  candidate_sets = c(body_adjusted_four = 4L, original_three = 10L),
  group_weight_sums = as.list(tapply(rankings$w, interaction(rankings$panel, rankings$candidate_set), sum)),
  maximum_full24_CICc_reproduction_difference = max(abs(reproduction$CICc_difference)),
  maximum_full24_weight_reproduction_difference = max(abs(reproduction$weight_difference)),
  coefficient_prediction_and_t_interval_checks = TRUE, n_warning_occurrences = nrow(warnings),
  R_version = R.version.string, phylolm_version = as.character(packageVersion("phylolm")),
  phylopath_version = as.character(packageVersion("phylopath")))
write_json(verification, file.path(OUT, "verification.json"), pretty = TRUE, auto_unbox = TRUE)
capture.output(sessionInfo(), file = file.path(OUT, "session_info.txt"))
writeLines(c(
  "# Body size and cellular traits: isolated exploratory analysis", "",
  "All new files are confined to this manuscript run. Existing figures, data and notes remain separate.", "",
  "## Data and transformations", "",
  "The 24 species come unchanged from the run's frozen Supplementary Data S04. Exact tree_tip matches join the canonical organismal_traits.csv species key; no taxon alias is newly inferred. All 24 have positive complete genome size, cell area, nucleus area and maximum SVL. Sources and SHA256 hashes are copied into inputs/body_size before fitting. Canonical and frozen traits and trees are checked for equality.",
  "SVL is collaborator-provided species maximum snout-vent length (mm), not the body length of the animals supplying microscopy. Genome size is the existing IOD-derived, D. fuscus assembly-anchored estimate. Cell and nucleus size are the existing species medians among retained upper-tail paired erythrocyte objects (areas in square micrometres), not unbiased species-wide means. No specimen-level body-size relationship is inferred.",
  "Every trait is log10-transformed and standardized using the full24 mean and sample SD. The subset22 retains exactly that scaling, so displayed coefficients use a common reference. scaling.csv preserves the constants. The comparison excludes exactly D. aeneus and D. wrighti, the two smallest sampled species and the two direct developers. This is a joint influence/life-history sensitivity, not an independent small-body group, not a test that distinguishes body size from development, and not an outlier cleanup. No threshold or species search was performed.", "",
  "## Regressions", "",
  "b = body size; g = genome size; c = cell area; n = nucleus area. Regressions in each panel are g~b, c~b, n~b; c~g, n~c, n~g; c~g+b, n~c+b, n~g+b; n~g+c; and n~g+c+b. Pagel-lambda regressions use phylolm with lambda constrained to [1e-7,1], maximum likelihood, no measurement-error model. Lambda bounds and warnings are retained. Confidence intervals are beta +/- qt(.975,n-p)*SE, with p the number of regression coefficients including the intercept, matching phylolm's t-test degrees of freedom. These intervals are conditional approximations and do not integrate uncertainty in estimated lambda.",
  "Pairwise curves use the fixed-effect regression mean on the standardized log scale, with variance X V X'. Curves and intervals are back-transformed without bias correction, representing the fitted conditional median on the original trait scale, not individual prediction intervals. Their range is restricted to the SVL range observed in each panel. coefficient_covariance.csv retains the full coefficient covariance matrices.", "",
  "## Four prespecified hypotheses", "",
  "All four hypotheses impose b->g, b->c and b->n. body_only adds no cellular edge; cell_bridge adds g->c and c->n; nucleus_bridge adds g->n and n->c; genome_bridge adds g->c and g->n. These are a finite exploratory comparison conditional on treating body size as a common predictor. The set is neither exhaustive over four-node DAGs nor a test proving body size exogenous. Each graph represents a Markov-equivalence class; the directions are hypotheses, not uniquely identified causal effects. A saturated graph is not included because it has no testable d-separation claim.",
  "phylopath::phylo_path and its summary implement d-separation, Fisher C, q and CICc. The package-selected d-separation regression orientation and estimated lambda are exported. k is the number of tested independence claims; q is the package's parameter count. Weights normalize within each candidate set and species panel, not across them. Absolute CICc values must not be compared between full24 and subset22 or between the four- and three-trait sets. Good global fit does not validate causal direction.",
  "Component path coefficients are fit with phylolm from each candidate's structural equations to preserve common full24 standardization. phylopath::est_DAG is intentionally not used for coefficients because it re-standardizes numeric variables in each subset. All candidate edges are exported, including weak or uncertain edges. No model-averaged causal effect is reported.",
  "The original ten testable three-trait equivalence classes are also rerun in each panel to distinguish changed species composition from body adjustment. The full24 ranking is numerically checked against frozen S30 before verification passes.", "",
  "## Limits and verification", "",
  "This is an exploratory, species-level analysis with n=24 or22. Genome-size calibration uncertainty, morphology measurement uncertainty, SVL uncertainty, tree uncertainty, and within-species sampling error are not propagated here. Direct-developer removal simultaneously changes body-size range, development composition and taxon sampling. It cannot attribute any change uniquely to one cause. Predictions and t intervals, exact membership, positive complete traits, common scaling, source preservation, tree matching, Fisher C, weights and full24 reproduction are checked by assertions in the analysis script. Warnings are printed and saved without user-code suppression.", "",
  "Reproduce from the run directory:", "",
  "    PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH Rscript --no-environ scripts/body_size_analysis.R", "",
  "CSV schemas: regression_coefficients.csv includes intercepts and uses standardized_coefficient, standard_error, ci_low, ci_high, p_value, fitted_phylogenetic_parameter, n_species and residual_df. path_edge_estimates.csv adds candidate_set, model, parent and child. path_model_comparison.csv identifies both candidate_set and panel. pairwise_predictions.csv contains max_svl_mm, fitted, ci_low and ci_high on raw scales."
), file.path(OUT, "README.md"))
print(rankings %>% select(panel, candidate_set, rank, model, k, q, p, CICc, delta_CICc, w))
message("PASS: body-size regressions, 4-hypothesis comparison and original 10-class sensitivity; all assertions passed.")

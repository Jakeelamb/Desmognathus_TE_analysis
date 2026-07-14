#!/usr/bin/env Rscript

# Parametric calibration for the corrected final-18 path-model candidate sets.
#
# Traits are simulated on the published time tree with independent Brownian
# residual processes. Signal scenarios use the standardized coefficients from
# the observed anchor fit; null scenarios generate mutually independent traits.
# This diagnoses model-selection behavior at n=18. It does not validate the
# biological meaning of relative IOD or turn observational paths into causes.

Sys.setenv(DESMO_PATH_LIBRARY_ONLY = "1")
source("scripts/processing/audit_corrected_path_models.R")

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
summarize_existing <- "--summarize-existing" %in% args
n_null <- if (quick) 5L else 200L
n_signal <- if (quick) 5L else 100L
seed <- 20260709L
effect_scales <- c(0.5, 1.0, 1.5)
simulation_families <- c("te_iod", "iod_morphology", "integrated")

edge_path <- file.path(
  project_root,
  "results/data/corrected/path_analysis/corrected_path_data_sensitivity_best_edges_analysis18_v1.csv"
)
tree_path <- file.path(
  project_root,
  "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk"
)
output_replicates <- file.path(
  output_dir,
  "corrected_path_simulation_replicates_analysis18_v1.csv"
)
output_summary <- file.path(
  output_dir,
  "corrected_path_simulation_calibration_analysis18_v1.csv"
)
output_coefficients <- file.path(
  output_dir,
  "corrected_path_simulation_coefficients_analysis18_v1.csv"
)
output_edge_replicates <- file.path(
  output_dir,
  "corrected_path_simulation_edge_replicates_analysis18_v1.csv"
)
output_edge_summary <- file.path(
  output_dir,
  "corrected_path_simulation_edge_calibration_analysis18_v1.csv"
)
output_manifest <- file.path(
  output_dir,
  "corrected_path_simulation_analysis18_v1.manifest.json"
)

observed_edges <- read_csv(edge_path, show_col_types = FALSE) %>%
  filter(
    .data$tree_id == "published_main",
    .data$morphology_estimator == anchor_morphology,
    .data$iod_subset == anchor_iod,
    .data$family %in% simulation_families
  ) %>%
  select("family", "best_model", "parent", "child", "coefficient")

expected_models <- c(
  te_iod = "mediated_evenness",
  iod_morphology = "iod_nucleus_cell_chain",
  integrated = "te_evenness_chain"
)

for (family in simulation_families) {
  family_edges <- observed_edges %>% filter(.data$family == !!family)
  if (!nrow(family_edges)) stop("No observed anchor edges for simulation family ", family)
  if (length(unique(family_edges$best_model)) != 1L || unique(family_edges$best_model) != expected_models[[family]]) {
    stop("Observed anchor winner changed for simulation family ", family)
  }
}

tree <- read.tree(tree_path)
tree$tip.label <- sub("^D\\.\\s*", "", tree$tip.label)
if (length(tree$tip.label) != 18L) stop("Calibration requires the published final-18 tree")

bm_trait <- function(tree) {
  values <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  values <- values[tree$tip.label]
  zscore(values)
}

scaled_beta <- function(beta, effect_scale) {
  pmax(-0.98, pmin(0.98, beta * effect_scale))
}

child_trait <- function(parent, beta, residual) {
  residual_sd <- sqrt(max(1 - beta^2, 0.01))
  zscore(beta * parent + residual_sd * residual)
}

edge_beta <- function(family, parent, child, effect_scale) {
  selected <- observed_edges %>%
    filter(.data$family == !!family, .data$parent == !!parent, .data$child == !!child)
  if (nrow(selected) != 1L) {
    stop("Expected exactly one anchor coefficient for ", family, ": ", parent, " -> ", child)
  }
  scaled_beta(selected$coefficient[[1]], effect_scale)
}

simulate_family <- function(tree, family, scenario, effect_scale) {
  traits <- switch(
    family,
    te_iod = c("ltr", "even", "iod"),
    iod_morphology = c("iod", "ns", "cs"),
    integrated = c("ltr", "even", "iod", "ns", "cs")
  )
  errors <- setNames(lapply(traits, function(x) bm_trait(tree)), traits)
  values <- errors
  if (scenario == "observed_chain") {
    if (family %in% c("te_iod", "integrated")) {
      values$ltr <- errors$ltr
      values$even <- child_trait(
        values$ltr, edge_beta(family, "ltr", "even", effect_scale), errors$even
      )
      values$iod <- child_trait(
        values$even, edge_beta(family, "even", "iod", effect_scale), errors$iod
      )
    }
    if (family == "iod_morphology") values$iod <- errors$iod
    if (family %in% c("iod_morphology", "integrated")) {
      values$ns <- child_trait(
        values$iod, edge_beta(family, "iod", "ns", effect_scale), errors$ns
      )
      values$cs <- child_trait(
        values$ns, edge_beta(family, "ns", "cs", effect_scale), errors$cs
      )
    }
  }
  data <- as.data.frame(values[traits], check.names = FALSE)
  rownames(data) <- tree$tip.label
  data[tree$tip.label, , drop = FALSE]
}

fit_simulation <- function(tree, family, scenario, effect_scale, replicate_id) {
  data <- simulate_family(tree, family, scenario, effect_scale)
  attempt <- tryCatch(
    {
      fit <- suppressWarnings(phylo_path(model_set_for_family(family), data = data, tree = tree))
      ranking <- as.data.frame(summary(fit), stringsAsFactors = FALSE)
      gate <- release_gate(ranking)
      top <- ranking[1, , drop = FALSE]
      expected <- if (scenario == "independent_null") {
        if (family == "iod_morphology") "morphology_null" else if (family == "integrated") "integrated_null" else "proxy_null"
      } else {
        expected_models[[family]]
      }
      selection <- tibble(
        family = family,
        scenario = scenario,
        effect_scale = effect_scale,
        replicate = replicate_id,
        seed = seed,
        n_species = nrow(data),
        expected_model = expected,
        top_model = as.character(top$model[[1]]),
        top_global_p = as.numeric(top$p[[1]]),
        top_cicc = as.numeric(top$CICc[[1]]),
        top_weight = as.numeric(top$w[[1]]),
        release_gate_status = gate,
        top_model_is_expected = .data$top_model == .data$expected_model,
        expected_unique_supported = .data$top_model_is_expected & gate == "supported_unique_top_model",
        false_supported_nonnull = scenario == "independent_null" &
          .data$top_model != .data$expected_model & is.finite(.data$top_global_p) & .data$top_global_p >= alpha,
        false_unique_nonnull = scenario == "independent_null" &
          .data$top_model != .data$expected_model & gate == "supported_unique_top_model",
        fit_success = TRUE,
        error = NA_character_
      )
      edge_rows <- tibble()
      if (scenario == "observed_chain") {
        fitted_true_dag <- est_DAG(
          fit$model_set[[expected]], fit$data, fit$tree, fit$model, fit$method
        )
        family_edges <- observed_edges %>% filter(.data$family == !!family)
        edge_rows <- bind_rows(lapply(seq_len(nrow(family_edges)), function(index) {
          parent <- family_edges$parent[[index]]
          child <- family_edges$child[[index]]
          truth <- edge_beta(family, parent, child, effect_scale)
          estimate <- as.numeric(fitted_true_dag$coef[parent, child])
          standard_error <- as.numeric(fitted_true_dag$se[parent, child])
          lower <- estimate - 1.96 * standard_error
          upper <- estimate + 1.96 * standard_error
          tibble(
            family = family,
            scenario = scenario,
            effect_scale = effect_scale,
            replicate = replicate_id,
            parent = parent,
            child = child,
            true_coefficient = truth,
            estimated_coefficient = estimate,
            standard_error = standard_error,
            coefficient_bias = estimate - truth,
            approx_ci_low = lower,
            approx_ci_high = upper,
            interval_covers_truth = is.finite(lower) & is.finite(upper) & lower <= truth & upper >= truth
          )
        }))
      }
      list(selection = selection, edges = edge_rows)
    },
    error = function(e) {
      selection <- tibble(
        family = family,
        scenario = scenario,
        effect_scale = effect_scale,
        replicate = replicate_id,
        seed = seed,
        n_species = length(tree$tip.label),
        expected_model = if (scenario == "independent_null") "null" else expected_models[[family]],
        top_model = NA_character_,
        top_global_p = NA_real_,
        top_cicc = NA_real_,
        top_weight = NA_real_,
        release_gate_status = "fit_failure",
        top_model_is_expected = FALSE,
        expected_unique_supported = FALSE,
        false_supported_nonnull = FALSE,
        false_unique_nonnull = FALSE,
        fit_success = FALSE,
        error = conditionMessage(e)
      )
      list(selection = selection, edges = tibble())
    }
  )
  attempt
}

wilson_interval <- function(successes, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- successes / n
  denominator <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denominator
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
  c(max(0, centre - half), min(1, centre + half))
}

if (summarize_existing) {
  if (!file.exists(output_replicates) || !file.exists(output_edge_replicates)) {
    stop("--summarize-existing requires existing replicate outputs")
  }
  replicates <- read_csv(output_replicates, show_col_types = FALSE) %>%
    mutate(
      false_unique_nonnull = .data$scenario == "independent_null" &
        .data$top_model != .data$expected_model &
        .data$release_gate_status == "supported_unique_top_model"
    )
  edge_replicates <- read_csv(output_edge_replicates, show_col_types = FALSE)
  null_counts <- replicates %>%
    filter(.data$scenario == "independent_null") %>%
    count(.data$family) %>%
    pull(.data$n)
  signal_counts <- replicates %>%
    filter(.data$scenario == "observed_chain") %>%
    count(.data$family, .data$effect_scale) %>%
    pull(.data$n)
  if (length(unique(null_counts)) != 1L || length(unique(signal_counts)) != 1L) {
    stop("Existing simulation design is unbalanced")
  }
  n_null <- as.integer(unique(null_counts))
  n_signal <- as.integer(unique(signal_counts))
} else {
  set.seed(seed)
  results <- list()
  for (family in simulation_families) {
    cat("Calibrating", family, "under independent null...\n")
    for (replicate_id in seq_len(n_null)) {
      results[[length(results) + 1L]] <- fit_simulation(
        tree, family, "independent_null", 0, replicate_id
      )
    }
    for (effect_scale in effect_scales) {
      cat("Calibrating", family, "at effect scale", effect_scale, "...\n")
      for (replicate_id in seq_len(n_signal)) {
        results[[length(results) + 1L]] <- fit_simulation(
          tree, family, "observed_chain", effect_scale, replicate_id
        )
      }
    }
  }
  replicates <- bind_rows(lapply(results, `[[`, "selection"))
  edge_replicates <- bind_rows(lapply(results, `[[`, "edges"))
}

summaries <- replicates %>%
  group_by(.data$family, .data$scenario, .data$effect_scale) %>%
  group_modify(function(data, key) {
    successful <- data %>% filter(.data$fit_success)
    n_success <- nrow(successful)
    exact <- sum(successful$top_model_is_expected)
    unique <- sum(successful$expected_unique_supported)
    false_nonnull <- sum(successful$false_supported_nonnull)
    false_unique <- sum(successful$false_unique_nonnull)
    target_count <- if (key$scenario[[1]] == "independent_null") false_unique else unique
    interval <- wilson_interval(target_count, n_success)
    top_counts <- sort(table(successful$top_model), decreasing = TRUE)
    tibble(
      n_attempted = nrow(data),
      n_successful = n_success,
      n_failures = sum(!data$fit_success),
      modal_top_model = if (length(top_counts)) names(top_counts)[[1]] else NA_character_,
      modal_top_model_rate = if (length(top_counts)) as.numeric(top_counts[[1]]) / n_success else NA_real_,
      expected_model_top_rate = if (n_success) exact / n_success else NA_real_,
      expected_unique_supported_rate = if (n_success) unique / n_success else NA_real_,
      false_nonnull_top_rate = if (n_success) false_nonnull / n_success else NA_real_,
      false_supported_nonnull_rate = if (n_success) false_nonnull / n_success else NA_real_,
      false_unique_nonnull_rate = if (n_success) false_unique / n_success else NA_real_,
      primary_calibration_rate = if (n_success) target_count / n_success else NA_real_,
      primary_rate_ci_low = interval[[1]],
      primary_rate_ci_high = interval[[2]],
      mean_top_weight = if (n_success) mean(successful$top_weight) else NA_real_,
      top_global_fit_rate = if (n_success) mean(successful$top_global_p >= alpha) else NA_real_
    )
  }) %>%
  ungroup()

simulation_coefficients <- tidyr::expand_grid(
  observed_edges,
  effect_scale = effect_scales
) %>%
  mutate(simulated_coefficient = scaled_beta(.data$coefficient, .data$effect_scale))

edge_summaries <- edge_replicates %>%
  group_by(
    .data$family, .data$scenario, .data$effect_scale,
    .data$parent, .data$child, .data$true_coefficient
  ) %>%
  summarise(
    n_successful = n(),
    mean_estimated_coefficient = mean(.data$estimated_coefficient),
    mean_bias = mean(.data$coefficient_bias),
    rmse = sqrt(mean((.data$estimated_coefficient - .data$true_coefficient)^2)),
    interval_coverage_rate = mean(.data$interval_covers_truth),
    mean_interval_width = mean(.data$approx_ci_high - .data$approx_ci_low),
    .groups = "drop"
  )

write_csv(replicates, output_replicates)
write_csv(summaries, output_summary)
write_csv(simulation_coefficients, output_coefficients)
write_csv(edge_replicates, output_edge_replicates)
write_csv(edge_summaries, output_edge_summary)
write_json(
  list(
    analysis_id = "corrected_path_simulation_analysis18_v1",
    seed = seed,
    n_null_replicates_per_family = n_null,
    n_signal_replicates_per_family_and_scale = n_signal,
    strict_null_calibration_metric = "unique supported non-null winner under independent null",
    effect_scales = effect_scales,
    families = simulation_families,
    simulation_process = "independent Brownian residuals on published final-18 time tree",
    observed_coefficients_from = sub(paste0(project_root, "/"), "", edge_path),
    absolute_genome_size_used = FALSE,
    causal_validation_claim_allowed = FALSE,
    outputs = lapply(
      c(
        output_replicates, output_summary, output_coefficients,
        output_edge_replicates, output_edge_summary
      ),
      function(path) sub(paste0(project_root, "/"), "", path)
    )
  ),
  output_manifest,
  pretty = TRUE,
  auto_unbox = TRUE
)
cat("Completed", nrow(replicates), "path-calibration simulations.\n")

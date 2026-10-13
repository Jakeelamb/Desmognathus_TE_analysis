from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUN_R = ROOT / "scripts/run_r.sh"

pytestmark = pytest.mark.skipif(
    shutil.which("conda") is None, reason="Requires the project R environment"
)

SETUP = r'''
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
})
source("analyses/01_transposable_elements/ltr_genome_pgls.R")
dataset_dir <- "Publication/datasets"
prepared <- prepare_ltr_genome_data(dataset_dir)
assert_close <- function(observed, expected, tolerance = 1e-9) {
  stopifnot(isTRUE(all.equal(as.numeric(observed), as.numeric(expected),
                            tolerance = tolerance, check.attributes = FALSE)))
}

# Reconstruct generalized least squares with ordinary matrix operations. This
# does not use the production fitter, its summary, or phylolm's covariance code.
independent_gls <- function(data, tree, lambda = 1, raw = FALSE) {
  covariance <- ape::vcv.phylo(tree)[data$species, data$species, drop = FALSE]
  covariance[row(covariance) != col(covariance)] <-
    covariance[row(covariance) != col(covariance)] * lambda
  predictor <- if (raw) data$genome_size_gbp else log10(data$genome_size_gbp)
  response <- if (raw) data$ratio_median else log10(data$ratio_median)
  design <- cbind(1, predictor)
  inverse_information <- solve(crossprod(design, solve(covariance, design)))
  coefficient <- inverse_information %*% crossprod(design, solve(covariance, response))
  residual <- response - drop(design %*% coefficient)
  residual_ss <- drop(crossprod(residual, solve(covariance, residual)))
  residual_df <- nrow(data) - ncol(design)
  coefficient_covariance <- inverse_information * residual_ss / residual_df
  standard_error <- sqrt(coefficient_covariance[2, 2])
  slope <- unname(coefficient[2])
  list(
    intercept = unname(coefficient[1]), slope = slope,
    standard_error = standard_error, t_value = slope / standard_error,
    p_value = 2 * pt(-abs(slope / standard_error), df = residual_df),
    ci_low = slope - qt(.975, df = residual_df) * standard_error,
    ci_high = slope + qt(.975, df = residual_df) * standard_error,
    covariance = coefficient_covariance,
    profile_nll = .5 * (as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
                       nrow(data) * log(residual_ss / nrow(data)))
  )
}
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", SETUP + code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=120,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_species_join_preserves_released_medians_and_measurement_intervals() -> None:
    run_r_assertions(
        r'''
summaries <- read.csv(file.path(dataset_dir,
  "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv"))
primary <- summaries[summaries$analysis_branch == "primary_iqr_filtered_coverage_ge_80pct", ]
traits <- read.csv(file.path(dataset_dir,
  "Supplementary_Data_S04_path_analysis_traits24.csv"))
elements <- read.csv(file.path(dataset_dir,
  "Supplementary_Data_S14_ltr_terminal_internal_elements.csv"))
elements <- elements[tolower(elements$primary_ltr_element) == "true", ]
data <- prepared$data
expected_species <- intersect(primary$species, traits$tree_tip)
expected_excluded <- setdiff(primary$species, traits$tree_tip)
stopifnot(nrow(data) == 19L, !anyDuplicated(data$species),
          setequal(data$species, expected_species),
          nrow(prepared$excluded) == 11L, !anyDuplicated(prepared$excluded$species),
          setequal(prepared$excluded$species, expected_excluded),
          setequal(expected_excluded, c("abditus", "adatsihi", "aureatus", "balsameus",
            "carolinensis", "cheaha", "conanti", "organi", "pascagoula", "santeetlah",
            "valentinei")))
summary_rows <- match(data$species, primary$species)
trait_rows <- match(data$species, traits$tree_tip)
assert_close(data$ratio_median, primary$ratio_median[summary_rows], 1e-14)
assert_close(data$median_bootstrap_ci_low,
             primary$median_bootstrap_ci_low[summary_rows], 1e-14)
assert_close(data$median_bootstrap_ci_high,
             primary$median_bootstrap_ci_high[summary_rows], 1e-14)
assert_close(data$genome_size_gbp,
             traits$genome_size_pg_fuscus_anchored[trait_rows] * .978, 1e-14)
assert_close(data$genome_size_ci_low_gbp,
             traits$genome_size_pg_ci_low[trait_rows] * .978, 1e-14)
assert_close(data$genome_size_ci_high_gbp,
             traits$genome_size_pg_ci_high[trait_rows] * .978, 1e-14)
assert_close(data$log10_ratio, log10(data$ratio_median), 1e-14)
assert_close(data$log10_genome, log10(data$genome_size_gbp), 1e-14)
raw_medians <- tapply(elements$ratio_terminal_internal_all_positions,
                     elements$species, median)[data$species]
log_medians <- tapply(log10(elements$ratio_terminal_internal_all_positions),
                     elements$species, median)[data$species]
assert_close(data$ratio_median, raw_medians, 1e-14)
stopifnot(any(abs(data$log10_ratio - log_medians) > 1e-6),
          all(is.finite(data$log10_ratio)), all(is.finite(data$log10_genome)),
          all(data$median_bootstrap_ci_low <= data$ratio_median),
          all(data$median_bootstrap_ci_high >= data$ratio_median))
'''
    )


def test_pruned_tree_retains_source_branch_lengths_and_canonical_species_order() -> None:
    run_r_assertions(
        r'''
source_tree <- ape::ladderize(ape::read.tree(file.path(dataset_dir,
  "trees/desmo900dated_source46.tre")), right = FALSE)
tip_height <- ape::node.height(source_tree)[seq_len(ape::Ntip(source_tree))]
source_order <- source_tree$tip.label[order(tip_height)]
expected_order <- source_order[source_order %in% prepared$data$species]
expected_tree <- ape::keep.tip(source_tree, prepared$data$species)
stopifnot(ape::Ntip(prepared$tree) == 19L,
          setequal(prepared$tree$tip.label, prepared$data$species),
          identical(as.character(prepared$canonical_order), source_order),
          identical(as.character(prepared$data$species), expected_order),
          ape::is.rooted(prepared$tree), ape::is.ultrametric(prepared$tree, tol = 1e-6))
assert_close(ape::vcv.phylo(prepared$tree)[expected_order, expected_order],
             ape::vcv.phylo(expected_tree)[expected_order, expected_order], 1e-12)
assert_close(ape::cophenetic.phylo(prepared$tree)[expected_order, expected_order],
             ape::cophenetic.phylo(expected_tree)[expected_order, expected_order], 1e-12)
'''
    )


def test_primary_and_sensitivity_inference_match_independent_gls() -> None:
    run_r_assertions(
        r'''
result <- fit_ltr_genome_pgls(prepared)
stopifnot(inherits(result$primary, "phylolm"), result$primary$n == 19L,
          nrow(result$summary) == 4L,
          setequal(result$summary$model,
            c("log_log_lambda", "log_log_BM", "log_log_OLS", "raw_raw_lambda")),
          all(result$summary$n_species == 19L))
primary <- result$summary[result$summary$model == "log_log_lambda", ]
assert_close(primary$lambda, result$primary$optpar)
stopifnot(primary$lambda >= 0, primary$lambda <= 1,
          identical(as.logical(primary$lambda_boundary),
                    primary$lambda < 1e-6 || primary$lambda > 1 - 1e-6))
for (model in c("log_log_lambda", "log_log_BM", "raw_raw_lambda")) {
  row <- result$summary[result$summary$model == model, ]
  lambda <- if (model == "log_log_BM") 1 else row$lambda
  expected <- independent_gls(prepared$data, prepared$tree, lambda,
                              raw = model == "raw_raw_lambda")
  for (column in c("intercept", "slope", "standard_error", "t_value",
                   "p_value", "ci_low", "ci_high")) {
    assert_close(row[[column]], expected[[column]], 1e-7)
  }
  if (model == "log_log_lambda") {
    assert_close(result$primary$vcov, expected$covariance, 1e-7)
  }
}
ols <- lm(log10_ratio ~ log10_genome, data = prepared$data)
row <- result$summary[result$summary$model == "log_log_OLS", ]
assert_close(row$slope, coef(ols)[[2]])
assert_close(row$standard_error, summary(ols)$coefficients[2, 2])
assert_close(row$p_value, summary(ols)$coefficients[2, 4])
assert_close(c(row$ci_low, row$ci_high), confint(ols)[2, ])

# Independently profile the lambda likelihood; a boundary fit must actually be
# supported by the likelihood rather than being a discarded optimizer warning.
criterion <- function(lambda) independent_gls(prepared$data, prepared$tree, lambda)$profile_nll
profile <- optimize(criterion, c(0, 1), tol = 1e-9)
candidates <- c(0, profile$minimum, 1)
best_lambda <- candidates[which.min(vapply(candidates, criterion, numeric(1)))]
assert_close(primary$lambda, best_lambda, 2e-4)
stopifnot(criterion(primary$lambda) <= min(vapply(candidates, criterion, numeric(1))) + 1e-6)
'''
    )


def test_leave_one_species_out_uses_all_nineteen_distinct_eighteen_species_panels() -> None:
    run_r_assertions(
        r'''
result <- fit_ltr_genome_pgls(prepared)
stopifnot(nrow(result$loo) == 19L, !anyDuplicated(result$loo$omitted_species),
          setequal(result$loo$omitted_species, prepared$data$species),
          all(result$loo$n_species == 18L))
for (index in seq_len(nrow(result$loo))) {
  row <- result$loo[index, ]
  data <- prepared$data[prepared$data$species != row$omitted_species, , drop = FALSE]
  tree <- ape::keep.tip(prepared$tree, data$species)
  expected <- independent_gls(data, tree, row$lambda)
  stopifnot(nrow(data) == 18L, !row$omitted_species %in% tree$tip.label,
            row$lambda >= 0, row$lambda <= 1)
  for (column in c("slope", "ci_low", "ci_high", "p_value")) {
    assert_close(row[[column]], expected[[column]], 1e-7)
  }
}
'''
    )


def test_prediction_is_gls_geometric_mean_with_coefficient_confidence_ribbon() -> None:
    run_r_assertions(
        r'''
result <- fit_ltr_genome_pgls(prepared)
prediction <- result$prediction
expected <- independent_gls(prepared$data, prepared$tree, result$primary$optpar)
log_x <- seq(log10(min(prepared$data$genome_size_gbp)),
             log10(max(prepared$data$genome_size_gbp)), length.out = 200L)
design <- cbind(1, log_x)
log_fit <- expected$intercept + expected$slope * log_x
# The diagonal is Var(X beta), without an added individual residual variance.
log_se <- sqrt(diag(design %*% expected$covariance %*% t(design)))
critical <- qt(.975, df = nrow(prepared$data) - ncol(design))
stopifnot(nrow(prediction) == 200L,
          all(diff(prediction$genome_size_gbp) > 0),
          all(is.finite(as.matrix(prediction))),
          all(prediction$ratio_ci_low < prediction$ratio_fit),
          all(prediction$ratio_ci_high > prediction$ratio_fit))
assert_close(prediction$genome_size_gbp, 10^log_x, 1e-14)
assert_close(range(prediction$genome_size_gbp),
             range(prepared$data$genome_size_gbp), 1e-14)
assert_close(prediction$ratio_fit, 10^log_fit, 1e-7)
assert_close(prediction$ratio_ci_low, 10^(log_fit - critical * log_se), 1e-7)
assert_close(prediction$ratio_ci_high, 10^(log_fit + critical * log_se), 1e-7)
assert_close(prediction$ratio_fit,
             sqrt(prediction$ratio_ci_low * prediction$ratio_ci_high), 1e-12)

# The released fit has lambda approximately zero, so an OLS substitution would
# be invisible there. A deterministic tree-structured response makes the same
# prediction contract sensitive to that error while retaining the real design.
structured <- prepared
observed_design <- cbind(1, structured$data$log10_genome)
residual_projection <- diag(nrow(structured$data)) - observed_design %*%
  solve(crossprod(observed_design)) %*% t(observed_design)
tree_covariance <- ape::vcv.phylo(structured$tree)[
  structured$data$species, structured$data$species]
phylogenetic_residual <- eigen(
  residual_projection %*% tree_covariance %*% residual_projection,
  symmetric = TRUE)$vectors[, 1]
structured$data$log10_ratio <- -.8 + .4 * structured$data$log10_genome +
  .15 * as.numeric(scale(phylogenetic_residual))
structured$data$ratio_median <- 10^structured$data$log10_ratio
structured_result <- fit_ltr_genome_pgls(structured)
structured_expected <- independent_gls(structured$data, structured$tree,
                                       structured_result$primary$optpar)
structured_log_fit <- structured_expected$intercept + structured_expected$slope * log_x
structured_log_se <- sqrt(diag(design %*% structured_expected$covariance %*% t(design)))
assert_close(structured_result$prediction$ratio_fit, 10^structured_log_fit, 1e-7)
assert_close(structured_result$prediction$ratio_ci_low,
             10^(structured_log_fit - critical * structured_log_se), 1e-7)
assert_close(structured_result$prediction$ratio_ci_high,
             10^(structured_log_fit + critical * structured_log_se), 1e-7)
ols_fit <- lm(log10_ratio ~ log10_genome, data = structured$data)
ordinary <- 10^predict(ols_fit, newdata = data.frame(log10_genome = log_x))
stopifnot(structured_result$primary$optpar > .5,
          max(abs(structured_result$prediction$ratio_fit - ordinary)) > .01)
'''
    )

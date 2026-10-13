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
suppressPackageStartupMessages({library(ape); library(dplyr); library(readr)})
source("scripts/publication/gbe_theme.R")
source("scripts/publication/figure6_panels.R")
source("scripts/publication/figure7_panels.R")
ROOT <- normalizePath(".")
DATA <- file.path(ROOT, "Publication/datasets")
expressions <- parse(file = "scripts/publication/build_gbe_figures.R")
evaluate_assignment <- function(name) {
  candidates <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name(name))
  }, expressions)
  stopifnot(length(candidates) == 1L)
  eval(candidates[[1]], envir = .GlobalEnv)
}
for (name in c("read_supp", "species_label", "source_tree_path", "source_tree",
               "source_tree_tip_y", "canonical_phylogenetic_order",
               "phylogenetic_species_order", "phylogenetic_species_levels",
               "path_traits", "morphology", "genome_size", "pair_data", "pgls",
               "pgls_panel", "genome_a", "genome_b", "genome_c",
               "path_models", "simulation", "path_a", "path_b")) {
  evaluate_assignment(name)
}
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", SETUP + code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_pgls_marks_intervals_and_back_transform_match_released_traits() -> None:
    run_r_assertions(
        r'''
stopifnot(nrow(pair_data) == 24L, n_distinct(pair_data$species) == 24L,
          sum(pair_data$genome_size_support == "One image and animal") == 3L,
          sum(pair_data$genome_size_support == "Quality-balance caution") == 4L,
          nrow(pgls) == 3L, all(pgls$pagel_lambda_at_lower_boundary))
panels <- list(genome_a, genome_b, genome_c)
stopifnot(identical(ggplot_build(genome_b)$layout$panel_params[[1]]$y$limits, c(180, 370)),
          identical(ggplot_build(genome_c)$layout$panel_params[[1]]$y$limits, c(180, 370)))
for (i in seq_along(panels)) {
  plot <- panels[[i]]
  stat <- pgls[i, ]
  data <- plot$data
  built <- ggplot_build(plot)
  observed <- built$data[[5]]
  stopifnot(nrow(observed) == 24L,
            identical(observed$x, data[[stat$left_metric]]),
            identical(observed$y, data[[stat$right_metric]]),
            identical(built$data[[3]]$ymin, data[[paste0(stat$right_metric, "_ci_low")]]),
            identical(built$data[[3]]$ymax, data[[paste0(stat$right_metric, "_ci_high")]]))
  x_prefix <- if (i <= 2) "genome_size_pg" else "nucleus_area_um2"
  stopifnot(identical(built$data[[4]]$xmin, data[[paste0(x_prefix, "_ci_low")]]),
            identical(built$data[[4]]$xmax, data[[paste0(x_prefix, "_ci_high")]]))
  fit <- plot$layers[[2]]$data
  z <- (log10(fit$x) - stat$predictor_log10_mean) / stat$predictor_log10_sd
  mean_z <- stat$pgls_intercept_z + stat$pgls_standardized_beta * z
  variance <- stat$pgls_var_intercept + 2 * z * stat$pgls_cov_intercept_slope +
    z^2 * stat$pgls_var_slope
  critical <- qt(.975, stat$pgls_degrees_freedom)
  expected <- 10^(stat$outcome_log10_mean + stat$outcome_log10_sd * mean_z)
  expected_low <- 10^(stat$outcome_log10_mean + stat$outcome_log10_sd *
                       (mean_z - critical * sqrt(variance)))
  expected_high <- 10^(stat$outcome_log10_mean + stat$outcome_log10_sd *
                        (mean_z + critical * sqrt(variance)))
  # Released transformation coefficients are rounded to eight decimals.
  stopifnot(max(abs(fit$fit / expected - 1)) < 1e-6,
            max(abs(fit$low / expected_low - 1)) < 1e-6,
            max(abs(fit$high / expected_high - 1)) < 1e-6,
            all(fit$low <= fit$fit), all(fit$fit <= fit$high),
            identical(plot$layers[[6]]$aes_params$family, "Arial"))
  for (layer in built$data[1:5]) {
    for (column in intersect(c("x", "y", "xmin", "xmax", "ymin", "ymax"), names(layer))) {
      stopifnot(all(is.finite(layer[[column]])))
    }
  }
}
'''
    )


def test_path_ranking_global_fit_and_simulations_preserve_frozen_results() -> None:
    run_r_assertions(
        r'''
released_models <- read_supp("Supplementary_Data_S30_path_model_comparison.csv") %>%
  filter(analysis_type == "primary") %>% arrange(delta_CICc)
released_simulation <- read_supp("Supplementary_Data_S36_path_simulation_calibration.csv")
stopifnot(nrow(path_models) == 10L, identical(path_models$model, released_models$model),
          identical(path_models$delta_CICc, released_models$delta_CICc),
          identical(path_models$globally_supported, released_models$globally_supported),
          path_models$model[[1]] == "cell_bridge", sum(path_models$globally_supported) == 1L)
a <- ggplot_build(path_a)
stopifnot(identical(a$data[[2]]$x, released_models$delta_CICc),
          identical(a$data[[2]]$y, structure(as.numeric(10:1), class = c("mapped_discrete", "numeric"))),
          identical(a$data[[1]]$xintercept, 2),
          identical(a$data[[2]]$colour,
                    ifelse(path_models$globally_supported, GBE_COLORS[["blue"]], GBE_COLORS[["vermillion"]])))
b <- ggplot_build(path_b)
stopifnot(nrow(simulation) == 8L, all(simulation$n_simulations == 100L),
          identical(b$data[[1]]$y, released_simulation$true_class_top_rate),
          identical(b$layout$panel_params[[1]]$x$limits, c(0, 1)),
          all(is.finite(b$data[[1]]$x)), all(is.finite(b$data[[1]]$y)))
'''
    )


def test_standalone_sensitivity_panels_preserve_data_scales_and_restore_legend() -> None:
    run_r_assertions(
        r'''
axis_state <- function(built, axis) {
  scale <- built$layout$panel_params[[1]][[axis]]
  list(limits = scale$limits, breaks = scale$get_breaks(), labels = scale$get_labels())
}
groups <- list(list(a = genome_a, b = genome_b, c = genome_c),
               list(a = path_a, b = path_b))
helpers <- list(figure6_review_panels, figure7_review_panels)
for (i in seq_along(groups)) {
  original <- groups[[i]]
  clones <- helpers[[i]](original)
  for (tag in names(original)) {
    before <- ggplot_build(original[[tag]])
    after <- ggplot_build(clones[[tag]])
    stopifnot(identical(before$data, after$data),
              identical(axis_state(before, "x"), axis_state(after, "x")),
              identical(axis_state(before, "y"), axis_state(after, "y")),
              identical(clones[[tag]]$labels$tag, tag),
              identical(clones[[tag]]$theme$plot.title$family, "Arial"))
  }
  for (bad in list(unname(original), rev(original), original[-1])) {
    stopifnot(inherits(try(helpers[[i]](bad), silent = TRUE), "try-error"))
  }
}
stopifnot(inherits(figure6_review_panels(groups[[1]])$b$guides$guides$shape, "GuideLegend"),
          identical(genome_b$guides$guides$shape, "none"))
'''
    )

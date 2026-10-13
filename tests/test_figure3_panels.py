from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_R = ROOT / "scripts/run_r.sh"


# Evaluate the actual preparation and ggplot expressions without rendering PDFs
# or executing other figures. The source data and species tree are the release.
FIGURE3_SETUP = r'''
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(readr)
})
source("scripts/publication/gbe_theme.R")
ROOT <- normalizePath(".")
DATA <- file.path(ROOT, "Publication/datasets")
expressions <- parse(file = "scripts/publication/build_gbe_figures.R")
assignment_index <- function(name) {
  indices <- which(vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name(name))
  }, logical(1)))
  stopifnot(length(indices) == 1L)
  indices
}
evaluate_assignment <- function(name) {
  eval(expressions[[assignment_index(name)]], envir = .GlobalEnv)
}
for (name in c("read_supp", "species_label", "source_tree_path", "source_tree",
               "source_tree_tip_y", "canonical_phylogenetic_order",
               "phylogenetic_species_order", "phylogenetic_species_levels",
               "te_strategy_levels", "te_strategy_colors", "te_strategy_shapes",
               "te_strategy_labels", "te_annotations")) {
  evaluate_assignment(name)
}
for (index in seq.int(assignment_index("landscape"),
                     assignment_index("Figure_3") - 1L)) {
  eval(expressions[[index]], envir = .GlobalEnv)
}
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", FIGURE3_SETUP + code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_profiles_preserve_all_hit_bp_denominator_zeros_and_equal_species_weights() -> None:
    run_r_assertions(
        r'''
raw <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
key <- c("species", "order", "divergence_bin_start_pct")
stopifnot(n_distinct(raw$species) == 34L,
          !anyDuplicated(raw[key]),
          max(abs(raw$percent_species_hit_bp -
                  100 * raw$hit_bp / raw$species_total_hit_bp)) < 1e-10)
totals <- raw %>% group_by(species) %>% summarise(
  total = sum(hit_bp), denominator = first(species_total_hit_bp), .groups = "drop"
)
stopifnot(identical(totals$total, totals$denominator))

stopifnot(nrow(landscape_profiles) == 34L * 3L * 50L,
          !anyDuplicated(landscape_profiles[key]),
          identical(sort(unique(landscape_profiles$divergence_bin_start_pct)),
                    as.numeric(0:49)),
          all(is.finite(landscape_profiles$percent_species_hit_bp)))
observed <- raw %>% filter(order %in% c("LTR", "LINE", "TIR"))
paired <- inner_join(landscape_profiles, observed, by = key,
                     suffix = c("_plotted", "_source"))
stopifnot(nrow(paired) == nrow(observed),
          max(abs(paired$percent_species_hit_bp_plotted -
                  paired$percent_species_hit_bp_source)) < 1e-12)
absent <- anti_join(landscape_profiles, observed, by = key)
stopifnot(nrow(absent) > 0L, all(absent$percent_species_hit_bp == 0))

group_sizes <- setNames(c(34L, 7L, 24L, 3L), profile_group_levels)
strategy_by_group <- setNames(c(NA_character_, te_strategy_levels), profile_group_levels)
for (group in profile_group_levels) {
  members <- landscape_profiles
  if (!is.na(strategy_by_group[[group]])) {
    members <- filter(members, reproductive_strategy == strategy_by_group[[group]])
  }
  expected <- members %>% group_by(order, divergence_bin_start_pct) %>%
    summarise(expected = sum(percent_species_hit_bp) / group_sizes[[group]],
              .groups = "drop")
  actual <- filter(profile_means, profile_group == group)
  comparison <- inner_join(actual, expected, by = c("order", "divergence_bin_start_pct"))
  stopifnot(nrow(comparison) == 3L * 50L,
            all(comparison$n_species == group_sizes[[group]]),
            max(abs(comparison$mean_percent - comparison$expected)) < 1e-12)
}

# Low-divergence values sum bins 0-4 against the same all-order denominator.
expected_young <- raw %>% filter(order %in% c("LTR", "LINE", "TIR"),
                               divergence_bin_start_pct < 5) %>%
  group_by(species, order) %>%
  summarise(expected = 100 * sum(hit_bp) / first(species_total_hit_bp),
            .groups = "drop")
comparison <- left_join(young_repeat_matrix, expected_young, by = c("species", "order"))
stopifnot(nrow(comparison) == 34L * 3L,
          !anyDuplicated(comparison[c("species", "order")]),
          max(abs(comparison$young_repeat_percent - replace_na(comparison$expected, 0))) < 1e-10)
'''
    )


def test_real_panels_keep_shared_axes_all_values_and_canonical_species_order() -> None:
    run_r_assertions(
        r'''
profiles <- ggplot_build(repeat_a)
stopifnot(nrow(profiles$layout$layout) == 3L,
          identical(as.character(profiles$layout$layout$order), c("LTR", "LINE", "TIR")),
          identical(repeat_a$scales$get_scales("x")$limits, c(0, 50)),
          identical(repeat_a$scales$get_scales("y")$limits, c(0, 3)),
          nrow(profiles$data[[1]]) == 34L * 3L * 50L,
          nrow(profiles$data[[2]]) == 3L * 3L * 50L)
for (axis in profiles$layout$panel_params) {
  stopifnot(identical(axis$x$limits, c(0, 50)),
            identical(axis$y$limits, c(0, 3)))
}
for (layer in profiles$data) {
  stopifnot(all(is.finite(layer$x)), all(is.finite(layer$y)))
}
for (panel in 1:3) {
  individuals <- filter(profiles$data[[1]], as.integer(PANEL) == panel)
  means <- filter(profiles$data[[2]], as.integer(PANEL) == panel)
  stopifnot(n_distinct(individuals$group) == 34L,
            n_distinct(means$group) == 3L,
            n_distinct(means$colour) == 3L,
            n_distinct(means$linetype) == 3L)
  expected <- filter(profile_means, as.character(order) == profile_orders[[panel]],
                     as.character(profile_group) != "All species (n = 34)")
  stopifnot(isTRUE(all.equal(sort(means$y), sort(expected$mean_percent))))
  for (group in profile_group_levels[-1]) {
    curve <- filter(means, colour == profile_group_colors[[.env$group]])
    group_expected <- filter(expected, as.character(profile_group) == group)
    stopifnot(nrow(curve) == 50L,
              all(curve$linetype == profile_group_linetypes[[group]]),
              isTRUE(all.equal(sort(curve$y), sort(group_expected$mean_percent))))
  }
}

low_divergence <- ggplot_build(repeat_b)
expected_species <- paste("D.", canonical_phylogenetic_order[
  canonical_phylogenetic_order %in% unique(landscape$species)
])
stopifnot(nrow(low_divergence$layout$layout) == 3L,
          identical(levels(young_repeat_matrix$species_label), expected_species),
          identical(repeat_b$scales$get_scales("x")$limits, c(0, 6)))
for (axis in low_divergence$layout$panel_params) {
  stopifnot(identical(axis$x$limits, c(0, 6)),
            identical(axis$y$get_labels(), expected_species))
}
point_index <- which(vapply(repeat_b$layers,
                           function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
stopifnot(length(point_index) == 1L)
points <- low_divergence$data[[point_index]]
stopifnot(nrow(points) == 102L,
          isTRUE(all.equal(sort(points$x), sort(young_repeat_matrix$young_repeat_percent))))
for (panel in 1:3) {
  species_points <- filter(points, as.integer(PANEL) == panel)
  species_points <- species_points[order(as.numeric(species_points$y)), ]
  expected <- filter(young_repeat_matrix, as.character(order) == profile_orders[[panel]])
  expected <- expected[order(as.numeric(expected$species_label)), ]
  strategies <- as.character(expected$reproductive_strategy)
  stopifnot(nrow(species_points) == 34L,
            identical(sort(as.numeric(species_points$y)), as.numeric(1:34)),
            identical(species_points$x, expected$young_repeat_percent),
            identical(species_points$colour, unname(te_strategy_colors[strategies])),
            identical(species_points$fill, unname(te_strategy_colors[strategies])),
            identical(species_points$shape, unname(te_strategy_shapes[strategies])))
}
for (layer in low_divergence$data) {
  for (position in intersect(c("x", "xend", "y", "yend"), names(layer))) {
    stopifnot(all(is.finite(layer[[position]])))
  }
}
text_style <- ggplot2::calc_element("axis.text.y", repeat_b$theme)
stopifnot(identical(text_style$family, "Arial"), identical(text_style$face, "italic"))
'''
    )


def test_isolated_real_panels_preserve_layers_scales_and_source_plots() -> None:
    run_r_assertions(
        r'''
source("scripts/publication/figure3_panels.R")
panels <- list(a = repeat_a, b = repeat_b)
source_state <- function(plot) list(data = plot$data, labels = plot$labels, theme = plot$theme)
axis_state <- function(built) lapply(built$layout$panel_params, function(panel) {
  lapply(panel[c("x", "y")], function(axis) {
    list(limits = axis$limits, breaks = axis$get_breaks(), labels = axis$get_labels())
  })
})
before <- lapply(panels, ggplot_build)
states <- lapply(panels, source_state)
isolated <- figure3_review_panels(panels)
stopifnot(identical(names(isolated), c("a", "b")))
for (tag in names(panels)) {
  after <- ggplot_build(isolated[[tag]])
  stopifnot(identical(before[[tag]]$data, after$data),
            identical(axis_state(before[[tag]]), axis_state(after)),
            identical(source_state(panels[[tag]]), states[[tag]]),
            identical(isolated[[tag]]$labels$tag, tag))
  for (aesthetic in c("colour", "fill", "shape", "linetype")) {
    original_scale <- before[[tag]]$plot$scales$get_scales(aesthetic)
    isolated_scale <- after$plot$scales$get_scales(aesthetic)
    if (!is.null(original_scale)) {
      stopifnot(identical(original_scale$map(original_scale$get_limits()),
                          isolated_scale$map(isolated_scale$get_limits())))
    }
  }
}
isolated$a <- isolated$a + labs(title = "Standalone review only")
stopifnot(identical(source_state(panels$a), states$a))
for (invalid in list(panels[2:1], panels[1], unname(panels))) {
  error <- tryCatch({figure3_review_panels(invalid); NULL}, error = conditionMessage)
  stopifnot(!is.null(error))
}
'''
    )

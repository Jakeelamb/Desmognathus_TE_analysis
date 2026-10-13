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

# Evaluate the production data preparation and ggplots against released inputs,
# without running another figure or creating exports.
FIGURE4_SETUP = r'''
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
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
  stopifnot(length(indices) >= 1L)
  indices[[1]]
}
evaluate_assignment <- function(name) {
  eval(expressions[[assignment_index(name)]], envir = .GlobalEnv)
}
for (name in c("read_supp", "species_label", "source_tree_path", "source_tree",
               "source_tree_tip_y", "canonical_phylogenetic_order",
               "phylogenetic_species_order", "prune_source_tree")) {
  evaluate_assignment(name)
}
for (index in seq.int(assignment_index("ltr"),
                     assignment_index("Figure_4") - 1L)) {
  eval(expressions[[index]], envir = .GlobalEnv)
}
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", FIGURE4_SETUP + code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_real_elements_and_summaries_keep_values_on_canonical_tree_rows() -> None:
    run_r_assertions(
        r'''
released_elements <- read_supp("Supplementary_Data_S14_ltr_terminal_internal_elements.csv") %>%
  filter(primary_ltr_element)
released_summaries <- read_supp("Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv") %>%
  filter(analysis_branch == "primary_iqr_filtered_coverage_ge_80pct")
stopifnot(nrow(ltr_elements) == 380L, n_distinct(ltr_elements$species) == 30L,
          nrow(ltr) == 30L, n_distinct(ltr$species) == 30L,
          identical(ltr_elements$element_id, released_elements$element_id),
          identical(ltr_elements$estimate,
                    released_elements$log2_ratio_terminal_internal_all_positions),
          identical(ltr$species, released_summaries$species),
          identical(ltr$estimate, log2(released_summaries$ratio_median)),
          identical(ltr$ci_low, log2(released_summaries$median_bootstrap_ci_low)),
          identical(ltr$ci_high, log2(released_summaries$median_bootstrap_ci_high)))

expected_order <- canonical_phylogenetic_order[
  canonical_phylogenetic_order %in% released_summaries$species
]
stopifnot(identical(ltr_tree_tips_sorted$species, expected_order),
          identical(ltr_tree_tips_sorted$y, as.numeric(1:30)),
          identical(ltr_tree_tips_sorted$species_label, paste("D.", expected_order)),
          identical(ltr_tree_tips_sorted$y,
                    as.numeric(node.height(ltr_tree)[match(expected_order, ltr_tree$tip.label)])),
          nrow(ltr_tree_horizontal) == 58L, nrow(ltr_tree_vertical) == 29L)
expected_y <- setNames(as.numeric(1:30), expected_order)
stopifnot(identical(ltr_elements$y, unname(expected_y[ltr_elements$species])),
          identical(ltr$y, unname(expected_y[ltr$species])))
element_counts <- table(ltr_elements$species)
stopifnot(identical(as.integer(element_counts[ltr$species]), as.integer(ltr$n_elements)))
'''
    )


def test_real_marks_are_unclipped_with_vertical_only_jitter_and_aligned_counts() -> None:
    run_r_assertions(
        r'''
built <- ggplot_build(ltr_distribution)
raw_index <- which(vapply(ltr_distribution$layers, function(layer) {
  inherits(layer$geom, "GeomPoint")
}, logical(1)))
median_index <- which(vapply(ltr_distribution$layers, function(layer) {
  inherits(layer$geom, "GeomSegment") &&
    identical(rlang::as_label(layer$mapping$x), "estimate") &&
    identical(rlang::as_label(layer$mapping$xend), "estimate")
}, logical(1)))
interval_index <- which(vapply(ltr_distribution$layers, function(layer) {
  inherits(layer$geom, "GeomSegment") &&
    identical(rlang::as_label(layer$mapping$x), "ci_low") &&
    identical(rlang::as_label(layer$mapping$xend), "ci_high")
}, logical(1)))
stopifnot(length(raw_index) == 1L, length(median_index) == 1L,
          length(interval_index) == 1L)
raw <- built$data[[raw_index]]
medians <- built$data[[median_index]]
intervals <- built$data[[interval_index]]
stopifnot(nrow(raw) == 380L, nrow(medians) == 30L, nrow(intervals) == 30L,
          identical(raw$x, ltr_elements$estimate),
          all(abs(raw$y - ltr_elements$y) <= .15),
          any(raw$y != ltr_elements$y),
          identical(medians$x, ltr$estimate), identical(medians$xend, ltr$estimate),
          identical(medians$y, ltr$y - .25), identical(medians$yend, ltr$y + .25),
          identical(intervals$x, ltr$ci_low), identical(intervals$xend, ltr$ci_high),
          identical(intervals$y, ltr$y), identical(intervals$yend, ltr$y),
          all(medians$colour == GBE_COLORS[["blue"]]),
          all(intervals$colour == GBE_COLORS[["blue"]]),
          all(intervals$linewidth < medians$linewidth))
for (layer in list(raw, medians, intervals)) {
  for (position in intersect(c("x", "xend", "y", "yend"), names(layer))) {
    stopifnot(all(is.finite(layer[[position]])))
  }
}
stopifnot(identical(raw, ggplot_build(ltr_distribution)$data[[raw_index]]))

axes <- built$layout$panel_params[[1]]
tree_axes <- ggplot_build(ltr_tree_left)$layout$panel_params[[1]]
stopifnot(identical(axes$x$limits, log2(c(.025, 1.6))),
          identical(axes$x$get_breaks(), log2(c(.025, .05, .1, .25, .5, 1, 1.5))),
          identical(axes$x$get_labels(), c("0.025", "0.05", "0.1", "0.25", "0.5", "1", "1.5")),
          identical(axes$y$limits, c(.5, 30.5)),
          identical(axes$y$limits, tree_axes$y$limits),
          identical(axes$y$get_breaks(), tree_axes$y$get_breaks()),
          identical(tree_axes$y$get_labels(), ltr_tree_tips_sorted$species_label),
          identical(axes$y.sec$get_breaks(), as.numeric(1:30)),
          identical(as.character(axes$y.sec$get_labels()),
                    ltr_genome_labels$label),
          grepl("log2", ltr_distribution$labels$x, fixed = TRUE),
          grepl("1 = equal depth", ltr_distribution$labels$x, fixed = TRUE))
right_style <- ggplot2::calc_element("axis.text.y.right", ltr_distribution$theme)
tree_style <- ggplot2::calc_element("axis.text.y.right", ltr_tree_left$theme)
stopifnot(inherits(right_style, "element_text"), identical(right_style$family, "Arial"),
          identical(tree_style$family, "Arial"), identical(tree_style$face, "italic"))
'''
    )


def test_genome_labels_join_released_traits_without_imputing_missing_species() -> None:
    run_r_assertions(
        r'''
released_traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
source_rows <- match(ltr_tree_tips_sorted$species, released_traits$tree_tip)
expected_gbp <- released_traits$genome_size_pg_fuscus_anchored[source_rows] * .978
expected_counts <- ltr$n_elements[match(ltr_tree_tips_sorted$species, ltr$species)]
expected_labels <- paste0(
  ifelse(is.na(expected_gbp), "n/a", sprintf("%.1f", expected_gbp)),
  " / ", expected_counts
)
stopifnot(nrow(ltr_genome_labels) == 30L,
          identical(ltr_genome_labels$species, ltr_tree_tips_sorted$species),
          identical(ltr_genome_labels$y, ltr_tree_tips_sorted$y),
          identical(ltr_genome_labels$genome_size_gbp, expected_gbp),
          identical(ltr_genome_labels$n_elements, expected_counts),
          sum(is.finite(ltr_genome_labels$genome_size_gbp)) == 19L,
          sum(is.na(ltr_genome_labels$genome_size_gbp)) == 11L,
          identical(is.na(ltr_genome_labels$genome_size_gbp), is.na(source_rows)),
          identical(ltr_genome_labels$label, expected_labels),
          sum(startsWith(ltr_genome_labels$label, "n/a / ")) == 11L)
axes <- ggplot_build(ltr_distribution)$layout$panel_params[[1]]
stopifnot(identical(as.character(axes$y.sec$get_labels()), expected_labels),
          identical(axes$y.sec$get_breaks(), ltr_genome_labels$y))
axis_title <- ltr_distribution$scales$get_scales("y")$secondary.axis$name
stopifnot(grepl("Gbp*", axis_title, fixed = TRUE),
          grepl("Estimated genome size", axis_title, fixed = TRUE),
          grepl("elements (n)", axis_title, fixed = TRUE),
          grepl("assembly-anchored estimate", ltr_distribution$labels$caption, fixed = TRUE),
          grepl("n/a", ltr_distribution$labels$caption, fixed = TRUE),
          grepl("unavailable", ltr_distribution$labels$caption, fixed = TRUE))
'''
    )


def test_horizontal_violins_use_one_bandwidth_and_stay_within_species_ranges() -> None:
    run_r_assertions(
        r'''
violin_index <- which(vapply(ltr_distribution$layers, function(layer) {
  inherits(layer$geom, "GeomViolin")
}, logical(1)))
reference_index <- which(vapply(ltr_distribution$layers, function(layer) {
  inherits(layer$geom, "GeomVline")
}, logical(1)))
stopifnot(length(violin_index) == 1L, length(reference_index) == 1L,
          !any(vapply(ltr_distribution$layers, function(layer) {
            inherits(layer$geom, "GeomRect")
          }, logical(1))))
violin <- ltr_distribution$layers[[violin_index]]
stopifnot(identical(rlang::as_label(violin$mapping$group), "species"),
          identical(violin$stat_params$orientation, "y"),
          identical(violin$stat_params$trim, TRUE),
          identical(violin$stat_params$scale, "width"),
          is.numeric(violin$stat_params$bw), length(violin$stat_params$bw) == 1L,
          identical(violin$stat_params$bw, stats::bw.nrd0(ltr_elements$estimate)))

built <- ggplot_build(ltr_distribution)
density <- built$data[[violin_index]]
stopifnot(n_distinct(density$group) == 30L, all(density$flipped_aes),
          all(is.finite(density$x)), all(is.finite(density$y)),
          all(is.finite(density$density)), all(density$density >= 0),
          all(density$violinwidth >= 0), all(density$violinwidth <= 1),
          all(density$ymax - density$ymin <= .8 + 1e-12))
for (group in split(density, density$group)) {
  row_y <- unique(group$y)
  stopifnot(length(row_y) == 1L, row_y %in% ltr_tree_tips_sorted$y)
  observations <- ltr_elements$estimate[ltr_elements$y == row_y]
  stopifnot(all(group$x >= min(observations) - 1e-12),
            all(group$x <= max(observations) + 1e-12),
            abs(min(group$x) - min(observations)) < 1e-12,
            abs(max(group$x) - max(observations)) < 1e-12,
            all(group$n == length(observations)),
            abs(max(group$violinwidth) - 1) < 1e-12,
            all(abs(group$ymin - (row_y - .4)) < 1e-12),
            all(abs(group$ymax - (row_y + .4)) < 1e-12))
}
stopifnot(setequal(unique(density$y), ltr_tree_tips_sorted$y))
reference <- built$data[[reference_index]]
stopifnot(nrow(reference) == 1L, identical(reference$xintercept, 0),
          reference$linetype %in% c(3, "dotted"))
'''
    )


def test_isolated_panels_preserve_marks_and_axes_and_restore_species_labels() -> None:
    run_r_assertions(
        r'''
source("scripts/publication/figure4_panels.R")
panels <- list(a = ltr_tree_left, b = ltr_distribution)
source_state <- function(plot) list(data = plot$data, labels = plot$labels, theme = plot$theme)
axis_state <- function(built, axis) {
  position <- built$layout$panel_params[[1]][[axis]]
  list(limits = position$limits, breaks = position$get_breaks(),
       labels = position$get_labels(), range = position$continuous_range)
}
before <- lapply(panels, ggplot_build)
states <- lapply(panels, source_state)
isolated <- figure4_review_panels(panels, ltr_tree_tips_sorted)
stopifnot(identical(names(isolated), c("a", "b")))
for (tag in names(panels)) {
  after <- ggplot_build(isolated[[tag]])
  stopifnot(identical(before[[tag]]$data, after$data),
            identical(axis_state(before[[tag]], "x"), axis_state(after, "x")),
            identical(source_state(panels[[tag]]), states[[tag]]),
            identical(isolated[[tag]]$labels$tag, tag))
  old_y <- axis_state(before[[tag]], "y")
  new_y <- axis_state(after, "y")
  stopifnot(identical(old_y[names(old_y) != "labels"], new_y[names(new_y) != "labels"]))
  if (tag == "a") {
    stopifnot(identical(old_y, new_y))
  } else {
    stopifnot(identical(new_y$labels, ltr_tree_tips_sorted$species_label),
              identical(axis_state(before[[tag]], "y.sec"), axis_state(after, "y.sec")))
  }
}
left_style <- ggplot2::calc_element("axis.text.y.left", isolated$b$theme)
right_style <- ggplot2::calc_element("axis.text.y.right", isolated$b$theme)
stopifnot(identical(left_style$family, "Arial"), identical(left_style$face, "italic"),
          inherits(right_style, "element_text"), identical(right_style$family, "Arial"))
isolated$a <- isolated$a + labs(title = "Standalone review only")
stopifnot(identical(source_state(panels$a), states$a))
for (invalid in list(panels[2:1], panels[1], unname(panels))) {
  error <- tryCatch({figure4_review_panels(invalid, ltr_tree_tips_sorted); NULL},
                    error = conditionMessage)
  stopifnot(!is.null(error))
}
duplicate_tips <- ltr_tree_tips_sorted
duplicate_tips$y[[2]] <- duplicate_tips$y[[1]]
nonfinite_tips <- ltr_tree_tips_sorted
nonfinite_tips$y[[1]] <- NA_real_
for (invalid in list(select(ltr_tree_tips_sorted, -species_label),
                     select(ltr_tree_tips_sorted, -y), duplicate_tips, nonfinite_tips)) {
  error <- tryCatch({figure4_review_panels(panels, invalid); NULL},
                    error = conditionMessage)
  stopifnot(!is.null(error))
}
'''
    )

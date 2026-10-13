from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_R = ROOT / "scripts/run_r.sh"


# Evaluate the production preparation and plots against the frozen release;
# no exports or unrelated figures are executed here.
FIGURE5_SETUP = r'''
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
for (name in c("read_supp", "species_label", "source_tree_path", "source_tree",
               "source_tree_tip_y", "canonical_phylogenetic_order",
               "phylogenetic_species_order", "phylogenetic_species_levels")) {
  eval(expressions[[assignment_index(name)]], envir = .GlobalEnv)
}
for (index in seq.int(assignment_index("path_traits"),
                     assignment_index("audit_tip_labels") - 1L)) {
  eval(expressions[[index]], envir = .GlobalEnv)
}
for (index in seq.int(assignment_index("morphology"),
                     assignment_index("Figure_5") - 1L)) {
  eval(expressions[[index]], envir = .GlobalEnv)
}
panels <- list(a = morph_tree_panel, b = morph_genome_panel,
               c = morph_nucleus_panel, d = morph_cell_panel)
tips <- arrange(morph_context, y)
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", FIGURE5_SETUP + code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_frozen_values_support_and_phylogenetic_rows_are_preserved() -> None:
    run_r_assertions(
        r'''
source_traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
source_morphology <- read_supp("Supplementary_Data_S19_cell_nucleus_species_estimates.csv")
source_objects <- read_supp("Supplementary_Data_S18_frozen_cell_nucleus_objects.csv")
source_quality <- read_supp("Supplementary_Data_S24_nuclear_iod_quality_balance.csv")
source_traits <- source_traits[match(morph_context$species, source_traits$species), ]
source_morphology <- source_morphology[
  match(morph_context$species, source_morphology$species), ]
for (column in c("genome_size_pg_fuscus_anchored", "genome_size_pg_ci_low",
                 "genome_size_pg_ci_high")) {
  stopifnot(identical(morph_context[[column]], source_traits[[column]]))
}
for (column in c("cell_area_um2", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
                 "nucleus_area_um2", "nucleus_area_um2_ci_low",
                 "nucleus_area_um2_ci_high", "n_size_specimens")) {
  stopifnot(identical(morph_context[[column]], source_morphology[[column]]))
}
stopifnot(nrow(morphology_objects) == 1152L,
          identical(morphology_objects$cell_area_um2, source_objects$cell_area_um2),
          identical(morphology_objects$nucleus_area_um2, source_objects$nuc_area_um2),
          identical(tips$tree_tip, phylogenetic_species_order(source_traits$species)),
          identical(tips$y, as.numeric(1:24)),
          nrow(path_topology_horizontal) == 46L,
          nrow(path_topology_vertical) == 23L,
          identical(sort(morph_context$tree_tip[
            morph_context$development_mode == "direct_development"]),
            c("aeneus", "wrighti")),
          sum(morph_context$development_mode == "aquatic_larva") == 22L,
          sum(morph_context$morphology_support == "One animal") == 7L,
          sum(morph_context$genome_support == "One image and animal") == 3L,
          sum(morph_context$genome_support == "Quality-balance caution") == 4L)
expected_cautions <- source_quality$species[
  !source_quality$passes_prespecified_balance_thresholds]
stopifnot(setequal(morph_context$species[
  morph_context$genome_support == "Quality-balance caution"], expected_cautions))
'''
    )


def test_fixed_axes_preserve_every_interval_and_support_symbol() -> None:
    run_r_assertions(
        r'''
specifications <- list(
  b = list(estimate = "genome_size_pg_fuscus_anchored", low = "genome_size_pg_ci_low",
           high = "genome_size_pg_ci_high", limits = c(5, 30),
           support = "genome_support",
           shapes = c(`Balanced multi-image support` = 21,
                      `One image and animal` = 24, `Quality-balance caution` = 4)),
  c = list(estimate = "nucleus_area_um2", low = "nucleus_area_um2_ci_low",
           high = "nucleus_area_um2_ci_high", limits = c(30, 60),
           support = "morphology_support",
           shapes = c(`Two or more animals` = 21, `One animal` = 24)),
  d = list(estimate = "cell_area_um2", low = "cell_area_um2_ci_low",
           high = "cell_area_um2_ci_high", limits = c(180, 340),
           support = "morphology_support",
           shapes = c(`Two or more animals` = 21, `One animal` = 24))
)
for (tag in names(specifications)) {
  specification <- specifications[[tag]]
  built <- ggplot_build(panels[[tag]])
  intervals <- built$data[[1]]
  points <- built$data[[2]]
  stopifnot(nrow(intervals) == 24L, nrow(points) == 24L,
            identical(intervals$x, morph_context[[specification$low]]),
            identical(intervals$xend, morph_context[[specification$high]]),
            identical(intervals$y, morph_context$y),
            identical(intervals$yend, morph_context$y),
            identical(points$x, morph_context[[specification$estimate]]),
            identical(points$y, morph_context$y),
            identical(points$shape, unname(specification$shapes[
              as.character(morph_context[[specification$support]])])),
            all(is.finite(intervals$x)), all(is.finite(intervals$xend)),
            all(is.finite(points$x)),
            identical(built$layout$panel_params[[1]]$x$limits, specification$limits),
            identical(built$layout$panel_params[[1]]$y$limits, c(.5, 24.5)))
}
tree <- ggplot_build(panels$a)
labels <- tree$data[[4]]
stopifnot(identical(labels$label, morph_context$display_species),
          identical(labels$y, morph_context$y), all(labels$family == "Arial"),
          all(labels$fontface == "italic"))
'''
    )


def test_isolated_panels_preserve_marks_and_restore_species_and_legends() -> None:
    run_r_assertions(
        r'''
source("scripts/publication/figure5_panels.R")
state <- function(plot) list(data = plot$data, labels = plot$labels, theme = plot$theme)
before <- lapply(panels, ggplot_build)
source_states <- lapply(panels, state)
isolated <- figure5_review_panels(panels, tips)
stopifnot(identical(names(isolated), letters[1:4]))
for (tag in names(panels)) {
  after <- ggplot_build(isolated[[tag]])
  old_axis <- before[[tag]]$layout$panel_params[[1]]
  new_axis <- after$layout$panel_params[[1]]
  stopifnot(identical(before[[tag]]$data, after$data),
            identical(old_axis$x$limits, new_axis$x$limits),
            identical(old_axis$x$get_breaks(), new_axis$x$get_breaks()),
            identical(old_axis$x$get_labels(), new_axis$x$get_labels()),
            identical(old_axis$y$limits, new_axis$y$limits),
            identical(state(panels[[tag]]), source_states[[tag]]),
            identical(isolated[[tag]]$labels$tag, tag))
  if (tag != "a") {
    stopifnot(identical(new_axis$y$get_breaks(), tips$y),
              identical(new_axis$y$get_labels(), tips$display_species),
              identical(isolated[[tag]]$theme$legend.position, "bottom"),
              inherits(isolated[[tag]]$scales$get_scales("shape")$guide, "GuideLegend"))
    label_style <- ggplot2::calc_element("axis.text.y.left", isolated[[tag]]$theme)
    stopifnot(identical(label_style$family, "Arial"),
              identical(label_style$face, "italic"))
  }
}
for (invalid in list(panels[4:1], panels[1], unname(panels))) {
  error <- tryCatch({figure5_review_panels(invalid, tips); NULL},
                    error = conditionMessage)
  stopifnot(!is.null(error))
}
duplicated_tips <- tips
duplicated_tips$display_species[[2]] <- duplicated_tips$display_species[[1]]
for (invalid in list(tips[24:1, ], tips[-1, ], select(tips, -display_species),
                     select(tips, -y), duplicated_tips)) {
  error <- tryCatch({figure5_review_panels(panels, invalid); NULL},
                    error = conditionMessage)
  stopifnot(!is.null(error))
}
'''
    )

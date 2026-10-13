from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


# Run the real Figure 8 preparation and plot expressions against the released
# data, without exporting figures or executing the other candidate blocks.
FIGURE8_SETUP = r'''
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})
source("scripts/publication/gbe_theme.R")
source("scripts/publication/auxiliary_panels.R")
ROOT <- normalizePath(".")
DATA <- file.path(ROOT, "Publication/datasets")
source_lines <- readLines("scripts/publication/build_gbe_figures.R")
expressions <- parse(text=source_lines)
evaluate_assignment <- function(name) {
  indices <- which(vapply(expressions, function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name(name))
  }, logical(1)))
  stopifnot(length(indices) == 1L)
  eval(expressions[[indices]], envir=.GlobalEnv)
}
for (name in c("read_supp", "species_label", "source_tree_path", "source_tree",
               "source_tree_tip_y", "canonical_phylogenetic_order",
               "phylogenetic_species_order", "phylogenetic_species_levels",
               "prune_source_tree", "te_strategy_levels", "te_strategy_colors",
               "te_strategy_shapes", "te_strategy_labels", "path_traits",
               "te_diversity", "te_annotations", "landscape")) {
  evaluate_assignment(name)
}
start <- match("## ---- wang-style-integrated-te-history-candidate", source_lines)
end <- match("## ---- wang-style-diversity-genome-size-candidate", source_lines)
stopifnot(!is.na(start), !is.na(end), end > start)
candidate_expressions <- parse(text=source_lines[seq.int(start+1L,end-1L)])
for (expression in candidate_expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name("record_auxiliary_figure"))) break
  eval(expression, envir=.GlobalEnv)
}
'''


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(ROOT / "scripts/run_r.sh"), "Rscript", "-e", FIGURE8_SETUP + code],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_all_species_keep_source_order_and_unavailable_genome_sizes() -> None:
    run_r_assertions(
        r'''
expected_species <- canonical_phylogenetic_order[
  canonical_phylogenetic_order %in% unique(landscape$species)
]
tip_rows <- integrated_tree_tips[order(integrated_tree_tips$y), ]
stopifnot(length(expected_species) == 34L,
          identical(integrated_species, expected_species),
          identical(tip_rows$species, expected_species),
          Ntip(integrated_tree) == 34L, Nnode(integrated_tree) == 33L,
          is.rooted(integrated_tree), is.binary(integrated_tree),
          !anyDuplicated(tip_rows$y), !anyDuplicated(tip_rows$species))
expected_genome <- path_traits %>% transmute(
  species=sub("^D\\. ", "", species), expected=genome_size_pg_fuscus_anchored
)
paired <- left_join(tip_rows, expected_genome, by="species")
stopifnot(nrow(paired) == 34L,
          sum(is.finite(paired$genome_size_pg_fuscus_anchored)) == 21L,
          sum(is.na(paired$genome_size_pg_fuscus_anchored)) == 13L,
          identical(is.na(paired$genome_size_pg_fuscus_anchored), is.na(paired$expected)),
          isTRUE(all.equal(paired$genome_size_pg_fuscus_anchored, paired$expected)))
expected_missing <- c("abditus", "adatsihi", "aureatus", "balsameus", "carolinensis",
  "catahoula", "cheaha", "conanti", "lycos", "organi", "pascagoula", "santeetlah", "valentinei")
stopifnot(setequal(paired$species[is.na(paired$expected)], expected_missing))
for (data in list(integrated_composition, integrated_fraction_data, integrated_profiles)) {
  rows <- distinct(data, species, y) %>% arrange(y)
  stopifnot(identical(rows$species, expected_species), identical(rows$y, tip_rows$y))
}
# The genome track must visibly retain the unavailable rows, not drop or impute them.
built <- ggplot_build(integrated_genome_panel)
text_layers <- which(vapply(integrated_genome_panel$layers,
  function(layer) inherits(layer$geom, "GeomText"), logical(1)))
labels <- bind_rows(built$data[text_layers])
stopifnot(nrow(labels) == 13L, all(labels$label == "n/a"),
          setequal(labels$y, tip_rows$y[is.na(tip_rows$genome_size_pg_fuscus_anchored)]))
point_index <- which(vapply(integrated_genome_panel$layers,
  function(layer) inherits(layer$geom,"GeomPoint"), logical(1)))
segment_index <- which(vapply(integrated_genome_panel$layers,
  function(layer) inherits(layer$geom,"GeomSegment"), logical(1)))
stopifnot(length(point_index) == 1L, length(segment_index) == 1L)
points <- built$data[[point_index]]
intervals <- built$data[[segment_index]]
expected <- tip_rows[is.finite(tip_rows$genome_size_pg_fuscus_anchored),]
points <- points[order(points$y),]; intervals <- intervals[order(intervals$y),]
stopifnot(nrow(points) == 21L, nrow(intervals) == 21L,
          identical(points$y,expected$y), identical(intervals$y,expected$y),
          identical(points$x,expected$genome_size_pg_fuscus_anchored),
          identical(intervals$x,expected$genome_size_pg_ci_low),
          identical(intervals$xend,expected$genome_size_pg_ci_high),
          setequal(c(points$y,labels$y),tip_rows$y))
'''
    )


def test_composition_and_profiles_preserve_the_full_hit_bp_denominator() -> None:
    run_r_assertions(
        r'''
raw <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
named_orders <- c("LTR", "LINE", "TIR", "DIRS", "SINE", "Maverick")
totals <- raw %>% group_by(species) %>% summarise(
  total=sum(hit_bp), denominator=first(species_total_hit_bp), .groups="drop"
)
stopifnot(identical(totals$total, totals$denominator))
expected_composition <- raw %>% mutate(
  repeat_category=if_else(order %in% named_orders, order, "Other")
) %>% group_by(species, repeat_category) %>% summarise(
  expected_bp=sum(hit_bp), denominator=first(species_total_hit_bp), .groups="drop"
)
paired <- inner_join(integrated_composition, expected_composition,
                     by=c("species", "repeat_category"))
stopifnot(nrow(integrated_composition) == 34L*7L, nrow(paired) == 34L*7L,
          !anyDuplicated(integrated_composition[c("species", "repeat_category")]),
          identical(paired$hit_bp, paired$expected_bp),
          max(abs(paired$repeat_fraction-paired$expected_bp/paired$denominator)) < 1e-12)
closure <- integrated_composition %>% group_by(species) %>% summarise(total=sum(repeat_fraction))
stopifnot(max(abs(closure$total-1)) < 1e-12)
fraction_expected <- expected_composition %>% filter(repeat_category %in% c("LTR", "LINE", "TIR"))
paired <- inner_join(integrated_fraction_data, fraction_expected,
                     by=c("species", "repeat_category"))
stopifnot(nrow(integrated_fraction_data) == 34L*3L, nrow(paired) == 34L*3L,
          max(abs(paired$fraction_percent-100*paired$expected_bp/paired$denominator)) < 1e-10)

key <- c("species", "repeat_category", "divergence_bin")
observed <- raw %>% filter(order %in% c("LTR", "LINE", "TIR")) %>%
  rename(repeat_category=order, divergence_bin=divergence_bin_start_pct)
paired <- inner_join(integrated_profiles, observed, by=key, suffix=c("_plot", "_source"))
stopifnot(nrow(integrated_profiles) == 34L*3L*50L,
          !anyDuplicated(integrated_profiles[key]), nrow(paired) == nrow(observed),
          identical(sort(as.integer(unique(integrated_profiles$divergence_bin))),0:49),
          max(abs(paired$percent_species_hit_bp_plot-paired$percent_species_hit_bp_source)) < 1e-12)
absent <- anti_join(integrated_profiles, observed, by=key)
stopifnot(nrow(absent) > 0L, all(absent$percent_species_hit_bp == 0),
          all(is.finite(integrated_profiles$percent_species_hit_bp)),
          max(abs(integrated_profiles$profile_bottom-(integrated_profiles$y-.34))) < 1e-12,
          max(abs((integrated_profiles$profile_top-integrated_profiles$profile_bottom)/.68*3-
                    integrated_profiles$percent_species_hit_bp)) < 1e-12)
'''
    )


def test_shared_axes_show_every_value_and_isolated_panels_retain_all_rows() -> None:
    run_r_assertions(
        r'''
stopifnot(identical(integrated_fraction_max,60), identical(integrated_profile_max,3),
          all(integrated_fraction_data$fraction_percent >= 0 &
              integrated_fraction_data$fraction_percent <= 60),
          all(integrated_profiles$percent_species_hit_bp >= 0 &
              integrated_profiles$percent_species_hit_bp <= 3))
for (spec in list(list(plot=integrated_fraction_panel, limits=c(0,60)),
                 list(plot=integrated_profile_panel, limits=c(0,50)))) {
  built <- ggplot_build(spec$plot)
  stopifnot(nrow(built$layout$layout) == 3L,
            identical(as.character(built$layout$layout$repeat_category), c("LTR","LINE","TIR")),
            identical(spec$plot$scales$get_scales("x")$limits, spec$limits))
  for (axis in built$layout$panel_params) {
    stopifnot(identical(axis$x$limits,spec$limits), identical(axis$y$limits,c(.5,34.5)))
  }
}
point_index <- which(vapply(integrated_fraction_panel$layers,
  function(layer) inherits(layer$geom,"GeomPoint"), logical(1)))
stopifnot(length(point_index) == 1L)
points <- ggplot_build(integrated_fraction_panel)$data[[point_index]]
stopifnot(nrow(points) == 102L, all(is.finite(points$x)), all(is.finite(points$y)),
          isTRUE(all.equal(sort(points$x),sort(integrated_fraction_data$fraction_percent))))
spec <- auxiliary_review_spec("Auxiliary_Wang_style_integrated_te_history_candidate")
originals <- list(integrated_tree_panel, integrated_genome_panel, integrated_composition_panel,
                  integrated_fraction_panel, integrated_profile_panel)
stopifnot(length(spec$plots) == 5L)
for (index in 2:5) {
  stopifnot(identical(ggplot_build(spec$plots[[index]])$data,
                      ggplot_build(originals[[index]])$data),
            setequal(spec$plots[[index]]$scales$get_scales("y")$labels,
                     integrated_tree_tips$display_species),
            length(spec$plots[[index]]$scales$get_scales("y")$labels) == 34L)
}
'''
    )

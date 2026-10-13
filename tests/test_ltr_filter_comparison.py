from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Publication/datasets"
ELEMENTS = DATA / "Supplementary_Data_S14_ltr_terminal_internal_elements.csv"
SUMMARIES = DATA / "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv"
RATIO = "ratio_terminal_internal_all_positions"
PRIMARY = "primary_iqr_filtered_coverage_ge_80pct"



def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, float_precision="round_trip", low_memory=False)


def run_r(code: str, *args: str) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [str(ROOT / "scripts/run_r.sh"), "Rscript", "-e", code, *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=120,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    return completed


@pytest.fixture(scope="module")
def comparison(tmp_path_factory: pytest.TempPathFactory) -> dict[str, pd.DataFrame]:
    output = tmp_path_factory.mktemp("ltr_filter_comparison")
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (ELEMENTS, SUMMARIES)}
    run_r(
        r"""
source("analyses/01_transposable_elements/ltr_filter_comparison.R")
prepared <- prepare_ltr_filter_comparison("Publication/datasets")
output <- commandArgs(trailingOnly = TRUE)[[1]]
for (name in c("elements", "summaries", "decisions", "counts", "paired_species_summary")) {
  readr::write_csv(prepared[[name]], file.path(output, paste0(name, ".csv")))
}
""",
        str(output),
    )
    assert before == {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before}
    return {path.stem: read_csv(path) for path in output.glob("*.csv")}


def independent_four_domain_decisions(elements: pd.DataFrame) -> pd.DataFrame:
    """Rebuild cohort selection with pandas, independently of the R implementation."""
    domain_count = elements["Domains"].str.split().str.len()
    eligible = (
        domain_count.ge(4)
        & elements["Order"].eq("LTR")
        & elements["Superfamily"].eq("Gypsy")
        & elements["element length"].ge(3000)
        & elements["tesorter_annotation_found"]
    )
    quartiles = elements.loc[eligible].groupby("species")[RATIO].quantile([0.25, 0.75]).unstack()
    q1 = elements["species"].map(quartiles[0.25]).where(eligible)
    q3 = elements["species"].map(quartiles[0.75]).where(eligible)
    lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
    iqr_keep = eligible & elements[RATIO].between(lower, upper, inclusive="both")
    coverage = (
        elements[
            [
                "left_ltr_positive_coverage_fraction",
                "right_ltr_positive_coverage_fraction",
                "internal_positive_coverage_fraction",
            ]
        ]
        .ge(0.8)
        .all(axis=1)
    )
    return pd.DataFrame(
        {
            "four_domain_count": domain_count,
            "four_domain_eligible": eligible,
            "four_domain_q1": q1,
            "four_domain_q3": q3,
            "four_domain_iqr": q3 - q1,
            "four_domain_lower_bound": lower,
            "four_domain_upper_bound": upper,
            "four_domain_iqr_keep": iqr_keep,
            "four_domain_coverage_keep": eligible & coverage,
            "four_domain_keep": iqr_keep & coverage,
        }
    )


def test_primary_branch_preserves_all_released_values_and_intervals(comparison: dict) -> None:
    released_elements = read_csv(ELEMENTS)
    released_elements = released_elements.loc[released_elements["primary_ltr_element"]]
    released_summaries = read_csv(SUMMARIES)
    released_summaries = released_summaries.loc[released_summaries["analysis_branch"].eq(PRIMARY)]
    elements = comparison["elements"].query("branch_id == 'five_domains'")
    summaries = comparison["summaries"].query("branch_id == 'five_domains'")
    assert len(elements) == 380
    assert elements["species"].nunique() == len(summaries) == 30
    assert set(elements["element_id"]) == set(released_elements["element_id"])
    expected = released_elements.set_index("element_id").sort_index()
    actual = elements.set_index("element_id").sort_index()
    for column in [RATIO, "log2_ratio_terminal_internal_all_positions"]:
        np.testing.assert_array_equal(actual[column], expected[column])
    expected_summary = released_summaries.set_index("species").sort_index()
    actual_summary = summaries.set_index("species").sort_index()
    for column in [
        "n_elements",
        "ratio_median",
        "median_bootstrap_ci_low",
        "median_bootstrap_ci_high",
    ]:
        np.testing.assert_array_equal(actual_summary[column], expected_summary[column])


def test_four_domain_eligibility_and_fences_match_independent_reconstruction(
    comparison: dict,
) -> None:
    elements = read_csv(ELEMENTS)
    expected = independent_four_domain_decisions(elements)
    actual = comparison["decisions"].set_index("element_id").loc[elements["element_id"]]
    assert len(actual) == 1086
    for column in expected:
        if expected[column].dtype.kind == "f":
            np.testing.assert_allclose(actual[column], expected[column], rtol=1e-13, atol=1e-14)
        else:
            np.testing.assert_array_equal(actual[column], expected[column])
    assert expected["four_domain_count"].value_counts().to_dict() == {4: 677, 5: 409}
    assert expected["four_domain_eligible"].sum() == 1083
    assert expected["four_domain_coverage_keep"].sum() == 1080
    assert expected["four_domain_keep"].sum() == 989
    # The legacy pipe count would classify every four-domain row as five domains.
    assert not np.array_equal(
        actual["four_domain_count"], elements["legacy_pipe_split_domain_count"]
    )
    assert not actual.loc["JAUBNV010740999.1_De_636_6770", "four_domain_eligible"]
    # Coverage is intersected AFTER fences are derived, as in the primary release.
    all_eligible = elements.loc[expected["four_domain_eligible"]]
    coverage_only = all_eligible.loc[expected["four_domain_coverage_keep"]]
    pre = all_eligible.groupby("species")[RATIO].quantile(0.75)
    post = coverage_only.groupby("species")[RATIO].quantile(0.75)
    assert not np.array_equal(pre, post)


def test_fresh_fences_determine_exact_membership_instead_of_forcing_a_superset(
    comparison: dict,
) -> None:
    released = read_csv(ELEMENTS)
    expected = independent_four_domain_decisions(released)
    elements = comparison["elements"]
    assert set(elements["branch_id"]) == {"five_domains", "four_domains"}
    assert elements["branch_retained"].all()
    assert not elements.duplicated(["branch_id", "element_id"]).any()
    five = set(elements.loc[elements["branch_id"].eq("five_domains"), "element_id"])
    four = set(elements.loc[elements["branch_id"].eq("four_domains"), "element_id"])
    assert four == set(released.loc[expected["four_domain_keep"], "element_id"])
    assert five - four == {"JAUBNT010480148.1_De_785_6985", "JAUBNT010598382.1_De_653_6905"}
    assert len(four - five) == 611
    assert len(four & five) == 378
    observed = elements.query("branch_id == 'four_domains'").set_index("element_id")
    np.testing.assert_array_equal(
        observed[RATIO], released.set_index("element_id").loc[observed.index, RATIO]
    )
    assert observed[RATIO].min() == pytest.approx(0.0297603723556847)
    assert observed[RATIO].max() == pytest.approx(2.203591072231693)


def test_species_summaries_use_raw_scale_medians_and_honest_paired_counts(comparison: dict) -> None:
    elements = comparison["elements"]
    summaries = comparison["summaries"].set_index(["branch_id", "species"]).sort_index()
    grouped = elements.groupby(["branch_id", "species"])[RATIO]
    expected_median = grouped.median().sort_index()
    expected_count = grouped.size().sort_index()
    np.testing.assert_allclose(summaries["ratio_median"], expected_median, rtol=1e-14, atol=1e-15)
    np.testing.assert_array_equal(summaries["n_elements"], expected_count)
    median_of_logs = (
        elements.assign(log_ratio=np.log2(elements[RATIO]))
        .groupby(["branch_id", "species"])["log_ratio"]
        .median()
        .sort_index()
    )
    assert np.any(np.abs(np.log2(summaries["ratio_median"]) - median_of_logs) > 1e-6)
    assert summaries["ratio_median"].lt(1).all()
    assert (summaries["median_bootstrap_ci_low"] <= summaries["ratio_median"]).all()
    assert (summaries["median_bootstrap_ci_high"] >= summaries["ratio_median"]).all()
    pairs = comparison["paired_species_summary"].set_index("species").sort_index()
    assert len(pairs) == 30
    for branch in ["five_domains", "four_domains"]:
        for value in ["n_elements", "ratio_median"]:
            np.testing.assert_array_equal(pairs[f"{value}_{branch}"], summaries.loc[branch, value])
    np.testing.assert_array_equal(
        pairs["n_elements_change"],
        pairs["n_elements_four_domains"] - pairs["n_elements_five_domains"],
    )
    np.testing.assert_allclose(
        pairs["ratio_median_change"],
        pairs["ratio_median_four_domains"] - pairs["ratio_median_five_domains"],
        atol=1e-15,
    )
    np.testing.assert_allclose(
        pairs["ratio_median_change_pct"],
        100 * (pairs["ratio_median_four_domains"] / pairs["ratio_median_five_domains"] - 1),
        atol=1e-12,
    )
    assert pairs.loc[
        ["pascagoula", "apalachicolae", "marmoratus"], "n_elements_five_domains"
    ].tolist() == [3, 4, 4]
    assert pairs.loc[
        ["pascagoula", "apalachicolae", "marmoratus"], "n_elements_four_domains"
    ].tolist() == [16, 12, 11]
    counts = comparison["counts"].set_index("branch_id")
    for branch, numbers in {
        "five_domains": [1086, 408, 407, 381, 380, 30],
        "four_domains": [1086, 1083, 1080, 992, 989, 30],
    }.items():
        assert (
            counts.loc[
                branch,
                [
                    "historical_usable",
                    "eligible",
                    "coverage_qualified",
                    "iqr_retained",
                    "retained",
                    "n_species",
                ],
            ].tolist()
            == numbers
        )


def test_bootstrap_intervals_are_deterministic_without_replacing_primary_intervals() -> None:
    run_r(r"""
source("analyses/01_transposable_elements/ltr_filter_comparison.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(71)
before_kind <- RNGkind()
before_seed <- .Random.seed
first <- prepare_ltr_filter_comparison("Publication/datasets")
stopifnot(identical(RNGkind(), before_kind), identical(.Random.seed, before_seed))
set.seed(912)
second <- prepare_ltr_filter_comparison("Publication/datasets")
stopifnot(identical(first, second))
rm(".Random.seed", envir = .GlobalEnv)
third <- prepare_ltr_filter_comparison("Publication/datasets")
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
          identical(RNGkind(), before_kind), identical(first, third))
""")


def test_broader_bootstrap_follows_species_identity_when_input_rows_are_reordered(
    comparison: dict, tmp_path: Path
) -> None:
    shuffled = read_csv(ELEMENTS).sample(frac=1, random_state=20260924)
    shuffled.to_csv(tmp_path / ELEMENTS.name, index=False)
    shutil.copyfile(SUMMARIES, tmp_path / SUMMARIES.name)
    run_r(
        r"""
source("analyses/01_transposable_elements/ltr_filter_comparison.R")
scratch <- commandArgs(trailingOnly = TRUE)[[1]]
prepared <- prepare_ltr_filter_comparison(scratch)
readr::write_csv(prepared$summaries, file.path(scratch, "shuffled_summaries.csv"))
""",
        str(tmp_path),
    )
    actual = read_csv(tmp_path / "shuffled_summaries.csv")
    expected = comparison["summaries"]
    pd.testing.assert_frame_equal(actual, expected, check_exact=True)
    broad = actual.loc[actual["branch_id"].eq("four_domains")]
    assert broad["bootstrap_replicates"].eq(2000).all()
    assert broad["bootstrap_seed"].notna().all()
    assert broad["bootstrap_seed"].nunique() == 30
    assert broad["bootstrap_method"].eq("R_Mersenne_Twister_Rejection_percentile_type7").all()
    assert (
        actual.loc[actual["branch_id"].eq("five_domains"), "bootstrap_method"]
        .eq("released_S15_NumPy_PCG64_percentile_linear")
        .all()
    )


PLOT_SETUP = r'''
suppressPackageStartupMessages(library(ggplot2))
source("scripts/publication/gbe_theme.R")
source("analyses/01_transposable_elements/ltr_filter_comparison.R")
comparison <- prepare_ltr_filter_comparison("Publication/datasets")
assert_close <- function(actual, expected, tolerance = 1e-12) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected),
                            tolerance = tolerance, check.attributes = FALSE)))
}
layer_index <- function(plot, class) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, class), logical(1)))
}
'''


def test_matched_violins_preserve_every_point_and_share_axes_species_order_and_smoothing() -> None:
    run_r(PLOT_SETUP + r'''
source("scripts/publication/ltr_filter_violin_panels.R")
views <- make_ltr_filter_violin_panels(comparison, "Publication/datasets")
tree <- ape::ladderize(ape::read.tree("Publication/datasets/trees/desmo900dated_source46.tre"),
                       right = FALSE)
canonical <- tree$tip.label[order(ape::node.height(tree)[seq_len(ape::Ntip(tree))])]
canonical <- canonical[canonical %in% comparison$summaries$species]
stopifnot(identical(views$rows$species, canonical),
          identical(views$rows$y, seq_len(30)),
          sum(!is.na(views$rows$genome_size_gbp)) == 19L,
          sum(is.na(views$rows$genome_size_gbp)) == 11L)
assert_close(views$bandwidth, stats::bw.nrd0(log2(comparison$elements[
  comparison$elements$branch_id == "five_domains", "ratio_terminal_internal_all_positions"])))
for (branch in c("five_domains", "four_domains")) {
  plot <- views$plots[[branch]]
  built <- ggplot_build(plot)
  compact <- ggplot_build(views$collage_panels[[branch]])
  stopifnot(identical(built$data, compact$data))
  data <- views$elements[views$elements$branch_id == branch, ]
  summary <- views$summaries[views$summaries$branch_id == branch, ]
  points <- built$data[[layer_index(plot, "GeomPoint")]]
  medians <- built$data[[5L]]
  intervals <- built$data[[4L]]
  stopifnot(nrow(points) == if (branch == "five_domains") 380L else 989L,
            nrow(medians) == 30L, nrow(intervals) == 30L,
            all(abs(points$y - data$y) <= .15))
  assert_close(points$x, log2(data$ratio_terminal_internal_all_positions))
  assert_close(medians$x, log2(summary$ratio_median))
  assert_close(medians$xend, medians$x)
  assert_close(medians$y, summary$y - .25)
  assert_close(medians$yend, summary$y + .25)
  assert_close(intervals$x, log2(summary$median_bootstrap_ci_low))
  assert_close(intervals$xend, log2(summary$median_bootstrap_ci_high))
  violin <- plot$layers[[layer_index(plot, "GeomViolin")]]
  stopifnot(identical(violin$stat_params$trim, TRUE),
            identical(violin$stat_params$scale, "width"),
            identical(violin$stat_params$orientation, "y"))
  assert_close(violin$stat_params$bw, views$bandwidth)
  density <- built$data[[layer_index(plot, "GeomViolin")]]
  stopifnot(length(unique(density$group)) == 30L)
  for (group in split(density, density$group)) {
    observed <- data$estimate[data$y == unique(group$y)]
    assert_close(range(group$x), range(observed))
    stopifnot(all(group$n == length(observed)), all(is.finite(group$density)))
  }
  axis <- built$layout$panel_params[[1]]
  compact_axis <- compact$layout$panel_params[[1]]
  assert_close(axis$x$limits, log2(c(.025, 2.5)))
  assert_close(axis$x$get_breaks(), log2(c(.025, .1, .5, 1, 2.5)))
  assert_close(axis$y$limits, c(.5, 30.5))
  stopifnot(identical(axis$y$get_labels(), paste("D.", canonical)),
            identical(axis$x$limits, compact_axis$x$limits),
            identical(axis$y$limits, compact_axis$y$limits))
  count <- summary$n_elements[match(views$rows$species, summary$species)]
  genome <- ifelse(is.na(views$rows$genome_size_gbp), "n/a",
                   sprintf("%.1f", views$rows$genome_size_gbp))
  stopifnot(identical(as.character(axis$y.sec$get_labels()), paste0(genome, " / ", count)),
            identical(as.character(compact_axis$y.sec$get_labels()), as.character(count)),
            identical(ggplot2::calc_element("axis.text.y.left", plot$theme)$family, "Arial"))
  for (layer in built$data) for (position in intersect(c("x", "xend", "y", "yend"), names(layer))) {
    stopifnot(all(is.finite(layer[[position]])))
  }
}
''')


def test_both_pgls_panels_use_same_nineteen_species_and_independent_gls_intervals() -> None:
    run_r(PLOT_SETUP + r'''
source("analyses/01_transposable_elements/ltr_genome_pgls.R")
source("scripts/publication/ltr_filter_pgls_panels.R")
views <- make_ltr_filter_pgls_panels(comparison, "Publication/datasets")
reference <- prepare_ltr_genome_data("Publication/datasets")
stopifnot(nrow(views$model_summary) == 2L, nrow(views$species_data) == 38L,
          nrow(views$predictions) == 400L, nrow(views$excluded_species) == 11L,
          all(views$model_summary$n_species == 19L))
assert_close(ape::vcv.phylo(views$tree), ape::vcv.phylo(reference$tree))
previous_axes <- NULL
previous_labels <- NULL
for (branch in c("five_domains", "four_domains")) {
  data <- views$species_data[views$species_data$branch_id == branch, ]
  fitted <- views$model_summary[views$model_summary$branch_id == branch, ]
  prediction <- views$predictions[views$predictions$branch_id == branch, ]
  source_summary <- comparison$summaries[comparison$summaries$branch_id == branch, ]
  selected <- match(data$species, source_summary$species)
  stopifnot(identical(data$species, reference$data$species),
            !anyDuplicated(data$species), fitted$lambda >= 0, fitted$lambda <= 1)
  for (column in c("genome_size_gbp", "genome_size_ci_low_gbp", "genome_size_ci_high_gbp")) {
    assert_close(data[[column]], reference$data[[column]])
  }
  for (column in c("n_elements", "ratio_median", "median_bootstrap_ci_low", "median_bootstrap_ci_high")) {
    assert_close(data[[column]], source_summary[[column]][selected])
  }
  assert_close(data$log10_ratio, log10(data$ratio_median))
  # Reconstruct generalized least squares with the actual tree covariance.
  # Neither phylolm coefficients nor the production prediction helper are used.
  covariance <- ape::vcv.phylo(views$tree)[data$species, data$species]
  covariance[row(covariance) != col(covariance)] <-
    covariance[row(covariance) != col(covariance)] * fitted$lambda
  design <- cbind(1, log10(data$genome_size_gbp))
  response <- log10(data$ratio_median)
  information_inverse <- solve(crossprod(design, solve(covariance, design)))
  coefficient <- information_inverse %*% crossprod(design, solve(covariance, response))
  residual <- response - drop(design %*% coefficient)
  coefficient_covariance <- information_inverse *
    drop(crossprod(residual, solve(covariance, residual))) / 17
  standard_error <- sqrt(coefficient_covariance[2, 2])
  critical <- qt(.975, df = 17)
  assert_close(fitted$intercept, coefficient[1], 1e-7)
  assert_close(fitted$slope, coefficient[2], 1e-7)
  assert_close(fitted$standard_error, standard_error, 1e-7)
  assert_close(fitted$p_value, 2 * pt(-abs(coefficient[2] / standard_error), 17), 1e-7)
  assert_close(c(fitted$ci_low, fitted$ci_high), coefficient[2] + c(-1, 1) * critical * standard_error, 1e-7)
  expected_x <- seq(min(design[, 2]), max(design[, 2]), length.out = 200)
  prediction_design <- cbind(1, expected_x)
  mean_log <- drop(prediction_design %*% coefficient)
  mean_se <- sqrt(diag(prediction_design %*% coefficient_covariance %*% t(prediction_design)))
  assert_close(prediction$genome_size_gbp, 10^expected_x)
  assert_close(prediction$ratio_fit, 10^mean_log, 1e-7)
  assert_close(prediction$ratio_ci_low, 10^(mean_log - critical * mean_se), 1e-7)
  assert_close(prediction$ratio_ci_high, 10^(mean_log + critical * mean_se), 1e-7)
  plot <- views$plots[[branch]]
  built <- ggplot_build(plot)
  point <- built$data[[layer_index(plot, "GeomPoint")]]
  ribbon <- built$data[[layer_index(plot, "GeomRibbon")]]
  line <- built$data[[layer_index(plot, "GeomLine")]]
  stopifnot(nrow(point) == 19L, nrow(ribbon) == 200L, nrow(line) == 200L)
  assert_close(point$x, log10(data$genome_size_gbp))
  assert_close(point$y, log10(data$ratio_median))
  assert_close(line$y, log10(prediction$ratio_fit))
  assert_close(ribbon$ymin, log10(prediction$ratio_ci_low))
  assert_close(ribbon$ymax, log10(prediction$ratio_ci_high))
  label_layers <- built$data[layer_index(plot, "GeomText")]
  labels <- Filter(function(layer) nrow(layer) == 19L, label_layers)[[1]]
  stopifnot(identical(labels$label, paste("D.", data$species)),
            all(labels$family == "Arial"))
  if (!is.null(previous_labels)) {
    stopifnot(identical(labels[, c("label", "x", "y", "hjust")], previous_labels))
  }
  previous_labels <- labels[, c("label", "x", "y", "hjust")]
  # Names may move for legibility, but every tether must start at its own species.
  segment_layers <- built$data[layer_index(plot, "GeomSegment")]
  tethers <- Filter(function(layer) all(layer$colour == "#999999"), segment_layers)[[1]]
  assert_close(tethers$x, point$x)
  assert_close(tethers$y, point$y)
  axes <- lapply(built$layout$panel_params[[1]][c("x", "y")], function(axis) {
    list(limits = axis$limits, breaks = axis$get_breaks(), labels = axis$get_labels())
  })
  if (!is.null(previous_axes)) stopifnot(identical(axes, previous_axes))
  previous_axes <- axes
  for (layer in built$data) for (position in intersect(c("x", "xend", "y", "yend", "ymin", "ymax"), names(layer))) {
    stopifnot(all(is.finite(layer[[position]])))
  }
}
''')

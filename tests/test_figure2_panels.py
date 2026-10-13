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


def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", code],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_isolated_bonsai_panels_preserve_source_plots_values_and_scales() -> None:
    run_r_assertions(
        r"""
source("scripts/publication/gbe_theme.R")
source("scripts/publication/figure2_panels.R")

groups <- c("aquatic_eggs", "terrestrial_eggs_aquatic_larvae", "direct_development")
palette <- setNames(c("#0072B2", "#D55E00", "#7040A0"), groups)
observations <- data.frame(
  species = c("D. aquatic", "D. terrestrial", "D. direct"),
  y = c(3, 2, 1), strategy = factor(groups, levels = groups),
  x = c(.275, .01, .35), lower = c(.218, .003, .287),
  upper = c(.341, .036, .418)
)
base <- ggplot(observations, aes(x, y, fill = strategy, shape = strategy)) +
  geom_point(size = 2) +
  scale_fill_manual(values = palette, drop = FALSE) +
  scale_shape_manual(values = setNames(c(21, 22, 24), groups), drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1),
                     labels = scales::label_percent()) +
  scale_y_continuous(limits = c(.5, 3.5), breaks = 1:3) +
  theme_gbe() + theme(legend.position = "bottom")
panels <- list(
  a = base + geom_text(aes(label = species), family = "Arial", fontface = "italic"),
  b = base + geom_segment(aes(x = lower, xend = upper, yend = y)),
  c = suppressMessages(base +
    scale_x_continuous(limits = c(0, 11), breaks = c(0, 2, 4, 6, 8, 10)))
)

source_state <- function(plot) list(
  data = plot$data, labels = plot$labels, theme = plot$theme,
  layers = length(plot$layers),
  x_limits = plot$scales$get_scales("x")$limits,
  y_limits = plot$scales$get_scales("y")$limits
)
axis_state <- function(built, axis) {
  position <- built$layout$panel_params[[1]][[axis]]
  list(limits = position$limits, breaks = position$get_breaks(),
       labels = position$get_labels(), range = position$continuous_range)
}
before_states <- lapply(panels, source_state)
before_builds <- lapply(panels, ggplot_build)
isolated <- figure2_review_panels(panels)
stopifnot(identical(names(isolated), letters[1:3]))
for (tag in names(panels)) {
  before <- before_builds[[tag]]
  after <- ggplot_build(isolated[[tag]])
  stopifnot(identical(before$data, after$data),
            identical(axis_state(before, "x"), axis_state(after, "x")),
            identical(axis_state(before, "y"), axis_state(after, "y")))
  for (aesthetic in c("fill", "shape")) {
    stopifnot(identical(
      before$plot$scales$get_scales(aesthetic)$map(groups),
      after$plot$scales$get_scales(aesthetic)$map(groups)
    ))
  }
  stopifnot(identical(isolated[[tag]]$labels$tag, tag),
            identical(isolated[[tag]]$data$species, observations$species),
            identical(isolated[[tag]]$theme$legend.position, "bottom"))
  stopifnot(identical(source_state(panels[[tag]]), before_states[[tag]]))
}
# A second edit to an isolated plot must also leave the retained collage untouched.
isolated$a <- isolated$a + labs(title = "Isolated review only")
stopifnot(identical(source_state(panels$a), before_states$a))
for (invalid in list(panels[c("b", "a", "c")], panels[c("a", "b")], unname(panels))) {
  message <- tryCatch({figure2_review_panels(invalid); NULL}, error = conditionMessage)
  stopifnot(identical(message, "Figure 2 requires panels a-c"))
}
"""
    )


def test_released_recovery_values_and_labels_fit_full_percentage_axis() -> None:
    run_r_assertions(
        r"""
suppressPackageStartupMessages(library(dplyr))
source("scripts/publication/gbe_theme.R")

te_strategy_levels <- c("aquatic_eggs", "terrestrial_eggs_aquatic_larvae",
                        "direct_development")
te_strategy_colors <- setNames(c("#0072B2", "#D55E00", "#7040A0"), te_strategy_levels)
te_strategy_shapes <- setNames(c(21, 22, 24), te_strategy_levels)
te_validation <- read.csv(
  "Publication/datasets/Supplementary_Data_S43_te_bonsai_strategy_validation.csv"
)
validation_before <- te_validation
expressions <- parse(file = "scripts/publication/build_gbe_figures.R")
evaluate_assignment <- function(name) {
  assignment <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name(name))
  }, expressions)
  stopifnot(length(assignment) == 1L)
  eval(assignment[[1]], envir = .GlobalEnv)
}
evaluate_assignment("z_value")
evaluate_assignment("te_recovery")
stopifnot(identical(te_validation, validation_before),
          nrow(te_recovery) == 3L,
          all(te_recovery$feature_bootstrap_replicates == 200L))
expected_counts <- c(aquatic_eggs = 55L, terrestrial_eggs_aquatic_larvae = 2L,
                     direct_development = 70L)
for (i in seq_len(nrow(te_recovery))) {
  count <- unname(expected_counts[as.character(te_recovery$strategy[[i]])])
  stopifnot(te_recovery$feature_bootstrap_edge_delimited_count[[i]] == count,
            abs(te_recovery$feature_bootstrap_edge_delimited_fraction[[i]] - count / 200) < 1e-14,
            te_recovery$count_label[[i]] == sprintf("%d/200", count))
}
stopifnot(all(te_recovery$ci_low >= 0), all(te_recovery$ci_high <= 1),
          all(te_recovery$ci_low < te_recovery$feature_bootstrap_edge_delimited_fraction),
          all(te_recovery$ci_high > te_recovery$feature_bootstrap_edge_delimited_fraction),
          max(te_recovery$ci_high) > .41, max(te_recovery$ci_high) < .43)
recovery_before <- te_recovery
evaluate_assignment("te_c")
built <- ggplot_build(te_c)
stopifnot(identical(te_recovery, recovery_before),
          identical(te_c$scales$get_scales("x")$limits, c(0, 1)),
          identical(built$layout$panel_params[[1]]$x$get_labels(), c("0%", "50%", "100%")))
stopifnot(length(built$data) == 3L)
intervals <- built$data[[1]]
points <- built$data[[2]]
counts <- built$data[[3]]
for (layer in built$data) {
  stopifnot(nrow(layer) == 3L)
  for (position in intersect(c("x", "xend", "y", "yend"), names(layer))) {
    stopifnot(all(is.finite(layer[[position]])))
  }
}
stopifnot(identical(intervals$x, te_recovery$ci_low),
          identical(intervals$xend, te_recovery$ci_high),
          identical(points$x, te_recovery$feature_bootstrap_edge_delimited_fraction),
          identical(counts$label, te_recovery$count_label),
          all(counts$x > max(intervals$xend) + .2), all(counts$x <= 1),
          all(counts$family == "Arial"))
"""
    )

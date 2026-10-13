from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_R = ROOT / "scripts/run_r.sh"



def test_isolated_panels_preserve_data_scales_and_species_identity() -> None:
    result = subprocess.run(
        [
            str(RUN_R),
            "Rscript",
            "-e",
            r'''
source("scripts/publication/gbe_theme.R")
source("scripts/publication/figure1_panels.R")

tips <- data.frame(
  y = 1:21,
  display_species = sprintf("D. species_%02d", 1:21)
)
observations <- transform(
  tips, x = seq(7, 27, length.out = 21),
  group = factor(rep(c("LTR", "LINE", "TIR"), 7),
                 levels = c("LTR", "LINE", "TIR"))
)
base <- ggplot(observations, aes(x, y, fill = group)) +
  geom_point(shape = 21, size = 2) +
  scale_fill_manual(values = GBE_TE_ORDER_COLORS, drop = FALSE) +
  scale_x_continuous(limits = c(5, 30), breaks = c(10, 20, 30)) +
  scale_y_continuous(limits = c(.5, 21.5), breaks = NULL,
                     expand = c(0, 0)) +
  theme_gbe() + theme(legend.position = "right")
panels <- setNames(rep(list(base), 5), letters[1:5])
panels$a <- base + geom_text(aes(label = display_species),
                            family = "Arial", fontface = "italic")
panels$b <- base + geom_segment(aes(x = x - 1, xend = x + 1, yend = y))
panels$c <- suppressMessages(base +
  scale_x_continuous(limits = c(0, 30), breaks = c(0, 15, 30),
                     labels = c("0%", "50%", "100%")))
panels$d <- base + geom_rect(aes(xmin = x - .3, xmax = x + .3,
                                ymin = y - .38, ymax = y + .15))
panels$e <- suppressMessages(base +
  scale_x_continuous(limits = c(5, 30), breaks = c(10, 20, 30),
    sec.axis = dup_axis(breaks = c(10, 20, 30),
                       labels = c("LTR", "LINE", "TIR"))))

before <- lapply(panels, ggplot_build)
x_state <- function(built) {
  params <- built$layout$panel_params[[1]]
  list(
    range = params$x.range,
    breaks = params$x$get_breaks(), labels = params$x$get_labels(),
    secondary_breaks = params$x.sec$get_breaks(),
    secondary_labels = params$x.sec$get_labels()
  )
}
source_state <- lapply(panels, function(p) list(
  labels = p$labels,
  theme = p$theme,
  y_breaks = p$scales$get_scales("y")$breaks,
  y_labels = p$scales$get_scales("y")$labels,
  y_limits = p$scales$get_scales("y")$limits,
  x_limits = p$scales$get_scales("x")$limits
))

for (seed in c(1, 53, 204)) {
  set.seed(seed)
  shuffled_tips <- tips[sample(nrow(tips)), , drop = FALSE]
  isolated <- suppressMessages(figure1_review_panels(panels, shuffled_tips))
  stopifnot(identical(names(isolated), letters[1:5]))
  for (name in names(panels)) {
    p <- isolated[[name]]
    after <- ggplot_build(p)
    # No data transformation, renormalization, palette, or axis-range changes.
    stopifnot(identical(before[[name]]$data, after$data))
    stopifnot(identical(x_state(before[[name]]), x_state(after)))
    stopifnot(identical(
      before[[name]]$plot$scales$get_scales("fill")$map(c("LTR", "LINE", "TIR")),
      after$plot$scales$get_scales("fill")$map(c("LTR", "LINE", "TIR"))
    ))
    stopifnot(identical(p$labels$tag, name))
    stopifnot(identical(p$theme$legend.position, "bottom"))
    stopifnot(identical(p$scales$get_scales("y")$limits, c(.5, 21.5)))
    if (name == "a") {
      stopifnot(is.null(p$scales$get_scales("y")$breaks))
    } else {
      axis <- after$layout$panel_params[[1]]$y
      species_by_y <- setNames(axis$get_labels(), axis$get_breaks())
      stopifnot(identical(unname(species_by_y[as.character(tips$y)]),
                          tips$display_species))
      text_style <- ggplot2::calc_element("axis.text.y", p$theme)
      stopifnot(identical(text_style$face, "italic"),
                identical(text_style$family, "Arial"))
    }
    # Modifying a standalone panel must not mutate its retained source plot.
    current <- list(
      labels = panels[[name]]$labels,
      theme = panels[[name]]$theme,
      y_breaks = panels[[name]]$scales$get_scales("y")$breaks,
      y_labels = panels[[name]]$scales$get_scales("y")$labels,
      y_limits = panels[[name]]$scales$get_scales("y")$limits,
      x_limits = panels[[name]]$scales$get_scales("x")$limits
    )
    stopifnot(identical(current, source_state[[name]]))
  }
}
cat("All five panel identities, plotted values, scales, and source plots preserved.\n")
''',
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "All five panel identities" in result.stdout


def test_order_panel_keeps_final_segment_when_total_has_floating_roundoff() -> None:
    result = subprocess.run(
        [
            str(RUN_R),
            "Rscript",
            "-e",
            r'''
source("scripts/publication/gbe_theme.R")
epsilon <- .Machine$double.eps * 2
fractions <- c(.6, .2, .15, .05 + epsilon)
overview_order_composition <- data.frame(
  species = "example", y = 1,
  repeat_order = factor(c("LTR", "LINE", "TIR", "Other"),
                        levels = c("LTR", "LINE", "TIR", "Other")),
  fraction = fractions,
  xmin = c(0, head(cumsum(fractions), -1)),
  xmax = cumsum(fractions)
)
original_data <- overview_order_composition
overview_order_colors <- GBE_TE_ORDER_COLORS
stopifnot(tail(overview_order_composition$xmax, 1) > 1)

expressions <- parse(file = "scripts/publication/build_gbe_figures.R")
panel_assignment <- Filter(function(expression) {
  is.call(expression) && identical(expression[[1]], as.name("<-")) &&
    identical(expression[[2]], as.name("overview_order_panel"))
}, expressions)
stopifnot(length(panel_assignment) == 1L)
eval(panel_assignment[[1]])

rectangles <- ggplot_build(overview_order_panel)$data[[1]]
stopifnot(nrow(rectangles) == 4L,
          all(is.finite(rectangles$xmin)),
          all(is.finite(rectangles$xmax)),
          all(rectangles$xmax - rectangles$xmin > 0),
          identical(tail(rectangles$xmax, 1), 1))
stopifnot(identical(overview_order_composition, original_data),
          identical(overview_order_panel$data$fraction, fractions),
          tail(overview_order_panel$data$xmax, 1) > 1)
cat("All four order segments remain visible without changing source fractions.\n")
''',
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "All four order segments remain visible" in result.stdout


def test_genome_track_labels_uniform_points_without_losing_intervals() -> None:
    result = subprocess.run(
        [str(RUN_R), "Rscript", "-e", r'''
suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("scripts/publication/gbe_theme.R")
traits <- read.csv("Publication/datasets/Supplementary_Data_S04_path_analysis_traits24.csv")
te <- read.csv("Publication/datasets/Supplementary_Data_S42_te_bonsai_tip_annotations.csv")
overview_tips <- traits[sub("^D\\. ", "", traits$species) %in% te$species,]
overview_tips$y <- seq_len(nrow(overview_tips))
stopifnot(nrow(overview_tips) == 21L)
original <- overview_tips
expressions <- parse(file="scripts/publication/build_gbe_figures.R")
assignment <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("overview_genome_panel")), expressions)
stopifnot(length(assignment) == 1L)
eval(assignment[[1]])
built <- ggplot_build(overview_genome_panel)
stopifnot(identical(overview_tips, original),
  identical(built$data[[1]]$x, original$genome_size_pg_ci_low),
  identical(built$data[[1]]$xend, original$genome_size_pg_ci_high),
  identical(built$data[[2]]$x, original$genome_size_pg_fuscus_anchored),
  all(built$data[[2]]$shape == 21), nrow(built$data[[3]]) == 21L,
  identical(built$data[[3]]$label, sprintf("%.1f", original$genome_size_pg_fuscus_anchored)),
  identical(built$data[[3]]$y, original$y + .38),
  is.null(overview_genome_panel$scales$get_scales("shape")),
  all(is.finite(built$data[[1]]$x)), all(is.finite(built$data[[1]]$xend)))
limits <- overview_genome_panel$scales$get_scales("x")$limits
stopifnot(identical(limits, c(8,27)),
  min(original$genome_size_pg_ci_low) >= limits[1],
  max(original$genome_size_pg_ci_high) <= limits[2])
'''], cwd=ROOT, text=True, capture_output=True, timeout=60, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr

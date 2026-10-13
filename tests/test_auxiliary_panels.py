from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]



def run_r_assertions(code: str) -> None:
    result = subprocess.run(
        [
            str(ROOT / "scripts/run_r.sh"),
            "Rscript",
            "-e",
            'source("scripts/publication/gbe_theme.R"); '
            'source("scripts/publication/auxiliary_panels.R"); ' + code,
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_isolated_species_labels_preserve_marks_and_source_axis() -> None:
    run_r_assertions(
        r'''
tips <- data.frame(y = c(3, 1, 2), display_species = c("D. c", "D. a", "D. b"))
points <- data.frame(x = c(0.2, 0.4, 0.8), y = c(1, 2, 3))
plot <- ggplot(points, aes(x, y)) + geom_point() +
  scale_y_continuous(limits = c(0.5, 3.5), expand = c(0, 0)) +
  theme_gbe() + theme(axis.text.y = element_blank())
result <- auxiliary_label_rows(plot, tips, "display_species")
stopifnot(identical(ggplot_build(plot)$data, ggplot_build(result)$data),
          identical(plot$scales$get_scales("y")$limits,
                    result$scales$get_scales("y")$limits),
          identical(result$scales$get_scales("y")$breaks, tips$y),
          identical(result$scales$get_scales("y")$labels, tips$display_species),
          inherits(plot$theme$axis.text.y, "element_blank"),
          identical(result$theme$axis.text.y$family, "Arial"))
bad <- tips; bad$y[[1]] <- bad$y[[2]]
stopifnot(inherits(try(auxiliary_label_rows(plot, bad, "display_species"), silent=TRUE), "try-error"))
'''
    )


def test_isolated_dag_contains_only_selected_edges_without_mutating_collage() -> None:
    run_r_assertions(
        r'''
edges <- data.frame(model=c("one", "two"), panel_label=c("First", "Second"),
                    x=c(.1,.2), y=c(.3,.4), xend=c(.5,.6), yend=c(.7,.8))
nodes <- data.frame(model=rep(c("one", "two"), each=3),
                    panel_label=rep(c("First", "Second"), each=3),
                    x=rep(c(.15,.85,.5), 2), y=rep(c(.23,.23,.78), 2))
plot <- ggplot() + geom_segment(data=edges, aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data=nodes, aes(x=x,y=y)) + facet_wrap(vars(panel_label)) +
  scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1))
before <- lapply(plot$layers, function(layer) layer$data)
result <- auxiliary_dag_panel(plot, "two", "Second")
stopifnot(identical(before, lapply(plot$layers, function(layer) layer$data)),
          identical(result$layers[[1]]$data, edges[2,,drop=FALSE]),
          identical(result$layers[[2]]$data, nodes[4:6,,drop=FALSE]),
          inherits(result$facet, "FacetNull"), inherits(plot$facet, "FacetWrap"),
          identical(result$scales$get_scales("x")$limits, c(0,1)),
          identical(result$scales$get_scales("y")$limits, c(0,1)))
'''
    )


def test_auxiliary_catalog_exports_every_collage_track_and_class() -> None:
    run_r_assertions(
        r'''
p <- ggplot(data.frame(x=1,y=1), aes(x,y)) + geom_point() +
     scale_y_continuous(limits=c(.5,1.5), expand=c(0,0))
integrated_tree_tips <- data.frame(y=1,display_species="D. example")
path24_development_traits <- data.frame(y=1,species="D. example")
for (name in c("integrated_tree_panel", "integrated_genome_panel", "integrated_composition_panel",
               "integrated_fraction_panel", "integrated_profile_panel", "wang_diversity_shannon",
               "wang_diversity_gini", "path24_tree_reference",
               "path24_genome_reference", "path24_nucleus_reference", "path24_cell_reference",
               "Path_DAG_reference", "wang_a", "wang_b", "wang_c", "wang_d")) assign(name,p)
representative_audit <- data.frame(model_order=1:10, model=paste0("m",1:10),
   display_name=paste0("Model ",1:10), panel_label=paste0("Model ",1:10), class_members=rep(1L,10))
ids <- c("Auxiliary_Wang_style_integrated_te_history_candidate",
         "Auxiliary_Wang_style_diversity_genome_size_candidate",
         "Auxiliary_Path24_development_traits_reference", "Auxiliary_Path_DAG_reference",
         "Auxiliary_Wang2025_cross_study_te_diversity_reference")
specs <- lapply(ids,auxiliary_review_spec)
stopifnot(identical(vapply(specs, function(x) x$number, integer(1)), 8:12),
          identical(vapply(specs, function(x) length(x$plots), integer(1)), c(5L,2L,4L,10L,4L)),
          identical(specs[[3]]$panel_ids, c("a", "c", "d", "e")))
for (spec in specs) {
  stopifnot(all(vapply(spec[c("height","titles","slugs","legends")], length, integer(1)) == length(spec$plots)),
            !anyDuplicated(spec$slugs), all(nzchar(spec$legends)))
}
stopifnot(inherits(try(auxiliary_review_spec("unknown"), silent=TRUE), "try-error"))
'''
    )


def test_path24_shared_axes_retain_every_released_measurement_interval() -> None:
    run_r_assertions(
        r'''
traits <- read.csv("Publication/datasets/Supplementary_Data_S04_path_analysis_traits24.csv")
morphology <- read.csv("Publication/datasets/Supplementary_Data_S19_cell_nucleus_species_estimates.csv")
expressions <- parse(file="scripts/publication/build_gbe_figures.R")
path24_development_traits <- merge(traits[,c("species","genome_size_pg_fuscus_anchored",
  "genome_size_pg_ci_low","genome_size_pg_ci_high")], morphology, by="species")
path24_development_traits$y <- seq_len(nrow(path24_development_traits))
tip_count <- 24L
assignment <- function(name) {
  matches <- which(vapply(expressions,function(expression) is.call(expression) &&
    identical(expression[[1]],as.name("<-")) && identical(expression[[2]],as.name(name)), logical(1)))
  expressions[[tail(matches,1)]]
}
eval(assignment("aligned_trait_panel"), envir=.GlobalEnv)
specs <- list(
  list(name="path24_genome_reference", estimate="genome_size_pg_fuscus_anchored",
       low="genome_size_pg_ci_low", high="genome_size_pg_ci_high", limits=c(8,27), tag="b"),
  list(name="path24_nucleus_reference", estimate="nucleus_area_um2",
       low="nucleus_area_um2_ci_low", high="nucleus_area_um2_ci_high", limits=c(35,55), tag="c"),
  list(name="path24_cell_reference", estimate="cell_area_um2",
       low="cell_area_um2_ci_low", high="cell_area_um2_ci_high", limits=c(190,335), tag="d"))
for (spec in specs) {
  data <- path24_development_traits
  eval(assignment(spec$name), envir=.GlobalEnv)
  result <- get(spec$name)
  built <- ggplot_build(result)
  stopifnot(identical(result$scales$get_scales("x")$limits, spec$limits),
            length(result$scales$get_scales("x")$breaks) >= 5L,
            identical(built$data[[1]]$xmin,data[[spec$low]]),
            identical(built$data[[1]]$xmax,data[[spec$high]]),
            identical(built$data[[2]]$x,data[[spec$estimate]]),
            identical(built$data[[2]]$y,as.numeric(data$y)),
            identical(built$data[[3]]$label,sprintf("%.1f",data[[spec$estimate]])),
            isTRUE(all.equal(built$data[[3]]$y,data$y + 0.38)),
            all(built$data[[3]]$family == "Arial"),
            all(built$data[[3]]$hjust[data[[spec$estimate]] >
              spec$limits[[2]] - 0.15 * diff(spec$limits)] == 1),
            all(built$data[[3]]$hjust[data[[spec$estimate]] <
              spec$limits[[1]] + 0.15 * diff(spec$limits)] == 0),
            identical(result$coordinates$clip,"off"),
            identical(result$scales$get_scales("x")$expand, expansion(mult=0.04)),
            all(complete.cases(built$data[[1]][,c("xmin","xmax","y")])),
            identical(result$labels$tag,spec$tag), nrow(data)==24L)
}
stopifnot(inherits(try(aligned_trait_panel(path24_development_traits,
  "cell_area_um2", "cell_area_um2_ci_low", "cell_area_um2_ci_high", "Cell", "blue",
  c(200,330), c(200,250,300)), silent=TRUE), "try-error"))
'''
    )


def test_path24_exports_keep_note_identifiers_with_contiguous_display_letters() -> None:
    run_r_assertions(
        r'''
library(dplyr)
library(readr)
p <- ggplot(data.frame(x=1,y=1), aes(x,y)) + geom_point() +
  scale_y_continuous(limits=c(.5,1.5), expand=c(0,0))
path24_development_traits <- data.frame(y=1,species="D. example")
path24_tree_reference <- path24_genome_reference <- path24_nucleus_reference <-
  path24_cell_reference <- p
saved <- list()
save_gbe_figure <- function(plot, basename, width_mm, height_mm, dpi) {
  saved[[length(saved)+1L]] <<- list(basename=basename, tag=plot$labels$tag)
}
resolve_gbe_font <- function() "Arial"
out <- tempfile(); dir.create(file.path(out,"panels"), recursive=TRUE)
manifest <- export_auxiliary_panels("Auxiliary_Path24_development_traits_reference",out)
stopifnot(identical(manifest$panel_id,c("a","c","d","e")),
  identical(manifest$display_panel_id,c("a","b","c","d")),
  identical(vapply(saved,function(x) x$tag, character(1)),c("a","b","c","d")),
  identical(manifest$pdf,c("Figure_10a_species_phylogeny.pdf",
    "Figure_10c_conditional_genome_size.pdf", "Figure_10d_nucleus_area.pdf",
    "Figure_10e_cell_area.pdf")))
legends <- readLines(file.path(out,"panels","PANEL_LEGENDS_AND_ALT_TEXT.md"))
stopifnot("## b - Estimated genome size" %in% legends,
  !any(grepl("Development mode$",legends)))
# The optional identifier mapping does not alter other auxiliary manifests.
wang_diversity_shannon <- wang_diversity_gini <- p
other <- export_auxiliary_panels("Auxiliary_Wang_style_diversity_genome_size_candidate",out)
stopifnot(identical(other$panel_id,c("a","b")), !"display_panel_id" %in% names(other),
  identical(other$pdf,c("Figure_9a_Shannon_diversity.pdf","Figure_9b_Gini_Simpson_diversity.pdf")))
unlink(out,recursive=TRUE)
'''
    )

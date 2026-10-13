from __future__ import annotations

import csv
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "analyses/03_phylogenetic_path/data/candidate_dags.csv"
COMPARISON = ROOT / "analyses/03_phylogenetic_path/data/path_model_comparison.csv"
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
FIGURE_DIR = ROOT / "Publication/figures"
FIGURE_MANIFEST = FIGURE_DIR / "FIGURE_MANIFEST.csv"
AUXILIARY_MANIFEST = FIGURE_DIR / "AUXILIARY_FIGURE_MANIFEST.csv"
BASENAME = "Path_DAG_equivalence_class_reference"

REPRESENTATIVES = {
    "independent": ("no_edges", "independent", 1),
    "genome_nucleus_only": ("GS->NS", "genome_nucleus_only", 2),
    "nucleus_cell_only": ("NS->CS", "nucleus_cell_only", 2),
    "genome_cell_only": ("GS->CS", "genome_cell_only", 2),
    "nucleus_bridge": ("GS->NS;NS->CS", "nucleus_bridge", 3),
    "nucleus_collider": ("GS->NS;CS->NS", "nucleus_collider", 1),
    "cell_bridge": ("GS->CS;CS->NS", "cell_bridge", 3),
    "cell_collider": ("GS->CS;NS->CS", "cell_collider", 1),
    "genome_bridge": ("GS->CS;NS->GS", "genome_bridge", 3),
    "genome_collider": ("NS->GS;CS->GS", "genome_collider", 1),
}


def test_representative_dags_match_the_frozen_candidate_classes() -> None:
    candidates = pd.read_csv(CANDIDATES).set_index("dag_id")
    comparison = pd.read_csv(COMPARISON)
    primary = comparison.loc[comparison["analysis_type"].eq("primary")]

    assert set(primary["model"]) == set(REPRESENTATIVES)
    assert sorted(primary["rank"].tolist()) == list(range(1, 11))
    assert len(candidates) == 25
    assert candidates["equivalence_class"].nunique() == 11
    assert candidates.loc[candidates["testable_by_dsep"]].equivalence_class.nunique() == 10

    for model, (dag_id, equivalence_class, member_count) in REPRESENTATIVES.items():
        row = candidates.loc[dag_id]
        assert row["equivalence_class"] == equivalence_class, model
        assert int(row["class_member_count"]) == member_count, model
        assert bool(row["testable_by_dsep"]), model


def test_dag_reference_is_built_as_an_unnumbered_auxiliary_figure() -> None:
    source = FIGURE_SOURCE.read_text(encoding="utf-8")
    block = source.split("## ---- path-dag-reference", maxsplit=1)[1].split(
        "## ---- figure-manifest", maxsplit=1
    )[0]

    assert "candidate_dags.csv" in block
    assert "Supplementary_Data_S30_path_model_comparison.csv" in block
    assert BASENAME in block
    assert "record_figure(" not in block
    assert "record_auxiliary_figure(" in block
    assert 'node_label = c("Genome\\nsize", "Nucleus\\narea", "Cell\\narea")' in block
    assert '"Relative\\nIOD"' not in block
    for model in REPRESENTATIVES:
        assert f'"{model}"' in block

    with FIGURE_MANIFEST.open(newline="", encoding="utf-8") as stream:
        manifest = list(csv.DictReader(stream))
    assert all(row["basename"] != BASENAME for row in manifest)

    with AUXILIARY_MANIFEST.open(newline="", encoding="utf-8") as stream:
        auxiliary_manifest = list(csv.DictReader(stream))
    dag_rows = [row for row in auxiliary_manifest if row["basename"] == BASENAME]
    assert len(dag_rows) == 1
    assert dag_rows[0]["auxiliary_id"] == "Auxiliary_Path_DAG_reference"


def test_dag_reference_assets_exist_in_all_release_formats() -> None:
    for suffix in ("pdf", "png", "tif"):
        path = FIGURE_DIR / f"{BASENAME}.{suffix}"
        assert path.is_file()
        assert path.stat().st_size > 10_000


def test_independent_dag_keeps_three_nodes_without_arrows_or_extra_annotation() -> None:
    code = r'''
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
source("scripts/publication/gbe_theme.R")
source("scripts/publication/auxiliary_panels.R")
ROOT <- normalizePath(".")
read_supp <- function(name) read_csv(file.path(ROOT, "Publication/datasets", name),
                                   show_col_types=FALSE, progress=FALSE)
lines <- readLines("scripts/publication/build_gbe_figures.R")
start <- match("## ---- path-dag-reference", lines)
end <- match("## ---- wang2025-cross-study-diversity-reference", lines)
stopifnot(!is.na(start), !is.na(end), end > start)
expressions <- parse(text=lines[seq.int(start+1L, end-1L)])
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name("record_auxiliary_figure"))) break
  eval(expression, envir=.GlobalEnv)
}
built <- ggplot_build(Path_DAG_reference)
independent <- built$layout$layout$PANEL[
  as.character(built$layout$layout$panel_label) == "Independent\n(1 equivalent DAG)"
]
stopifnot(length(independent) == 1L, nrow(built$layout$layout) == 10L)
text_indices <- which(vapply(Path_DAG_reference$layers,
  function(layer) inherits(layer$geom, "GeomText"), logical(1)))
labels <- bind_rows(built$data[text_indices])
stopifnot(nrow(labels) == 30L,
          setequal(labels$label, c("Genome\nsize", "Nucleus\narea", "Cell\narea")))
for (plot in list(Path_DAG_reference,
                 auxiliary_dag_panel(Path_DAG_reference, "independent", "Independent"))) {
  data <- ggplot_build(plot)$data
  points <- bind_rows(data[which(vapply(plot$layers,
    function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))])
  edges <- bind_rows(data[which(vapply(plot$layers,
    function(layer) inherits(layer$geom, "GeomSegment"), logical(1)))])
  if (inherits(plot$facet, "FacetWrap")) {
    points <- points[points$PANEL == independent, ]
    edges <- edges[edges$PANEL == independent, ]
  }
  stopifnot(nrow(points) == 3L, nrow(edges) == 0L,
            identical(sort(points$x), c(.15, .5, .85)),
            identical(sort(points$y), c(.23, .23, .78)))
}
caption <- tolower(Path_DAG_reference$labels$caption)
stopifnot(grepl("d. fuscus", caption, fixed=TRUE),
          grepl("assembly-anchored", caption, fixed=TRUE),
          grepl("pg/1c", caption, fixed=TRUE),
          grepl("not a direct species c-value", caption, fixed=TRUE),
          grepl("arrows are not unique within classes", caption, fixed=TRUE),
          grepl("saturated class (6 members) unscored", caption, fixed=TRUE))
traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
stopifnot(nrow(traits) == 24L,
          all(grepl("derived from nuclear IOD", traits$genome_size_definition, fixed=TRUE)),
          all(grepl("conditional D. fuscus assembly-anchored", traits$genome_size_definition,
                    fixed=TRUE)))
'''
    result = subprocess.run(
        [str(ROOT / "scripts/run_r.sh"), "Rscript", "-e", code],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr

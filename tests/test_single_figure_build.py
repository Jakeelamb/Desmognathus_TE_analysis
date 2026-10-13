from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILD = "scripts/publication/build_gbe_figures.R"
RUN_R = ROOT / "scripts/run_r.sh"



def _run_r(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(RUN_R), "Rscript", *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=30,
        check=False,
    )


@pytest.mark.parametrize(
    ("args", "error"),
    [
        (["--figure", "Figure_13"], "Supported figures: Figure_1 through Figure_12"),
        (["--figure", "Figure_1"], "--figure requires --output-dir"),
        (["--figure"], "Missing value for --figure"),
        (["--output-dir"], "Missing value for --output-dir"),
        (["--output-dir", ""], "Missing value for --output-dir"),
        (["--figure", "--output-dir", "unused"], "Missing value for --figure"),
        (["--figure", "Figure_1", "--figure", "Figure_1"], "Duplicate option: --figure"),
        (["--unknown"], "Unknown option: --unknown"),
    ],
)
def test_invalid_selection_fails_before_rendering(args: list[str], error: str) -> None:
    result = _run_r(BUILD, *args)
    assert result.returncode != 0
    assert error in result.stderr


@pytest.mark.parametrize("use_symlink", [False, True])
def test_single_figure_cannot_replace_full_release(tmp_path: Path, use_symlink: bool) -> None:
    release_dir = ROOT / "Publication/figures"
    output = release_dir
    if use_symlink:
        output = tmp_path / "release-alias"
        output.symlink_to(release_dir, target_is_directory=True)
    protected = [release_dir / name for name in ("FIGURE_MANIFEST.csv", "FIGURE_LEGENDS_AND_ALT_TEXT.md")]
    before = [path.read_bytes() for path in protected]

    result = _run_r(BUILD, "--figure", "Figure_1", "--output-dir", str(output))

    assert result.returncode != 0
    assert "separate from Publication/figures" in result.stderr
    assert [path.read_bytes() for path in protected] == before


def _load_function(name: str) -> str:
    """Load a production helper without evaluating the figure builds."""
    return f"""expressions <- parse(file = {json.dumps(BUILD)})
for (expression in expressions) {{
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      identical(expression[[2]], as.name({json.dumps(name)}))) {{
    eval(expression)
  }}
}}
"""


def test_single_figure_cannot_replace_another_candidate(tmp_path: Path) -> None:
    manifest = tmp_path / "FIGURE_MANIFEST.csv"
    original = "figure_id,basename\nFigure_1,retained_overview\n"
    manifest.write_text(original)
    result = _run_r(BUILD, "--figure", "Figure_2", "--output-dir", str(tmp_path))
    assert result.returncode != 0
    assert "belongs to a different figure" in result.stderr
    assert manifest.read_text() == original


def test_selection_parser_preserves_full_build_default_and_accepts_candidate(tmp_path: Path) -> None:
    result = _run_r(
        "-e",
        _load_function("parse_figure_options")
        + f"""
stopifnot(identical(parse_figure_options(character()), list(figure = NULL, output_dir = NULL)))
selected <- parse_figure_options(c("--output-dir", {json.dumps(str(tmp_path))}, "--figure", "Figure_1"))
stopifnot(identical(selected$figure, "Figure_1"), identical(selected$output_dir, {json.dumps(str(tmp_path))}))
second <- parse_figure_options(c("--output-dir", {json.dumps(str(tmp_path))}, "--figure", "Figure_2"))
stopifnot(identical(second$figure, "Figure_2"))
third <- parse_figure_options(c("--output-dir", {json.dumps(str(tmp_path))}, "--figure", "Figure_3"))
stopifnot(identical(third$figure, "Figure_3"))
for (id in paste0("Figure_", 4:12)) {{
  selected <- parse_figure_options(c("--output-dir", {json.dumps(str(tmp_path))}, "--figure", id))
  stopifnot(identical(selected$figure, id))
}}
""",
    )
    assert result.returncode == 0, result.stderr


def test_metadata_writer_keeps_candidate_records_in_selected_directory(tmp_path: Path) -> None:
    result = _run_r(
        "-e",
        _load_function("write_numbered_figure_metadata")
        + f"""
suppressPackageStartupMessages({{ library(dplyr); library(readr) }})
OUT <- {json.dumps(str(tmp_path))}
manifest_rows <- list(tibble(
  figure_id = "Figure_1", basename = "candidate", legend = "Candidate legend.",
  alt_text = "Candidate alt text.", pdf = "candidate.pdf", png = "candidate.png", tiff = "candidate.tif"
))
result <- write_numbered_figure_metadata()
stopifnot(nrow(result) == 1L, identical(result$figure_id, "Figure_1"))
""",
    )
    assert result.returncode == 0, result.stderr
    with (tmp_path / "FIGURE_MANIFEST.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["figure_id"] == "Figure_1"
    assert rows[0]["pdf"] == "candidate.pdf"
    legend = (tmp_path / "FIGURE_LEGENDS_AND_ALT_TEXT.md").read_text()
    assert "## Figure_1\n\nCandidate legend." in legend
    assert "**Alt text:** Candidate alt text." in legend
    assert {path.name for path in tmp_path.iterdir()} == {
        "FIGURE_MANIFEST.csv", "FIGURE_LEGENDS_AND_ALT_TEXT.md"
    }


def test_auxiliary_selection_skips_unrequested_exports() -> None:
    result = _run_r(
        "-e",
        _load_function("record_auxiliary_figure")
        + '''
figure_options <- list(figure = "Figure_5")
save_gbe_figure <- function(...) stop("An unselected figure was rendered")
result <- record_auxiliary_figure(
  "Auxiliary_Path_DAG_reference", "Reference", NULL, "unused", 180, 220,
  "S30", "Boundary", "Legend", "Alt", "review"
)
stopifnot(is.null(result))
''',
    )
    assert result.returncode == 0, result.stderr


def test_auxiliary_selector_writes_review_id_and_preserves_boundaries(tmp_path: Path) -> None:
    result = _run_r(
        "-e",
        _load_function("record_auxiliary_figure")
        + _load_function("write_numbered_figure_metadata")
        + f'''
suppressPackageStartupMessages({{library(dplyr); library(readr)}})
figure_options <- list(figure = "Figure_11")
ROOT <- "."
OUT <- {json.dumps(str(tmp_path))}
GBE_MIN_TEXT_PT <- 7
resolve_gbe_font <- function() "Arial"
manifest_rows <- list()
save_gbe_figure <- function(...) invisible(NULL)
source <- function(...) invisible(NULL)
exported_id <- NULL
export_auxiliary_panels <- function(id, output_dir) {{
  stopifnot(identical(output_dir, OUT))
  exported_id <<- id
}}
quit <- function(...) stop("selected export complete")
finished <- tryCatch(record_auxiliary_figure(
  "Auxiliary_Path_DAG_reference", "Path reference", NULL, "dag", 180, 220,
  "S30", "No causal orientation.", "Ten equivalence classes.", "Ten diagrams.", "review"
), error = conditionMessage)
stopifnot(identical(finished, "selected export complete"),
          identical(exported_id, "Auxiliary_Path_DAG_reference"))
''',
    )
    assert result.returncode == 0, result.stderr
    with (tmp_path / "FIGURE_MANIFEST.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["figure_id"] == "Figure_11"
    assert rows[0]["title"] == "Path reference"
    assert rows[0]["legend"] == "Ten equivalence classes. No causal orientation."
    assert rows[0]["pdf"] == "dag.pdf"

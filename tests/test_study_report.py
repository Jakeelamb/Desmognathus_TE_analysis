from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Publication/Desmognathus_study_data.Rmd"
OUTPUT = ROOT / "Publication/Desmognathus_study_data.html"
PATH_README = ROOT / "analyses/03_phylogenetic_path/README.md"


def test_study_report_replaces_notebook_surface() -> None:
    assert SOURCE.exists()
    assert OUTPUT.exists()
    assert not list((ROOT / "analyses").glob("*/explore.ipynb"))
    assert not (ROOT / "scripts/build_notebooks.py").exists()

    project = (ROOT / "pyproject.toml").read_text()
    assert "jupyterlab" not in project
    assert "nbformat" not in project

    environment = (ROOT / "environment.yml").read_text()
    assert "r-rmarkdown=2.31" in environment

    assert not any(line.rstrip("\r\n") != line.rstrip() for line in OUTPUT.open())


def test_study_report_exposes_current_release_contract() -> None:
    source = SOURCE.read_text()
    rendered = OUTPUT.read_text()

    required_source_phrases = {
        "nrow(te_diversity) == 68L",
        "Unresolved aligned-base",
        "S07",
        "primary_iqr_filtered_coverage_ge_80pct",
        "relative_iod_index",
        "0.0144",
        "Q1 - 1.5",
        "Gini-Simpson",
        "scripts/publication/build_gbe_figures.R",
    }
    assert all(phrase in source for phrase in required_source_phrases)
    assert "mass_aware_unresolved_bin" not in source
    assert "mass_aware_unresolved_bin" not in rendered
    assert "composition_mode" not in source
    assert "write_csv(" not in source
    assert "write.csv(" not in source
    assert "ggsave(" not in source

    source_sha256 = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert source_sha256 in rendered
    for expected in [
        "TE34",
        "Path24",
        "overlap21",
        "380",
        "1,152",
        "805",
        "Figure 2. TE diversity and compositional ordination.",
        "Figure 4. LTR terminal:internal mapping/deletion-footprint proxy.",
        "Figure 5. Reviewed upper-tail cell and nucleus morphology.",
        "Figure 6. Relative IOD, nucleus area, and cell area sensitivities.",
    ]:
        assert expected in rendered
    assert rendered.count("data:image/png;base64") >= 8

    for error in ["Quitting from lines", "Execution halted", "Error in "]:
        assert error not in rendered


def test_study_report_sources_include_unnumbered_dag_reference_contract() -> None:
    source = SOURCE.read_text()
    readme = PATH_README.read_text()

    required_source_phrases = {
        "Auxiliary DAG reference (unnumbered)",
        "10\ntestable Markov-equivalence classes",
        "The saturated 6-member class is omitted",
        "not absolute genome size",
        "arrows shown here are representative only, not uniquely oriented causal claims",
        'figure_path("Path_DAG_equivalence_class_reference.png")',
    }
    assert all(phrase in source for phrase in required_source_phrases)
    assert "dag_reference <- tibble::tribble(" not in source
    assert 'title = "Auxiliary DAG reference for Figure 7"' not in source
    assert "Representative of 3 DAGs; ΔCICc = 0.000; Passed" not in source

    required_readme_phrases = {
        "25 possible three-node DAG members",
        "10 testable",
        "6-member saturated class",
        "one representative DAG per scored",
        "relative nuclear IOD is",
        "not absolute genome size",
    }
    assert all(phrase in readme for phrase in required_readme_phrases)

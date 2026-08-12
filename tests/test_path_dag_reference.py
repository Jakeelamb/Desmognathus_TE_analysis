from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "analyses/03_phylogenetic_path/data/candidate_dags.csv"
COMPARISON = ROOT / "analyses/03_phylogenetic_path/data/path_model_comparison.csv"
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
FIGURE_DIR = ROOT / "Publication/figures"
FIGURE_MANIFEST = FIGURE_DIR / "FIGURE_MANIFEST.csv"
BASENAME = "Path_DAG_equivalence_class_reference"

REPRESENTATIVES = {
    "independent": ("no_edges", "independent", 1),
    "iod_nucleus_only": ("GS->NS", "genome_nucleus_only", 2),
    "nucleus_cell_only": ("NS->CS", "nucleus_cell_only", 2),
    "iod_cell_only": ("GS->CS", "genome_cell_only", 2),
    "nucleus_bridge": ("GS->NS;NS->CS", "nucleus_bridge", 3),
    "nucleus_collider": ("GS->NS;CS->NS", "nucleus_collider", 1),
    "cell_bridge": ("GS->CS;CS->NS", "cell_bridge", 3),
    "cell_collider": ("GS->CS;NS->CS", "cell_collider", 1),
    "iod_bridge": ("GS->CS;NS->GS", "genome_bridge", 3),
    "iod_collider": ("NS->GS;CS->GS", "genome_collider", 1),
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

    for model, (dag_id, legacy_class, member_count) in REPRESENTATIVES.items():
        row = candidates.loc[dag_id]
        assert row["equivalence_class"] == legacy_class, model
        assert int(row["class_member_count"]) == member_count, model
        assert bool(row["testable_by_dsep"]), model


def test_dag_reference_is_built_as_an_unnumbered_auxiliary_figure() -> None:
    source = FIGURE_SOURCE.read_text(encoding="utf-8")
    block = source.split("## ---- path-dag-reference", maxsplit=1)[1].split(
        "## ---- figure-s1", maxsplit=1
    )[0]

    assert "candidate_dags.csv" in block
    assert "Supplementary_Data_S30_path_model_comparison.csv" in block
    assert BASENAME in block
    assert "record_figure(" not in block
    assert "save_gbe_figure(" in block
    assert 'node_label = c("Relative\\nIOD", "Nucleus\\narea", "Cell\\narea")' in block
    assert "standardized log10 traits" in block
    assert "arrows are not unique within classes" in block
    assert "Saturated class (6 members) unscored" in block
    assert "genome size" not in block.lower()
    for model in REPRESENTATIVES:
        assert f'"{model}"' in block

    with FIGURE_MANIFEST.open(newline="", encoding="utf-8") as stream:
        manifest = list(csv.DictReader(stream))
    assert all(row["basename"] != BASENAME for row in manifest)


def test_dag_reference_assets_exist_in_all_release_formats() -> None:
    for suffix in ("pdf", "png", "tif"):
        path = FIGURE_DIR / f"{BASENAME}.{suffix}"
        assert path.is_file()
        assert path.stat().st_size > 10_000

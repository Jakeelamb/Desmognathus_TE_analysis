from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METHODS = ROOT / "Publication/methods"


def test_methods_evidence_matrix_uses_only_active_sources() -> None:
    evidence = pd.read_csv(METHODS / "METHODS_EVIDENCE_MATRIX.csv")
    assert evidence["section_id"].tolist() == [f"M{number:02d}" for number in range(1, 21)]

    forbidden = ("input_data/", "path_analysis/", "results/", "plans/")
    for source_list in evidence["primary_sources"].astype(str):
        for source in (item.strip() for item in source_list.split(";")):
            assert not source.startswith(forbidden), source
            if source == "none" or source.startswith("doi:"):
                continue
            assert (ROOT / source).is_file(), source


def test_methods_scaffold_records_finalized_release_boundaries() -> None:
    scaffold = (METHODS / "METHODS_BULLET_SCAFFOLD.md").read_text(encoding="utf-8")
    required_phrases = {
        "136 rows total",
        "classified-conditional order-level stratum is the primary paper-facing diversity analysis",
        "Natural-log Shannon entropy was `H′ = −Σpᵢ ln(pᵢ)`",
        "explicitly named Gini-Simpson index was `1 − Σpᵢ²`",
        "These are the only paper-facing diversity indices",
        "Observed richness is retained as an audit/support count",
        "Rare features excluded from the shared-feature PCA remain in these diversity calculations",
        "left LTR, right LTR, and internal region",
        "inclusive, within-species, two-sided Tukey filter",
        "380 primary elements",
        "whitespace-delimited `DOMAIN|MODEL` annotations",
        "0.12 µm/pixel",
        "Seven species are represented by one animal",
        "`nuclear_iod_by_image.csv`",
        "fail this gate and remain visibly flagged",
    }
    assert all(phrase in scaffold for phrase in required_phrases)
    assert "Hill q1" not in scaffold
    assert "Hill q2" not in scaffold
    assert "Pielou" not in scaffold

def test_unresolved_methods_are_recorded_without_recreating_old_paths() -> None:
    required = pd.read_csv(METHODS / "METHODS_REQUIRED_AUTHOR_INPUT.csv").set_index(
        "item_id"
    )
    assert required.loc["W12", "status"] == "partial"
    assert required.loc["W16", "status"] == "partial"
    assert "multimapper policy" in required.loc["W12", "required_detail"]
    assert "checkpoint-to-mask-release mapping" in required.loc["W16", "required_detail"]

    old_image_table = (
        ROOT / "analyses/02_morphology/data/relative_nuclear_iod_by_image.csv"
    )
    new_image_table = ROOT / "analyses/02_morphology/data/nuclear_iod_by_image.csv"
    assert not old_image_table.exists()
    assert new_image_table.is_file()

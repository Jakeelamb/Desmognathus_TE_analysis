#!/usr/bin/env python3
"""Tests for the corrected all-position LTR30 sensitivity branch."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.processing.build_corrected_ectopic_ltr30 import (
    ELEMENT_OUTPUT,
    EXCLUDED_OUTPUT,
    EXPECTED_CORRUPT_DEPTH_ELEMENT_IDS,
    HISTORICAL_LTR30,
    MANIFEST_OUTPUT,
    RESOURCE_OUTPUT,
    SPECIES_OUTPUT,
    calculate_zero_aware_depth_metrics,
    inspect_depth_source,
    merge_elements_with_tesorter,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_exact_tesorter_join_uses_sequence_and_element_coordinates() -> None:
    elements = pd.DataFrame(
        {
            "sequence": ["contig", "contig"],
            "element start": [10, 100],
            "element end": [50, 150],
        }
    )
    tesorter = pd.DataFrame(
        {
            "#TE": ["contig_100_150", "contig_10_50"],
            "Domains": ["GAG|RT|RH", "GAG|RT"],
        }
    )

    merged = merge_elements_with_tesorter(elements, tesorter)

    assert merged["element_id"].tolist() == ["contig_10_50", "contig_100_150"]
    assert merged["#TE"].tolist() == ["contig_10_50", "contig_100_150"]
    assert merged["Domains"].tolist() == ["GAG|RT", "GAG|RT|RH"]


def test_zero_aware_metrics_retain_unreported_and_explicit_zero_positions() -> None:
    depth = pd.DataFrame(
        {
            "position": [1, 2, 3, 5, 6, 8, 9, 10],
            "depth": [2, 0, 1, 0, 1, 1, 0, 2],
        }
    )

    metrics = calculate_zero_aware_depth_metrics(
        depth=depth,
        element_length=10,
        left_ltr_length=2,
        right_ltr_length=2,
    )

    assert metrics["depth_positions_expected"] == 10
    assert metrics["depth_positions_reported"] == 8
    assert metrics["depth_positions_missing"] == 2
    assert metrics["depth_positions_explicit_zero"] == 3
    assert metrics["mean_depth_terminal_all_positions"] == 1.0
    assert metrics["mean_depth_internal_all_positions"] == 0.5
    assert metrics["ratio_terminal_internal_all_positions"] == 2.0


def test_depth_audit_rejects_a_partial_final_line_without_modifying_source(
    tmp_path: Path,
) -> None:
    depth_path = tmp_path / "element.fa.depth.txt"
    original = b"element\t1\t4\nelement\t2\t0\nelement\t3"
    depth_path.write_bytes(original)

    audit = inspect_depth_source(
        path=depth_path,
        element_id="element",
        element_length=4,
    )

    assert audit["usable"] is False
    assert audit["partial_final_line"] is True
    assert audit["expected_terminal_position"] == 4
    assert audit["reported_terminal_position"] == 2
    assert audit["exclusion_reason"] == "truncated_depth_file_partial_final_line"
    assert depth_path.read_bytes() == original


def test_ltr30_outputs_are_complete_identity_safe_and_non_destructive() -> None:
    for path in (
        ELEMENT_OUTPUT,
        SPECIES_OUTPUT,
        EXCLUDED_OUTPUT,
        RESOURCE_OUTPUT,
        MANIFEST_OUTPUT,
    ):
        assert path.is_file(), path

    manifest = json.loads(MANIFEST_OUTPUT.read_text(encoding="utf-8"))
    assert manifest["analysis_scope"] == "audited_te_resource_panel34_ltr30"
    assert manifest["n_panel_species"] == 34
    assert manifest["n_species_with_usable_elements"] == 30
    assert manifest["n_selected_elements"] == 1088
    assert manifest["n_source_corrupt_exclusions"] == 2
    assert manifest["n_usable_elements"] == 1086
    assert manifest["zero_depth_positions_retained"] is True
    assert manifest["exact_tesorter_join"] is True
    assert manifest["expensive_upstream_tools_executed"] is False
    assert manifest["historical_ltr30_source_sha256"] == sha256(HISTORICAL_LTR30)

    elements = pd.read_csv(ELEMENT_OUTPUT)
    assert len(elements) == 1086
    assert elements["species"].nunique() == 30
    assert elements["element_id"].is_unique
    assert elements["element_id"].equals(elements["#TE"])
    assert elements["tesorter_annotation_found"].all()
    assert elements["ratio_terminal_internal_all_positions"].notna().all()
    assert elements["ratio_terminal_internal_all_positions"].gt(0).all()
    assert elements["nonzero_ratio_reproduction_abs_error"].max() < 1e-10

    excluded = pd.read_csv(EXCLUDED_OUTPUT)
    assert len(excluded) == 2
    assert set(excluded["element_id"]) == EXPECTED_CORRUPT_DEPTH_ELEMENT_IDS
    assert excluded["partial_final_line"].all()
    assert excluded["depth_file_size_bytes"].eq(245760).all()
    assert excluded["exclusion_reason"].eq(
        "truncated_depth_file_partial_final_line"
    ).all()
    observed_positions = dict(
        zip(excluded["element_id"], excluded["reported_terminal_position"])
    )
    assert observed_positions["JAUEJG010675316.1_De_4598_10648"] == 5945
    assert observed_positions["JASANM010280042.1_De_179_6669"] == 6396

    resources = pd.read_csv(RESOURCE_OUTPUT)
    assert len(resources) == 34
    assert resources["species"].is_unique
    assert resources["n_selected_5_6_domain_elements"].sum() == 1088
    assert resources["n_source_corrupt_exclusions"].sum() == 2
    assert resources["n_usable_elements"].sum() == 1086
    assert set(resources.loc[resources["resource_status"].eq("missing_tabout"), "species"]) == {
        "catahoula",
        "kanawha",
        "valtos",
    }
    assert set(
        resources.loc[resources["resource_status"].eq("no_5_6_domain_elements"), "species"]
    ) == {"lycos"}

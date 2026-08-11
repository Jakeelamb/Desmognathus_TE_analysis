from __future__ import annotations

import hashlib
import json
import runpy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
TE_DATA = ROOT / "analyses/01_transposable_elements/data"
PROVENANCE = (
    ROOT
    / "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_provenance() -> dict[str, object]:
    return json.loads(PROVENANCE.read_text(encoding="utf-8"))


def test_ltr_release_files_match_provenance() -> None:
    provenance = read_provenance()
    expected_files = {
        "S14": "ltr_element_metrics.csv",
        "S15": "ltr_species_robustness.csv",
        "S16": "ltr_resource_coverage.csv",
        "S17": "ltr_excluded_elements.csv",
    }

    assert set(provenance["canonical_files"]) == set(expected_files)
    for supplement_id, filename in expected_files.items():
        record = provenance["canonical_files"][supplement_id]
        path = ROOT / record["path"]
        frame = pd.read_csv(path, low_memory=False)
        assert path == TE_DATA / filename
        assert sha256(path) == record["sha256"]
        assert [len(frame), len(frame.columns)] == [record["rows"], record["columns"]]

    expected_archive = {
        "corrected_release_builder": (
            (
                "../Desmognathus_TE_archive/scripts/processing/"
                "build_corrected_ectopic_ltr30.py"
            ),
            "8f0a3cf408da25472c035ed4be49190c5dfb7e0945aa937c6bf2badf46d4536b",
        ),
        "corrected_release_manifest": (
            (
                "../Desmognathus_TE_archive/results/data/corrected/"
                "ectopic_ltr30/ectopic_ltr30_v1.manifest.json"
            ),
            "38ed29e570da8b58ae7ac19e58a13a687e6e0e1567ee3b5b0aba07fe745954d1",
        ),
        "historical_screen_table": (
            (
                "../Desmognathus_TE_archive/results/data/"
                "research_review/ltr_mapping_artifact_screen_ltr30_v1.csv"
            ),
            "28fb235f0e9c4a9f757b06dab299eff82d2ff40b64afb9b075a52987034c8536",
        ),
        "historical_screen_rule_source": (
            (
                "../Desmognathus_TE_archive/scripts/processing/"
                "build_audited_historical_style_figures.R"
            ),
            "64713710880b18fd7521e370b733dfbcd82838e9d20051d096bf220bbb5aef4a",
        ),
    }
    for evidence_id, (expected_path, expected_hash) in expected_archive.items():
        record = provenance["archived_evidence"][evidence_id]
        assert record == {"path": expected_path, "sha256": expected_hash}
        archive_path = ROOT / expected_path
        if archive_path.is_file():
            assert sha256(archive_path) == expected_hash

    replay = provenance["active_replay"]
    replay_path = ROOT / replay["path"]
    assert replay["default_mode"] == "read_only_verification"
    assert replay["promotion_flag"] == "--write"
    assert sha256(replay_path) == replay["sha256"]

    availability = provenance["availability_projection"]
    availability_path = ROOT / availability["path"]
    availability_frame = pd.read_csv(availability_path)
    assert sha256(availability_path) == availability["sha256"]
    assert [len(availability_frame), len(availability_frame.columns)] == [
        availability["rows"],
        availability["columns"],
    ]

    original = provenance["original_iqr_method_evidence"]
    original_path = ROOT / original["path"]
    assert original["notebook_cell_index_zero_based"] == 7
    assert original["method"] == "within_species_two_sided_tukey_1.5_iqr"
    if original_path.is_file():
        assert sha256(original_path) == original["sha256"]


def test_ltr_element_formulas_and_positive_depth_gate() -> None:
    provenance = read_provenance()
    elements = pd.read_csv(TE_DATA / "ltr_element_metrics.csv", low_memory=False)

    assert (len(elements), elements["species"].nunique()) == (1086, 30)
    assert elements["domain_count"].value_counts().sort_index().to_dict() == {
        4: 677,
        5: 409,
    }
    assert set(elements["legacy_pipe_split_domain_count"]) == {5, 6}
    assert elements["element length"].ge(3000).all()
    assert elements["tesorter_annotation_found"].all()
    assert elements["source_depth_status"].eq("usable").all()
    assert (
        elements["depth_positions_reported"] + elements["depth_positions_missing"]
        == elements["depth_positions_expected"]
    ).all()
    assert elements["depth_positions_expected"].eq(elements["element length"]).all()

    terminal_weighted_mean = (
        elements["mean_depth_left_ltr_all_positions"] * elements["lLTR length"]
        + elements["mean_depth_right_ltr_all_positions"] * elements["rLTR length"]
    ) / (elements["lLTR length"] + elements["rLTR length"])
    np.testing.assert_allclose(
        elements["mean_depth_terminal_all_positions"], terminal_weighted_mean
    )
    np.testing.assert_allclose(
        elements["ratio_terminal_internal_all_positions"],
        elements["mean_depth_terminal_all_positions"]
        / elements["mean_depth_internal_all_positions"],
    )
    np.testing.assert_allclose(
        elements["log2_ratio_terminal_internal_all_positions"],
        np.log2(elements["ratio_terminal_internal_all_positions"]),
    )

    expected_gate = (
        elements["left_ltr_positive_coverage_fraction"].ge(0.8)
        & elements["right_ltr_positive_coverage_fraction"].ge(0.8)
        & elements["internal_positive_coverage_fraction"].ge(0.8)
    )
    assert expected_gate.equals(elements["coverage_ge_80pct"])
    assert int(expected_gate.sum()) == 1083
    assert int((elements["domain_gate_keep"] & expected_gate).sum()) == 407
    assert elements["reported_position_fraction"].eq(1.0).all()
    assert provenance["coverage_gate"] == {
        "metric": "regional_positive_depth_position_fraction",
        "threshold": 0.8,
        "operator": "greater_than_or_equal",
        "required_regions": ["left_ltr", "right_ltr", "internal"],
        "all_regions_required": True,
        "not_a_reported_position_fraction_gate": True,
    }


def test_ltr_domain_gate_uses_real_tesorter_token_grammar() -> None:
    elements = pd.read_csv(TE_DATA / "ltr_element_metrics.csv", low_memory=False)

    annotation_tokens = elements["Domains"].astype(str).str.split()
    actual_domain_count = annotation_tokens.str.len()
    tokens_have_tesorter_shape = annotation_tokens.map(
        lambda tokens: all(token.count("|") == 1 for token in tokens)
    )
    eligible = (
        elements["Order"].eq("LTR")
        & elements["Superfamily"].eq("Gypsy")
        & actual_domain_count.ge(5)
    )

    assert tokens_have_tesorter_shape.all()
    assert elements["domain_count"].equals(actual_domain_count)
    assert elements["domain_gate_keep"].equals(eligible)
    assert int(eligible.sum()) == 408
    assert elements.loc[eligible, "species"].nunique() == 30
    assert elements.loc[eligible, "domain_count"].eq(5).all()
    assert elements.loc[
        actual_domain_count.ge(5) & ~eligible,
        ["element_id", "Superfamily"],
    ].to_dict(orient="records") == [
        {
            "element_id": "JAUBNV010740999.1_De_636_6770",
            "Superfamily": "Copia",
        }
    ]
    assert read_provenance()["selection"] == {
        "minimum_element_length_bp": 3000,
        "required_order": "LTR",
        "required_superfamily": "Gypsy",
        "minimum_domain_annotations": 5,
        "domain_field_grammar": "whitespace-delimited DOMAIN|MODEL annotations",
        "domain_count_method": (
            "count validated whitespace-delimited annotation tokens"
        ),
        "annotation_join": "exact_sequence_and_coordinate_tesorter_join",
        "legacy_error": (
            "the archived parser split the Domains field on the within-token pipe "
            "delimiter, so four annotations were stored as five and five annotations "
            "were stored as six"
        ),
        "historical_rows_preserved": True,
    }


def test_tesorter_domain_parser_rejects_malformed_fields() -> None:
    module = runpy.run_path(
        str(TE_DATA.parent / "recompute_ltr_release.py"),
        run_name="ltr_release_test",
    )
    count = module["_tesorter_domain_count"]
    five_gypsy = (
        "GAG|Ty3_gypsy PROT|Ty3_gypsy INT|Ty3_gypsy "
        "RT|Ty3_gypsy RH|Ty3_gypsy"
    )

    assert count(five_gypsy) == 5
    assert count(" ".join(five_gypsy.split()[:-1])) == 4
    assert count(f"{five_gypsy} ENV|Ty3_gypsy") == 6
    assert count("none") == 0
    for malformed in (pd.NA, "", "GAG|", "GAG|Ty3_gypsy BROKEN"):
        with pytest.raises(ValueError, match="TEsorter"):
            count(malformed)


def test_ltr_counts_exclusions_and_mapping_boundary() -> None:
    provenance = read_provenance()
    elements = pd.read_csv(TE_DATA / "ltr_element_metrics.csv", low_memory=False)
    species = pd.read_csv(TE_DATA / "ltr_species_robustness.csv")
    resources = pd.read_csv(TE_DATA / "ltr_resource_coverage.csv")
    excluded = pd.read_csv(TE_DATA / "ltr_excluded_elements.csv")
    availability = pd.read_csv(ROOT / "data/identity/analysis_availability.csv")

    assert provenance["release_counts"] == {
        "te_panel_species": 34,
        "species_with_eligible_elements": 30,
        "historical_pipe_split_selected_elements": 1088,
        "historical_source_corrupt_elements": 2,
        "historical_usable_audit_elements": 1086,
        "selected_ltr_gypsy_ge5_domain_elements": 408,
        "selected_source_corrupt_exclusions": 0,
        "iqr_retained_elements": 381,
        "iqr_outliers_removed": 27,
        "iqr_lower_outliers_removed": 1,
        "iqr_upper_outliers_removed": 26,
        "positive_depth_gate_elements": 407,
        "primary_elements": 380,
    }
    assert len(elements) == provenance["release_counts"][
        "historical_usable_audit_elements"
    ]
    assert resources["n_historical_pipe_split_5_6_domain_elements"].sum() == 1088
    assert resources["n_historical_source_corrupt_elements"].sum() == 2
    assert resources["n_selected_ltr_gypsy_ge5_domain_elements"].sum() == 408
    assert resources["n_usable_elements"].sum() == 408
    assert resources["n_selected_source_corrupt_exclusions"].sum() == 0
    assert resources["n_coverage_ge_80pct_elements"].sum() == 407
    assert resources["n_iqr_retained_elements"].sum() == 381
    assert resources["n_primary_iqr_filtered_elements"].sum() == 380
    historical_rows = pd.concat(
        [elements[["species"]], excluded[["species"]]], ignore_index=True
    )["species"].value_counts()
    eligible_rows = elements.loc[elements["domain_gate_keep"], "species"].value_counts()
    assert resources["species"].map(historical_rows).fillna(0).astype(int).equals(
        resources["n_historical_pipe_split_5_6_domain_elements"]
    )
    assert resources["species"].map(eligible_rows).fillna(0).astype(int).equals(
        resources["n_selected_ltr_gypsy_ge5_domain_elements"]
    )
    assert 1088 == 408 + int((~elements["domain_gate_keep"]).sum()) + len(excluded)
    branch_counts = species.groupby("analysis_branch")["n_elements"].sum().to_dict()
    assert branch_counts == {
        "all_ltr_gypsy_ge5_domain_elements": 408,
        "coverage_ge_80pct": 407,
        "coverage_ge_80pct_and_complete_yes": 237,
        "primary_iqr_filtered_coverage_ge_80pct": 380,
        "tesorter_complete_yes": 238,
    }

    assert len(excluded) == 2
    assert excluded["partial_final_line"].all()
    assert ~excluded["source_ends_with_newline"].any()
    assert ~excluded["usable"].any()
    assert ~excluded["input_modified"].any()
    assert excluded["domain_count"].eq(4).all()
    assert excluded["legacy_pipe_split_domain_count"].eq(5).all()
    assert ~excluded["domain_gate_keep"].any()
    assert excluded["exclusion_reason"].eq(
        "truncated_depth_file_partial_final_line"
    ).all()
    observed_exclusions = excluded[
        ["species", "element_id", "depth_file_sha256"]
    ].to_dict(orient="records")
    assert observed_exclusions == provenance["corrupt_exclusions"]

    projected = availability.merge(
        resources,
        on="species",
        how="inner",
        validate="one_to_one",
    )
    assert projected["step04_ltr_resource_status"].equals(
        projected["resource_status"]
    )
    for availability_column, resource_column in {
        "step04_ltr_tabout_elements": "n_tabout_elements",
        "step04_ltr_selected_elements": (
            "n_selected_ltr_gypsy_ge5_domain_elements"
        ),
        "step04_ltr_usable_elements": "n_usable_elements",
        "step04_ltr_source_corrupt_exclusions": (
            "n_selected_source_corrupt_exclusions"
        ),
    }.items():
        assert projected[availability_column].equals(
            projected[resource_column].astype(float)
        )

    mapping = provenance["mapping_workflow"]
    assert mapping["input_reads"] == "paired_trimmed_mitochondrial_filtered_reads"
    assert mapping["reference"] == "per_species_ltr_contigs"
    assert mapping["alignment"] == {
        "tool": "Bowtie2",
        "known_options": ["-L", "20", "--very-sensitive-local"],
    }
    assert mapping["post_alignment"] == {
        "tool": "SAMtools",
        "known_operations": ["BAM_conversion", "sort", "index", "depth"],
    }
    assert set(mapping["unresolved"]) == {
        "read_trimming_and_mitochondrial_filtering_commands_and_versions",
        "bowtie2_version",
        "samtools_version",
        "mapq_threshold",
        "multimapper_policy",
        "secondary_and_supplementary_alignment_policy",
        "duplicate_policy",
    }
    assert provenance["claim_boundary"] == {
        "allowed": "mapping_and_deletion_footprint_proxy",
        "not_direct_measures": [
            "ectopic_recombination_rate",
            "solo_ltr_rate",
            "dna_loss_rate",
        ],
        "release_status": "sensitivity_only",
    }


def test_original_two_sided_iqr_filter_is_primary_and_row_auditable() -> None:
    provenance = read_provenance()
    elements = pd.read_csv(TE_DATA / "ltr_element_metrics.csv", low_memory=False)
    ratio = "ratio_terminal_internal_all_positions"

    eligible = elements["domain_gate_keep"]
    q1 = pd.Series(np.nan, index=elements.index)
    q3 = pd.Series(np.nan, index=elements.index)
    grouped_ratio = elements.loc[eligible].groupby("species")[ratio]
    q1.loc[eligible] = grouped_ratio.transform(lambda values: values.quantile(0.25))
    q3.loc[eligible] = grouped_ratio.transform(lambda values: values.quantile(0.75))
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    below = eligible & elements[ratio].lt(lower)
    above = eligible & elements[ratio].gt(upper)
    keep = eligible & ~(below | above)
    primary = keep & elements["coverage_ge_80pct"]

    np.testing.assert_allclose(elements["species_ratio_q1"], q1)
    np.testing.assert_allclose(elements["species_ratio_q3"], q3)
    np.testing.assert_allclose(elements["species_ratio_iqr"], iqr)
    np.testing.assert_allclose(elements["species_ratio_iqr_lower_bound"], lower)
    np.testing.assert_allclose(elements["species_ratio_iqr_upper_bound"], upper)
    assert elements["iqr_filter_keep"].equals(keep)
    assert elements["primary_ltr_element"].equals(primary)
    assert elements["iqr_filter_rule_id"].eq(
        "within_species_two_sided_tukey_1.5_iqr_v1"
    ).all()
    assert int(eligible.sum()) == 408
    assert int(keep.sum()) == 381
    assert int(below.sum()) == 1
    assert int(above.sum()) == 26
    assert int(primary.sum()) == 380
    assert elements.loc[primary, "species"].nunique() == 30
    assert elements.loc[~eligible, "species_ratio_q1"].isna().all()
    assert elements.loc[~eligible, "iqr_filter_direction"].eq("not_eligible").all()

    indexed = elements.set_index("element_id")
    coverage_failure = indexed.loc["JAWLIA010002421.1_De_459_5558"]
    assert coverage_failure["iqr_filter_keep"]
    assert not coverage_failure["primary_ltr_element"]
    order_sentinel = indexed.loc["JAWLIA010001284.1_De_507_5688"]
    assert order_sentinel["iqr_filter_keep"]
    assert order_sentinel["primary_ltr_element"]

    tilleyi = elements.loc[
        elements["species"].eq("tilleyi") & elements["domain_gate_keep"]
    ].iloc[0]
    assert tilleyi["species_ratio_q1"] == pytest.approx(0.45162764736004235)
    assert tilleyi["species_ratio_q3"] == pytest.approx(0.6993587832941568)
    assert tilleyi["species_ratio_iqr_lower_bound"] == pytest.approx(
        0.08003094345887069
    )
    assert tilleyi["species_ratio_iqr_upper_bound"] == pytest.approx(
        1.0709554871953284
    )

    rule = provenance["primary_iqr_filter"]
    assert rule["role"] == "required_primary_element_filter"
    assert rule["applied_to_primary_release"] is True
    assert rule["lower_bound"] == "Q1 - 1.5 * IQR"
    assert rule["upper_bound"] == "Q3 + 1.5 * IQR"
    assert rule["n_input_elements"] == 408
    assert rule["n_retained_by_iqr"] == 381
    assert rule["n_removed_by_iqr"] == 27
    assert rule["n_primary_elements"] == 380
    assert rule["primary_analysis_branch"] == (
        "primary_iqr_filtered_coverage_ge_80pct"
    )


def test_later_compound_screen_is_explicitly_retired() -> None:
    provenance = read_provenance()
    elements = pd.read_csv(TE_DATA / "ltr_element_metrics.csv", low_memory=False)
    ratio = "ratio_terminal_internal_all_positions"
    grouped_ratio = elements.groupby("species")[ratio]
    q1 = grouped_ratio.transform(lambda values: values.quantile(0.25))
    q3 = grouped_ratio.transform(lambda values: values.quantile(0.75))
    ratio_sum = grouped_ratio.transform("sum")
    flagged = elements.loc[
        (elements[ratio] > q3 + 1.5 * (q3 - q1))
        & elements["left_right_terminal_log2_imbalance"].abs().gt(6)
        & (elements[ratio] / ratio_sum).gt(0.5)
    ]

    sensitivity = provenance["retired_compound_screen_reconstruction"]
    assert sensitivity["rule_id"] == "one_sided_terminal_pileup_v1"
    assert sensitivity["applied_to_primary_release"] is False
    assert sensitivity["role"] == (
        "later_audit_reconstruction_not_the_original_iqr_method"
    )
    assert sensitivity["n_input_elements"] == 1086
    assert sensitivity["n_flagged_elements"] == 1
    assert sensitivity["n_retained_elements"] == 1085
    assert sensitivity["conditions"] == [
        "ratio above the within-species Q3 + 1.5 * IQR upper fence",
        "absolute left:right terminal log2 imbalance > 6 (>64-fold)",
        "element contributes >50% of its species ratio sum",
    ]
    assert flagged[["species", "element_id"]].to_dict(orient="records") == [
        sensitivity["flagged_element"]
    ]
    assert sensitivity["flagged_element"]["element_id"] in set(elements["element_id"])


def test_ltr_release_replay_is_read_only_and_passes() -> None:
    script = ROOT / "analyses/01_transposable_elements/recompute_ltr_release.py"
    completed = subprocess.run(
        [sys.executable, str(script)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert "1,086 historical usable -> 408 LTR/Gypsy >=5-domain eligible" in (
        completed.stdout
    )
    assert "381 IQR-retained; 407 coverage-qualified -> 380 primary" in (
        completed.stdout
    )


def test_ltr_promotion_serialization_is_byte_stable(tmp_path: Path) -> None:
    module = runpy.run_path(
        str(TE_DATA.parent / "recompute_ltr_release.py"),
        run_name="ltr_release_idempotence_test",
    )
    tables = module["build_release"]()
    canonical_paths = (
        TE_DATA / "ltr_element_metrics.csv",
        TE_DATA / "ltr_species_robustness.csv",
        TE_DATA / "ltr_resource_coverage.csv",
        TE_DATA / "ltr_excluded_elements.csv",
        ROOT / "data/identity/analysis_availability.csv",
    )

    for table, canonical in zip(tables, canonical_paths, strict=True):
        rendered = tmp_path / canonical.name
        table.to_csv(rendered, index=False)
        assert rendered.read_bytes() == canonical.read_bytes()

#!/usr/bin/env python3
"""Build a corrected all-position LTR-depth branch for the available TE34 resources.

This branch reuses frozen LTRharvest/TEsorter/depth intermediates. It never maps
reads, discovers LTRs, or modifies the historical 1,088-element table. Two depth
files that are demonstrably truncated mid-line at 245,760 bytes are excluded with
an explicit source audit; every other selected element is recalculated with
unreported and explicit-zero positions retained in regional means.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.processing.build_corrected_ectopic_recombination import (  # noqa: E402
    DOMAIN_COUNTS,
    MIN_ELEMENT_LENGTH,
    _canonical_species,
    _count_domains,
    _element_ids,
    _sha256,
    calculate_zero_aware_depth_metrics,
    merge_elements_with_tesorter,
    summarize_species,
)


INPUT_DIR = PROJECT_ROOT / "input_data" / "ectopic_recombination"
LOOKUP = PROJECT_ROOT / "input_data" / "lookup_table.txt"
PANEL = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "derived"
    / "panels"
    / "study_te_resource_panel34_v1.csv"
)
TESORTER = INPUT_DIR / "combined_sequences.fasta.rexdb-metazoa.cls.tsv"
HISTORICAL_MASTER = PROJECT_ROOT / "results" / "data" / "ectopic_recombination_master.csv"
HISTORICAL_LTR30 = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"
)
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "corrected" / "ectopic_ltr30"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "corrected" / "ectopic_ltr30"
ELEMENT_OUTPUT = OUTPUT_DIR / "ectopic_element_metrics_ltr30_v1.csv"
SPECIES_OUTPUT = OUTPUT_DIR / "ectopic_species_robustness_ltr30_v1.csv"
EXCLUDED_OUTPUT = OUTPUT_DIR / "ectopic_excluded_elements_ltr30_v1.csv"
RESOURCE_OUTPUT = OUTPUT_DIR / "ectopic_resource_coverage_ltr30_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "ectopic_ltr30_v1.manifest.json"
QC_PNG = FIGURE_DIR / "ectopic_ltr30_coverage_qc_v1.png"
QC_PDF = FIGURE_DIR / "ectopic_ltr30_coverage_qc_v1.pdf"

EXPECTED_CORRUPT_DEPTH_ELEMENT_IDS = {
    "JAUEJG010675316.1_De_4598_10648",
    "JASANM010280042.1_De_179_6669",
}


def portable(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def inspect_depth_source(
    path: Path,
    element_id: str,
    element_length: int,
) -> dict[str, Any]:
    """Audit a depth source without repairing, truncating, or imputing it."""

    expected = int(element_length)
    audit: dict[str, Any] = {
        "element_id": element_id,
        "depth_file": portable(path) if path.is_absolute() and PROJECT_ROOT in path.parents else str(path),
        "depth_file_exists": path.is_file(),
        "depth_file_size_bytes": path.stat().st_size if path.is_file() else 0,
        "depth_file_sha256": _sha256(path) if path.is_file() else "",
        "expected_terminal_position": expected,
        "reported_terminal_position": 0,
        "complete_depth_rows": 0,
        "missing_terminal_suffix_positions": expected,
        "reported_position_fraction": 0.0,
        "partial_final_line": False,
        "source_ends_with_newline": False,
        "usable": False,
        "exclusion_reason": "missing_depth_file",
        "input_modified": False,
    }
    if not path.is_file():
        return audit
    payload = path.read_bytes()
    audit["source_ends_with_newline"] = payload.endswith(b"\n")
    if not payload:
        audit["exclusion_reason"] = "empty_depth_file"
        return audit

    lines = payload.splitlines()
    positions: list[int] = []
    seen_positions: set[int] = set()
    malformed_reason = ""
    for line_number, raw_line in enumerate(lines, start=1):
        fields = raw_line.split(b"\t")
        if len(fields) != 3:
            is_partial_tail = line_number == len(lines) and not payload.endswith(b"\n")
            audit["partial_final_line"] = is_partial_tail
            malformed_reason = (
                "truncated_depth_file_partial_final_line"
                if is_partial_tail
                else "malformed_depth_row"
            )
            break
        try:
            observed_id = fields[0].decode("ascii")
            position = int(fields[1])
            depth = float(fields[2])
        except (UnicodeDecodeError, ValueError):
            is_partial_tail = line_number == len(lines) and not payload.endswith(b"\n")
            audit["partial_final_line"] = is_partial_tail
            malformed_reason = (
                "truncated_depth_file_partial_final_line"
                if is_partial_tail
                else "non_numeric_or_non_ascii_depth_row"
            )
            break
        if observed_id != element_id:
            malformed_reason = "depth_sequence_identity_mismatch"
            break
        if position < 1 or position > expected or position in seen_positions:
            malformed_reason = "invalid_or_duplicate_depth_position"
            break
        if not np.isfinite(depth) or depth < 0:
            malformed_reason = "invalid_depth_value"
            break
        positions.append(position)
        seen_positions.add(position)

    terminal = max(positions, default=0)
    audit["reported_terminal_position"] = terminal
    audit["complete_depth_rows"] = len(positions)
    audit["missing_terminal_suffix_positions"] = max(expected - terminal, 0)
    audit["reported_position_fraction"] = len(positions) / expected if expected else 0.0
    if malformed_reason:
        audit["exclusion_reason"] = malformed_reason
        return audit

    audit["usable"] = True
    audit["exclusion_reason"] = ""
    return audit


def load_ltr30_candidates() -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    """Load current TE34 resources and exact-join 5/6-domain LTR candidates."""

    panel = pd.read_csv(PANEL)
    required_panel = {"species", "te_sra_accession", "te_assembly_accession"}
    missing_panel = sorted(required_panel.difference(panel.columns))
    if missing_panel:
        raise ValueError(f"TE34 panel lacks columns: {missing_panel}")
    panel["species"] = panel["species"].map(_canonical_species)
    if len(panel) != 34 or panel["species"].duplicated().any():
        raise ValueError("TE34 panel must contain 34 unique species")

    lookup = pd.read_csv(LOOKUP, sep="\t")
    lookup["species"] = lookup["Species"].map(_canonical_species)
    lookup = lookup.loc[lookup["species"].isin(set(panel["species"]))].copy()
    if len(lookup) != 34 or lookup["species"].duplicated().any():
        raise ValueError("Active lookup must map TE34 one-to-one")
    crosswalk = panel[
        ["species", "te_sra_accession", "te_assembly_accession"]
    ].merge(
        lookup[["species", "SRA_Accension", "Genome_Accension"]],
        on="species",
        validate="one_to_one",
    )
    if not crosswalk["te_sra_accession"].eq(crosswalk["SRA_Accension"]).all():
        raise ValueError("Panel and lookup disagree on SRA accessions")
    if not crosswalk["te_assembly_accession"].eq(crosswalk["Genome_Accension"]).all():
        raise ValueError("Panel and lookup disagree on assembly accessions")

    frames: list[pd.DataFrame] = []
    raw_paths: list[Path] = []
    resource_rows: list[dict[str, Any]] = []
    numeric_columns = [
        "element start",
        "element end",
        "element length",
        "lLTR length",
        "rLTR length",
    ]
    for row in panel.sort_values("species").itertuples(index=False):
        path = INPUT_DIR / f"{row.te_assembly_accession}_tabout.csv"
        resource_record: dict[str, Any] = {
            "species": row.species,
            "te_sra_accession": row.te_sra_accession,
            "te_assembly_accession": row.te_assembly_accession,
            "tabout_path": portable(path),
            "tabout_exists": path.is_file(),
            "tabout_sha256": _sha256(path) if path.is_file() else "",
            "n_tabout_elements": 0,
            "n_elements_ge_3000bp": 0,
        }
        if not path.is_file():
            resource_rows.append(resource_record)
            continue
        frame = pd.read_csv(path, sep="\t", index_col=False)
        missing_columns = sorted(set(numeric_columns + ["sequence"]).difference(frame.columns))
        if missing_columns:
            raise ValueError(f"{path} lacks required columns: {missing_columns}")
        converted = frame[numeric_columns].apply(pd.to_numeric, errors="coerce")
        if converted.isna().any().any():
            raise ValueError(f"{path} has non-numeric required LTR coordinates")
        frame[numeric_columns] = converted.astype(int)
        resource_record["n_tabout_elements"] = len(frame)
        resource_record["n_elements_ge_3000bp"] = int(
            frame["element length"].ge(MIN_ELEMENT_LENGTH).sum()
        )
        frame = frame.loc[frame["element length"].ge(MIN_ELEMENT_LENGTH)].copy()
        frame["species"] = row.species
        frame["te_sra_accession"] = row.te_sra_accession
        frame["te_assembly_accession"] = row.te_assembly_accession
        frames.append(frame)
        raw_paths.append(path)
        resource_rows.append(resource_record)

    if not frames:
        raise RuntimeError("No TE34 tabout resources are available")
    candidates = pd.concat(frames, ignore_index=True)
    tesorter = pd.read_csv(TESORTER, sep="\t")
    merged = merge_elements_with_tesorter(candidates, tesorter)
    merged["domain_count"] = merged["Domains"].apply(_count_domains)
    selected = merged.loc[merged["domain_count"].isin(DOMAIN_COUNTS)].copy()
    if len(selected) != 1088 or selected["species"].nunique() != 30:
        raise ValueError(
            "Expected 1,088 exact 5/6-domain elements across 30 current species"
        )

    resources = pd.DataFrame(resource_rows)
    annotation_counts = (
        merged.groupby("species")["tesorter_annotation_found"].sum().astype(int)
    )
    selected_counts = selected.groupby("species").size()
    resources["n_exact_tesorter_annotations"] = (
        resources["species"].map(annotation_counts).fillna(0).astype(int)
    )
    resources["n_selected_5_6_domain_elements"] = (
        resources["species"].map(selected_counts).fillna(0).astype(int)
    )
    return selected, resources, raw_paths


def add_zero_aware_metrics_and_audit(
    selected: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate metrics for usable files and fail on unapproved exclusions."""

    usable_rows: list[pd.Series] = []
    metric_rows: list[dict[str, Any]] = []
    excluded_rows: list[dict[str, Any]] = []
    for _, row in selected.sort_values(["species", "element_id"]).iterrows():
        path = INPUT_DIR / f"{row['element_id']}.fa.depth.txt"
        audit = inspect_depth_source(path, row["element_id"], row["element length"])
        if not audit["usable"]:
            excluded_rows.append(
                {
                    "species": row["species"],
                    "te_sra_accession": row["te_sra_accession"],
                    "te_assembly_accession": row["te_assembly_accession"],
                    "element_id": row["element_id"],
                    "sequence": row["sequence"],
                    "element_start": int(row["element start"]),
                    "element_end": int(row["element end"]),
                    "element_length": int(row["element length"]),
                    "domain_count": int(row["domain_count"]),
                    **{key: value for key, value in audit.items() if key != "element_id"},
                }
            )
            continue
        depth = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["depth_sequence_id", "position", "depth"],
        )
        try:
            metric = calculate_zero_aware_depth_metrics(
                depth,
                row["element length"],
                row["lLTR length"],
                row["rLTR length"],
            )
        except Exception as error:  # fail-explicit audit below names any new case
            audit["usable"] = False
            audit["exclusion_reason"] = f"depth_metric_validation_failed:{type(error).__name__}"
            excluded_rows.append(
                {
                    "species": row["species"],
                    "te_sra_accession": row["te_sra_accession"],
                    "te_assembly_accession": row["te_assembly_accession"],
                    "element_id": row["element_id"],
                    "sequence": row["sequence"],
                    "element_start": int(row["element start"]),
                    "element_end": int(row["element end"]),
                    "element_length": int(row["element length"]),
                    "domain_count": int(row["domain_count"]),
                    **{key: value for key, value in audit.items() if key != "element_id"},
                }
            )
            continue
        metric["depth_file"] = portable(path)
        metric["depth_file_sha256"] = audit["depth_file_sha256"]
        metric["depth_file_size_bytes"] = audit["depth_file_size_bytes"]
        metric["source_depth_status"] = "usable"
        metric["reported_position_fraction"] = audit["reported_position_fraction"]
        usable_rows.append(row)
        metric_rows.append(metric)

    excluded = pd.DataFrame(excluded_rows).sort_values("element_id").reset_index(drop=True)
    observed_exclusions = set(excluded["element_id"]) if not excluded.empty else set()
    if observed_exclusions != EXPECTED_CORRUPT_DEPTH_ELEMENT_IDS:
        raise RuntimeError(
            "LTR30 source exclusions changed. Expected "
            f"{sorted(EXPECTED_CORRUPT_DEPTH_ELEMENT_IDS)}, observed "
            f"{sorted(observed_exclusions)}"
        )
    if not excluded["exclusion_reason"].eq(
        "truncated_depth_file_partial_final_line"
    ).all():
        raise RuntimeError("Approved LTR30 exclusions no longer have the expected truncation signature")

    usable = pd.DataFrame(usable_rows).reset_index(drop=True)
    metrics = pd.DataFrame(metric_rows).reset_index(drop=True)
    elements = pd.concat([usable, metrics], axis=1)
    elements["coverage_ge_80pct"] = (
        elements["left_ltr_positive_coverage_fraction"].ge(0.8)
        & elements["right_ltr_positive_coverage_fraction"].ge(0.8)
        & elements["internal_positive_coverage_fraction"].ge(0.8)
    )
    if len(elements) != 1086 or elements["species"].nunique() != 30:
        raise RuntimeError("Corrected LTR30 must retain 1,086 usable elements across 30 species")
    return elements.sort_values(["species", "element_id"]).reset_index(drop=True), excluded


def add_historical_comparison(elements: pd.DataFrame) -> pd.DataFrame:
    historical = pd.read_csv(HISTORICAL_MASTER, sep="\t")
    historical["element_id"] = _element_ids(historical)
    historical = historical.loc[
        historical["domain_count"].isin(DOMAIN_COUNTS),
        ["element_id", "ratio_terminal_internal"],
    ].copy()
    if historical["element_id"].duplicated().any():
        raise ValueError("Historical selected branch has duplicate exact element IDs")
    historical = historical.rename(
        columns={"ratio_terminal_internal": "historical_nonzero_only_ratio"}
    )
    result = elements.merge(historical, on="element_id", how="left", validate="one_to_one")
    if result["historical_nonzero_only_ratio"].isna().any():
        raise ValueError("Corrected LTR30 element is absent from the preserved historical table")
    result["nonzero_ratio_reproduction_abs_error"] = (
        result["ratio_terminal_internal_nonzero_only"]
        - result["historical_nonzero_only_ratio"]
    ).abs()
    if result["nonzero_ratio_reproduction_abs_error"].max() >= 1e-10:
        raise ValueError("Nonzero-only reconstruction no longer reproduces historical ratios")
    return result


def finalize_resource_audit(
    resources: pd.DataFrame,
    elements: pd.DataFrame,
    excluded: pd.DataFrame,
) -> pd.DataFrame:
    usable_counts = elements.groupby("species").size()
    excluded_counts = excluded.groupby("species").size()
    result = resources.copy()
    result["n_usable_elements"] = result["species"].map(usable_counts).fillna(0).astype(int)
    result["n_source_corrupt_exclusions"] = (
        result["species"].map(excluded_counts).fillna(0).astype(int)
    )
    result["resource_status"] = np.select(
        [
            ~result["tabout_exists"],
            result["n_selected_5_6_domain_elements"].eq(0),
            result["n_source_corrupt_exclusions"].gt(0),
        ],
        [
            "missing_tabout",
            "no_5_6_domain_elements",
            "usable_with_source_corrupt_exclusions",
        ],
        default="usable",
    )
    return result.sort_values("species").reset_index(drop=True)


def write_qc_figure(elements: pd.DataFrame, excluded: pd.DataFrame) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    axes[0].scatter(
        elements["internal_positive_coverage_fraction"],
        elements["terminal_positive_coverage_fraction"],
        s=14,
        alpha=0.35,
        color="#2563eb",
        edgecolors="none",
    )
    axes[0].axhline(0.8, color="#b91c1c", linestyle="--", linewidth=1)
    axes[0].axvline(0.8, color="#b91c1c", linestyle="--", linewidth=1)
    axes[0].set(
        xlim=(0.75, 1.005),
        ylim=(0.75, 1.005),
        xlabel="Internal positive-position fraction",
        ylabel="Combined terminal positive-position fraction",
        title="Regional positive-position coverage",
    )

    axes[1].scatter(
        elements["depth_positions_expected"],
        elements["reported_position_fraction"],
        s=14,
        alpha=0.3,
        color="#4b5563",
        edgecolors="none",
        label="1,086 usable elements",
    )
    axes[1].scatter(
        excluded["expected_terminal_position"],
        excluded["reported_position_fraction"],
        marker="x",
        s=90,
        linewidths=2,
        color="#dc2626",
        label="2 truncated exclusions",
        zorder=3,
    )
    annotation_offsets = {
        "aureatus": (-8, -32, "right"),
        "santeetlah": (8, 12, "left"),
    }
    for row in excluded.itertuples(index=False):
        x_offset, y_offset, horizontal_alignment = annotation_offsets[row.species]
        axes[1].annotate(
            f"{row.species}\n{row.reported_terminal_position}/{row.expected_terminal_position}",
            (row.expected_terminal_position, row.reported_position_fraction),
            xytext=(x_offset, y_offset),
            textcoords="offset points",
            fontsize=8,
            color="#991b1b",
            ha=horizontal_alignment,
            arrowprops={"arrowstyle": "-", "color": "#991b1b", "lw": 0.6},
        )
    axes[1].axhline(1, color="#111827", linestyle="--", linewidth=1)
    axes[1].set(
        ylim=(0.88, 1.01),
        xlabel="Expected element positions",
        ylabel="Reported rows / expected positions",
        title="Depth-source completeness",
    )
    axes[1].legend(frameon=False, loc="lower left")

    fig.suptitle(
        "Corrected LTR30 coverage audit: 1,086 usable elements / 30 species",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.01,
        "Two 245,760-byte files ended mid-line and were excluded; historical 1,088-row input is unchanged.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(QC_PNG, dpi=300, facecolor="white")
    fig.savefig(QC_PDF, facecolor="white")
    plt.close(fig)
    return [QC_PNG, QC_PDF]


def main() -> None:
    all_outputs = [
        ELEMENT_OUTPUT,
        SPECIES_OUTPUT,
        EXCLUDED_OUTPUT,
        RESOURCE_OUTPUT,
        MANIFEST_OUTPUT,
        QC_PNG,
        QC_PDF,
    ]
    occupied = [path for path in all_outputs if path.exists()]
    if occupied:
        raise FileExistsError(
            "Non-destructive corrected LTR30 build requires unused outputs: "
            + ", ".join(str(path) for path in occupied)
        )

    selected, resources, tabout_paths = load_ltr30_candidates()
    elements, excluded = add_zero_aware_metrics_and_audit(selected)
    elements = add_historical_comparison(elements)
    species = summarize_species(elements)
    resources = finalize_resource_audit(resources, elements, excluded)

    ratio = elements["ratio_terminal_internal_all_positions"]
    if not (np.isfinite(ratio).all() and ratio.gt(0).all()):
        raise ValueError("Corrected LTR30 contains a non-finite or non-positive ratio")
    if resources["n_selected_5_6_domain_elements"].sum() != 1088:
        raise ValueError("Resource audit lost selected LTR30 elements")
    if resources["n_usable_elements"].sum() != 1086:
        raise ValueError("Resource audit lost usable LTR30 elements")

    historical_hash_before = _sha256(HISTORICAL_LTR30)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    elements.to_csv(ELEMENT_OUTPUT, index=False)
    species.to_csv(SPECIES_OUTPUT, index=False)
    excluded.to_csv(EXCLUDED_OUTPUT, index=False)
    resources.to_csv(RESOURCE_OUTPUT, index=False)
    figure_paths = write_qc_figure(elements, excluded)
    historical_hash_after = _sha256(HISTORICAL_LTR30)
    if historical_hash_before != historical_hash_after:
        raise RuntimeError("Historical LTR30 source changed during corrected build")

    data_outputs = {
        "element_metrics": {
            "path": portable(ELEMENT_OUTPUT),
            "sha256": _sha256(ELEMENT_OUTPUT),
            "rows": len(elements),
        },
        "species_robustness": {
            "path": portable(SPECIES_OUTPUT),
            "sha256": _sha256(SPECIES_OUTPUT),
            "rows": len(species),
        },
        "excluded_elements": {
            "path": portable(EXCLUDED_OUTPUT),
            "sha256": _sha256(EXCLUDED_OUTPUT),
            "rows": len(excluded),
        },
        "resource_coverage": {
            "path": portable(RESOURCE_OUTPUT),
            "sha256": _sha256(RESOURCE_OUTPUT),
            "rows": len(resources),
        },
    }
    manifest = {
        "analysis_id": "ectopic_ltr30_v1",
        "analysis_scope": "audited_te_resource_panel34_ltr30",
        "analysis_status": "sensitivity_only_not_validated_ectopic_rate",
        "eligible_for_confirmatory_path_analysis": False,
        "n_panel_species": 34,
        "n_tabout_resources": int(resources["tabout_exists"].sum()),
        "n_species_with_selected_elements": int(selected["species"].nunique()),
        "n_species_with_usable_elements": int(elements["species"].nunique()),
        "n_selected_elements": len(selected),
        "n_source_corrupt_exclusions": len(excluded),
        "n_usable_elements": len(elements),
        "n_elements_coverage_ge_80pct": int(elements["coverage_ge_80pct"].sum()),
        "exact_tesorter_join": True,
        "zero_depth_positions_retained": True,
        "expensive_upstream_tools_executed": False,
        "historical_ltr30_source_unchanged": historical_hash_before == historical_hash_after,
        "historical_ltr30_source": portable(HISTORICAL_LTR30),
        "historical_ltr30_source_sha256": historical_hash_after,
        "source_corrupt_exclusion_policy": (
            "exclude only depth files proven truncated by a partial final line; "
            "fail on any unapproved exclusion"
        ),
        "source_corrupt_exclusions": excluded[
            [
                "species",
                "element_id",
                "depth_file",
                "depth_file_sha256",
                "depth_file_size_bytes",
                "expected_terminal_position",
                "reported_terminal_position",
                "missing_terminal_suffix_positions",
                "exclusion_reason",
            ]
        ].to_dict(orient="records"),
        "missing_tabout_resources": resources.loc[
            resources["resource_status"].eq("missing_tabout"),
            ["species", "te_sra_accession", "te_assembly_accession", "tabout_path"],
        ].to_dict(orient="records"),
        "resources_without_selected_elements": resources.loc[
            resources["resource_status"].eq("no_5_6_domain_elements"), "species"
        ].tolist(),
        "panel": portable(PANEL),
        "panel_sha256": _sha256(PANEL),
        "lookup": portable(LOOKUP),
        "lookup_sha256": _sha256(LOOKUP),
        "tesorter": portable(TESORTER),
        "tesorter_sha256": _sha256(TESORTER),
        "historical_master": portable(HISTORICAL_MASTER),
        "historical_master_sha256": _sha256(HISTORICAL_MASTER),
        "tabout_inputs": [
            {"path": portable(path), "sha256": _sha256(path)} for path in tabout_paths
        ],
        "selected_depth_file_count": len(selected),
        "usable_depth_file_hashes_embedded_in_element_output": True,
        "excluded_depth_file_hashes_embedded_in_exclusion_output": True,
        "data_outputs": data_outputs,
        "figures": [
            {"path": portable(path), "sha256": _sha256(path)} for path in figure_paths
        ],
        "bootstrap_replicates": 2000,
        "coverage_sensitivity_threshold": 0.8,
        "blocking_provenance": [
            "read_mapping_command_and_reference",
            "multimapper_policy",
            "mapq_threshold",
            "secondary_and_supplementary_alignment_policy",
            "duplicate_policy",
            "structural_solo_ltr_validation",
        ],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {ELEMENT_OUTPUT}")
    print(f"Wrote {SPECIES_OUTPUT}")
    print(f"Wrote {EXCLUDED_OUTPUT}")
    print(f"Wrote {RESOURCE_OUTPUT}")
    print(f"Wrote {MANIFEST_OUTPUT}")
    for path in figure_paths:
        print(f"Wrote {path}")
    print("Corrected LTR30: 1,086 usable elements / 30 species; 2 source-corrupt exclusions")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summarize publication-facing readiness for the merged path-analysis workflow."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PATH_ANALYSIS_ROOT = PROJECT_ROOT / "path_analysis"
DERIVED_DIR = PATH_ANALYSIS_ROOT / "data" / "derived"
EXTERNAL_DERIVED_DIR = PATH_ANALYSIS_ROOT / "data" / "external" / "derived"
RESULTS_DIR = PATH_ANALYSIS_ROOT / "results"

SOURCE_GAPS_PATH = DERIVED_DIR / "source_traceability_gaps.csv"
CELLPROFILER_BUNDLE_PATH = EXTERNAL_DERIVED_DIR / "cellprofiler_final_species_results.csv"
CELLPROFILER_STATE_SUMMARY_PATH = EXTERNAL_DERIVED_DIR / "cellprofiler_genome_state_summary.csv"
CELLPROFILER_SOURCE_DISCOVERY_PATH = EXTERNAL_DERIVED_DIR / "cellprofiler_source_discovery.json"
CELLPROFILER_TRACE_GAPS_PATH = EXTERNAL_DERIVED_DIR / "cellprofiler_traceability_audit_gaps.csv"
CELLPROFILER_TRACE_SUMMARY_PATH = EXTERNAL_DERIVED_DIR / "cellprofiler_traceability_audit_summary.csv"
PANEL_SUMMARY_PATH = DERIVED_DIR / "analysis_panel_summary.csv"
PAPER_FREEZE_MANIFEST_PATH = PROJECT_ROOT / "PAPER_FREEZE_MANIFEST.md"

OUTPUT_CHECKS = DERIVED_DIR / "publication_readiness_checks.csv"
OUTPUT_SUMMARY = DERIVED_DIR / "publication_readiness_summary.json"

PRIMARY_MODEL_FILES = [
    "te_genome_primary_mediumplus_model_ranking.csv",
    "te_genome_primary_strict_body_model_ranking.csv",
    "te_genome_ectopic_primary_mediumplus_model_ranking.csv",
    "te_genome_ectopic_primary_strict_body_model_ranking.csv",
    "te_genome_organismal_primary_mediumplus_model_ranking.csv",
    "te_genome_organismal_primary_strict_body_model_ranking.csv",
    "te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv",
    "te_genome_ectopic_organismal_primary_strict_body_model_ranking.csv",
]

WARNING_GAP_TYPES = {
    "missing_object_mask_trace",
    "missing_raw_imagej_results",
    "missing_roi_zip_trace",
    "missing_source_image_files",
    "metadata_species_missing_from_run",
    "upstream_bundle_missing_reconstructed_from_linked_runs",
}

FAIL_GAP_TYPES = {
    "missing_linked_genome_bundle",
    "missing_linked_genome_trace",
    "empty_linked_genome_trace",
    "missing_linked_genome_state_summary",
    "missing_linked_reference_state",
}


def add_check(rows: list[dict[str, object]], check_id: str, status: str, details: str) -> None:
    rows.append({"check_id": check_id, "status": status, "details": details})


def read_csv_if_present(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def read_json_if_present(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def summarize_primary_winners() -> list[dict[str, object]]:
    winners: list[dict[str, object]] = []
    for filename in PRIMARY_MODEL_FILES:
        path = RESULTS_DIR / filename
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if df.empty:
            continue
        first = df.iloc[0]
        winners.append(
            {
                "file": filename,
                "winner": first.get("model"),
                "weight": first.get("w"),
                "delta_CICc": first.get("delta_CICc"),
            }
        )
    return winners


def main() -> None:
    check_rows: list[dict[str, object]] = []

    source_gaps = read_csv_if_present(SOURCE_GAPS_PATH)
    if not SOURCE_GAPS_PATH.exists():
        add_check(check_rows, "source_traceability_gaps_file", "fail", "Missing source_traceability_gaps.csv")
    elif source_gaps.empty:
        add_check(check_rows, "source_traceability_gaps_empty", "pass", "No source-traceability gaps are currently recorded.")
    else:
        add_check(
            check_rows,
            "source_traceability_gaps_empty",
            "fail",
            f"Found {len(source_gaps)} source-traceability gaps in {SOURCE_GAPS_PATH.name}.",
        )

    bundle = read_csv_if_present(CELLPROFILER_BUNDLE_PATH)
    if not CELLPROFILER_BUNDLE_PATH.exists():
        add_check(check_rows, "cellprofiler_bundle_present", "fail", "Missing CellProfiler genome bundle snapshot.")
    else:
        origins = sorted(set(bundle.get("bundle_origin", pd.Series(dtype="object")).dropna().astype(str)))
        if "recovered_from_path_analysis_tables" in origins:
            add_check(
                check_rows,
                "cellprofiler_bundle_provenance",
                "fail",
                "CellProfiler genome bundle still falls back to downstream path-analysis recovery.",
            )
        elif "reconstructed_from_upstream_raw_runs" in origins:
            add_check(
                check_rows,
                "cellprofiler_bundle_provenance",
                "fail",
                "CellProfiler genome bundle is still reconstructed from standalone raw nucleus-IOD runs rather than linked YOLO nuclei.",
            )
        elif "reconstructed_from_upstream_linked_runs" in origins:
            add_check(
                check_rows,
                "cellprofiler_bundle_provenance",
                "pass",
                "CellProfiler genome bundle is reconstructed from linked YOLO nucleus measurements after cell linkage.",
            )
        else:
            add_check(
                check_rows,
                "cellprofiler_bundle_provenance",
                "pass",
                f"CellProfiler genome bundle imports directly from upstream origins: {origins}",
            )

        status_counts = bundle.get("result_status", pd.Series(dtype="object")).value_counts(dropna=False).to_dict()
        add_check(
            check_rows,
            "cellprofiler_bundle_status_distribution",
            "warn" if status_counts.get("sensitivity_limited", 0) else "pass",
            f"Genome status counts: {status_counts}",
        )

    if CELLPROFILER_STATE_SUMMARY_PATH.exists():
        state_summary = pd.read_csv(CELLPROFILER_STATE_SUMMARY_PATH)
        add_check(
            check_rows,
            "cellprofiler_state_summary_present",
            "pass",
            f"Genome state summary is present with {len(state_summary)} species/state rows.",
        )
    else:
        state_summary = pd.DataFrame()
        add_check(check_rows, "cellprofiler_state_summary_present", "fail", "Missing cellprofiler_genome_state_summary.csv")

    source_discovery = read_json_if_present(CELLPROFILER_SOURCE_DISCOVERY_PATH)
    if CELLPROFILER_SOURCE_DISCOVERY_PATH.exists():
        active_run_tags = source_discovery.get("active_run_tags", {})
        add_check(
            check_rows,
            "cellprofiler_source_discovery_present",
            "pass",
            f"CellProfiler source discovery is present with active_run_tags={active_run_tags}",
        )
    else:
        add_check(
            check_rows,
            "cellprofiler_source_discovery_present",
            "warn",
            "Missing cellprofiler_source_discovery.json",
        )

    trace_gaps = read_csv_if_present(CELLPROFILER_TRACE_GAPS_PATH)
    if not CELLPROFILER_TRACE_GAPS_PATH.exists():
        add_check(check_rows, "cellprofiler_trace_gaps_file", "fail", "Missing CellProfiler traceability gap report.")
    elif trace_gaps.empty:
        add_check(check_rows, "cellprofiler_trace_gaps", "pass", "No CellProfiler traceability gaps are currently recorded.")
    else:
        for gap_type, gap_df in trace_gaps.groupby("gap_type"):
            if gap_type in FAIL_GAP_TYPES:
                status = "fail"
            elif gap_type in WARNING_GAP_TYPES:
                status = "warn"
            else:
                status = "warn"
            species = sorted(gap_df["species"].dropna().astype(str).unique().tolist())
            details = f"{gap_type}: {len(gap_df)} rows"
            if species:
                details += f"; species={species}"
            add_check(check_rows, f"cellprofiler_gap_{gap_type}", status, details)

    missing_model_files = [name for name in PRIMARY_MODEL_FILES if not (RESULTS_DIR / name).exists()]
    if missing_model_files:
        add_check(
            check_rows,
            "primary_model_rankings_present",
            "fail",
            f"Missing primary model ranking files: {missing_model_files}",
        )
    else:
        add_check(
            check_rows,
            "primary_model_rankings_present",
            "pass",
            f"All {len(PRIMARY_MODEL_FILES)} primary model ranking files are present.",
        )

    panel_summary = read_csv_if_present(PANEL_SUMMARY_PATH)
    if PANEL_SUMMARY_PATH.exists() and not panel_summary.empty:
        primary_panels = panel_summary[panel_summary["panel_name"].str.contains("primary", na=False)]
        add_check(
            check_rows,
            "analysis_panel_summary_present",
            "pass",
            f"Panel summary is present with {len(primary_panels)} primary/sensitivity panel rows.",
        )
    else:
        add_check(check_rows, "analysis_panel_summary_present", "fail", "Missing or empty analysis_panel_summary.csv")

    add_check(
        check_rows,
        "paper_freeze_manifest_present",
        "pass" if PAPER_FREEZE_MANIFEST_PATH.exists() else "warn",
        "Paper freeze manifest is present." if PAPER_FREEZE_MANIFEST_PATH.exists() else "Paper freeze manifest is missing.",
    )

    checks = pd.DataFrame(check_rows)
    checks.to_csv(OUTPUT_CHECKS, index=False)

    status_priority = {"fail": 2, "warn": 1, "pass": 0}
    overall_status = "pass"
    if not checks.empty:
        max_priority = max(status_priority.get(status, 1) for status in checks["status"])
        overall_status = {0: "pass", 1: "warn", 2: "fail"}[max_priority]

    summary = {
        "overall_status": overall_status,
        "generated_from": str(Path(__file__).resolve()),
        "checks_path": str(OUTPUT_CHECKS.resolve()),
        "source_traceability_gap_count": int(len(source_gaps)),
        "cellprofiler_trace_gap_count": int(len(trace_gaps)),
        "cellprofiler_source_discovery_path": (
            str(CELLPROFILER_SOURCE_DISCOVERY_PATH.resolve()) if CELLPROFILER_SOURCE_DISCOVERY_PATH.exists() else None
        ),
        "cellprofiler_bundle_path": str(CELLPROFILER_BUNDLE_PATH.resolve()) if CELLPROFILER_BUNDLE_PATH.exists() else None,
        "cellprofiler_state_summary_path": (
            str(CELLPROFILER_STATE_SUMMARY_PATH.resolve()) if CELLPROFILER_STATE_SUMMARY_PATH.exists() else None
        ),
        "active_run_tags": source_discovery.get("active_run_tags", {}),
        "primary_model_winners": summarize_primary_winners(),
    }
    OUTPUT_SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"Wrote {OUTPUT_CHECKS}")
    print(f"Wrote {OUTPUT_SUMMARY}")
    print(f"Overall status: {overall_status}")


if __name__ == "__main__":
    main()

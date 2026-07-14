#!/usr/bin/env python3
"""Rebuild canonical TE diversity tables into scratch outputs.

This script is intentionally non-destructive. It reconstructs candidate
`diversity_order_stats.csv` and `diversity_superfamily_stats.csv` tables from
the repo-local TE breakdown tables plus the threshold `0.0` metrics written
by `scripts/processing/diversity_stats.py`, then audits the candidates against
the current `results/data/` snapshots.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DATA = PROJECT_ROOT / "results" / "data"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"

TOLERANCE = 1e-12


@dataclass(frozen=True)
class DiversityTableConfig:
    dataset: str
    breakdown_path: Path
    comparison_path: Path
    canonical_path: Path
    candidate_path: Path


CONFIGS = [
    DiversityTableConfig(
        dataset="order",
        breakdown_path=RESULTS_DATA / "dnaPipeTE_order_breakdown.csv",
        comparison_path=RESULTS_DATA / "comparison_diversity_order_stats_granular_0_5pct.csv",
        canonical_path=RESULTS_DATA / "diversity_order_stats.csv",
        candidate_path=DERIVED_DIR / "diversity_order_stats_candidate.csv",
    ),
    DiversityTableConfig(
        dataset="superfamily",
        breakdown_path=RESULTS_DATA / "dnaPipeTE_superfamily_breakdown.csv",
        comparison_path=RESULTS_DATA / "comparison_diversity_superfamily_stats_granular_0_5pct.csv",
        canonical_path=RESULTS_DATA / "diversity_superfamily_stats.csv",
        candidate_path=DERIVED_DIR / "diversity_superfamily_stats_candidate.csv",
    ),
]

AUDIT_PATH = DERIVED_DIR / "diversity_canonicalization_audit.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def calculate_simpson_diversity(row: pd.Series) -> float:
    """Return standard Gini-Simpson diversity, 1 - sum(p_i^2).

    The source rows are relative abundances (percentages), not integer counts,
    so a finite-count correction is neither defined nor scale invariant here.
    """

    row = pd.to_numeric(row, errors="coerce")
    row = row[np.isfinite(row) & (row > 0)]
    total = row.sum()
    if row.empty or total <= 0 or pd.isna(total):
        return 0.0
    proportions = row / total
    return float(1.0 - np.sum(proportions * proportions))


def calculate_shannon_diversity(row: pd.Series) -> float:
    row = row[row > 0]
    if row.empty or pd.isna(row.sum()) or row.sum() == 0:
        return 0.0
    proportions = row / row.sum()
    return float(-np.sum(proportions * np.log(proportions)))


def calculate_pielou_evenness(row: pd.Series) -> float:
    row = row[row > 0]
    num_categories = len(row)
    if num_categories <= 1:
        return 0.0
    shannon_index = calculate_shannon_diversity(row)
    if pd.isna(shannon_index):
        return 0.0
    return float(shannon_index / np.log(num_categories))


def build_candidate_table(breakdown_df: pd.DataFrame, comparison_df: pd.DataFrame) -> pd.DataFrame:
    species_col = breakdown_df.columns[0]
    comparison_species_col = comparison_df.columns[0]
    metrics = comparison_df[
        [comparison_species_col, "Simpson_0.0", "Shannon_0.0", "Pielou_0.0"]
    ].rename(
        columns={
            comparison_species_col: species_col,
            "Simpson_0.0": "Simpson_Diversity",
            "Shannon_0.0": "Shannon_Diversity",
            "Pielou_0.0": "Pielou_Evenness",
        }
    )
    return breakdown_df.merge(metrics, on=species_col, how="left", validate="one_to_one")


def recompute_metrics_from_breakdown(breakdown_df: pd.DataFrame) -> pd.DataFrame:
    numeric = breakdown_df.select_dtypes(include=np.number)
    recomputed = breakdown_df[[breakdown_df.columns[0]]].copy()
    recomputed["Simpson_Diversity"] = numeric.apply(calculate_simpson_diversity, axis=1)
    recomputed["Shannon_Diversity"] = numeric.apply(calculate_shannon_diversity, axis=1)
    recomputed["Pielou_Evenness"] = numeric.apply(calculate_pielou_evenness, axis=1)
    return recomputed


def max_abs_numeric_diff(left: pd.DataFrame, right: pd.DataFrame) -> float:
    numeric_cols = [
        col
        for col in left.columns
        if col in right.columns
        and pd.api.types.is_numeric_dtype(left[col])
        and pd.api.types.is_numeric_dtype(right[col])
    ]
    if not numeric_cols:
        return 0.0
    max_diff = 0.0
    for col in numeric_cols:
        diff = (left[col] - right[col]).abs().max()
        if pd.notna(diff):
            max_diff = max(max_diff, float(diff))
    return max_diff


def audit_dataset(config: DiversityTableConfig) -> dict[str, object]:
    breakdown_df = pd.read_csv(config.breakdown_path)
    comparison_df = pd.read_csv(config.comparison_path)
    canonical_df = pd.read_csv(config.canonical_path)

    candidate_df = build_candidate_table(breakdown_df, comparison_df)
    formula_df = recompute_metrics_from_breakdown(breakdown_df)

    candidate_df.to_csv(config.candidate_path, index=False)

    metric_cols = ["Simpson_Diversity", "Shannon_Diversity", "Pielou_Evenness"]
    candidate_metrics = candidate_df[[candidate_df.columns[0]] + metric_cols]
    canonical_metrics = canonical_df[[canonical_df.columns[0]] + metric_cols]

    formula_merged = candidate_metrics.merge(
        formula_df,
        on=candidate_df.columns[0],
        suffixes=("_candidate", "_formula"),
        how="left",
        validate="one_to_one",
    )

    metric_formula_diffs = []
    for metric in metric_cols:
        diff = (
            formula_merged[f"{metric}_candidate"] - formula_merged[f"{metric}_formula"]
        ).abs().max()
        metric_formula_diffs.append(float(diff) if pd.notna(diff) else 0.0)

    candidate_csv_text = config.candidate_path.read_text()
    canonical_csv_text = config.canonical_path.read_text()

    return {
        "dataset": config.dataset,
        "rows": len(candidate_df),
        "columns": len(candidate_df.columns),
        "breakdown_path": str(config.breakdown_path.relative_to(PROJECT_ROOT)),
        "comparison_path": str(config.comparison_path.relative_to(PROJECT_ROOT)),
        "canonical_path": str(config.canonical_path.relative_to(PROJECT_ROOT)),
        "candidate_path": str(config.candidate_path.relative_to(PROJECT_ROOT)),
        "breakdown_sha256": sha256_for_file(config.breakdown_path),
        "comparison_sha256": sha256_for_file(config.comparison_path),
        "canonical_sha256": sha256_for_file(config.canonical_path),
        "row_count_match": len(candidate_df) == len(canonical_df),
        "column_names_match": list(candidate_df.columns) == list(canonical_df.columns),
        "species_order_match": candidate_df.iloc[:, 0].tolist() == canonical_df.iloc[:, 0].tolist(),
        "candidate_missing_metric_rows": int(candidate_df[metric_cols].isna().any(axis=1).sum()),
        "candidate_vs_canonical_all_numeric_max_abs_diff": max_abs_numeric_diff(
            candidate_df,
            canonical_df,
        ),
        "candidate_vs_canonical_metrics_max_abs_diff": max_abs_numeric_diff(
            candidate_metrics,
            canonical_metrics,
        ),
        "candidate_vs_formula_metrics_max_abs_diff": max(metric_formula_diffs, default=0.0),
        "candidate_numeric_match_within_tolerance": max_abs_numeric_diff(candidate_df, canonical_df)
        <= TOLERANCE,
        "formula_metrics_match_within_tolerance": max(metric_formula_diffs, default=0.0)
        <= TOLERANCE,
        "candidate_csv_text_exact_match": candidate_csv_text == canonical_csv_text,
        "promotion_ready_non_destructive": (
            len(candidate_df) == len(canonical_df)
            and list(candidate_df.columns) == list(canonical_df.columns)
            and candidate_df.iloc[:, 0].tolist() == canonical_df.iloc[:, 0].tolist()
            and int(candidate_df[metric_cols].isna().any(axis=1).sum()) == 0
            and max_abs_numeric_diff(candidate_df, canonical_df) <= TOLERANCE
            and max(metric_formula_diffs, default=0.0) <= TOLERANCE
        ),
    }


def main() -> None:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    audit_rows = [audit_dataset(config) for config in CONFIGS]
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(AUDIT_PATH, index=False)

    for row in audit_rows:
        print(
            f"{row['dataset']}: promotion_ready_non_destructive="
            f"{row['promotion_ready_non_destructive']} "
            f"(candidate_vs_canonical_all_numeric_max_abs_diff="
            f"{row['candidate_vs_canonical_all_numeric_max_abs_diff']:.3e}, "
            f"candidate_vs_formula_metrics_max_abs_diff="
            f"{row['candidate_vs_formula_metrics_max_abs_diff']:.3e})"
        )
    print(f"Wrote {AUDIT_PATH}")


if __name__ == "__main__":
    main()

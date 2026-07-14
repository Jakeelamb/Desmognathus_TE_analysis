#!/usr/bin/env python3
"""Build corrected RepeatMasker divergence/gap summaries without replacing history."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "corrected"
    / "repeatmasker_detailed_classification_hit_level_analysis18_v1.csv"
)
DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "corrected"
    / "divergence"
    / "divergence_summary_statistics_by_species_analysis18_v1.csv"
)
HISTORICAL_OUTPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "divergence"
    / "divergence_summary_statistics_by_species.csv"
)
THRESHOLDS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
STAT_COLUMNS = (
    "score",
    "percent_divergence",
    "percent_deletions",
    "percent_insertions",
)
GROUP_COLUMNS = {
    "class": "Class",
    "order": "Order",
    "superfamily": "Superfamily",
}
SUMMARY_COLUMNS = [
    *GROUP_COLUMNS.values(),
    "Desmognathus_Species",
    *STAT_COLUMNS,
    "hitlength_contiglength",
]
INPUT_COLUMNS = [
    *SUMMARY_COLUMNS,
    "query_start",
    "query_end",
]


def _portable(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_summary(
    hits: pd.DataFrame, thresholds: Sequence[float] = THRESHOLDS
) -> pd.DataFrame:
    """Summarize hit-level metrics with the historical threshold design."""

    required = set(SUMMARY_COLUMNS)
    missing = sorted(required.difference(hits.columns))
    if missing:
        raise ValueError(f"Corrected divergence input is missing columns: {missing}")

    frame = hits[SUMMARY_COLUMNS].copy()
    for column in [*STAT_COLUMNS, "hitlength_contiglength"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    outputs = []
    for group_level, group_column in GROUP_COLUMNS.items():
        level = frame.dropna(subset=[group_column, "Desmognathus_Species"]).copy()
        support = level[[group_column, "Desmognathus_Species"]].drop_duplicates()

        def append_summary(
            filtered: pd.DataFrame, threshold: float, threshold_basis: str
        ) -> None:
            grouped = filtered.groupby(
                [group_column, "Desmognathus_Species"], dropna=False
            )
            count = grouped.size().rename("count")
            statistics = grouped[list(STAT_COLUMNS)].agg(["mean", "median", "std"])
            statistics.columns = [
                f"{metric}_{statistic}" for metric, statistic in statistics.columns
            ]
            statistics = pd.concat([count, statistics], axis=1).reset_index()

            result = support.merge(
                statistics,
                on=[group_column, "Desmognathus_Species"],
                how="left",
                validate="one_to_one",
            )
            result["count"] = result["count"].fillna(0).astype(int)
            result = result.rename(columns={group_column: "group_name"})
            result.insert(0, "group_level", group_level)
            result["threshold"] = float(threshold)
            result["threshold_basis"] = threshold_basis
            result["classification_basis"] = "repeatmasker_native_hit_repeat_class"
            result["gap_metric_semantics"] = "alignment_gap_relative_to_repeat_consensus"
            outputs.append(result)

        append_summary(level, -1.0, "none_all_hits")
        for threshold in thresholds:
            filtered = level.loc[level["hitlength_contiglength"].ge(float(threshold))]
            append_summary(
                filtered,
                float(threshold),
                "dnapipete_contig_hitlength_ratio",
            )

    summary = pd.concat(outputs, ignore_index=True)
    ordered = [
        "group_level",
        "group_name",
        "Desmognathus_Species",
        "threshold",
        "threshold_basis",
        "classification_basis",
        "gap_metric_semantics",
        "count",
    ] + [
        f"{metric}_{statistic}"
        for metric in STAT_COLUMNS
        for statistic in ("mean", "median", "std")
    ]
    return summary[ordered].sort_values(
        ["group_level", "group_name", "Desmognathus_Species", "threshold"]
    ).reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def build_threshold_coverage(hits: pd.DataFrame) -> pd.DataFrame:
    """Report the hit/bp mass eligible for legacy dnaPipeTE thresholds."""

    frame = hits[
        [
            "Desmognathus_Species",
            "hitlength_contiglength",
            "query_start",
            "query_end",
        ]
    ].copy()
    start = pd.to_numeric(frame["query_start"], errors="coerce")
    end = pd.to_numeric(frame["query_end"], errors="coerce")
    frame["aligned_bp"] = (end - start).abs() + 1
    if frame["aligned_bp"].isna().any() or frame["aligned_bp"].le(0).any():
        raise ValueError("Corrected divergence input has invalid coordinates")
    frame["has_threshold_context"] = pd.to_numeric(
        frame["hitlength_contiglength"], errors="coerce"
    ).notna()

    rows = []
    for species, group in frame.groupby("Desmognathus_Species", dropna=False):
        has_context = group["has_threshold_context"]
        total_bp = int(group["aligned_bp"].sum())
        context_bp = int(group.loc[has_context, "aligned_bp"].sum())
        rows.append(
            {
                "Desmognathus_Species": species,
                "n_hits_total": len(group),
                "n_hits_with_threshold_context": int(has_context.sum()),
                "n_hits_missing_threshold_context": int((~has_context).sum()),
                "threshold_context_hit_fraction": float(has_context.mean()),
                "aligned_bp_total": total_bp,
                "aligned_bp_with_threshold_context": context_bp,
                "aligned_bp_missing_threshold_context": total_bp - context_bp,
                "threshold_context_bp_fraction": context_bp / total_bp,
            }
        )
    return pd.DataFrame(rows).sort_values("Desmognathus_Species").reset_index(drop=True)


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    coverage_path = output.with_name(
        output.name.replace(
            "divergence_summary_statistics_by_species",
            "divergence_threshold_context_coverage_by_species",
        )
    )
    if output == HISTORICAL_OUTPUT.resolve():
        raise ValueError("Refusing to replace the historical divergence summary")
    occupied = [path for path in (output, manifest_path, coverage_path) if path.exists()]
    if occupied:
        raise FileExistsError(
            "Non-destructive summary requires unused output paths: "
            + ", ".join(str(path) for path in occupied)
        )

    hits = pd.read_csv(args.input, usecols=INPUT_COLUMNS, low_memory=False)
    summary = build_summary(hits)
    coverage = build_threshold_coverage(hits)
    n_hits_with_context = int(coverage["n_hits_with_threshold_context"].sum())
    for group_level in GROUP_COLUMNS:
        threshold_zero = summary.loc[
            summary["group_level"].eq(group_level) & summary["threshold"].eq(0.0)
        ]
        if int(threshold_zero["count"].sum()) != n_hits_with_context:
            raise RuntimeError(
                f"Threshold-zero {group_level} counts do not conserve context-eligible hits"
            )
        all_hits = summary.loc[
            summary["group_level"].eq(group_level) & summary["threshold"].eq(-1.0)
        ]
        if int(all_hits["count"].sum()) != len(hits):
            raise RuntimeError(f"All-hit {group_level} counts do not conserve input hits")

    output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output, index=False)
    coverage.to_csv(coverage_path, index=False)
    manifest = {
        "output_path": _portable(output),
        "output_sha256": _sha256(output),
        "threshold_context_coverage_path": _portable(coverage_path),
        "threshold_context_coverage_sha256": _sha256(coverage_path),
        "input_path": _portable(args.input),
        "input_sha256": _sha256(args.input),
        "historical_output_preserved": _portable(HISTORICAL_OUTPUT),
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "n_input_hits": len(hits),
        "n_input_species": int(hits["Desmognathus_Species"].nunique()),
        "n_summary_rows": len(summary),
        "n_hits_with_threshold_context": n_hits_with_context,
        "n_hits_missing_threshold_context": int(
            coverage["n_hits_missing_threshold_context"].sum()
        ),
        "threshold_context_hit_fraction": n_hits_with_context / len(hits),
        "aligned_bp_with_threshold_context": int(
            coverage["aligned_bp_with_threshold_context"].sum()
        ),
        "aligned_bp_missing_threshold_context": int(
            coverage["aligned_bp_missing_threshold_context"].sum()
        ),
        "thresholds": list(THRESHOLDS),
        "threshold_basis": "dnapipete_contig_hitlength_ratio",
        "classification_basis": "repeatmasker_native_hit_repeat_class",
        "percent_deletions_semantics": "alignment_gap_relative_to_repeat_consensus",
        "percent_insertions_semantics": "alignment_gap_relative_to_repeat_consensus",
        "threshold_zero_counts_by_level": {
            group_level: int(
                summary.loc[
                    summary["group_level"].eq(group_level)
                    & summary["threshold"].eq(0.0),
                    "count",
                ].sum()
            )
            for group_level in GROUP_COLUMNS
        },
        "all_hit_counts_by_level": {
            group_level: int(
                summary.loc[
                    summary["group_level"].eq(group_level)
                    & summary["threshold"].eq(-1.0),
                    "count",
                ].sum()
            )
            for group_level in GROUP_COLUMNS
        },
        "eligible_for_path_sensitivity": True,
        "eligible_as_direct_dna_loss_measure": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Wrote {coverage_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()

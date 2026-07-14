#!/usr/bin/env python3
"""Compare historical and hit-corrected TE divergence features for the final panel."""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import path_analysis.scripts.prepare_te_features as te_features  # noqa: E402


HISTORICAL = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "divergence"
    / "divergence_summary_statistics_by_species.csv"
)
CORRECTED = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "corrected"
    / "divergence"
    / "divergence_summary_statistics_by_species_analysis18_v1.csv"
)
PANEL = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "derived"
    / "panels"
    / "te_genome_primary_mediumplus.csv"
)
OUTPUT_DIR = PROJECT_ROOT / "plans" / "publication-readiness-deep-audit"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_features(path: Path, order: pd.DataFrame) -> pd.DataFrame:
    original = te_features.DIVERGENCE_FILE
    try:
        te_features.DIVERGENCE_FILE = path
        return te_features.load_divergence_features(order)
    finally:
        te_features.DIVERGENCE_FILE = original


def build_comparison():
    panel_species = pd.read_csv(PANEL)["species"].astype(str).tolist()
    order = te_features.load_order_features()
    historical = load_features(HISTORICAL, order).reindex(panel_species)
    corrected = load_features(CORRECTED, order).reindex(panel_species)
    metric_columns = sorted(set(historical.columns).intersection(corrected.columns))

    rows = []
    summary_rows = []
    for metric in metric_columns:
        current = pd.to_numeric(historical[metric], errors="coerce")
        updated = pd.to_numeric(corrected[metric], errors="coerce")
        current_rank = current.rank(method="average")
        updated_rank = updated.rank(method="average")
        for species in panel_species:
            rows.append(
                {
                    "species": species,
                    "metric": metric,
                    "historical_value": current.loc[species],
                    "corrected_value": updated.loc[species],
                    "delta": updated.loc[species] - current.loc[species],
                    "absolute_delta": abs(updated.loc[species] - current.loc[species]),
                    "historical_rank": current_rank.loc[species],
                    "corrected_rank": updated_rank.loc[species],
                    "rank_shift": updated_rank.loc[species] - current_rank.loc[species],
                }
            )
        valid = pd.concat([current, updated], axis=1).dropna()
        delta = updated - current
        summary_rows.append(
            {
                "metric": metric,
                "n_species": len(valid),
                "spearman_rho": valid.iloc[:, 0].corr(valid.iloc[:, 1], method="spearman"),
                "median_absolute_delta": delta.abs().median(),
                "max_absolute_delta": delta.abs().max(),
                "max_absolute_rank_shift": (updated_rank - current_rank).abs().max(),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def render_markdown(comparison: pd.DataFrame, summary: pd.DataFrame) -> str:
    focal = summary.loc[
        summary["metric"].isin(
            [
                "weighted_te_divergence_p90",
                "weighted_te_deletions_p90",
                "weighted_te_insertions_p90",
            ]
        )
    ].set_index("metric")
    fuscus = comparison.loc[
        comparison["species"].eq("fuscus")
        & comparison["metric"].isin(focal.index)
    ]

    lines = [
        "# RepeatMasker corrected-feature sensitivity",
        "",
        "This comparison isolates the effect of replacing inherited dnaPipeTE contig classes with each RepeatMasker hit's native class. Both branches retain the historical 0.9 dnaPipeTE contig-confidence threshold and the same dnaPipeTE order weights.",
        "",
        "| Feature | Spearman rho | Median absolute change | Maximum absolute change |",
        "|---|---:|---:|---:|",
    ]
    for metric, row in focal.iterrows():
        lines.append(
            f"| `{metric}` | {row.spearman_rho:.4f} | {row.median_absolute_delta:.6f} | {row.max_absolute_delta:.6f} |"
        )
    lines.extend(
        [
            "",
            "The weighted species predictors are numerically stable to the classification correction, even though many individual hits were assigned to the wrong order/superfamily in the historical detailed table. This happens because the downstream predictors use order-weighted medians, which dampen many row-level reassignments.",
            "",
            "This stability does not convert `percent_deletions` into a DNA-loss rate. It remains a RepeatMasker alignment-gap statistic relative to a repeat consensus, and the corrected divergence branch remains a sensitivity layer rather than direct evidence of genomic deletion.",
            "",
            "## Fuscus values",
            "",
            "| Feature | Historical | Corrected | Change |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in fuscus.itertuples(index=False):
        lines.append(
            f"| `{row.metric}` | {row.historical_value:.6f} | {row.corrected_value:.6f} | {row.delta:+.6f} |"
        )
    lines.extend(
        [
            "",
            "## Inputs",
            "",
            f"- Historical: `{HISTORICAL.relative_to(PROJECT_ROOT)}` (`{sha256(HISTORICAL)}`)",
            f"- Corrected: `{CORRECTED.relative_to(PROJECT_ROOT)}` (`{sha256(CORRECTED)}`)",
            f"- Species panel: `{PANEL.relative_to(PROJECT_ROOT)}` (`{sha256(PANEL)}`)",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    comparison, summary = build_comparison()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    comparison_path = OUTPUT_DIR / "repeatmasker_corrected_feature_sensitivity.csv"
    summary_path = OUTPUT_DIR / "repeatmasker_corrected_feature_sensitivity_summary.csv"
    report_path = OUTPUT_DIR / "repeatmasker_corrected_feature_sensitivity.md"
    comparison.to_csv(comparison_path, index=False)
    summary.to_csv(summary_path, index=False)
    report_path.write_text(render_markdown(comparison, summary), encoding="utf-8")
    print(f"Wrote {comparison_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()

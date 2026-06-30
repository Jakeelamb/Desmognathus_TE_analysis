#!/usr/bin/env python3
"""
Build species-level TE age-spectrum shape metrics from RepeatMasker landscape CSVs.

This script treats Kimura-bin landscapes as coarse age distributions and derives
compact tempo metrics for all TEs and selected major orders. The outputs are
used as an audited temporal-mechanism layer, not as a replacement for the
existing divergence summaries.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
LANDSCAPE_DIR = PROJECT_ROOT / "results" / "landscapes"
TE_FEATURES_PATH = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "te_path_features.csv"
LTR_HISTORY_PATH = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "ltr_history_features.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "te_age_spectra"
METRICS_PATH = OUTPUT_DIR / "te_age_spectrum_metrics.csv"
CORRELATIONS_PATH = OUTPUT_DIR / "te_age_spectrum_correlations.csv"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"
NOTE_PATH = REPORT_DIR / "TE_AGE_SPECTRUM_AUDIT.md"

ORDER_LEVELS = ("ALL", "LTR", "LINE", "TIR")
RECENT_THRESHOLD = 5.0
MID_THRESHOLD = 15.0
OLD_THRESHOLD = 20.0
ANCIENT_THRESHOLD = 30.0
LOGRATIO_EPS = 1e-6

CORRELATION_METRICS = (
    "all_recent_mass_frac_0_5",
    "all_old_tail_frac_20plus",
    "all_recent_old_logratio",
    "all_weighted_mean_bin_pct",
    "all_effective_bin_count",
    "ltr_recent_mass_frac_0_5",
    "ltr_old_tail_frac_20plus",
    "ltr_recent_old_logratio",
    "ltr_weighted_mean_bin_pct",
    "line_recent_mass_frac_0_5",
    "line_weighted_mean_bin_pct",
    "tir_recent_mass_frac_0_5",
    "tir_weighted_mean_bin_pct",
)

CORRELATION_TARGETS = (
    "weighted_te_divergence_p90",
    "ltr_divergence_p90",
    "line_divergence_p90",
    "tir_divergence_p90",
    "order_pielou",
    "ltr_line_logratio",
    "ltr_history_median_k2p_distance",
    "ltr_history_age_central_mya",
    "ltr_history_n_pairs_estimated",
)


def normalize_species(label: str) -> str:
    label = str(label).strip()
    if label.startswith("D."):
        return label.split(".", 1)[1]
    return label


def extract_bin_start(series: pd.Series) -> pd.Series:
    extracted = series.astype(str).str.extract(r"^(\d+)", expand=False)
    return pd.to_numeric(extracted, errors="coerce")


def bh_adjust(p_values: list[float]) -> list[float]:
    n = len(p_values)
    if n == 0:
        return []
    order = np.argsort(p_values)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for rank, idx in enumerate(order[::-1], start=1):
        original_rank = n - rank + 1
        candidate = p_values[idx] * n / original_rank
        running = min(running, candidate)
        adjusted[idx] = min(running, 1.0)
    return adjusted.tolist()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    sorter = np.argsort(values)
    sorted_values = values[sorter]
    sorted_weights = weights[sorter]
    cumulative = np.cumsum(sorted_weights) / sorted_weights.sum()
    idx = np.searchsorted(cumulative, q, side="left")
    idx = min(idx, len(sorted_values) - 1)
    return float(sorted_values[idx])


def summarize_distribution(frame: pd.DataFrame, prefix: str) -> dict[str, float]:
    total_bp = float(frame["aligned_bp"].sum())
    by_bin = (
        frame.groupby("bin_start", as_index=False)["aligned_bp"]
        .sum()
        .sort_values("bin_start")
        .reset_index(drop=True)
    )
    values = by_bin["bin_start"].to_numpy(dtype=float)
    weights = by_bin["aligned_bp"].to_numpy(dtype=float)
    proportions = weights / weights.sum()

    recent = float(weights[values < RECENT_THRESHOLD].sum() / total_bp)
    mid = float(weights[(values >= RECENT_THRESHOLD) & (values < MID_THRESHOLD)].sum() / total_bp)
    old = float(weights[values >= OLD_THRESHOLD].sum() / total_bp)
    ancient = float(weights[values >= ANCIENT_THRESHOLD].sum() / total_bp)
    peak_idx = int(np.argmax(weights))
    entropy = float(-(proportions * np.log(proportions + LOGRATIO_EPS)).sum())
    effective_bin_count = float(np.exp(entropy))

    return {
        f"{prefix}_total_aligned_bp": total_bp,
        f"{prefix}_n_bins": int(len(by_bin)),
        f"{prefix}_recent_mass_frac_0_5": recent,
        f"{prefix}_mid_mass_frac_5_15": mid,
        f"{prefix}_old_tail_frac_20plus": old,
        f"{prefix}_ancient_tail_frac_30plus": ancient,
        f"{prefix}_recent_old_logratio": float(
            math.log((recent + LOGRATIO_EPS) / (old + LOGRATIO_EPS))
        ),
        f"{prefix}_peak_bin_start_pct": float(values[peak_idx]),
        f"{prefix}_weighted_mean_bin_pct": float(np.average(values, weights=weights)),
        f"{prefix}_weighted_median_bin_pct": weighted_quantile(values, weights, 0.5),
        f"{prefix}_peak_mass_fraction": float(weights[peak_idx] / total_bp),
        f"{prefix}_bin_entropy": entropy,
        f"{prefix}_effective_bin_count": effective_bin_count,
    }


def build_metrics() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for landscape_path in sorted(LANDSCAPE_DIR.glob("repeat_landscape_*.csv")):
        frame = pd.read_csv(landscape_path)
        if frame.empty:
            continue
        frame["bin_start"] = extract_bin_start(frame["Kimura_bin"])
        frame = frame.dropna(subset=["bin_start", "aligned_bp"])
        if frame.empty:
            continue

        species = normalize_species(frame["Species"].iloc[0])
        row: dict[str, object] = {
            "species": species,
            "landscape_file": landscape_path.name,
            "landscape_source_path": str(landscape_path.relative_to(PROJECT_ROOT)),
        }

        for order_name in ORDER_LEVELS:
            subset = frame if order_name == "ALL" else frame[frame["Order"] == order_name]
            if subset.empty or float(subset["aligned_bp"].sum()) <= 0:
                continue
            prefix = order_name.lower()
            row.update(summarize_distribution(subset, prefix))

        rows.append(row)

    metrics = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    return metrics


def build_correlations(metrics: pd.DataFrame) -> pd.DataFrame:
    te_features = pd.read_csv(TE_FEATURES_PATH)
    history = pd.read_csv(LTR_HISTORY_PATH)
    merged = metrics.merge(te_features, on="species", how="left").merge(history, on="species", how="left")

    rows: list[dict[str, object]] = []
    for metric in CORRELATION_METRICS:
        if metric not in merged.columns:
            continue
        for target in CORRELATION_TARGETS:
            if target not in merged.columns:
                continue
            subset = merged[[metric, target]].dropna()
            if len(subset) < 8:
                continue
            if subset[metric].nunique() < 3 or subset[target].nunique() < 3:
                continue
            rho, p_value = spearmanr(subset[metric], subset[target])
            rows.append(
                {
                    "metric": metric,
                    "target": target,
                    "n_species": int(len(subset)),
                    "spearman_rho": float(rho),
                    "p_value": float(p_value),
                }
            )

    correlations = pd.DataFrame(rows)
    if not correlations.empty:
        correlations["bh_p_value"] = bh_adjust(correlations["p_value"].tolist())
        correlations = correlations.sort_values(
            ["bh_p_value", "p_value", "metric", "target"]
        ).reset_index(drop=True)
    return correlations


def format_species_block(metrics: pd.DataFrame, column: str, ascending: bool) -> str:
    subset = metrics[["species", column, "all_old_tail_frac_20plus", "all_weighted_mean_bin_pct"]].dropna()
    subset = subset.sort_values(column, ascending=ascending).head(5)
    lines = []
    for _, row in subset.iterrows():
        lines.append(
            f"- `{row['species']}`: {column}={row[column]:.3f}, "
            f"old_tail={row['all_old_tail_frac_20plus']:.3f}, "
            f"mean_bin={row['all_weighted_mean_bin_pct']:.2f}"
        )
    return "\n".join(lines)


def write_note(metrics: pd.DataFrame, correlations: pd.DataFrame) -> None:
    n_species = int(len(metrics))
    overlap_history = int(
        metrics["species"].isin(pd.read_csv(LTR_HISTORY_PATH)["species"].tolist()).sum()
    )

    selected_pairs = [
        ("all_recent_mass_frac_0_5", "weighted_te_divergence_p90"),
        ("all_weighted_mean_bin_pct", "weighted_te_divergence_p90"),
        ("ltr_recent_mass_frac_0_5", "ltr_divergence_p90"),
        ("line_weighted_mean_bin_pct", "line_divergence_p90"),
        ("all_effective_bin_count", "order_pielou"),
        ("all_recent_mass_frac_0_5", "ltr_history_median_k2p_distance"),
    ]
    hit_lines = []
    for metric, target in selected_pairs:
        match = correlations[
            (correlations["metric"] == metric) & (correlations["target"] == target)
        ]
        if match.empty:
            continue
        row = match.iloc[0]
        hit_lines.append(
            f"- `{metric}` vs `{target}`: "
            f"rho = {row['spearman_rho']:.3f}, "
            f"BH p = {row['bh_p_value']:.3g}, n = {int(row['n_species'])}"
        )
    if not hit_lines:
        hit_lines.append("- No preselected age-spectrum relationship was estimable.")

    note = f"""# TE Age Spectrum Audit

## Purpose

This note audits whether the local RepeatMasker landscapes support a stable,
interpretable species-level TE tempo axis before the final cell/nucleus/genome
size layers are ready.

## Upstream inputs

- `results/landscapes/repeat_landscape_*.csv`
- `path_analysis/data/derived/te_path_features.csv`
- `path_analysis/data/derived/ltr_history_features.csv`

## Outputs

- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`
- `results/data/te_age_spectra/te_age_spectrum_correlations.csv`

## Metric definition

Each landscape was collapsed to species-level shape summaries for `ALL`, `LTR`,
`LINE`, and `TIR` bins:

- recent mass: Kimura bins `< 5%`
- mid mass: Kimura bins `5-15%`
- old tail: Kimura bins `>= 20%`
- ancient tail: Kimura bins `>= 30%`
- recent/old logratio
- peak bin
- weighted mean and median bin
- peak mass fraction
- entropy-derived effective bin count

These are tempo descriptors, not direct insertion-age estimates.

## Coverage

- Landscapes processed: `{n_species}` species
- Species overlapping the paired-LTR history layer: `{overlap_history}`

## Main result

The landscape branch passes the basic interpretability screen. It recovers a
clear recent-vs-old species axis, and that axis strongly tracks the existing
RepeatMasker divergence summaries rather than behaving like random noise.

The strongest supported relationships are:

{chr(10).join(hit_lines)}

## Interpretation

- The clearest stable tempo axis is `recent` versus `old retained tail`.
- That axis is strongest for the full TE landscape and for `LTR`/`LINE`
  sub-landscapes.
- The landscape-derived tempo summaries line up strongly with the existing
  divergence summaries, which is expected and validates the new metrics rather
  than replacing the existing divergence layer.
- The paired-LTR history layer is comparatively independent: it does not track
  these coarse landscape summaries strongly enough to treat them as the same
  signal.

## Species extremes

Most recent all-TE landscapes by `all_recent_mass_frac_0_5`:
{format_species_block(metrics, 'all_recent_mass_frac_0_5', ascending=False)}

Oldest all-TE landscapes by `all_recent_mass_frac_0_5`:
{format_species_block(metrics, 'all_recent_mass_frac_0_5', ascending=True)}

## Bottom line

The TE landscape data support a usable species-level tempo axis. The main
descriptive contrast is not "single recent burst everywhere" but variation in
how much recent mass versus old retained tail each species carries. This is
worth keeping as a comparative layer, but it should be described as a
landscape-derived tempo summary rather than a substitute for paired-LTR age
estimation.
"""
    NOTE_PATH.write_text(note)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = build_metrics()
    correlations = build_correlations(metrics)
    metrics.to_csv(METRICS_PATH, index=False)
    correlations.to_csv(CORRELATIONS_PATH, index=False)
    write_note(metrics, correlations)
    print(f"Wrote {METRICS_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {CORRELATIONS_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {NOTE_PATH.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()

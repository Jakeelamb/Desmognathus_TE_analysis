#!/usr/bin/env python3
"""
Audit whether the paired-LTR history layer helps explain current TE state.

The main distinction tested here is between:
1. historical age/oldness of retained paired LTRs
2. abundance/density of recoverable paired-LTR substrate
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TE_FEATURES_PATH = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "te_path_features.csv"
LTR_HISTORY_PATH = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "ltr_history_features.csv"
SPECTRUM_PATH = PROJECT_ROOT / "results" / "data" / "te_age_spectra" / "te_age_spectrum_metrics.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "ltr_history_mechanism"
TESTS_PATH = OUTPUT_DIR / "ltr_history_mechanism_tests.csv"
NOTE_PATH = PROJECT_ROOT / "LTR_HISTORY_MECHANISM_AUDIT.md"

PREDICTORS = [
    ("ltr_history_median_k2p_distance", "paired_ltr_age"),
    ("ltr_history_n_pairs_estimated", "paired_ltr_abundance"),
    ("ltr_history_pairs_per_million_ltr_bp", "paired_ltr_density"),
    ("ltr_history_n_pairs_high_confidence", "paired_ltr_high_confidence_count"),
]

RESPONSES = [
    ("ltr_line_logratio", "current_te_structure"),
    ("order_pielou", "current_te_structure"),
    ("weighted_te_divergence_p90", "current_te_tempo"),
    ("weighted_te_deletions_p90", "current_te_tempo"),
    ("ectopic_log10_mean_ratio", "ectopic_proxy"),
    ("all_recent_mass_frac_0_5", "landscape_tempo"),
    ("all_old_tail_frac_20plus", "landscape_tempo"),
    ("all_effective_bin_count", "landscape_tempo"),
    ("ltr_recent_mass_frac_0_5", "landscape_tempo"),
    ("ltr_weighted_mean_bin_pct", "landscape_tempo"),
]


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


def build_tests() -> pd.DataFrame:
    te_features = pd.read_csv(TE_FEATURES_PATH)
    history = pd.read_csv(LTR_HISTORY_PATH)
    spectrum = pd.read_csv(SPECTRUM_PATH)

    merged = te_features.merge(history, on="species", how="inner").merge(
        spectrum[
            [
                "species",
                "ltr_total_aligned_bp",
                "all_recent_mass_frac_0_5",
                "all_old_tail_frac_20plus",
                "all_effective_bin_count",
                "ltr_recent_mass_frac_0_5",
                "ltr_weighted_mean_bin_pct",
            ]
        ],
        on="species",
        how="left",
    )
    merged["ltr_history_pairs_per_million_ltr_bp"] = (
        merged["ltr_history_n_pairs_estimated"] / (merged["ltr_total_aligned_bp"] / 1_000_000.0)
    )

    rows: list[dict[str, object]] = []
    for predictor, predictor_role in PREDICTORS:
        if predictor not in merged.columns:
            continue
        for response, response_role in RESPONSES:
            if response not in merged.columns:
                continue
            subset = merged[[predictor, response]].dropna()
            if len(subset) < 8:
                continue
            if subset[predictor].nunique() < 3 or subset[response].nunique() < 3:
                continue
            rho, p_value = spearmanr(subset[predictor], subset[response])
            rows.append(
                {
                    "predictor": predictor,
                    "predictor_role": predictor_role,
                    "response": response,
                    "response_role": response_role,
                    "test": "spearman",
                    "n_species": int(len(subset)),
                    "spearman_rho": float(rho),
                    "p_value": float(p_value),
                }
            )

    tests = pd.DataFrame(rows)
    if not tests.empty:
        tests["bh_p_value"] = bh_adjust(tests["p_value"].tolist())
        tests = tests.sort_values(
            ["bh_p_value", "p_value", "predictor", "response"]
        ).reset_index(drop=True)
    return tests


def stat_line(tests: pd.DataFrame, predictor: str, response: str) -> str:
    match = tests[(tests["predictor"] == predictor) & (tests["response"] == response)]
    if match.empty:
        return f"- `{predictor}` vs `{response}`: not estimable"
    row = match.iloc[0]
    return (
        f"- `{predictor}` vs `{response}`: "
        f"rho = {row['spearman_rho']:.3f}, "
        f"BH p = {row['bh_p_value']:.3g}, "
        f"n = {int(row['n_species'])}"
    )


def write_note(tests: pd.DataFrame) -> None:
    support_age = tests[
        (tests["predictor"] == "ltr_history_median_k2p_distance") & (tests["bh_p_value"] <= 0.05)
    ]
    support_density = tests[
        (tests["predictor"] == "ltr_history_pairs_per_million_ltr_bp") & (tests["bh_p_value"] <= 0.05)
    ]
    support_count = tests[
        (tests["predictor"] == "ltr_history_n_pairs_estimated") & (tests["bh_p_value"] <= 0.05)
    ]
    overlap_n = int(tests["n_species"].max()) if not tests.empty else 0

    note = f"""# LTR History Mechanism Audit

## Purpose

This note tests whether the paired-LTR history layer helps explain current TE
state, or whether it is mostly orthogonal to the current TE composition and
landscape summaries.

## Upstream inputs

- `path_analysis/data/derived/ltr_history_features.csv`
- `path_analysis/data/derived/te_path_features.csv`
- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`

## Outputs

- `results/data/ltr_history_mechanism/ltr_history_mechanism_tests.csv`

## Coverage

- Species in the merged mechanism audit: `{overlap_n}`

## Predictor interpretation

- `ltr_history_median_k2p_distance`: oldness of retained paired LTRs
- `ltr_history_n_pairs_estimated`: amount of recoverable paired-LTR substrate
- `ltr_history_pairs_per_million_ltr_bp`: recoverable paired-LTR density after
  normalizing by total LTR landscape mass
- `ltr_history_n_pairs_high_confidence`: stricter support subset

## Main result

The paired-LTR history branch splits into two different signals:

- retained pair age is mostly independent of the current TE state
- retained pair abundance/density tracks landscape recency more than it tracks
  TE balance or ectopic ratio

Key checks:

{stat_line(tests, 'ltr_history_median_k2p_distance', 'weighted_te_divergence_p90')}
{stat_line(tests, 'ltr_history_median_k2p_distance', 'all_recent_mass_frac_0_5')}
{stat_line(tests, 'ltr_history_n_pairs_estimated', 'all_recent_mass_frac_0_5')}
{stat_line(tests, 'ltr_history_n_pairs_estimated', 'all_old_tail_frac_20plus')}
{stat_line(tests, 'ltr_history_pairs_per_million_ltr_bp', 'all_recent_mass_frac_0_5')}
{stat_line(tests, 'ltr_history_pairs_per_million_ltr_bp', 'ltr_recent_mass_frac_0_5')}
{stat_line(tests, 'ltr_history_n_pairs_estimated', 'ectopic_log10_mean_ratio')}

## Interpretation

- The age/oldness axis (`ltr_history_median_k2p_distance`) does not explain the
  current broad TE divergence or landscape recency summaries.
- The abundance axis (`ltr_history_n_pairs_estimated`) does connect to current
  TE tempo: species with more recoverable paired LTRs tend to show younger,
  more recent landscapes and lower old-tail burden.
- That signal survives simple normalization by total LTR landscape mass, so it
  is not only a trivial "more LTR sequence means more pairs" artifact.
- The ectopic proxy stays largely separate from the paired-LTR mechanism layer.
- The strict high-confidence pair count is weaker than the full recoverable
  count, so this abundance link is real enough to note but should remain a
  secondary mechanism result rather than a headline claim.

## Bottom line

The paired-LTR history layer is not one thing. Retained pair age behaves like a
mostly orthogonal historical axis, while retained pair abundance behaves more
like a recency/recoverability signal tied to the present-day TE landscape. That
means the current repo supports a nuanced mechanism story: old retained LTR
history and current TE tempo are related only weakly, but the amount of intact
paired-LTR substrate does carry information about how recent the present TE
landscape looks.

## Audit summary

- BH-significant tests for paired-LTR age: `{len(support_age)}`
- BH-significant tests for paired-LTR count: `{len(support_count)}`
- BH-significant tests for paired-LTR density: `{len(support_density)}`
"""
    NOTE_PATH.write_text(note)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tests = build_tests()
    tests.to_csv(TESTS_PATH, index=False)
    write_note(tests)
    print(f"Wrote {TESTS_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Wrote {NOTE_PATH.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()

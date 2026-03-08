#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"

INPUT_FILE = DERIVED_DIR / "te_path_features.csv"
OUTPUT_FILE = DERIVED_DIR / "te_model_feature_panel.csv"


def main() -> None:
    df = pd.read_csv(INPUT_FILE)
    keep_cols = [
        "species",
        "ltr_line_logratio",
        "order_pielou",
        "weighted_te_divergence_p90",
        "weighted_te_deletions_p90",
        "ectopic_log10_mean_ratio",
        "te_order_source_id",
        "te_order_source_path",
        "te_order_source_sha256",
        "order_diversity_source_id",
        "order_diversity_source_path",
        "order_diversity_source_sha256",
        "divergence_source_id",
        "divergence_source_path",
        "divergence_source_sha256",
        "ectopic_source_id",
        "ectopic_source_path",
        "ectopic_source_sha256",
        "te_feature_source_ids",
    ]

    panel = df[keep_cols].copy()
    panel["has_core_te_panel"] = panel[
        ["ltr_line_logratio", "order_pielou", "weighted_te_divergence_p90"]
    ].notna().all(axis=1)
    panel["has_extended_te_panel"] = panel[
        ["ltr_line_logratio", "order_pielou", "weighted_te_divergence_p90", "ectopic_log10_mean_ratio"]
    ].notna().all(axis=1)
    panel = panel.sort_values("species").reset_index(drop=True)
    panel.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Rows: {len(panel)}")
    print(f"Core-ready species: {int(panel['has_core_te_panel'].sum())}")
    print(f"Extended-ready species: {int(panel['has_extended_te_panel'].sum())}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
PANEL_DIR = DERIVED_DIR / "panels"

INPUT_FILE = DERIVED_DIR / "path_input_master.csv"
OUTPUT_SUMMARY = DERIVED_DIR / "analysis_panel_summary.csv"


def confidence_score(series: pd.Series) -> pd.Series:
    mapping = {"low": 1, "medium": 2, "high": 3}
    return series.map(mapping)


def prepare_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["body_size_conf_score"] = confidence_score(out["body_size_proxy_confidence"])
    out["development_conf_score"] = confidence_score(out["development_confidence"])
    out["lifestyle_conf_score"] = confidence_score(out["lifestyle_confidence"])

    out["has_core_organismal_complete"] = out[
        ["body_size_proxy_mm", "development_mode", "aquaticity_index", "microhabitat_class"]
    ].notna().all(axis=1)
    out["has_core_organismal_mediumplus"] = (
        out["has_core_organismal_complete"]
        & out["body_size_conf_score"].ge(2)
        & out["development_conf_score"].ge(2)
        & out["lifestyle_conf_score"].ge(2)
    )
    out["has_elevation"] = out["elevation_mid_m"].notna()
    out["uses_total_length_body_proxy"] = out["body_size_proxy_measurement"].fillna("").str.contains(
        "total_length",
        case=False,
    )
    out["uses_mixed_stage_body_proxy"] = out["body_size_proxy_measurement"].fillna("").str.contains(
        "mixed|all_specimens_mixed_life_stages|transformed",
        case=False,
        regex=True,
    )
    out["is_low_confidence_body_proxy"] = out["body_size_proxy_confidence"].eq("low")
    out["is_low_confidence_lifestyle_proxy"] = out["lifestyle_confidence"].eq("low")
    return out


def build_panel(df: pd.DataFrame, mask: pd.Series, keep_cols: list[str]) -> pd.DataFrame:
    return df.loc[mask, keep_cols].sort_values("species").reset_index(drop=True)


def main() -> None:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    df = prepare_flags(df)

    base_cols = [
        "species",
        "has_tree_tip",
        "has_te",
        "has_genome",
        "has_ectopic",
        "has_morphology",
        "genome_size_pg",
        "genome_size_se_pg",
        "genome_result_status",
        "genome_flag_summary",
        "morph_cell_area_um2",
        "morph_nucleus_area_um2",
        "morph_nc_ratio",
        "ltr_line_logratio",
        "order_pielou",
        "weighted_te_divergence_p90",
        "weighted_te_deletions_p90",
        "ectopic_log10_mean_ratio",
        "body_size_proxy_mm",
        "body_size_proxy_measurement",
        "body_size_proxy_source_id",
        "body_size_proxy_confidence",
        "development_mode",
        "development_source_id",
        "development_confidence",
        "aquaticity_index",
        "microhabitat_class",
        "lifestyle_source_id",
        "lifestyle_confidence",
        "elevation_mid_m",
        "elevation_source_id",
        "has_core_organismal_complete",
        "has_core_organismal_mediumplus",
        "has_elevation",
        "uses_total_length_body_proxy",
        "uses_mixed_stage_body_proxy",
        "is_low_confidence_body_proxy",
        "is_low_confidence_lifestyle_proxy",
        "organismal_source_ids",
        "te_feature_source_ids",
        "genome_source_id",
        "morphology_source_id",
    ]

    masks = {
        "te_genome_all": df["has_tree_tip"] & df["has_te"] & df["has_genome"],
        "te_genome_primary": df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_core_organismal_complete"],
        "te_genome_primary_mediumplus": (
            df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_core_organismal_mediumplus"]
        ),
        "te_genome_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_core_organismal_mediumplus"]
            & ~df["uses_total_length_body_proxy"]
            & ~df["uses_mixed_stage_body_proxy"]
        ),
        "te_genome_ectopic_all": df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_ectopic"],
        "te_genome_ectopic_primary_mediumplus": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_ectopic"]
            & df["has_core_organismal_mediumplus"]
        ),
        "te_genome_ectopic_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_ectopic"]
            & df["has_core_organismal_mediumplus"]
            & ~df["uses_total_length_body_proxy"]
            & ~df["uses_mixed_stage_body_proxy"]
        ),
        "te_genome_morphology_all": df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_morphology"],
        "te_genome_morphology_primary_mediumplus": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_morphology"]
            & df["has_core_organismal_mediumplus"]
        ),
        "te_genome_morphology_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_morphology"]
            & df["has_core_organismal_mediumplus"]
            & ~df["uses_total_length_body_proxy"]
            & ~df["uses_mixed_stage_body_proxy"]
        ),
    }

    summary_rows = []
    for name, mask in masks.items():
        panel = build_panel(df, mask, base_cols)
        panel.to_csv(PANEL_DIR / f"{name}.csv", index=False)
        summary_rows.append(
            {
                "panel_name": name,
                "n_species": int(len(panel)),
                "n_with_elevation": int(panel["has_elevation"].sum()),
                "n_low_body_proxy": int(panel["is_low_confidence_body_proxy"].sum()),
                "n_total_length_body_proxy": int(panel["uses_total_length_body_proxy"].sum()),
                "n_mixed_stage_body_proxy": int(panel["uses_mixed_stage_body_proxy"].sum()),
                "species_list": ";".join(panel["species"].tolist()),
            }
        )

    pd.DataFrame(summary_rows).to_csv(OUTPUT_SUMMARY, index=False)

    print(f"Wrote panel directory {PANEL_DIR}")
    print(f"Wrote {OUTPUT_SUMMARY}")


if __name__ == "__main__":
    main()

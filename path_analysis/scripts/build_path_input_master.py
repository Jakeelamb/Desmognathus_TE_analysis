#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
EXTERNAL_DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

MASTER_FILE = DERIVED_DIR / "master_species_table.csv"
TE_FEATURE_FILE = DERIVED_DIR / "te_path_features.csv"
AMPHIBIO_FILE = EXTERNAL_DERIVED_DIR / "amphibio_desmognathus_traits.csv"
CONANTI_FILE = EXTERNAL_DERIVED_DIR / "conanti_fuscus_morphometrics_summary.csv"
ORGANISMAL_FILE = DERIVED_DIR / "organismal_traits_curated.csv"

OUTPUT_MASTER = DERIVED_DIR / "path_input_master.csv"
OUTPUT_COVERAGE = DERIVED_DIR / "path_input_coverage_summary.csv"
OUTPUT_PRIORITY = DERIVED_DIR / "manual_trait_priority_species.csv"


def with_prefix(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    rename_map = {c: f"{prefix}{c}" for c in df.columns if c != "species"}
    return df.rename(columns=rename_map)


def summarize_coverage(df: pd.DataFrame) -> pd.DataFrame:
    datasets = {
        "all_species": df["species"].notna(),
        "te_genome_overlap": df["has_tree_tip"] & df["has_te"] & df["has_genome"],
        "te_genome_morph_overlap": df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_morphology"],
    }

    traits = [
        "ltr_line_logratio",
        "order_pielou",
        "weighted_te_divergence_p90",
        "weighted_te_deletions_p90",
        "ectopic_log10_mean_ratio",
        "body_size_proxy_mm",
        "development_mode",
        "aquaticity_index",
        "microhabitat_class",
        "elevation_mid_m",
        "amphibio_body_size_mm",
        "amphibio_body_mass_g",
        "amphibio_size_at_maturity_min_mm",
        "amphibio_size_at_maturity_max_mm",
        "amphibio_dir",
        "amphibio_lar",
        "amphibio_aqu",
        "amphibio_ter",
        "cf_svl_q90_all_mm",
        "cf_svl_max_mm",
    ]

    rows = []
    for trait in traits:
        if trait not in df.columns:
            continue
        row = {"trait_name": trait}
        for dataset_name, mask in datasets.items():
            row[f"n_{dataset_name}"] = int(df.loc[mask, trait].notna().sum())
        rows.append(row)

    return pd.DataFrame(rows)


def summarize_priority_species(df: pd.DataFrame) -> pd.DataFrame:
    summary = df[["species", "has_tree_tip", "has_te", "has_genome", "has_morphology"]].copy()
    summary["te_genome"] = summary["has_tree_tip"] & summary["has_te"] & summary["has_genome"]
    summary["te_genome_morph"] = summary["te_genome"] & summary["has_morphology"]
    summary["has_body_size_proxy"] = df["body_size_proxy_mm"].notna()
    summary["has_development_proxy"] = df["development_mode"].notna()
    summary["has_lifestyle_proxy"] = df["aquaticity_index"].notna()
    summary["missing_any_priority_proxy"] = ~(
        summary["has_body_size_proxy"] & summary["has_development_proxy"] & summary["has_lifestyle_proxy"]
    )
    summary = summary[summary["te_genome"] & summary["missing_any_priority_proxy"]].copy()
    summary["priority_group"] = summary["te_genome_morph"].map(
        {True: "integrated_subset_gap", False: "te_genome_subset_gap"}
    )
    keep_cols = [
        "species",
        "priority_group",
        "te_genome",
        "te_genome_morph",
        "has_body_size_proxy",
        "has_development_proxy",
        "has_lifestyle_proxy",
    ]
    return summary[keep_cols].sort_values(["priority_group", "species"])


def main() -> None:
    master = pd.read_csv(MASTER_FILE)
    te = pd.read_csv(TE_FEATURE_FILE)
    amphibio = pd.read_csv(AMPHIBIO_FILE) if AMPHIBIO_FILE.exists() else pd.DataFrame(columns=["species"])
    conanti = pd.read_csv(CONANTI_FILE) if CONANTI_FILE.exists() else pd.DataFrame(columns=["species"])
    organismal = pd.read_csv(ORGANISMAL_FILE) if ORGANISMAL_FILE.exists() else pd.DataFrame(columns=["species"])

    merged = master.merge(te, on="species", how="left", suffixes=("", "_te"))
    merged = merged.merge(amphibio, on="species", how="left")
    merged = merged.merge(with_prefix(conanti, "cf_"), on="species", how="left")
    merged = merged.merge(organismal, on="species", how="left")

    merged["has_te_feature_table"] = merged["ltr_line_logratio"].notna()
    merged["has_amphibio_traits"] = merged["amphibio_source_id"].notna()
    merged["has_conanti_fuscus_morphometrics"] = merged["cf_source_id"].notna()
    merged["has_curated_organismal_traits"] = merged["organismal_source_ids"].notna()

    coverage = summarize_coverage(merged)
    priority_species = summarize_priority_species(merged)
    merged.to_csv(OUTPUT_MASTER, index=False)
    coverage.to_csv(OUTPUT_COVERAGE, index=False)
    priority_species.to_csv(OUTPUT_PRIORITY, index=False)

    print(f"Wrote {OUTPUT_MASTER}")
    print(f"Wrote {OUTPUT_COVERAGE}")
    print(f"Wrote {OUTPUT_PRIORITY}")


if __name__ == "__main__":
    main()

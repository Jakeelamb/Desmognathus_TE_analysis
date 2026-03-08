#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
EXTERNAL_DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

MASTER_FILE = DERIVED_DIR / "master_species_table.csv"
AMPHIBIO_FILE = EXTERNAL_DERIVED_DIR / "amphibio_desmognathus_traits.csv"
CONANTI_FILE = EXTERNAL_DERIVED_DIR / "conanti_fuscus_morphometrics_summary.csv"
NC_FILE = EXTERNAL_DERIVED_DIR / "nc_biodiversity_desmognathus_traits.csv"
MANUAL_EXTRACTION_FILE = PROJECT_ROOT / "path_analysis" / "data" / "templates" / "literature_trait_extraction.csv"

OUTPUT_FILE = DERIVED_DIR / "organismal_traits_curated.csv"
OUTPUT_COVERAGE = DERIVED_DIR / "organismal_traits_coverage_summary.csv"


def load_optional_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["species"])
    return pd.read_csv(path)


def load_manual_trait_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["species"])
    raw = pd.read_csv(path)
    if raw.empty:
        return pd.DataFrame(columns=["species"])

    raw = raw.rename(columns={"current_species": "species"})
    raw["value_num"] = pd.to_numeric(raw["value"], errors="coerce")
    raw["curation_flag"] = raw["curation_flag"].fillna("")
    raw = raw[~raw["curation_flag"].str.contains("exclude", case=False)].copy()

    def join_unique(series: pd.Series) -> str | float:
        values = sorted({str(v) for v in series if pd.notna(v) and str(v).strip()})
        return ";".join(values) if values else np.nan

    def numeric_trait_block(traits: list[str], rename_map: dict[str, str], prefix: str) -> pd.DataFrame:
        block = raw[raw["trait_name"].isin(traits)].copy()
        if block.empty:
            return pd.DataFrame(columns=["species"])
        pivot = (
            block.pivot_table(index="species", columns="trait_name", values="value_num", aggfunc="first")
            .reset_index()
            .rename(columns=rename_map)
        )
        for column in rename_map.values():
            if column not in pivot.columns:
                pivot[column] = np.nan

        if prefix == "manual_body_size":
            if "manual_adult_svl_mid_mm" not in pivot.columns:
                pivot["manual_adult_svl_mid_mm"] = np.nan
            if "manual_adult_tl_mid_mm" not in pivot.columns:
                pivot["manual_adult_tl_mid_mm"] = np.nan

            svl_mid_missing = (
                pivot["manual_adult_svl_mid_mm"].isna()
                & pivot["manual_adult_svl_min_mm"].notna()
                & pivot["manual_adult_svl_max_mm"].notna()
            )
            pivot.loc[svl_mid_missing, "manual_adult_svl_mid_mm"] = (
                pivot.loc[svl_mid_missing, "manual_adult_svl_min_mm"]
                + pivot.loc[svl_mid_missing, "manual_adult_svl_max_mm"]
            ) / 2.0

            tl_mid_missing = (
                pivot["manual_adult_tl_mid_mm"].isna()
                & pivot["manual_adult_tl_min_mm"].notna()
                & pivot["manual_adult_tl_max_mm"].notna()
            )
            pivot.loc[tl_mid_missing, "manual_adult_tl_mid_mm"] = (
                pivot.loc[tl_mid_missing, "manual_adult_tl_min_mm"]
                + pivot.loc[tl_mid_missing, "manual_adult_tl_max_mm"]
            ) / 2.0

        if prefix == "manual_elevation":
            if "manual_elevation_mid_m" not in pivot.columns:
                pivot["manual_elevation_mid_m"] = np.nan
            elev_mid_missing = (
                pivot["manual_elevation_mid_m"].isna()
                & pivot["manual_elevation_min_m"].notna()
                & pivot["manual_elevation_max_m"].notna()
            )
            pivot.loc[elev_mid_missing, "manual_elevation_mid_m"] = (
                pivot.loc[elev_mid_missing, "manual_elevation_min_m"]
                + pivot.loc[elev_mid_missing, "manual_elevation_max_m"]
            ) / 2.0

        provenance = (
            block.groupby("species")
            .agg(
                **{
                    f"{prefix}_source_id": ("source_id", join_unique),
                    f"{prefix}_source_locator": ("source_locator", join_unique),
                    f"{prefix}_source_citation": ("source_citation", join_unique),
                }
            )
            .reset_index()
        )
        return pivot.merge(provenance, on="species", how="left")

    def text_trait_block(traits: list[str], rename_map: dict[str, str], prefix: str) -> pd.DataFrame:
        block = raw[raw["trait_name"].isin(traits)].copy()
        if block.empty:
            return pd.DataFrame(columns=["species"])
        pivot = (
            block.pivot_table(index="species", columns="trait_name", values="value", aggfunc="first")
            .reset_index()
            .rename(columns=rename_map)
        )
        for column in rename_map.values():
            if column not in pivot.columns:
                pivot[column] = np.nan
        provenance = (
            block.groupby("species")
            .agg(
                **{
                    f"{prefix}_source_id": ("source_id", join_unique),
                    f"{prefix}_source_locator": ("source_locator", join_unique),
                    f"{prefix}_source_citation": ("source_citation", join_unique),
                }
            )
            .reset_index()
        )
        return pivot.merge(provenance, on="species", how="left")

    manual_size = numeric_trait_block(
        [
            "adult_svl_min_mm",
            "adult_svl_max_mm",
            "adult_svl_mid_mm",
            "adult_total_length_min_mm",
            "adult_total_length_max_mm",
            "adult_total_length_mid_mm",
        ],
        {
            "adult_svl_min_mm": "manual_adult_svl_min_mm",
            "adult_svl_max_mm": "manual_adult_svl_max_mm",
            "adult_svl_mid_mm": "manual_adult_svl_mid_mm",
            "adult_total_length_min_mm": "manual_adult_tl_min_mm",
            "adult_total_length_max_mm": "manual_adult_tl_max_mm",
            "adult_total_length_mid_mm": "manual_adult_tl_mid_mm",
        },
        "manual_body_size",
    )
    manual_development = text_trait_block(
        ["development_mode"],
        {"development_mode": "manual_development_mode"},
        "manual_development",
    )
    manual_lifestyle = (
        numeric_trait_block(
            ["aquaticity_index"],
            {"aquaticity_index": "manual_aquaticity_index"},
            "manual_lifestyle",
        )
        .merge(
            text_trait_block(
                ["microhabitat_class"],
                {"microhabitat_class": "manual_microhabitat_class"},
                "manual_lifestyle",
            ),
            on="species",
            how="outer",
            suffixes=("", "_text"),
        )
        .rename(
            columns={
                "manual_lifestyle_source_id_text": "manual_lifestyle_source_id_text",
                "manual_lifestyle_source_locator_text": "manual_lifestyle_source_locator_text",
                "manual_lifestyle_source_citation_text": "manual_lifestyle_source_citation_text",
            }
        )
    )
    if not manual_lifestyle.empty:
        for base_col, text_col in [
            ("manual_lifestyle_source_id", "manual_lifestyle_source_id_text"),
            ("manual_lifestyle_source_locator", "manual_lifestyle_source_locator_text"),
            ("manual_lifestyle_source_citation", "manual_lifestyle_source_citation_text"),
        ]:
            if text_col in manual_lifestyle.columns:
                manual_lifestyle[base_col] = manual_lifestyle[[base_col, text_col]].apply(
                    lambda row: join_unique(pd.Series(row.values)),
                    axis=1,
                )
                manual_lifestyle = manual_lifestyle.drop(columns=[text_col])
    manual_elevation = numeric_trait_block(
        ["elevation_min_m", "elevation_max_m", "elevation_mid_m"],
        {
            "elevation_min_m": "manual_elevation_min_m",
            "elevation_max_m": "manual_elevation_max_m",
            "elevation_mid_m": "manual_elevation_mid_m",
        },
        "manual_elevation",
    )

    merged = pd.DataFrame({"species": sorted(raw["species"].dropna().unique())})
    for block in [manual_size, manual_development, manual_lifestyle, manual_elevation]:
        if not block.empty:
            merged = merged.merge(block, on="species", how="left")
    return merged


def derive_lifestyle_from_amphibio(row: pd.Series) -> tuple[float | None, str | None, str | None]:
    aqu = row.get("amphibio_aqu")
    ter = row.get("amphibio_ter")
    development = row.get("amphibio_development_mode")

    if pd.isna(aqu) and pd.isna(ter):
        return (None, None, None)
    if aqu == 1 and ter != 1:
        return (2.0, "stream_aquatic", "amphibio_aquatic_flag")
    if ter == 1 and aqu == 1:
        return (1.0, "streamside", "amphibio_mixed_aquatic_terrestrial_flags")
    if ter == 1 and aqu != 1:
        if development == "direct_development":
            return (0.0, "mountain_woodland", "amphibio_terrestrial_flag_direct_development")
        return (1.0, "streamside", "amphibio_terrestrial_flag")
    return (None, None, None)


def choose_body_size(row: pd.Series) -> tuple[float | None, str | None, str | None, str | None]:
    nc_phrase_type = row.get("nc_size_phrase_type")
    has_conanti_proxy = pd.notna(row.get("svl_q90_all_mm"))

    if pd.notna(row.get("manual_adult_svl_mid_mm")):
        return (
            float(row["manual_adult_svl_mid_mm"]),
            "manual_adult_svl_range_mid_mm",
            row.get("manual_body_size_source_id"),
            "high",
        )
    if pd.notna(row.get("manual_adult_svl_max_mm")):
        return (
            float(row["manual_adult_svl_max_mm"]),
            "manual_adult_svl_max_mm",
            row.get("manual_body_size_source_id"),
            "medium",
        )
    if pd.notna(row.get("nc_adult_svl_mid_mm")) and nc_phrase_type == "adult_range":
        return (
            float(row["nc_adult_svl_mid_mm"]),
            "adult_svl_range_mid_mm",
            "nc_biodiversity_amphibians",
            "high",
        )
    if pd.notna(row.get("nc_adult_svl_mid_mm")) and nc_phrase_type == "adult_sex_specific_range":
        return (
            float(row["nc_adult_svl_mid_mm"]),
            "adult_svl_sex_specific_range_mid_mm",
            "nc_biodiversity_amphibians",
            "high",
        )
    if pd.notna(row.get("nc_adult_svl_max_mm")) and nc_phrase_type == "adult_max":
        return (
            float(row["nc_adult_svl_max_mm"]),
            "adult_svl_max_mm",
            "nc_biodiversity_amphibians",
            "medium",
        )
    if pd.notna(row.get("svl_q90_all_mm")):
        return (
            float(row["svl_q90_all_mm"]),
            "svl_q90_all_specimens_mixed_life_stages_mm",
            row.get("source_id"),
            "medium",
        )
    if pd.notna(row.get("manual_adult_tl_mid_mm")):
        return (
            float(row["manual_adult_tl_mid_mm"]),
            "manual_adult_total_length_range_mid_mm",
            row.get("manual_body_size_source_id"),
            "low",
        )
    if pd.notna(row.get("manual_adult_tl_max_mm")):
        return (
            float(row["manual_adult_tl_max_mm"]),
            "manual_adult_total_length_max_mm",
            row.get("manual_body_size_source_id"),
            "low",
        )
    if (
        pd.notna(row.get("nc_adult_tl_min_mm"))
        and pd.notna(row.get("nc_adult_tl_max_mm"))
        and pd.isna(row.get("nc_adult_svl_mid_mm"))
    ):
        return (
            float((row["nc_adult_tl_min_mm"] + row["nc_adult_tl_max_mm"]) / 2.0),
            "adult_total_length_range_mid_mm",
            "nc_biodiversity_amphibians",
            "low",
        )
    if pd.notna(row.get("nc_adult_tl_max_mm")) and pd.isna(row.get("nc_adult_svl_mid_mm")):
        return (
            float(row["nc_adult_tl_max_mm"]),
            "adult_total_length_max_mm",
            "nc_biodiversity_amphibians",
            "low",
        )
    if pd.notna(row.get("nc_adult_svl_mid_mm")) and nc_phrase_type == "description_range_unspecified":
        return (
            float(row["nc_adult_svl_mid_mm"]),
            "species_description_svl_range_mid_mm",
            "nc_biodiversity_amphibians",
            "medium",
        )
    if pd.notna(row.get("nc_adult_svl_mid_mm")) and nc_phrase_type == "metamorphosed_range" and not has_conanti_proxy:
        return (
            float(row["nc_adult_svl_mid_mm"]),
            "metamorphosed_svl_range_mid_mm",
            "nc_biodiversity_amphibians",
            "medium",
        )
    if pd.notna(row.get("amphibio_body_size_mm")):
        return (
            float(row["amphibio_body_size_mm"]),
            "amphibio_body_size_mm_measurement",
            row.get("amphibio_source_id"),
            "low",
        )
    if pd.notna(row.get("nc_adult_svl_mid_mm")) and nc_phrase_type == "mixed_stage_range" and not has_conanti_proxy:
        return (
            float(row["nc_adult_svl_mid_mm"]),
            "transformed_or_mixed_stage_svl_range_mid_mm",
            "nc_biodiversity_amphibians",
            "low",
        )
    return (None, None, None, None)


def choose_development(row: pd.Series) -> tuple[str | None, str | None, str | None]:
    value = row.get("manual_development_mode")
    if pd.notna(value):
        return (str(value), row.get("manual_development_source_id"), "high")
    value = row.get("nc_development_mode")
    if pd.notna(value):
        return (str(value), row.get("nc_source_id"), "medium")
    value = row.get("amphibio_development_mode")
    if pd.notna(value) and value != "mixed_or_unclear":
        return (str(value), row.get("amphibio_source_id"), "medium")
    return (None, None, None)


def choose_lifestyle(row: pd.Series) -> tuple[float | None, str | None, str | None, str | None]:
    if pd.notna(row.get("manual_aquaticity_index")) and pd.notna(row.get("manual_microhabitat_class")):
        return (
            float(row["manual_aquaticity_index"]),
            str(row["manual_microhabitat_class"]),
            row.get("manual_lifestyle_source_id"),
            "high",
        )
    if pd.notna(row.get("nc_aquaticity_index")) and pd.notna(row.get("nc_microhabitat_class")):
        return (
            float(row["nc_aquaticity_index"]),
            str(row["nc_microhabitat_class"]),
            row.get("nc_source_id"),
            "medium",
        )

    aquaticity_index, microhabitat_class, _ = derive_lifestyle_from_amphibio(row)
    if aquaticity_index is not None and microhabitat_class is not None:
        return (aquaticity_index, microhabitat_class, row.get("amphibio_source_id"), "low")

    return (None, None, None, None)


def choose_elevation(row: pd.Series) -> tuple[float | None, float | None, float | None, str | None]:
    if pd.notna(row.get("manual_elevation_min_m")) or pd.notna(row.get("manual_elevation_max_m")):
        low = row.get("manual_elevation_min_m")
        high = row.get("manual_elevation_max_m")
        midpoint = row.get("manual_elevation_mid_m")
        return (
            float(low) if pd.notna(low) else None,
            float(high) if pd.notna(high) else None,
            float(midpoint) if pd.notna(midpoint) else None,
            row.get("manual_elevation_source_id"),
        )
    if pd.notna(row.get("nc_elevation_min_m")) or pd.notna(row.get("nc_elevation_max_m")):
        low = row.get("nc_elevation_min_m")
        high = row.get("nc_elevation_max_m")
        midpoint = row.get("nc_elevation_mid_m")
        return (
            float(low) if pd.notna(low) else None,
            float(high) if pd.notna(high) else None,
            float(midpoint) if pd.notna(midpoint) else None,
            row.get("nc_source_id"),
        )
    return (None, None, None, None)


def summarize_coverage(df: pd.DataFrame) -> pd.DataFrame:
    datasets = {
        "all_species": df["species"].notna(),
        "te_genome_overlap": df["has_tree_tip"] & df["has_te"] & df["has_genome"],
        "te_genome_morph_overlap": df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_morphology"],
    }
    traits = [
        "body_size_proxy_mm",
        "adult_svl_mid_mm",
        "development_mode",
        "aquaticity_index",
        "microhabitat_class",
        "elevation_mid_m",
    ]

    rows = []
    for trait in traits:
        row = {"trait_name": trait}
        for dataset_name, mask in datasets.items():
            row[f"n_{dataset_name}"] = int(df.loc[mask, trait].notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    master = pd.read_csv(MASTER_FILE)
    amphibio = load_optional_csv(AMPHIBIO_FILE)
    conanti = load_optional_csv(CONANTI_FILE)
    nc = load_optional_csv(NC_FILE)
    manual = load_manual_trait_summary(MANUAL_EXTRACTION_FILE)

    merged = master[["species", "has_tree_tip", "has_te", "has_genome", "has_morphology"]].copy()
    merged = merged.merge(amphibio, on="species", how="left")
    merged = merged.merge(conanti, on="species", how="left")
    merged = merged.merge(nc, on="species", how="left")
    merged = merged.merge(manual, on="species", how="left")

    for column in [
        "source_id",
        "source_filename",
        "amphibio_source_id",
        "amphibio_source_filename",
        "nc_source_id",
        "nc_source_filename",
        "manual_body_size_source_id",
        "manual_body_size_source_locator",
        "manual_body_size_source_citation",
        "manual_development_mode",
        "manual_development_source_id",
        "manual_development_source_locator",
        "manual_development_source_citation",
        "manual_aquaticity_index",
        "manual_microhabitat_class",
        "manual_lifestyle_source_id",
        "manual_lifestyle_source_locator",
        "manual_lifestyle_source_citation",
        "manual_elevation_min_m",
        "manual_elevation_max_m",
        "manual_elevation_mid_m",
        "manual_elevation_source_id",
        "manual_elevation_source_locator",
        "manual_elevation_source_citation",
    ]:
        if column not in merged.columns:
            merged[column] = np.nan

    body_size = merged.apply(choose_body_size, axis=1, result_type="expand")
    body_size.columns = [
        "body_size_proxy_mm",
        "body_size_proxy_measurement",
        "body_size_proxy_source_id",
        "body_size_proxy_confidence",
    ]
    merged = pd.concat([merged, body_size], axis=1)

    development = merged.apply(choose_development, axis=1, result_type="expand")
    development.columns = ["development_mode", "development_source_id", "development_confidence"]
    merged = pd.concat([merged, development], axis=1)

    lifestyle = merged.apply(choose_lifestyle, axis=1, result_type="expand")
    lifestyle.columns = [
        "aquaticity_index",
        "microhabitat_class",
        "lifestyle_source_id",
        "lifestyle_confidence",
    ]
    merged = pd.concat([merged, lifestyle], axis=1)

    elevation = merged.apply(choose_elevation, axis=1, result_type="expand")
    elevation.columns = ["elevation_min_m", "elevation_max_m", "elevation_mid_m", "elevation_source_id"]
    merged = pd.concat([merged, elevation], axis=1)

    merged["adult_svl_min_mm"] = merged["nc_adult_svl_min_mm"]
    merged["adult_svl_max_mm"] = merged["nc_adult_svl_max_mm"]
    merged["adult_svl_mid_mm"] = merged["nc_adult_svl_mid_mm"]
    merged["adult_total_length_min_mm"] = merged["nc_adult_tl_min_mm"]
    merged["adult_total_length_max_mm"] = merged["nc_adult_tl_max_mm"]
    merged["adult_total_length_mid_mm"] = np.where(
        merged["nc_adult_tl_min_mm"].notna() & merged["nc_adult_tl_max_mm"].notna(),
        (merged["nc_adult_tl_min_mm"] + merged["nc_adult_tl_max_mm"]) / 2.0,
        np.nan,
    )
    manual_mask = merged["body_size_proxy_source_id"] == merged["manual_body_size_source_id"]
    merged.loc[manual_mask, "adult_svl_min_mm"] = merged.loc[manual_mask, "manual_adult_svl_min_mm"]
    merged.loc[manual_mask, "adult_svl_max_mm"] = merged.loc[manual_mask, "manual_adult_svl_max_mm"]
    merged.loc[manual_mask, "adult_svl_mid_mm"] = merged.loc[manual_mask, "manual_adult_svl_mid_mm"]
    merged.loc[manual_mask, "adult_total_length_min_mm"] = merged.loc[manual_mask, "manual_adult_tl_min_mm"]
    merged.loc[manual_mask, "adult_total_length_max_mm"] = merged.loc[manual_mask, "manual_adult_tl_max_mm"]
    merged.loc[manual_mask, "adult_total_length_mid_mm"] = merged.loc[manual_mask, "manual_adult_tl_mid_mm"]

    def resolve_body_size_detail(row: pd.Series) -> str | float:
        if row.get("body_size_proxy_source_id") == row.get("manual_body_size_source_id"):
            return row.get("manual_body_size_source_locator") or row.get("manual_body_size_source_citation")
        if row.get("body_size_proxy_source_id") == "nc_biodiversity_amphibians":
            return row.get("nc_source_filename")
        if row.get("body_size_proxy_source_id") == row.get("source_id"):
            return row.get("source_filename")
        if row.get("body_size_proxy_source_id") == row.get("amphibio_source_id"):
            return row.get("amphibio_source_filename")
        return np.nan

    merged["body_size_source_detail"] = merged.apply(resolve_body_size_detail, axis=1)
    merged["development_source_detail"] = np.where(
        merged["development_source_id"] == merged["manual_development_source_id"],
        merged["manual_development_source_locator"],
        np.where(
            merged["development_source_id"] == "nc_biodiversity_amphibians",
            merged["nc_source_filename"],
            np.where(
                merged["development_source_id"] == merged["amphibio_source_id"],
                merged["amphibio_source_filename"],
                np.nan,
            ),
        ),
    )
    merged["lifestyle_source_detail"] = np.where(
        merged["lifestyle_source_id"] == merged["manual_lifestyle_source_id"],
        merged["manual_lifestyle_source_locator"],
        np.where(
            merged["lifestyle_source_id"] == "nc_biodiversity_amphibians",
            merged["nc_source_filename"],
            np.where(
                merged["lifestyle_source_id"] == merged["amphibio_source_id"],
                merged["amphibio_source_filename"],
                np.nan,
            ),
        ),
    )
    merged["elevation_source_detail"] = np.where(
        merged["elevation_source_id"] == merged["manual_elevation_source_id"],
        merged["manual_elevation_source_locator"],
        np.where(
            merged["elevation_source_id"] == "nc_biodiversity_amphibians",
            merged["nc_source_filename"],
            np.nan,
        ),
    )

    def combine_source_ids(row: pd.Series) -> str | None:
        ids = [
            row.get("body_size_proxy_source_id"),
            row.get("development_source_id"),
            row.get("lifestyle_source_id"),
            row.get("elevation_source_id"),
        ]
        cleaned = sorted({str(value) for value in ids if pd.notna(value) and value})
        return ";".join(cleaned) if cleaned else None

    merged["organismal_source_ids"] = merged.apply(combine_source_ids, axis=1)

    keep_cols = [
        "species",
        "adult_svl_min_mm",
        "adult_svl_max_mm",
        "adult_svl_mid_mm",
        "adult_total_length_min_mm",
        "adult_total_length_max_mm",
        "adult_total_length_mid_mm",
        "body_size_proxy_mm",
        "body_size_proxy_measurement",
        "body_size_proxy_source_id",
        "body_size_source_detail",
        "body_size_proxy_confidence",
        "development_mode",
        "development_source_id",
        "development_source_detail",
        "development_confidence",
        "aquaticity_index",
        "microhabitat_class",
        "lifestyle_source_id",
        "lifestyle_source_detail",
        "lifestyle_confidence",
        "elevation_min_m",
        "elevation_max_m",
        "elevation_mid_m",
        "elevation_source_id",
        "elevation_source_detail",
        "organismal_source_ids",
    ]
    curated = merged[keep_cols].sort_values("species").reset_index(drop=True)
    coverage = summarize_coverage(merged)

    curated.to_csv(OUTPUT_FILE, index=False)
    coverage.to_csv(OUTPUT_COVERAGE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Wrote {OUTPUT_COVERAGE}")


if __name__ == "__main__":
    main()

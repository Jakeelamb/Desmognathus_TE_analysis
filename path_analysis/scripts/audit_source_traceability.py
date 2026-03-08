#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"
TEMPLATES_DIR = PROJECT_ROOT / "path_analysis" / "data" / "templates"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
EXTERNAL_RAW_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "raw"
EXTERNAL_DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

SOURCE_MANIFEST = TEMPLATES_DIR / "source_manifest.csv"
RAW_INVENTORY = EXTERNAL_RAW_DIR / "source_file_inventory.csv"

OUTPUT_USAGE = DERIVED_DIR / "source_usage_summary.csv"
OUTPUT_GAPS = DERIVED_DIR / "source_traceability_gaps.csv"
OUTPUT_FILE_REGISTRY = DERIVED_DIR / "source_file_registry.csv"

TABLES_TO_SCAN = [
    TEMPLATES_DIR / "literature_trait_extraction.csv",
    EXTERNAL_DERIVED_DIR / "amphibio_desmognathus_traits.csv",
    EXTERNAL_DERIVED_DIR / "conanti_fuscus_morphometrics_summary.csv",
    EXTERNAL_DERIVED_DIR / "nc_biodiversity_desmognathus_traits.csv",
    DERIVED_DIR / "master_species_table.csv",
    DERIVED_DIR / "te_path_features.csv",
    DERIVED_DIR / "te_model_feature_panel.csv",
    DERIVED_DIR / "organismal_traits_curated.csv",
    DERIVED_DIR / "path_input_master.csv",
]

TRACKED_FILES = [
    ("repo_dnapipete_order_breakdown", PROJECT_ROOT / "results" / "data" / "dnaPipeTE_order_breakdown.csv"),
    ("repo_dnapipete_superfamily_breakdown", PROJECT_ROOT / "results" / "data" / "dnaPipeTE_superfamily_breakdown.csv"),
    ("repo_diversity_order_stats", PROJECT_ROOT / "results" / "data" / "diversity_order_stats.csv"),
    ("repo_diversity_superfamily_stats", PROJECT_ROOT / "results" / "data" / "diversity_superfamily_stats.csv"),
    (
        "repo_divergence_summary_statistics_by_species",
        PROJECT_ROOT / "results" / "data" / "divergence" / "divergence_summary_statistics_by_species.csv",
    ),
    (
        "repo_ectopic_recombination_filtered_3000bp_5plusdomains",
        PROJECT_ROOT / "results" / "data" / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv",
    ),
    ("repo_desmognathus_phylogeny_tree", PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"),
    (
        "cellprofiler_final_species_results",
        CELLPROFILER_ROOT / "output" / "qc_report_blockbalanced" / "final_species_results.csv",
    ),
    (
        "cellprofiler_species_morphology_summary",
        CELLPROFILER_ROOT / "output" / "publication_analysis" / "species_morphology_summary.csv",
    ),
]


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_file_registry() -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    if RAW_INVENTORY.exists():
        raw_inventory = pd.read_csv(RAW_INVENTORY)
        for _, row in raw_inventory.iterrows():
            local_rel = row.get("local_filename")
            local_path = (EXTERNAL_RAW_DIR / str(local_rel)).resolve() if pd.notna(local_rel) else pd.NA
            rows.append(
                {
                    "source_id": row.get("source_id"),
                    "local_path": str(local_path) if pd.notna(local_path) else pd.NA,
                    "sha256": row.get("sha256"),
                    "exists": pd.notna(local_path) and Path(local_path).exists(),
                    "origin": "external_raw_inventory",
                }
            )

    for source_id, path in TRACKED_FILES:
        rows.append(
            {
                "source_id": source_id,
                "local_path": str(path.resolve()),
                "sha256": sha256_for_file(path) if path.exists() else pd.NA,
                "exists": path.exists(),
                "origin": "tracked_input",
            }
        )

    registry = pd.DataFrame(rows).drop_duplicates(subset=["source_id", "local_path"]).reset_index(drop=True)
    return registry


def extract_used_source_ids(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["table_name", "column_name", "species", "source_id"])
    df = pd.read_csv(path)
    source_cols = [
        col
        for col in df.columns
        if col == "source_id" or col.endswith("_source_id") or col.endswith("_source_ids")
    ]
    if not source_cols:
        return pd.DataFrame(columns=["table_name", "column_name", "species", "source_id"])

    species_col = "species" if "species" in df.columns else ("current_species" if "current_species" in df.columns else None)
    rows: list[dict[str, object]] = []
    for _, row in df.iterrows():
        for col in source_cols:
            value = row.get(col)
            if pd.isna(value):
                continue
            for source_id in [part.strip() for part in str(value).split(";") if part.strip()]:
                rows.append(
                    {
                        "table_name": path.name,
                        "column_name": col,
                        "species": row.get(species_col) if species_col else pd.NA,
                        "source_id": source_id,
                    }
                )
    return pd.DataFrame(rows)


def build_usage_summary(manifest: pd.DataFrame, registry: pd.DataFrame) -> pd.DataFrame:
    usage_frames = [extract_used_source_ids(path) for path in TABLES_TO_SCAN]
    used = pd.concat(usage_frames, ignore_index=True) if usage_frames else pd.DataFrame(columns=["source_id"])
    if used.empty:
        return pd.DataFrame(columns=["source_id"])

    summary = (
        used.groupby("source_id")
        .agg(
            n_occurrences=("source_id", "size"),
            n_tables=("table_name", "nunique"),
            tables=("table_name", lambda s: ";".join(sorted(set(s)))),
            columns=("column_name", lambda s: ";".join(sorted(set(s)))),
        )
        .reset_index()
    )
    summary["in_manifest"] = summary["source_id"].isin(set(manifest["source_id"]))
    summary["has_local_file_record"] = summary["source_id"].isin(set(registry["source_id"]))
    summary = summary.merge(
        manifest[["source_id", "short_citation", "full_reference", "url"]],
        on="source_id",
        how="left",
    )
    return summary.sort_values(["in_manifest", "source_id"], ascending=[True, True]).reset_index(drop=True)


def build_gap_report(
    manifest: pd.DataFrame,
    registry: pd.DataFrame,
    usage: pd.DataFrame,
) -> pd.DataFrame:
    gaps: list[dict[str, object]] = []
    manifest_by_id = manifest.set_index("source_id")

    for _, row in usage.iterrows():
        if not bool(row["in_manifest"]):
            gaps.append(
                {
                    "gap_type": "missing_manifest_entry",
                    "source_id": row["source_id"],
                    "table_name": row["tables"],
                    "species": pd.NA,
                    "details": "Source id is used in derived tables but absent from source_manifest.csv",
                }
            )
        source_type = (
            manifest_by_id.loc[row["source_id"], "source_type"]
            if row["source_id"] in manifest_by_id.index
            else pd.NA
        )
        local_file_expected = source_type in {"website", "dataset", "local_repo_table", "local_generated_table"}
        if local_file_expected and not bool(row["has_local_file_record"]):
            gaps.append(
                {
                    "gap_type": "missing_local_file_record",
                    "source_id": row["source_id"],
                    "table_name": row["tables"],
                    "species": pd.NA,
                    "details": "Source id is used in derived tables but has no tracked local file or registry entry",
                }
            )

    if (DERIVED_DIR / "organismal_traits_curated.csv").exists():
        df = pd.read_csv(DERIVED_DIR / "organismal_traits_curated.csv")
        block_checks = [
            ("body_size_proxy_mm", "body_size_proxy_source_id"),
            ("development_mode", "development_source_id"),
            ("aquaticity_index", "lifestyle_source_id"),
            ("microhabitat_class", "lifestyle_source_id"),
            ("elevation_mid_m", "elevation_source_id"),
        ]
        for value_col, source_col in block_checks:
            missing = df[df[value_col].notna() & df[source_col].isna()]
            for _, row in missing.iterrows():
                gaps.append(
                    {
                        "gap_type": "missing_block_source_id",
                        "source_id": pd.NA,
                        "table_name": "organismal_traits_curated.csv",
                        "species": row["species"],
                        "details": f"{value_col} is populated but {source_col} is missing",
                    }
                )

    if (DERIVED_DIR / "te_path_features.csv").exists():
        df = pd.read_csv(DERIVED_DIR / "te_path_features.csv")
        required = [
            "te_order_source_id",
            "order_diversity_source_id",
            "divergence_source_id",
            "ectopic_source_id",
            "te_feature_source_ids",
        ]
        for column in required:
            if column not in df.columns or df[column].isna().any():
                gaps.append(
                    {
                        "gap_type": "missing_te_provenance_column",
                        "source_id": pd.NA,
                        "table_name": "te_path_features.csv",
                        "species": pd.NA,
                        "details": f"Missing or incomplete TE provenance column: {column}",
                    }
                )

    if (DERIVED_DIR / "master_species_table.csv").exists():
        df = pd.read_csv(DERIVED_DIR / "master_species_table.csv")
        for value_col, source_col in [
            ("genome_size_pg", "genome_source_id"),
            ("morph_cell_area_um2", "morphology_source_id"),
            ("ectopic_n_elements", "ectopic_source_id"),
        ]:
            if source_col not in df.columns:
                gaps.append(
                    {
                        "gap_type": "missing_master_provenance_column",
                        "source_id": pd.NA,
                        "table_name": "master_species_table.csv",
                        "species": pd.NA,
                        "details": f"Missing master-table provenance column: {source_col}",
                    }
                )
                continue
            missing = df[df[value_col].notna() & df[source_col].isna()]
            for _, row in missing.iterrows():
                gaps.append(
                    {
                        "gap_type": "missing_master_block_source_id",
                        "source_id": pd.NA,
                        "table_name": "master_species_table.csv",
                        "species": row["species"],
                        "details": f"{value_col} is populated but {source_col} is missing",
                    }
                )

    columns = ["gap_type", "source_id", "table_name", "species", "details"]
    if not gaps:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(gaps, columns=columns)


def main() -> None:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(SOURCE_MANIFEST)
    registry = load_file_registry()
    usage = build_usage_summary(manifest, registry)
    gaps = build_gap_report(manifest, registry, usage)

    registry.to_csv(OUTPUT_FILE_REGISTRY, index=False)
    usage.to_csv(OUTPUT_USAGE, index=False)
    gaps.to_csv(OUTPUT_GAPS, index=False)

    print(f"Wrote {OUTPUT_FILE_REGISTRY}")
    print(f"Wrote {OUTPUT_USAGE}")
    print(f"Wrote {OUTPUT_GAPS}")


if __name__ == "__main__":
    main()

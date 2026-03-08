#!/usr/bin/env python3
"""Build consolidated path-analysis datasets for the Desmognathus chapter."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


DEFAULT_CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_project_root(start: Path | None = None) -> Path:
    current = (start or Path(__file__).resolve()).parent
    for _ in range(8):
        if (current / "paths.yaml").exists():
            return current
        if current.parent == current:
            break
        current = current.parent
    raise FileNotFoundError("Could not locate project root containing paths.yaml")


def canonical_species(value: object) -> str | pd.NA:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return pd.NA
    text = str(value).strip()
    if not text:
        return pd.NA
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def slugify_column(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9]+", "_", name.strip())
    name = re.sub(r"_+", "_", name).strip("_")
    return name.lower()


def parse_tree_tips(tree_path: Path) -> set[str]:
    text = tree_path.read_text()
    tips = re.findall(r"(?<=[(,])([^():,]+?)(?=[:),])", text)
    return {canonical_species(tip) for tip in tips if canonical_species(tip) is not pd.NA}


def load_prefixed_table(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    species_col = df.columns[0]
    df = df.rename(columns={species_col: "species"})
    df["species"] = df["species"].map(canonical_species)
    rename_map = {
        col: f"{prefix}_{slugify_column(col)}"
        for col in df.columns
        if col != "species"
    }
    return (
        df.rename(columns=rename_map)
        .dropna(subset=["species"])
        .drop_duplicates(subset=["species"])
        .reset_index(drop=True)
    )


def load_diversity_table(path: Path, prefix: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    species_col = df.columns[0]
    df = df.rename(columns={species_col: "species"})
    df["species"] = df["species"].map(canonical_species)

    keep_map = {
        "Simpson_Diversity": f"{prefix}_simpson",
        "Shannon_Diversity": f"{prefix}_shannon",
        "Pielou_Evenness": f"{prefix}_pielou",
    }
    present = [col for col in keep_map if col in df.columns]
    out = df[["species"] + present].rename(columns={col: keep_map[col] for col in present})
    return (
        out.dropna(subset=["species"])
        .drop_duplicates(subset=["species"])
        .reset_index(drop=True)
    )


def load_ectopic_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["species"] = df["species"].map(canonical_species)
    complete_state = df["Complete"].astype(str).str.strip().str.lower()
    df["complete_flag_known"] = complete_state.map({"yes": 1.0, "no": 0.0})
    numeric_cols = [
        "ratio_terminal_internal",
        "mean_depth_terminal",
        "mean_depth_internal",
        "domain_count",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    out = (
        df.groupby("species", dropna=True)
        .agg(
            ectopic_n_rows_total=("sequence", "size"),
            ectopic_n_elements=("ratio_terminal_internal", "count"),
            ectopic_mean_ratio=("ratio_terminal_internal", "mean"),
            ectopic_median_ratio=("ratio_terminal_internal", "median"),
            ectopic_mean_terminal_depth=("mean_depth_terminal", "mean"),
            ectopic_mean_internal_depth=("mean_depth_internal", "mean"),
            ectopic_mean_domain_count=("domain_count", "mean"),
            ectopic_n_complete_known=("complete_flag_known", "count"),
            ectopic_complete_fraction=("complete_flag_known", "mean"),
        )
        .reset_index()
    )
    return out


def load_genome_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["species"] = df["species"].map(canonical_species)
    keep = {
        "primary_genome_pg": "genome_size_pg",
        "primary_genome_se_pg": "genome_size_se_pg",
        "primary_genome_gb": "genome_size_gb",
        "primary_n_images": "genome_n_images",
        "primary_n_specimens": "genome_n_specimens",
        "primary_n_nuclei": "genome_n_nuclei",
        "primary_cv_area_pct": "genome_cv_area_pct",
        "result_status": "genome_result_status",
        "flag_summary": "genome_flag_summary",
    }
    present = [col for col in keep if col in df.columns]
    out = df[["species"] + present].rename(columns={col: keep[col] for col in present})
    return (
        out.dropna(subset=["species"])
        .drop_duplicates(subset=["species"])
        .reset_index(drop=True)
    )


def load_morphology_results(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["species"] = df["species"].map(canonical_species)
    keep = {
        "n_specimens_strict": "morph_n_specimens",
        "n_images_strict": "morph_n_images",
        "n_analysis_ready_images": "morph_n_ready_images",
        "n_pairs_strict": "morph_n_pairs",
        "species_median_cell_area_um2": "morph_cell_area_um2",
        "species_median_nuc_area_um2": "morph_nucleus_area_um2",
        "species_median_nc_ratio": "morph_nc_ratio",
        "species_median_cytoplasm_area_um2": "morph_cytoplasm_area_um2",
        "low_support_for_species_median": "morph_low_support",
    }
    present = [col for col in keep if col in df.columns]
    out = df[["species"] + present].rename(columns={col: keep[col] for col in present})
    return (
        out.dropna(subset=["species"])
        .drop_duplicates(subset=["species"])
        .reset_index(drop=True)
    )


def build_master_table(project_root: Path, cellprofiler_root: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    te_order_path = project_root / "results" / "data" / "dnaPipeTE_order_breakdown.csv"
    te_superfamily_path = project_root / "results" / "data" / "dnaPipeTE_superfamily_breakdown.csv"
    order_diversity_path = project_root / "results" / "data" / "diversity_order_stats.csv"
    superfamily_diversity_path = project_root / "results" / "data" / "diversity_superfamily_stats.csv"
    ectopic_path = project_root / "results" / "data" / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"
    genome_path = cellprofiler_root / "output" / "qc_report_blockbalanced" / "final_species_results.csv"
    morphology_path = cellprofiler_root / "output" / "publication_analysis" / "species_morphology_summary.csv"

    te_order = load_prefixed_table(te_order_path, "order")
    te_order["te_order_source_id"] = "repo_dnapipete_order_breakdown"
    te_order["te_order_source_path"] = str(te_order_path.resolve())
    te_order["te_order_source_sha256"] = sha256_for_file(te_order_path)

    te_superfamily = load_prefixed_table(
        te_superfamily_path,
        "superfamily",
    )
    te_superfamily["te_superfamily_source_id"] = "repo_dnapipete_superfamily_breakdown"
    te_superfamily["te_superfamily_source_path"] = str(te_superfamily_path.resolve())
    te_superfamily["te_superfamily_source_sha256"] = sha256_for_file(te_superfamily_path)

    div_order = load_diversity_table(order_diversity_path, "order")
    div_order["order_diversity_source_id"] = "repo_diversity_order_stats"
    div_order["order_diversity_source_path"] = str(order_diversity_path.resolve())
    div_order["order_diversity_source_sha256"] = sha256_for_file(order_diversity_path)

    div_superfamily = load_diversity_table(
        superfamily_diversity_path,
        "superfamily",
    )
    div_superfamily["superfamily_diversity_source_id"] = "repo_diversity_superfamily_stats"
    div_superfamily["superfamily_diversity_source_path"] = str(superfamily_diversity_path.resolve())
    div_superfamily["superfamily_diversity_source_sha256"] = sha256_for_file(superfamily_diversity_path)

    ectopic = load_ectopic_summary(ectopic_path)
    ectopic["ectopic_source_id"] = "repo_ectopic_recombination_filtered_3000bp_5plusdomains"
    ectopic["ectopic_source_path"] = str(ectopic_path.resolve())
    ectopic["ectopic_source_sha256"] = sha256_for_file(ectopic_path)

    genome = load_genome_results(genome_path)
    genome["genome_source_id"] = "cellprofiler_final_species_results"
    genome["genome_source_path"] = str(genome_path.resolve())
    genome["genome_source_sha256"] = sha256_for_file(genome_path)

    morphology = load_morphology_results(morphology_path)
    morphology["morphology_source_id"] = "cellprofiler_species_morphology_summary"
    morphology["morphology_source_path"] = str(morphology_path.resolve())
    morphology["morphology_source_sha256"] = sha256_for_file(morphology_path)

    tree_tips = parse_tree_tips(project_root / "input_data" / "phylogeny" / "desmo900dated_test.tre")

    sources = {
        "te_order": te_order,
        "te_superfamily": te_superfamily,
        "order_diversity": div_order,
        "superfamily_diversity": div_superfamily,
        "ectopic": ectopic,
        "genome": genome,
        "morphology": morphology,
    }

    all_species = sorted(
        {
            species
            for frame in sources.values()
            for species in frame["species"].dropna().tolist()
        }
    )
    master = pd.DataFrame({"species": all_species})
    for name, frame in sources.items():
        master = master.merge(frame, on="species", how="left", validate="one_to_one")
        master[f"has_{name}"] = master["species"].isin(set(frame["species"]))

    master["has_tree_tip"] = master["species"].isin(tree_tips)
    master["has_te"] = master["has_te_order"] & master["has_te_superfamily"]
    master["has_genome"] = master["genome_size_pg"].notna()
    master["has_ectopic"] = master["ectopic_mean_ratio"].notna()
    master["has_morphology"] = (
        master["morph_cell_area_um2"].notna() & master["morph_nucleus_area_um2"].notna()
    )

    return master.sort_values("species").reset_index(drop=True), sources


def build_dataset_summary(master: pd.DataFrame) -> pd.DataFrame:
    dataset_masks = {
        "dataset_te_genome": master["has_tree_tip"] & master["has_te"] & master["has_genome"],
        "dataset_te_genome_ectopic": master["has_tree_tip"] & master["has_te"] & master["has_genome"] & master["has_ectopic"],
        "dataset_genome_morphology": master["has_tree_tip"] & master["has_genome"] & master["has_morphology"],
        "dataset_te_genome_morphology": master["has_tree_tip"] & master["has_te"] & master["has_genome"] & master["has_morphology"],
        "dataset_te_genome_ectopic_morphology": (
            master["has_tree_tip"] & master["has_te"] & master["has_genome"] & master["has_ectopic"] & master["has_morphology"]
        ),
    }

    rows = []
    for name, mask in dataset_masks.items():
        species = master.loc[mask, "species"].tolist()
        rows.append(
            {
                "dataset": name,
                "n_species": len(species),
                "species_list": ";".join(species),
            }
        )
    return pd.DataFrame(rows), dataset_masks


def write_outputs(master: pd.DataFrame, summary: pd.DataFrame, masks: dict[str, pd.Series], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    master.to_csv(out_dir / "master_species_table.csv", index=False)
    summary.to_csv(out_dir / "dataset_overlap_summary.csv", index=False)

    dataset_columns = ["species"] + [col for col in master.columns if col != "species"]
    for name, mask in masks.items():
        master.loc[mask, dataset_columns].sort_values("species").to_csv(out_dir / f"{name}.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=find_project_root(),
        help="Project root for Desmognathus_TE",
    )
    parser.add_argument(
        "--cellprofiler-root",
        type=Path,
        default=DEFAULT_CELLPROFILER_ROOT,
        help="Root of the cellprofiler_test repository",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for derived path-analysis tables",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    out_dir = args.out_dir or (project_root / "path_analysis" / "data" / "derived")

    master, _ = build_master_table(project_root, args.cellprofiler_root.resolve())
    summary, masks = build_dataset_summary(master)
    write_outputs(master, summary, masks, out_dir)

    payload = {
        row["dataset"]: int(row["n_species"])
        for _, row in summary.iterrows()
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

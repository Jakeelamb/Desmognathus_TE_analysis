#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DATA = PROJECT_ROOT / "results" / "data"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"

ORDER_FILE = RESULTS_DATA / "dnaPipeTE_order_breakdown.csv"
ORDER_DIVERSITY_FILE = RESULTS_DATA / "diversity_order_stats.csv"
DIVERGENCE_FILE = RESULTS_DATA / "divergence" / "divergence_summary_statistics_by_species.csv"
ECTOPIC_FILE = RESULTS_DATA / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"

OUTPUT_FILE = DERIVED_DIR / "te_path_features.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def standardize_species(series: pd.Series) -> pd.Series:
    cleaned = series.astype(str).str.strip()
    cleaned = cleaned.str.replace(r"^D\.\s*", "", regex=True)
    return cleaned


def clr_transform(df: pd.DataFrame) -> pd.DataFrame:
    positive_values = df.where(df > 0).stack()
    pseudocount = float(positive_values.min() / 2) if not positive_values.empty else 1e-6
    adjusted = df + pseudocount
    geometric_mean = np.exp(np.log(adjusted).mean(axis=1))
    clr = np.log(adjusted).sub(np.log(geometric_mean), axis=0)
    return clr


def pca_scores_2d(df: pd.DataFrame) -> pd.DataFrame:
    centered = df - df.mean(axis=0)
    matrix = centered.to_numpy()
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    scores = u[:, :2] * s[:2]
    result = pd.DataFrame(index=df.index)
    result["major_order_clr_pc1"] = scores[:, 0]
    result["major_order_clr_pc2"] = scores[:, 1] if scores.shape[1] > 1 else 0.0
    return result


def load_order_features() -> pd.DataFrame:
    order = pd.read_csv(ORDER_FILE, index_col=0)
    order.index = standardize_species(order.index.to_series())
    order.columns = [f"order_{c.lower()}" for c in order.columns]

    positive_values = order.where(order > 0).stack()
    pseudocount = float(positive_values.min() / 2) if not positive_values.empty else 1e-6

    feature_df = order.copy()
    feature_df["ltr_line_logratio"] = np.log(
        (feature_df["order_ltr"] + pseudocount) / (feature_df["order_line"] + pseudocount)
    )
    feature_df["ltr_tir_logratio"] = np.log(
        (feature_df["order_ltr"] + pseudocount) / (feature_df["order_tir"] + pseudocount)
    )

    retro = feature_df[["order_ltr", "order_line", "order_sine", "order_ple"]].sum(axis=1)
    dna = feature_df[["order_tir", "order_maverick", "order_helitron"]].sum(axis=1)
    feature_df["retro_dna_logratio"] = np.log((retro + pseudocount) / (dna + pseudocount))

    major_order_cols = ["order_ltr", "order_line", "order_tir", "order_dirs", "order_sine"]
    clr = clr_transform(feature_df[major_order_cols])
    pcs = pca_scores_2d(clr)
    feature_df = feature_df.join(pcs)

    return feature_df


def load_order_diversity() -> pd.DataFrame:
    diversity = pd.read_csv(ORDER_DIVERSITY_FILE)
    diversity = diversity.rename(columns={diversity.columns[0]: "species"})
    diversity["species"] = standardize_species(diversity["species"])
    diversity = diversity.rename(
        columns={
            "Simpson_Diversity": "order_simpson",
            "Shannon_Diversity": "order_shannon",
            "Pielou_Evenness": "order_pielou",
        }
    )
    return diversity.set_index("species")[["order_simpson", "order_shannon", "order_pielou"]]


def load_divergence_features(order_weights: pd.DataFrame) -> pd.DataFrame:
    divergence = pd.read_csv(DIVERGENCE_FILE)
    divergence = divergence[
        (divergence["group_level"] == "order") & (divergence["threshold"].astype(float) == 0.9)
    ].copy()
    divergence["species"] = standardize_species(divergence["Desmognathus_Species"])
    divergence["group_name_lower"] = divergence["group_name"].str.lower()

    metrics = {}
    for metric in [
        "percent_divergence_median",
        "percent_deletions_median",
        "percent_insertions_median",
    ]:
        wide = divergence.pivot(
            index="species",
            columns="group_name_lower",
            values=metric,
        )
        wide = wide.add_prefix(f"{metric}_")
        metrics[metric] = wide

    merged = pd.concat(metrics.values(), axis=1)

    weight_cols = [
        "order_dirs",
        "order_helitron",
        "order_line",
        "order_ltr",
        "order_maverick",
        "order_ple",
        "order_sine",
        "order_tir",
        "order_yr",
    ]
    weights = order_weights[weight_cols].copy()
    weights.columns = [c.replace("order_", "") for c in weights.columns]

    feature_rows = []
    for species, row in weights.iterrows():
        values = {"species": species}
        aligned_weights = row / row.sum()
        for metric, out_name in [
            ("percent_divergence_median", "weighted_te_divergence_p90"),
            ("percent_deletions_median", "weighted_te_deletions_p90"),
            ("percent_insertions_median", "weighted_te_insertions_p90"),
        ]:
            metric_cols = [f"{metric}_{name}" for name in aligned_weights.index]
            if species not in merged.index:
                values[out_name] = np.nan
                continue
            metric_values = merged.loc[species, metric_cols]
            values[out_name] = float((aligned_weights.to_numpy() * metric_values.to_numpy()).sum())

        for order_name in ["ltr", "line", "tir", "dirs", "sine"]:
            for metric, suffix in [
                ("percent_divergence_median", "divergence_p90"),
                ("percent_deletions_median", "deletions_p90"),
            ]:
                col = f"{metric}_{order_name}"
                values[f"{order_name}_{suffix}"] = (
                    float(merged.loc[species, col]) if species in merged.index and col in merged.columns else np.nan
                )

        feature_rows.append(values)

    return pd.DataFrame(feature_rows).set_index("species")


def load_ectopic_features() -> pd.DataFrame:
    ectopic = pd.read_csv(ECTOPIC_FILE, sep="\t")
    ectopic["species"] = standardize_species(ectopic["species"])
    ectopic["complete_yes"] = (ectopic["Complete"].astype(str).str.lower() == "yes").astype(float)

    grouped = ectopic.groupby("species", dropna=False)
    summary = grouped.agg(
        ectopic_n_elements=("ratio_terminal_internal", "size"),
        ectopic_mean_ratio=("ratio_terminal_internal", "mean"),
        ectopic_median_ratio=("ratio_terminal_internal", "median"),
        ectopic_mean_terminal_depth=("mean_depth_terminal", "mean"),
        ectopic_mean_internal_depth=("mean_depth_internal", "mean"),
        ectopic_mean_domain_count=("domain_count", "mean"),
        ectopic_complete_fraction=("complete_yes", "mean"),
    )
    summary["ectopic_log10_mean_ratio"] = summary["ectopic_mean_ratio"].apply(
        lambda x: math.log10(x) if pd.notna(x) and x > 0 else np.nan
    )
    return summary


def main() -> None:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)

    order = load_order_features()
    diversity = load_order_diversity()
    divergence = load_divergence_features(order)
    ectopic = load_ectopic_features()

    output = order.join(diversity, how="outer")
    output = output.join(divergence, how="outer")
    output = output.join(ectopic, how="outer")
    output = output.reset_index().rename(columns={"index": "species"})
    output["te_order_source_id"] = "repo_dnapipete_order_breakdown"
    output["te_order_source_path"] = str(ORDER_FILE.relative_to(PROJECT_ROOT))
    output["te_order_source_sha256"] = sha256_for_file(ORDER_FILE)
    output["order_diversity_source_id"] = "repo_diversity_order_stats"
    output["order_diversity_source_path"] = str(ORDER_DIVERSITY_FILE.relative_to(PROJECT_ROOT))
    output["order_diversity_source_sha256"] = sha256_for_file(ORDER_DIVERSITY_FILE)
    output["divergence_source_id"] = "repo_divergence_summary_statistics_by_species"
    output["divergence_source_path"] = str(DIVERGENCE_FILE.relative_to(PROJECT_ROOT))
    output["divergence_source_sha256"] = sha256_for_file(DIVERGENCE_FILE)
    output["ectopic_source_id"] = "repo_ectopic_recombination_filtered_3000bp_5plusdomains"
    output["ectopic_source_path"] = str(ECTOPIC_FILE.relative_to(PROJECT_ROOT))
    output["ectopic_source_sha256"] = sha256_for_file(ECTOPIC_FILE)
    output["te_feature_source_ids"] = (
        "repo_divergence_summary_statistics_by_species;"
        "repo_diversity_order_stats;"
        "repo_dnapipete_order_breakdown;"
        "repo_ectopic_recombination_filtered_3000bp_5plusdomains"
    )

    output.to_csv(OUTPUT_FILE, index=False)
    print(f"Wrote {OUTPUT_FILE}")
    print(f"Rows: {len(output)}")
    print(f"Columns: {len(output.columns)}")


if __name__ == "__main__":
    main()

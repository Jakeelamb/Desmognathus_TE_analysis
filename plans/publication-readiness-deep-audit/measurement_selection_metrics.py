#!/usr/bin/env python3
"""Read-only metrics for the publication-readiness measurement audit."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"
REVIEW_ROOT = (
    CELLPROFILER_ROOT
    / "output"
    / "runs"
    / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
    / "top50_linked_pair_review"
)
CANDIDATES_CSV = REVIEW_ROOT / "pair_quality_model" / "quality_scored_candidates.csv"
SELECTED_CSV = (
    REVIEW_ROOT
    / "final_curated_top50_latest"
    / "final_curated_top50_linked_pairs.csv"
)
BALANCED_ROOT = (
    Path(__file__).resolve().parents[2]
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "balanced_genome_iod_sensitivity"
)
BALANCED_PAIRS_CSV = BALANCED_ROOT / "balanced_genome_iod_selected_pairs.csv.gz"
BALANCED_GENOME_COMPARISON_CSV = BALANCED_ROOT / "balanced_genome_iod_comparison.csv"

METRICS = [
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_mean_od",
    "nuc_iod",
    "top50_score",
    "pair_quality_probability",
    "quality_rank_score",
]


def read_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(CANDIDATES_CSV, low_memory=False)
    selected = pd.read_csv(SELECTED_CSV, low_memory=False)
    for frame in (candidates, selected):
        frame["species"] = frame["species"].astype(str)
        frame["review_key"] = frame["review_key"].astype(str)
        for metric in METRICS:
            frame[metric] = pd.to_numeric(frame[metric], errors="coerce")
    return candidates, selected


def species_metrics(candidates: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    selected_keys = set(selected["review_key"])
    candidates = candidates.copy()
    candidates["selected"] = candidates["review_key"].isin(selected_keys)
    for metric in METRICS:
        candidates[f"{metric}_percentile"] = candidates.groupby("species")[metric].rank(
            pct=True,
            method="average",
        )

    rows: list[dict[str, object]] = []
    for species, group in candidates.groupby("species", sort=True):
        chosen = group[group["selected"]].copy()
        largest_50 = set(
            group.nlargest(min(50, len(group)), "cell_area_um2")["review_key"]
        )
        chosen_keys = set(chosen["review_key"])
        rows.append(
            {
                "species": species,
                "candidate_rows": len(group),
                "selected_rows": len(chosen),
                "selected_fraction": len(chosen) / len(group),
                "selected_images": chosen["filename"].nunique(),
                "selected_specimens": chosen["specimen_id"].nunique(),
                "overlap_with_largest_50": len(chosen_keys & largest_50),
                "selected_cell_area_median_um2": chosen["cell_area_um2"].median(),
                "candidate_cell_area_median_um2": group["cell_area_um2"].median(),
                "selected_to_candidate_cell_median_ratio": (
                    chosen["cell_area_um2"].median() / group["cell_area_um2"].median()
                ),
                "selected_cell_area_percentile_median": chosen[
                    "cell_area_um2_percentile"
                ].median(),
                "selected_cell_area_percentile_min": chosen[
                    "cell_area_um2_percentile"
                ].min(),
                "selected_nucleus_area_percentile_median": chosen[
                    "nuc_area_um2_percentile"
                ].median(),
                "selected_nucleus_od_percentile_median": chosen[
                    "nuc_mean_od_percentile"
                ].median(),
                "selected_nucleus_iod_percentile_median": chosen[
                    "nuc_iod_percentile"
                ].median(),
                "selected_top50_score_percentile_median": chosen[
                    "top50_score_percentile"
                ].median(),
                "selected_quality_score_percentile_median": chosen[
                    "quality_rank_score_percentile"
                ].median(),
            }
        )
    return pd.DataFrame(rows)


def print_distribution(selected: pd.DataFrame, column: str) -> None:
    if column not in selected:
        print(f"{column}: MISSING")
        return
    counts = selected[column].fillna("<NA>").astype(str).value_counts(dropna=False)
    print(f"{column}: {counts.to_dict()}")


def analysis_weights(group: pd.DataFrame) -> pd.Series:
    weights = pd.Series(0.0, index=group.index, dtype=float)
    specimens = list(dict.fromkeys(group["specimen_group"].fillna("").astype(str)))
    for specimen in specimens:
        specimen_mask = group["specimen_group"].fillna("").astype(str).eq(specimen)
        specimen_rows = group.loc[specimen_mask]
        images = list(dict.fromkeys(specimen_rows["filename"].fillna("").astype(str)))
        for filename in images:
            image_mask = specimen_mask & group["filename"].fillna("").astype(str).eq(filename)
            weights.loc[image_mask] = (
                1.0
                / max(len(specimens), 1)
                / max(len(images), 1)
                / max(int(image_mask.sum()), 1)
            )
    return weights / weights.sum()


def weighted_quantile(values: pd.Series, weights: pd.Series, quantile: float) -> float:
    frame = pd.DataFrame(
        {
            "value": pd.to_numeric(values, errors="coerce"),
            "weight": pd.to_numeric(weights, errors="coerce"),
        }
    ).dropna()
    frame = frame.loc[frame["weight"].gt(0)].sort_values("value")
    if frame.empty:
        return np.nan
    cumulative = frame["weight"].cumsum() / frame["weight"].sum()
    return float(frame.loc[cumulative.ge(quantile), "value"].iloc[0])


def panel_morphology_comparison() -> pd.DataFrame:
    pairs = pd.read_csv(BALANCED_PAIRS_CSV, low_memory=False)
    panels = {"current_top50_all_selected", "balanced_qc_curated"}
    pairs = pairs.loc[pairs["panel"].isin(panels)].copy()
    rows: list[dict[str, object]] = []
    for (panel, species), group in pairs.groupby(["panel", "species"], sort=True):
        weights = analysis_weights(group)
        rows.append(
            {
                "panel": panel,
                "species": species,
                "cell_area_um2": weighted_quantile(group["cell_area_um2"], weights, 0.5),
                "nucleus_area_um2": weighted_quantile(group["nuc_area_um2"], weights, 0.5),
                "nc_area_ratio": weighted_quantile(group["nc_area_ratio"], weights, 0.5),
            }
        )
    summaries = pd.DataFrame(rows)
    current = summaries.loc[summaries["panel"].eq("current_top50_all_selected")].drop(
        columns="panel"
    )
    balanced = summaries.loc[summaries["panel"].eq("balanced_qc_curated")].drop(
        columns="panel"
    )
    comparison = current.merge(
        balanced,
        on="species",
        suffixes=("_top50", "_balanced"),
        validate="one_to_one",
    )
    for metric in ["cell_area_um2", "nucleus_area_um2", "nc_area_ratio"]:
        comparison[f"{metric}_balanced_vs_top50_pct"] = (
            comparison[f"{metric}_balanced"] - comparison[f"{metric}_top50"]
        ) / comparison[f"{metric}_top50"] * 100.0
    return comparison


def main() -> None:
    candidates, selected = read_inputs()
    audit = species_metrics(candidates, selected)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.width", 240)

    print("SPECIES METRICS")
    print(audit.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("\nCROSS-SPECIES SUMMARY")
    print(
        audit.drop(columns="species")
        .describe()
        .transpose()
        .to_string(float_format=lambda value: f"{value:.3f}")
    )

    print("\nSELECTION STATE COUNTS")
    for column in [
        "grid_effective_decision",
        "grid_selection_tier",
        "decision",
        "current_decision",
        "ranking_profile",
    ]:
        print_distribution(selected, column)

    print("\nSELECTED-ROW SPEARMAN CORRELATIONS")
    print(selected[METRICS].corr(method="spearman").round(3).to_string())

    print("\nCROSS-SPECIES SELECTION ASSOCIATIONS")
    association_columns = [
        "candidate_rows",
        "selected_fraction",
        "selected_cell_area_median_um2",
        "selected_to_candidate_cell_median_ratio",
        "selected_cell_area_percentile_median",
    ]
    print(audit[association_columns].corr(method="spearman").round(3).to_string())

    morphology_comparison = panel_morphology_comparison()
    shift_columns = [
        "cell_area_um2_balanced_vs_top50_pct",
        "nucleus_area_um2_balanced_vs_top50_pct",
        "nc_area_ratio_balanced_vs_top50_pct",
    ]
    print("\nBALANCED-QC VERSUS TOP50 MORPHOLOGY")
    print(
        morphology_comparison[["species", *shift_columns]]
        .sort_values("cell_area_um2_balanced_vs_top50_pct")
        .to_string(index=False, float_format=lambda value: f"{value:.3f}")
    )
    print("\nBALANCED-QC MORPHOLOGY SHIFT SUMMARY")
    print(
        morphology_comparison[shift_columns]
        .describe()
        .transpose()
        .to_string(float_format=lambda value: f"{value:.3f}")
    )
    print("\nTOP50 VERSUS BALANCED-QC SPECIES RANK CORRELATIONS")
    for metric in ["cell_area_um2", "nucleus_area_um2", "nc_area_ratio"]:
        rho = morphology_comparison[
            [f"{metric}_top50", f"{metric}_balanced"]
        ].corr(method="spearman").iloc[0, 1]
        print(f"{metric}: rho={rho:.3f}")

    genome_comparison = pd.read_csv(BALANCED_GENOME_COMPARISON_CSV)
    nonreference = genome_comparison.loc[genome_comparison["species"].ne("D. fuscus")]
    genome_shift = nonreference["balanced_vs_current_pct"]
    print("\nBALANCED-QC VERSUS CURRENT-PRIMARY GENOME SHIFT")
    print(genome_shift.describe().to_string(float_format=lambda value: f"{value:.3f}"))
    print(f"all_nonreference_positive={bool(genome_shift.gt(0).all())}")
    print(f"species_ge_40pct={int(genome_shift.ge(40).sum())}/{len(genome_shift)}")
    print(f"species_ge_50pct={int(genome_shift.ge(50).sum())}/{len(genome_shift)}")

    uncertain = selected["grid_effective_decision"].fillna("").eq("maybe").sum()
    stale = (
        selected["decision"].fillna("").astype(str)
        != selected["grid_effective_decision"].fillna("").astype(str)
    ).sum()
    assert len(selected) == 1_050
    assert selected["species"].nunique() == 21
    assert audit["selected_rows"].eq(50).all()
    print("\nINVARIANTS")
    print(f"selected_rows={len(selected)}")
    print(f"selected_species={selected['species'].nunique()}")
    print(f"grid_maybe_rows={uncertain}")
    print(f"decision_vs_grid_effective_disagreements={stale}")


if __name__ == "__main__":
    main()

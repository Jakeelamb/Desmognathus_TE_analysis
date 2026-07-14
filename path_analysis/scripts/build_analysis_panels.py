#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
PANEL_DIR = DERIVED_DIR / "panels"

INPUT_FILE = DERIVED_DIR / "path_input_master.csv"
OUTPUT_SUMMARY = DERIVED_DIR / "analysis_panel_summary.csv"
OUTPUT_READINESS = DERIVED_DIR / "analysis_species_readiness.csv"
OUTPUT_PHYLO_COMPARE = DERIVED_DIR / "phylofill_panel_comparison.csv"


def confidence_score(series: pd.Series) -> pd.Series:
    mapping = {"low": 1, "medium": 2, "high": 3}
    return series.map(mapping)


def annotate_phylofill_trait(
    observed: pd.Series,
    confidence_score_series: pd.Series,
    inferred: pd.Series,
    status: pd.Series,
) -> tuple[pd.Series, pd.Series, pd.Series]:
    effective = observed.copy()
    mode = pd.Series(pd.NA, index=observed.index, dtype="object")
    usable = pd.Series(False, index=observed.index, dtype="boolean")

    observed_mediumplus = observed.notna() & confidence_score_series.ge(2)
    mode.loc[observed_mediumplus] = "observed_mediumplus"
    usable.loc[observed_mediumplus] = True

    low_supported = observed.notna() & confidence_score_series.eq(1) & status.eq(
        "low_confidence_observed_supported_by_phylogeny"
    )
    mode.loc[low_supported] = "observed_low_conf_supported_by_phylogeny"
    usable.loc[low_supported] = True

    missing_inferred = observed.isna() & inferred.notna() & status.eq("phylo_inferred_for_missing")
    effective.loc[missing_inferred] = inferred.loc[missing_inferred]
    mode.loc[missing_inferred] = "phylo_inferred_for_missing"
    usable.loc[missing_inferred] = True

    observed_low_unusable = observed.notna() & confidence_score_series.eq(1) & mode.isna()
    mode.loc[observed_low_unusable] = "observed_low_conf_not_phylo_supported"

    observed_missing = observed.isna() & mode.isna()
    mode.loc[observed_missing] = "missing_no_phylofill"

    return effective, mode, usable


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
        "all_specimens|mixed|transformed",
        case=False,
        regex=True,
    )
    out["is_low_confidence_body_proxy"] = out["body_size_proxy_confidence"].eq("low")
    out["is_low_confidence_lifestyle_proxy"] = out["lifestyle_confidence"].eq("low")
    if "ltr_history_n_pairs_high_confidence" not in out.columns:
        out["ltr_history_n_pairs_high_confidence"] = 0
    out["ltr_history_n_pairs_high_confidence"] = pd.to_numeric(
        out["ltr_history_n_pairs_high_confidence"],
        errors="coerce",
    ).fillna(0)
    out["has_ltr_high_confidence"] = out["ltr_history_n_pairs_high_confidence"].gt(0)

    (
        out["body_size_proxy_mm_phylofill"],
        out["body_size_phylofill_mode"],
        out["body_size_phylofill_usable"],
    ) = annotate_phylofill_trait(
        out["body_size_proxy_mm"],
        out["body_size_conf_score"],
        out["phylo_body_size_proxy_mm"],
        out["phylo_body_size_proxy_mm_status"],
    )
    (
        out["development_mode_phylofill"],
        out["development_phylofill_mode"],
        out["development_phylofill_usable"],
    ) = annotate_phylofill_trait(
        out["development_mode"],
        out["development_conf_score"],
        out["phylo_development_mode"],
        out["phylo_development_mode_status"],
    )
    (
        out["aquaticity_index_phylofill"],
        out["aquaticity_phylofill_mode"],
        out["aquaticity_phylofill_usable"],
    ) = annotate_phylofill_trait(
        out["aquaticity_index"],
        out["lifestyle_conf_score"],
        out["phylo_aquaticity_index"],
        out["phylo_aquaticity_index_status"],
    )
    (
        out["microhabitat_class_phylofill"],
        out["microhabitat_phylofill_mode"],
        out["microhabitat_phylofill_usable"],
    ) = annotate_phylofill_trait(
        out["microhabitat_class"],
        out["lifestyle_conf_score"],
        out["phylo_microhabitat_class"],
        out["phylo_microhabitat_class_status"],
    )
    out["has_core_organismal_phylofill"] = (
        out["body_size_phylofill_usable"]
        & out["development_phylofill_usable"]
        & out["aquaticity_phylofill_usable"]
        & out["microhabitat_phylofill_usable"]
    )
    phylo_assist_modes = {
        "phylo_inferred_for_missing",
        "observed_low_conf_supported_by_phylogeny",
    }
    out["n_core_traits_with_phylo_assistance"] = (
        out["body_size_phylofill_mode"].isin(phylo_assist_modes).astype(int)
        + out["development_phylofill_mode"].isin(phylo_assist_modes).astype(int)
        + out["aquaticity_phylofill_mode"].isin(phylo_assist_modes).astype(int)
        + out["microhabitat_phylofill_mode"].isin(phylo_assist_modes).astype(int)
    )
    out["uses_any_phylo_assistance"] = out["n_core_traits_with_phylo_assistance"].gt(0)
    return out


def build_panel(df: pd.DataFrame, mask: pd.Series, keep_cols: list[str]) -> pd.DataFrame:
    return df.loc[mask, keep_cols].sort_values("species").reset_index(drop=True)


def build_exclusion_reasons(
    row: pd.Series,
    *,
    require_ectopic: bool = False,
    require_morphology: bool = False,
    require_ltr_history: bool = False,
    require_ltr_high_confidence: bool = False,
    require_core: bool = False,
    require_mediumplus: bool = False,
    require_phylofill: bool = False,
    require_strict_body: bool = False,
) -> str:
    reasons: list[str] = []

    if not bool(row.get("has_tree_tip")):
        reasons.append("missing_tree_tip")
    if not bool(row.get("has_te")):
        reasons.append("missing_te")
    if not bool(row.get("has_genome")):
        reasons.append("missing_genome")
    if require_ectopic and not bool(row.get("has_ectopic")):
        reasons.append("missing_ectopic")
    if require_morphology and not bool(row.get("has_morphology")):
        reasons.append("missing_morphology")
    if require_ltr_history and not bool(row.get("has_ltr_history")):
        reasons.append("missing_ltr_history")
    if require_ltr_high_confidence and not bool(row.get("has_ltr_high_confidence")):
        reasons.append("no_high_confidence_ltr_pairs")

    if require_core or require_mediumplus or require_strict_body:
        if pd.isna(row.get("body_size_proxy_mm")):
            reasons.append("missing_body_size")
        if pd.isna(row.get("development_mode")):
            reasons.append("missing_development")
        if pd.isna(row.get("aquaticity_index")):
            reasons.append("missing_aquaticity")
        if pd.isna(row.get("microhabitat_class")):
            reasons.append("missing_microhabitat")

    if require_phylofill:
        if not bool(row.get("body_size_phylofill_usable")):
            reasons.append("body_size_not_usable_after_phylofill")
        if not bool(row.get("development_phylofill_usable")):
            reasons.append("development_not_usable_after_phylofill")
        if not bool(row.get("aquaticity_phylofill_usable")):
            reasons.append("aquaticity_not_usable_after_phylofill")
        if not bool(row.get("microhabitat_phylofill_usable")):
            reasons.append("microhabitat_not_usable_after_phylofill")

    if require_mediumplus or require_strict_body:
        body_score = row.get("body_size_conf_score")
        if pd.isna(body_score) or body_score < 2:
            reasons.append("low_body_confidence")
        development_score = row.get("development_conf_score")
        if pd.isna(development_score) or development_score < 2:
            reasons.append("low_development_confidence")
        lifestyle_score = row.get("lifestyle_conf_score")
        if pd.isna(lifestyle_score) or lifestyle_score < 2:
            reasons.append("low_lifestyle_confidence")

    if require_strict_body:
        if bool(row.get("uses_total_length_body_proxy")):
            reasons.append("total_length_body_proxy")
        if bool(row.get("uses_mixed_stage_body_proxy")):
            reasons.append("mixed_stage_body_proxy")

    return ";".join(reasons)


def main() -> None:
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_FILE)
    df = prepare_flags(df)
    df_indexed = df.set_index("species")

    base_cols = [
        "species",
        "te_sra_accession",
        "te_assembly_accession",
        "te_resource_lookup_path",
        "te_resource_lookup_sha256",
        "has_tree_tip",
        "has_te",
        "has_genome",
        "has_ectopic",
        "has_morphology",
        "has_ltr_history",
        "has_ltr_high_confidence",
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
        "ltr_history_median_k2p_distance",
        "ltr_history_n_pairs_estimated",
        "ltr_history_n_pairs_high_confidence",
        "ltr_history_age_low_mya",
        "ltr_history_age_central_mya",
        "ltr_history_age_high_mya",
        "ltr_history_source_id",
        "ltr_history_calibration_source_id",
        "ltr_history_rate_source_id",
        "ltr_history_source_ids",
        "body_size_proxy_mm",
        "body_size_proxy_measurement",
        "body_size_proxy_source_id",
        "body_size_proxy_confidence",
        "body_size_proxy_mm_phylofill",
        "body_size_phylofill_mode",
        "body_size_phylofill_usable",
        "development_mode",
        "development_source_id",
        "development_confidence",
        "development_mode_phylofill",
        "development_phylofill_mode",
        "development_phylofill_usable",
        "aquaticity_index",
        "aquaticity_index_phylofill",
        "aquaticity_phylofill_mode",
        "aquaticity_phylofill_usable",
        "microhabitat_class",
        "microhabitat_class_phylofill",
        "microhabitat_phylofill_mode",
        "microhabitat_phylofill_usable",
        "lifestyle_source_id",
        "lifestyle_confidence",
        "elevation_mid_m",
        "elevation_source_id",
        "has_core_organismal_complete",
        "has_core_organismal_mediumplus",
        "has_core_organismal_phylofill",
        "has_elevation",
        "uses_total_length_body_proxy",
        "uses_mixed_stage_body_proxy",
        "is_low_confidence_body_proxy",
        "is_low_confidence_lifestyle_proxy",
        "n_core_traits_with_phylo_assistance",
        "uses_any_phylo_assistance",
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
        "te_genome_organismal_primary_mediumplus": (
            df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_core_organismal_mediumplus"]
        ),
        "te_genome_organismal_primary_phylofill": (
            df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_core_organismal_phylofill"]
        ),
        "te_genome_primary_phylofill": (
            df["has_tree_tip"] & df["has_te"] & df["has_genome"] & df["has_core_organismal_phylofill"]
        ),
        "te_genome_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_core_organismal_mediumplus"]
            & ~df["uses_total_length_body_proxy"]
            & ~df["uses_mixed_stage_body_proxy"]
        ),
        "te_genome_organismal_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_core_organismal_mediumplus"]
            & ~df["uses_total_length_body_proxy"]
            & ~df["uses_mixed_stage_body_proxy"]
        ),
        "te_genome_ltr_history_primary_mediumplus": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_core_organismal_mediumplus"]
            & df["has_ltr_history"]
            & df["has_ltr_high_confidence"]
        ),
        "te_genome_ltr_history_primary_strict_body": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_core_organismal_mediumplus"]
            & df["has_ltr_history"]
            & df["has_ltr_high_confidence"]
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
        "te_genome_ectopic_organismal_primary_mediumplus": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_ectopic"]
            & df["has_core_organismal_mediumplus"]
        ),
        "te_genome_ectopic_organismal_primary_phylofill": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_ectopic"]
            & df["has_core_organismal_phylofill"]
        ),
        "te_genome_ectopic_primary_phylofill": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_ectopic"]
            & df["has_core_organismal_phylofill"]
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
        "te_genome_ectopic_organismal_primary_strict_body": (
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
        "te_genome_morphology_primary_phylofill": (
            df["has_tree_tip"]
            & df["has_te"]
            & df["has_genome"]
            & df["has_morphology"]
            & df["has_core_organismal_phylofill"]
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
    panel_requirements = {
        "te_genome_all": {},
        "te_genome_primary": {"require_core": True},
        "te_genome_primary_mediumplus": {"require_core": True, "require_mediumplus": True},
        "te_genome_organismal_primary_mediumplus": {"require_core": True, "require_mediumplus": True},
        "te_genome_organismal_primary_phylofill": {"require_phylofill": True},
        "te_genome_primary_phylofill": {"require_phylofill": True},
        "te_genome_primary_strict_body": {
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
        "te_genome_organismal_primary_strict_body": {
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
        "te_genome_ltr_history_primary_mediumplus": {
            "require_ltr_history": True,
            "require_ltr_high_confidence": True,
            "require_core": True,
            "require_mediumplus": True,
        },
        "te_genome_ltr_history_primary_strict_body": {
            "require_ltr_history": True,
            "require_ltr_high_confidence": True,
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
        "te_genome_ectopic_all": {"require_ectopic": True},
        "te_genome_ectopic_primary_mediumplus": {
            "require_ectopic": True,
            "require_core": True,
            "require_mediumplus": True,
        },
        "te_genome_ectopic_organismal_primary_mediumplus": {
            "require_ectopic": True,
            "require_core": True,
            "require_mediumplus": True,
        },
        "te_genome_ectopic_organismal_primary_phylofill": {
            "require_ectopic": True,
            "require_phylofill": True,
        },
        "te_genome_ectopic_primary_phylofill": {
            "require_ectopic": True,
            "require_phylofill": True,
        },
        "te_genome_ectopic_primary_strict_body": {
            "require_ectopic": True,
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
        "te_genome_ectopic_organismal_primary_strict_body": {
            "require_ectopic": True,
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
        "te_genome_morphology_all": {"require_morphology": True},
        "te_genome_morphology_primary_mediumplus": {
            "require_morphology": True,
            "require_core": True,
            "require_mediumplus": True,
        },
        "te_genome_morphology_primary_phylofill": {
            "require_morphology": True,
            "require_phylofill": True,
        },
        "te_genome_morphology_primary_strict_body": {
            "require_morphology": True,
            "require_core": True,
            "require_mediumplus": True,
            "require_strict_body": True,
        },
    }

    summary_rows = []
    for name, mask in masks.items():
        panel = build_panel(df, mask, base_cols)
        panel.to_csv(PANEL_DIR / f"{name}.csv", index=False)
        genome_status_counts = panel["genome_result_status"].fillna("").value_counts()
        ltr_high_conf = pd.to_numeric(panel["ltr_history_n_pairs_high_confidence"], errors="coerce")
        summary_rows.append(
            {
                "panel_name": name,
                "n_species": int(len(panel)),
                "n_genome_stable": int(genome_status_counts.get("stable", 0)),
                "n_genome_minor_caution": int(genome_status_counts.get("minor_caution", 0)),
                "n_genome_caution": int(genome_status_counts.get("caution", 0)),
                "n_genome_sensitivity_limited": int(genome_status_counts.get("sensitivity_limited", 0)),
                "n_ltr_high_confidence_species": int(panel["has_ltr_high_confidence"].sum()),
                "min_ltr_high_confidence_pairs": (
                    int(ltr_high_conf.min()) if len(panel) and name.startswith("te_genome_ltr_history") else pd.NA
                ),
                "n_with_elevation": int(panel["has_elevation"].sum()),
                "n_low_body_proxy": int(panel["is_low_confidence_body_proxy"].sum()),
                "n_total_length_body_proxy": int(panel["uses_total_length_body_proxy"].sum()),
                "n_mixed_stage_body_proxy": int(panel["uses_mixed_stage_body_proxy"].sum()),
                "n_with_phylo_assistance": int(panel["uses_any_phylo_assistance"].sum()),
                "species_with_phylo_assistance": ";".join(
                    panel.loc[panel["uses_any_phylo_assistance"], "species"].tolist()
                ),
                "species_list": ";".join(panel["species"].tolist()),
            }
        )

    readiness_cols = [
        "species",
        "has_tree_tip",
        "has_te",
        "has_genome",
        "has_ectopic",
        "has_morphology",
        "has_ltr_history",
        "has_ltr_high_confidence",
        "ltr_history_n_pairs_high_confidence",
        "genome_result_status",
        "genome_flag_summary",
        "body_size_proxy_measurement",
        "body_size_proxy_confidence",
        "development_confidence",
        "lifestyle_confidence",
        "has_core_organismal_complete",
        "has_core_organismal_mediumplus",
        "has_core_organismal_phylofill",
        "body_size_phylofill_mode",
        "development_phylofill_mode",
        "aquaticity_phylofill_mode",
        "microhabitat_phylofill_mode",
        "n_core_traits_with_phylo_assistance",
        "uses_any_phylo_assistance",
        "uses_total_length_body_proxy",
        "uses_mixed_stage_body_proxy",
        "has_elevation",
        "organismal_source_ids",
        "te_feature_source_ids",
        "genome_source_id",
        "morphology_source_id",
    ]
    readiness = df[readiness_cols].sort_values("species").reset_index(drop=True)
    for name, requirements in panel_requirements.items():
        eligibility_map = df.assign(_eligible=masks[name]).set_index("species")["_eligible"]
        readiness[f"eligible_{name}"] = readiness["species"].map(eligibility_map).fillna(False)
        readiness[f"{name}_exclusion_reasons"] = readiness.apply(
            lambda row: (
                ""
                if bool(row[f"eligible_{name}"])
                else build_exclusion_reasons(
                    df_indexed.loc[row["species"]],
                    **requirements,
                )
            ),
            axis=1,
        )

    pd.DataFrame(summary_rows).to_csv(OUTPUT_SUMMARY, index=False)
    comparison_specs = [
        ("te_genome_primary_mediumplus", "te_genome_primary_phylofill"),
        ("te_genome_organismal_primary_mediumplus", "te_genome_organismal_primary_phylofill"),
        ("te_genome_ectopic_primary_mediumplus", "te_genome_ectopic_primary_phylofill"),
        ("te_genome_ectopic_organismal_primary_mediumplus", "te_genome_ectopic_organismal_primary_phylofill"),
        ("te_genome_morphology_primary_mediumplus", "te_genome_morphology_primary_phylofill"),
    ]
    comparison_rows = []
    for observed_panel, phylofill_panel in comparison_specs:
        observed_species = set(df.loc[masks[observed_panel], "species"])
        phylofill_species = set(df.loc[masks[phylofill_panel], "species"])
        added_species = sorted(phylofill_species - observed_species)
        removed_species = sorted(observed_species - phylofill_species)
        phylofill_subset = df.loc[masks[phylofill_panel]]
        species_with_assistance = sorted(
            phylofill_subset.loc[phylofill_subset["uses_any_phylo_assistance"], "species"].tolist()
        )
        comparison_rows.append(
            {
                "observed_panel": observed_panel,
                "phylofill_panel": phylofill_panel,
                "n_observed_panel": int(len(observed_species)),
                "n_phylofill_panel": int(len(phylofill_species)),
                "n_added_species": int(len(added_species)),
                "added_species": ";".join(added_species),
                "n_removed_species": int(len(removed_species)),
                "removed_species": ";".join(removed_species),
                "n_species_with_phylo_assistance": int(len(species_with_assistance)),
                "species_with_phylo_assistance": ";".join(species_with_assistance),
                "is_species_set_identical": observed_species == phylofill_species,
            }
        )
    pd.DataFrame(comparison_rows).to_csv(OUTPUT_PHYLO_COMPARE, index=False)
    readiness.to_csv(OUTPUT_READINESS, index=False)

    print(f"Wrote panel directory {PANEL_DIR}")
    print(f"Wrote {OUTPUT_SUMMARY}")
    print(f"Wrote {OUTPUT_PHYLO_COMPARE}")
    print(f"Wrote {OUTPUT_READINESS}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results" / "data"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"
PANELS_DIR = DERIVED_DIR / "panels"
TEMPLATES_DIR = PROJECT_ROOT / "path_analysis" / "data" / "templates"

ORDER_FILE = RESULTS_DIR / "dnaPipeTE_order_breakdown.csv"
SUPERFAMILY_FILE = RESULTS_DIR / "dnaPipeTE_superfamily_breakdown.csv"
ORDER_DIVERSITY_FILE = RESULTS_DIR / "diversity_order_stats.csv"
SUPERFAMILY_DIVERSITY_FILE = RESULTS_DIR / "diversity_superfamily_stats.csv"
DIVERGENCE_FILE = RESULTS_DIR / "divergence" / "divergence_summary_statistics_by_species.csv"
ECTOPIC_FILE = RESULTS_DIR / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"

TE_FEATURE_FILE = DERIVED_DIR / "te_path_features.csv"
TE_PANEL_FILE = DERIVED_DIR / "te_model_feature_panel.csv"
MASTER_INPUT_FILE = DERIVED_DIR / "path_input_master.csv"
PANEL_SUMMARY_FILE = DERIVED_DIR / "analysis_panel_summary.csv"
TRACEABILITY_GAPS_FILE = DERIVED_DIR / "source_traceability_gaps.csv"
SOURCE_MANIFEST_FILE = TEMPLATES_DIR / "source_manifest.csv"

SUMMARY_OUT = DERIVED_DIR / "te_audit_summary.csv"
EDGE_CASES_OUT = DERIVED_DIR / "te_audit_edge_cases.csv"

TE_SOURCE_IDS = {
    "repo_dnapipete_order_breakdown",
    "repo_dnapipete_superfamily_breakdown",
    "repo_diversity_order_stats",
    "repo_diversity_superfamily_stats",
    "repo_divergence_summary_statistics_by_species",
    "repo_ectopic_recombination_filtered_3000bp_5plusdomains",
}


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


def clr_transform(df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    positive_values = df.where(df > 0).stack()
    pseudocount = float(positive_values.min() / 2) if not positive_values.empty else 1e-6
    adjusted = df + pseudocount
    geometric_mean = np.exp(np.log(adjusted).mean(axis=1))
    clr = np.log(adjusted).sub(np.log(geometric_mean), axis=0)
    return clr, pseudocount


def pca_scores_2d(df: pd.DataFrame) -> pd.DataFrame:
    centered = df - df.mean(axis=0)
    matrix = centered.to_numpy()
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    scores = u[:, :2] * s[:2]
    result = pd.DataFrame(index=df.index)
    result["major_order_clr_pc1"] = scores[:, 0]
    result["major_order_clr_pc2"] = scores[:, 1] if scores.shape[1] > 1 else 0.0
    return result


def max_abs_diff(a: pd.Series, b: pd.Series) -> float:
    diff = (pd.to_numeric(a, errors="coerce") - pd.to_numeric(b, errors="coerce")).abs()
    if diff.dropna().empty:
        return 0.0
    return float(diff.max())


def max_abs_diff_with_sign_flip(a: pd.Series, b: pd.Series) -> tuple[float, bool]:
    direct = max_abs_diff(a, b)
    flipped = max_abs_diff(a, -pd.to_numeric(b, errors="coerce"))
    if flipped < direct:
        return flipped, True
    return direct, False


def compare_frames_exact(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    columns: list[str],
    allow_sign_flip: set[str] | None = None,
    tol: float = 1e-12,
) -> tuple[bool, list[str]]:
    allow_sign_flip = allow_sign_flip or set()
    merged = observed[["species"] + columns].merge(
        expected[["species"] + columns],
        on="species",
        how="outer",
        suffixes=("_obs", "_exp"),
        validate="one_to_one",
    )
    failures: list[str] = []
    for column in columns:
        obs = merged[f"{column}_obs"]
        exp = merged[f"{column}_exp"]
        missing_mismatch = int((obs.isna() != exp.isna()).sum())
        is_numeric = (
            (pd.api.types.is_numeric_dtype(obs) and not pd.api.types.is_bool_dtype(obs))
            or (pd.api.types.is_numeric_dtype(exp) and not pd.api.types.is_bool_dtype(exp))
        )
        if is_numeric:
            if column in allow_sign_flip:
                diff, sign_flipped = max_abs_diff_with_sign_flip(obs, exp)
                if diff > tol or missing_mismatch:
                    failures.append(
                        f"{column}: max_diff={diff:.3g}, missing_mismatch={missing_mismatch}, sign_flipped={sign_flipped}"
                    )
            else:
                diff = max_abs_diff(obs, exp)
                if diff > tol or missing_mismatch:
                    failures.append(f"{column}: max_diff={diff:.3g}, missing_mismatch={missing_mismatch}")
        else:
            mismatch = int((obs.fillna("__NA__") != exp.fillna("__NA__")).sum())
            if mismatch:
                failures.append(f"{column}: mismatch={mismatch}")
    return not failures, failures


def reconstruct_te_features() -> pd.DataFrame:
    order = pd.read_csv(ORDER_FILE, index_col=0)
    order.index = standardize_species(order.index.to_series())
    order.columns = [f"order_{c.lower()}" for c in order.columns]

    feature_df = order.copy()
    positive_values = order.where(order > 0).stack()
    pseudocount = float(positive_values.min() / 2) if not positive_values.empty else 1e-6

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
    clr, _ = clr_transform(feature_df[major_order_cols])
    feature_df = feature_df.join(pca_scores_2d(clr))

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
        metrics[metric] = divergence.pivot(
            index="species",
            columns="group_name_lower",
            values=metric,
        ).add_prefix(f"{metric}_")
    divergence_wide = pd.concat(metrics.values(), axis=1)

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
    weights = order[weight_cols].copy()
    weights.columns = [c.replace("order_", "") for c in weights.columns]

    rows = []
    for species, row in weights.iterrows():
        aligned_weights = row / row.sum()
        values: dict[str, float | str] = {"species": species}
        for metric, out_name in [
            ("percent_divergence_median", "weighted_te_divergence_p90"),
            ("percent_deletions_median", "weighted_te_deletions_p90"),
            ("percent_insertions_median", "weighted_te_insertions_p90"),
        ]:
            metric_cols = [f"{metric}_{name}" for name in aligned_weights.index]
            metric_values = divergence_wide.loc[species, metric_cols]
            values[out_name] = float((aligned_weights.to_numpy() * metric_values.to_numpy()).sum())
        for order_name in ["ltr", "line", "tir", "dirs", "sine"]:
            for metric, suffix in [
                ("percent_divergence_median", "divergence_p90"),
                ("percent_deletions_median", "deletions_p90"),
            ]:
                values[f"{order_name}_{suffix}"] = float(
                    divergence_wide.loc[species, f"{metric}_{order_name}"]
                )
        rows.append(values)
    divergence_features = pd.DataFrame(rows)

    ectopic = pd.read_csv(ECTOPIC_FILE, sep="\t")
    ectopic["species"] = standardize_species(ectopic["species"])
    ectopic["ratio_terminal_internal"] = pd.to_numeric(ectopic["ratio_terminal_internal"], errors="coerce")
    complete_state = ectopic["Complete"].astype(str).str.strip().str.lower()
    ectopic["complete_flag_known"] = complete_state.map({"yes": 1.0, "no": 0.0})
    ectopic_summary = (
        ectopic.groupby("species", dropna=False)
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
    ectopic_summary["ectopic_log10_mean_ratio"] = ectopic_summary["ectopic_mean_ratio"].apply(
        lambda x: math.log10(x) if pd.notna(x) and x > 0 else np.nan
    )

    merged = feature_df.reset_index().rename(columns={"index": "species"})
    merged = merged.merge(diversity, on="species", how="left", validate="one_to_one")
    merged = merged.merge(divergence_features, on="species", how="left", validate="one_to_one")
    merged = merged.merge(ectopic_summary, on="species", how="left", validate="one_to_one")
    return merged.sort_values("species").reset_index(drop=True)


def build_panel_expectation(te_features: pd.DataFrame) -> pd.DataFrame:
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
    panel = te_features[keep_cols].copy()
    panel["has_core_te_panel"] = panel[
        ["ltr_line_logratio", "order_pielou", "weighted_te_divergence_p90"]
    ].notna().all(axis=1)
    panel["has_extended_te_panel"] = panel[
        ["ltr_line_logratio", "order_pielou", "weighted_te_divergence_p90", "ectopic_log10_mean_ratio"]
    ].notna().all(axis=1)
    return panel.sort_values("species").reset_index(drop=True)


def record_check(rows: list[dict[str, object]], check_name: str, status: str, details: str, value: object = "") -> None:
    rows.append(
        {
            "check_name": check_name,
            "status": status,
            "value": value,
            "details": details,
        }
    )


def main() -> None:
    summary_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []

    te_features = pd.read_csv(TE_FEATURE_FILE)
    te_panel = pd.read_csv(TE_PANEL_FILE)
    path_input = pd.read_csv(MASTER_INPUT_FILE)
    panel_summary = pd.read_csv(PANEL_SUMMARY_FILE)
    traceability_gaps = pd.read_csv(TRACEABILITY_GAPS_FILE)
    source_manifest = pd.read_csv(SOURCE_MANIFEST_FILE)

    source_checks = [
        ("te_order", ORDER_FILE, "te_order_source_sha256"),
        ("order_diversity", ORDER_DIVERSITY_FILE, "order_diversity_source_sha256"),
        ("divergence", DIVERGENCE_FILE, "divergence_source_sha256"),
        ("ectopic", ECTOPIC_FILE, "ectopic_source_sha256"),
    ]
    for label, path, column in source_checks:
        actual = sha256_for_file(path)
        stored = te_panel[column].dropna().unique().tolist()
        status = "pass" if stored == [actual] else "fail"
        record_check(
            summary_rows,
            f"{label}_hash_matches_stored_value",
            status,
            f"stored_unique={stored}; actual={actual}",
            actual,
        )

    order = pd.read_csv(ORDER_FILE, index_col=0)
    row_sums = order.sum(axis=1)
    max_deviation = float((row_sums - 100).abs().max())
    negative_cells = int((order < 0).sum().sum())
    record_check(
        summary_rows,
        "order_percentages_sum_to_100",
        "pass" if max_deviation < 1e-9 else "fail",
        f"max_abs_deviation={max_deviation:.3g}",
        max_deviation,
    )
    record_check(
        summary_rows,
        "order_table_has_no_negative_values",
        "pass" if negative_cells == 0 else "fail",
        f"negative_cells={negative_cells}",
        negative_cells,
    )

    order_species = standardize_species(order.index.to_series())
    order.index = order_species
    diversity = pd.read_csv(ORDER_DIVERSITY_FILE, index_col=0)
    diversity.index = standardize_species(diversity.index.to_series())
    proportions = order.div(order.sum(axis=1), axis=0)
    shannon = -(proportions * np.log(proportions)).sum(axis=1)
    richness = (order > 0).sum(axis=1)
    pielou = shannon / np.log(richness)
    gini_simpson = 1 - (proportions**2).sum(axis=1)

    shannon_diff = float((shannon - diversity["Shannon_Diversity"]).abs().max())
    pielou_diff = float((pielou - diversity["Pielou_Evenness"]).abs().max())
    simpson_ratio = diversity["Simpson_Diversity"] / gini_simpson
    simpson_ratio_min = float(simpson_ratio.min())
    simpson_ratio_max = float(simpson_ratio.max())
    simpson_ratio_target = 100.0 / 99.0

    record_check(
        summary_rows,
        "order_shannon_matches_order_table",
        "pass" if shannon_diff < 1e-12 else "fail",
        f"max_abs_diff={shannon_diff:.3g}",
        shannon_diff,
    )
    record_check(
        summary_rows,
        "order_pielou_matches_order_table",
        "pass" if pielou_diff < 1e-12 else "fail",
        f"max_abs_diff={pielou_diff:.3g}",
        pielou_diff,
    )
    simpson_status = "warn"
    if abs(simpson_ratio_min - simpson_ratio_target) > 1e-12 or abs(simpson_ratio_max - simpson_ratio_target) > 1e-12:
        simpson_status = "fail"
    record_check(
        summary_rows,
        "order_simpson_definition_is_corrected_vs_raw_gini_simpson",
        simpson_status,
        (
            "stored Simpson equals raw_gini_simpson * 100/99 across all species; "
            "document the correction if this metric is reported"
        ),
        f"{simpson_ratio_min:.12f}..{simpson_ratio_max:.12f}",
    )

    divergence = pd.read_csv(DIVERGENCE_FILE)
    divergence_sub = divergence[
        (divergence["group_level"] == "order") & (divergence["threshold"].astype(float) == 0.9)
    ].copy()
    divergence_species = int(divergence_sub["Desmognathus_Species"].nunique())
    divergence_orders = int(divergence_sub["group_name"].nunique())
    divergence_dupes = int(divergence_sub.duplicated(["Desmognathus_Species", "group_name"]).sum())
    divergence_missing = int(
        divergence_sub[
            [
                "percent_divergence_median",
                "percent_deletions_median",
                "percent_insertions_median",
            ]
        ].isna().sum().sum()
    )
    record_check(
        summary_rows,
        "divergence_threshold_p90_order_table_is_complete",
        "pass" if divergence_species == 34 and divergence_orders == 9 and divergence_dupes == 0 and divergence_missing == 0 else "fail",
        (
            f"species={divergence_species}; orders={divergence_orders}; "
            f"duplicates={divergence_dupes}; missing_metric_cells={divergence_missing}"
        ),
        len(divergence_sub),
    )

    reconstructed = reconstruct_te_features()
    stored = te_features.sort_values("species").reset_index(drop=True)
    feature_columns = [
        "order_dirs",
        "order_helitron",
        "order_line",
        "order_ltr",
        "order_maverick",
        "order_ple",
        "order_sine",
        "order_tir",
        "order_yr",
        "ltr_line_logratio",
        "ltr_tir_logratio",
        "retro_dna_logratio",
        "major_order_clr_pc1",
        "major_order_clr_pc2",
        "order_simpson",
        "order_shannon",
        "order_pielou",
        "weighted_te_divergence_p90",
        "weighted_te_deletions_p90",
        "weighted_te_insertions_p90",
        "ltr_divergence_p90",
        "ltr_deletions_p90",
        "line_divergence_p90",
        "line_deletions_p90",
        "tir_divergence_p90",
        "tir_deletions_p90",
        "dirs_divergence_p90",
        "dirs_deletions_p90",
        "sine_divergence_p90",
        "sine_deletions_p90",
        "ectopic_n_rows_total",
        "ectopic_n_elements",
        "ectopic_mean_ratio",
        "ectopic_median_ratio",
        "ectopic_mean_terminal_depth",
        "ectopic_mean_internal_depth",
        "ectopic_mean_domain_count",
        "ectopic_n_complete_known",
        "ectopic_complete_fraction",
        "ectopic_log10_mean_ratio",
    ]
    features_match, feature_failures = compare_frames_exact(
        stored,
        reconstructed,
        feature_columns,
        allow_sign_flip={"major_order_clr_pc1", "major_order_clr_pc2"},
    )
    record_check(
        summary_rows,
        "te_path_features_reconstruct_from_repo_local_sources",
        "pass" if features_match else "fail",
        "; ".join(feature_failures) if feature_failures else "all compared TE feature columns matched exactly",
        len(feature_columns),
    )

    panel_expected = build_panel_expectation(te_features)
    panel_columns = [
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
        "has_core_te_panel",
        "has_extended_te_panel",
    ]
    panel_match, panel_failures = compare_frames_exact(te_panel, panel_expected, panel_columns)
    record_check(
        summary_rows,
        "te_model_feature_panel_matches_expected_subset",
        "pass" if panel_match else "fail",
        "; ".join(panel_failures) if panel_failures else "subset columns and booleans matched expected values",
        len(te_panel),
    )

    join_map: dict[str, str] = {}
    for column in te_features.columns:
        if column == "species":
            continue
        suffixed = f"{column}_te"
        if suffixed in path_input.columns:
            join_map[column] = suffixed
        elif column in path_input.columns:
            join_map[column] = column
    path_input_indexed = path_input.set_index("species")
    te_features_indexed = te_features.set_index("species")
    shared_species = sorted(set(path_input_indexed.index) & set(te_features_indexed.index))
    join_failures: list[str] = []
    for te_column, path_column in join_map.items():
        left = path_input_indexed.loc[shared_species, path_column]
        right = te_features_indexed.loc[shared_species, te_column]
        is_numeric = (
            (pd.api.types.is_numeric_dtype(left) and not pd.api.types.is_bool_dtype(left))
            or (pd.api.types.is_numeric_dtype(right) and not pd.api.types.is_bool_dtype(right))
        )
        if is_numeric:
            diff = max_abs_diff(left, right)
            missing_mismatch = int((left.isna() != right.isna()).sum())
            if diff > 1e-12 or missing_mismatch:
                join_failures.append(
                    f"{te_column}<-{path_column}: max_diff={diff:.3g}, missing_mismatch={missing_mismatch}"
                )
        else:
            mismatch = int((left.fillna("__NA__") != right.fillna("__NA__")).sum())
            if mismatch:
                join_failures.append(f"{te_column}<-{path_column}: mismatch={mismatch}")
    record_check(
        summary_rows,
        "path_input_master_preserves_te_feature_values",
        "pass" if not join_failures else "fail",
        (
            "; ".join(join_failures)
            if join_failures
            else "all TE feature columns matched te_path_features.csv after join; overlapping provenance columns are read from *_te fields"
        ),
        len(join_map),
    )

    panel_files = {
        "te_genome_all": PANELS_DIR / "te_genome_all.csv",
        "te_genome_primary_mediumplus": PANELS_DIR / "te_genome_primary_mediumplus.csv",
        "te_genome_primary_strict_body": PANELS_DIR / "te_genome_primary_strict_body.csv",
        "te_genome_ectopic_all": PANELS_DIR / "te_genome_ectopic_all.csv",
        "te_genome_ectopic_primary_mediumplus": PANELS_DIR / "te_genome_ectopic_primary_mediumplus.csv",
        "te_genome_ectopic_primary_strict_body": PANELS_DIR / "te_genome_ectopic_primary_strict_body.csv",
        "te_genome_morphology_all": PANELS_DIR / "te_genome_morphology_all.csv",
        "te_genome_morphology_primary_mediumplus": PANELS_DIR / "te_genome_morphology_primary_mediumplus.csv",
        "te_genome_morphology_primary_strict_body": PANELS_DIR / "te_genome_morphology_primary_strict_body.csv",
    }
    summary_lookup = dict(zip(panel_summary["panel_name"], panel_summary["n_species"]))
    count_failures = []
    for panel_name, path in panel_files.items():
        actual = len(pd.read_csv(path))
        expected = int(summary_lookup.get(panel_name, -1))
        if actual != expected:
            count_failures.append(f"{panel_name}: file_rows={actual}, summary_rows={expected}")
    record_check(
        summary_rows,
        "analysis_panel_summary_matches_panel_files",
        "pass" if not count_failures else "fail",
        "; ".join(count_failures) if count_failures else "panel file row counts matched analysis_panel_summary.csv",
        len(panel_files),
    )

    te_gaps = traceability_gaps[traceability_gaps["source_id"].isin(TE_SOURCE_IDS)]
    manifest_te = source_manifest[source_manifest["source_id"].isin(TE_SOURCE_IDS)]
    record_check(
        summary_rows,
        "te_source_manifest_entries_present",
        "pass" if len(manifest_te) == len(TE_SOURCE_IDS) else "fail",
        f"manifest_entries={len(manifest_te)}; expected={len(TE_SOURCE_IDS)}",
        len(manifest_te),
    )
    record_check(
        summary_rows,
        "te_traceability_gaps_absent",
        "pass" if te_gaps.empty else "fail",
        f"te_gap_rows={len(te_gaps)}",
        len(te_gaps),
    )

    ectopic = pd.read_csv(ECTOPIC_FILE, sep="\t")
    ectopic["species_std"] = standardize_species(ectopic["species"])
    missing_ratio_rows = ectopic[ectopic["ratio_terminal_internal"].isna()].copy()
    unknown_complete_rows = ectopic[ectopic["Complete"].astype(str).str.lower() == "unknown"].copy()
    record_check(
        summary_rows,
        "ectopic_missing_ratio_rows_present",
        "warn" if len(missing_ratio_rows) else "pass",
        "missing ratio rows are excluded from usable ectopic support counts and mean-ratio summaries",
        len(missing_ratio_rows),
    )
    record_check(
        summary_rows,
        "ectopic_unknown_complete_rows_present",
        "warn" if len(unknown_complete_rows) else "pass",
        "Complete=unknown rows are excluded from the ectopic_complete_fraction denominator",
        len(unknown_complete_rows),
    )

    for _, row in missing_ratio_rows.iterrows():
        edge_rows.append(
            {
                "edge_case_type": "ectopic_missing_ratio",
                "species": row["species_std"],
                "source_file": str(ECTOPIC_FILE.relative_to(PROJECT_ROOT)),
                "detail": f"sequence={row['sequence']}",
                "audit_implication": "excluded from ectopic_n_elements and mean_ratio summaries because ratio is missing",
            }
        )
    for _, row in unknown_complete_rows.iterrows():
        edge_rows.append(
            {
                "edge_case_type": "ectopic_unknown_complete_flag",
                "species": row["species_std"],
                "source_file": str(ECTOPIC_FILE.relative_to(PROJECT_ROOT)),
                "detail": f"sequence={row['sequence']}",
                "audit_implication": "excluded from ectopic_complete_fraction denominator because Complete is unknown",
            }
        )

    extended_missing = te_panel.loc[
        te_panel["has_core_te_panel"] & ~te_panel["has_extended_te_panel"],
        "species",
    ].tolist()
    for species in extended_missing:
        edge_rows.append(
            {
                "edge_case_type": "no_ectopic_panel_support",
                "species": species,
                "source_file": str(ECTOPIC_FILE.relative_to(PROJECT_ROOT)),
                "detail": "absent from stored ectopic source table",
                "audit_implication": "eligible for TE core panel but excluded from ectopic-extended panels",
            }
        )
    record_check(
        summary_rows,
        "extended_te_panel_missingness_is_explicit",
        "pass",
        f"species_missing_ectopic_support={','.join(extended_missing)}",
        len(extended_missing),
    )

    old_vs_panel_checks = [
        ("dataset_te_genome.csv", PANELS_DIR / "te_genome_primary_mediumplus.csv"),
        ("dataset_te_genome_ectopic.csv", PANELS_DIR / "te_genome_ectopic_primary_mediumplus.csv"),
        ("dataset_te_genome_morphology.csv", PANELS_DIR / "te_genome_morphology_primary_mediumplus.csv"),
    ]
    old_vs_panel_failures = []
    for old_name, panel_path in old_vs_panel_checks:
        old_species = set(pd.read_csv(DERIVED_DIR / old_name)["species"])
        panel_species = set(pd.read_csv(panel_path)["species"])
        if old_species != panel_species:
            old_vs_panel_failures.append(old_name)
    record_check(
        summary_rows,
        "legacy_overlap_datasets_match_current_mediumplus_species_sets",
        "pass" if not old_vs_panel_failures else "warn",
        (
            "legacy overlap dataset species sets matched current panel species sets"
            if not old_vs_panel_failures
            else f"mismatched datasets={','.join(old_vs_panel_failures)}"
        ),
        len(old_vs_panel_checks) - len(old_vs_panel_failures),
    )

    summary = pd.DataFrame(summary_rows)
    edges = pd.DataFrame(edge_rows)
    summary.to_csv(SUMMARY_OUT, index=False)
    edges.to_csv(EDGE_CASES_OUT, index=False)

    n_pass = int((summary["status"] == "pass").sum())
    n_warn = int((summary["status"] == "warn").sum())
    n_fail = int((summary["status"] == "fail").sum())
    print(f"Wrote {SUMMARY_OUT}")
    print(f"Wrote {EDGE_CASES_OUT}")
    print(f"Checks: pass={n_pass}, warn={n_warn}, fail={n_fail}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build a non-destructive, identity-safe ectopic-depth sensitivity branch."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = PROJECT_ROOT / "input_data" / "ectopic_recombination"
LOOKUP = PROJECT_ROOT / "input_data" / "lookup_table.txt"
PANEL = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "derived"
    / "panels"
    / "te_genome_primary_mediumplus.csv"
)
TESORTER = INPUT_DIR / "combined_sequences.fasta.rexdb-metazoa.cls.tsv"
HISTORICAL_MASTER = PROJECT_ROOT / "results" / "data" / "ectopic_recombination_master.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "corrected" / "ectopic"
FIGURE_DIR = PROJECT_ROOT / "results" / "figures" / "corrected" / "ectopic"
AUDIT_DIR = PROJECT_ROOT / "plans" / "publication-readiness-deep-audit"
ELEMENT_OUTPUT = OUTPUT_DIR / "ectopic_element_metrics_analysis18_v1.csv"
SPECIES_OUTPUT = OUTPUT_DIR / "ectopic_species_robustness_analysis18_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "ectopic_analysis18_v1.manifest.json"
REPORT_OUTPUT = AUDIT_DIR / "ectopic_recombination_corrected_analysis18_v1.md"
MIN_ELEMENT_LENGTH = 3000
DOMAIN_COUNTS = {5, 6}
BOOTSTRAP_REPLICATES = 2000


def _canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def _count_domains(value: object) -> int:
    if pd.isna(value) or not str(value).strip() or str(value).lower() in {"no", "unknown"}:
        return 0
    return len(str(value).split("|"))


def _element_ids(elements: pd.DataFrame) -> pd.Series:
    return (
        elements["sequence"].astype(str)
        + "_"
        + pd.to_numeric(elements["element start"], errors="raise").astype(int).astype(str)
        + "_"
        + pd.to_numeric(elements["element end"], errors="raise").astype(int).astype(str)
    )


def merge_elements_with_tesorter(
    elements: pd.DataFrame, tesorter: pd.DataFrame
) -> pd.DataFrame:
    """Join TEsorter annotations to exact sequence/start/end element IDs."""

    required_elements = {"sequence", "element start", "element end"}
    missing_elements = sorted(required_elements.difference(elements.columns))
    if missing_elements:
        raise ValueError(f"LTR element table is missing columns: {missing_elements}")
    if "#TE" not in tesorter:
        raise ValueError("TEsorter table is missing #TE")
    left = elements.copy()
    left["element_id"] = _element_ids(left)
    if left["element_id"].duplicated().any():
        raise ValueError("LTR input contains duplicate exact element IDs")
    right = tesorter.copy()
    if right["#TE"].duplicated().any():
        raise ValueError("TEsorter input contains duplicate exact element IDs")
    merged = left.merge(
        right,
        left_on="element_id",
        right_on="#TE",
        how="left",
        validate="one_to_one",
        indicator="tesorter_join_state",
    )
    merged["tesorter_annotation_found"] = merged["tesorter_join_state"].eq("both")
    return merged.drop(columns="tesorter_join_state")


def calculate_zero_aware_depth_metrics(
    depth: pd.DataFrame,
    element_length: int,
    left_ltr_length: int,
    right_ltr_length: int,
) -> Dict[str, object]:
    """Calculate full-position and historical nonzero-only regional depths."""

    element_length = int(element_length)
    left_ltr_length = int(left_ltr_length)
    right_ltr_length = int(right_ltr_length)
    internal_length = element_length - left_ltr_length - right_ltr_length
    if min(element_length, left_ltr_length, right_ltr_length, internal_length) <= 0:
        raise ValueError("LTR element has invalid regional lengths")
    required = {"position", "depth"}
    missing = sorted(required.difference(depth.columns))
    if missing:
        raise ValueError(f"Depth table is missing columns: {missing}")
    table = depth.loc[:, ["position", "depth"]].copy()
    table["position"] = pd.to_numeric(table["position"], errors="raise").astype(int)
    table["depth"] = pd.to_numeric(table["depth"], errors="raise")
    if table["position"].duplicated().any():
        raise ValueError("Depth file contains duplicate positions")
    if table["position"].lt(1).any() or table["position"].gt(element_length).any():
        raise ValueError("Depth file contains positions outside the LTR element")
    if table["depth"].lt(0).any() or not np.isfinite(table["depth"]).all():
        raise ValueError("Depth file contains invalid depth values")

    values = np.zeros(element_length, dtype=float)
    values[table["position"].to_numpy() - 1] = table["depth"].to_numpy()
    left = values[:left_ltr_length]
    internal = values[left_ltr_length : element_length - right_ltr_length]
    right = values[element_length - right_ltr_length :]
    terminal = np.concatenate([left, right])

    def nonzero_mean(region: np.ndarray) -> float:
        positive = region[region > 0]
        return float(positive.mean()) if len(positive) else 0.0

    mean_terminal = float(terminal.mean())
    mean_internal = float(internal.mean())
    mean_terminal_nonzero = nonzero_mean(terminal)
    mean_internal_nonzero = nonzero_mean(internal)
    ratio = mean_terminal / mean_internal if mean_internal > 0 else np.nan
    ratio_nonzero = (
        mean_terminal_nonzero / mean_internal_nonzero
        if mean_internal_nonzero > 0
        else np.nan
    )
    return {
        "depth_positions_expected": element_length,
        "depth_positions_reported": int(len(table)),
        "depth_positions_missing": int(element_length - len(table)),
        "depth_positions_explicit_zero": int(table["depth"].eq(0).sum()),
        "left_ltr_positive_coverage_fraction": float(np.mean(left > 0)),
        "right_ltr_positive_coverage_fraction": float(np.mean(right > 0)),
        "terminal_positive_coverage_fraction": float(np.mean(terminal > 0)),
        "internal_positive_coverage_fraction": float(np.mean(internal > 0)),
        "mean_depth_left_ltr_all_positions": float(left.mean()),
        "mean_depth_right_ltr_all_positions": float(right.mean()),
        "mean_depth_terminal_all_positions": mean_terminal,
        "mean_depth_internal_all_positions": mean_internal,
        "ratio_terminal_internal_all_positions": ratio,
        "log2_ratio_terminal_internal_all_positions": (
            float(np.log2(ratio)) if np.isfinite(ratio) and ratio > 0 else np.nan
        ),
        "mean_depth_terminal_nonzero_only": mean_terminal_nonzero,
        "mean_depth_internal_nonzero_only": mean_internal_nonzero,
        "ratio_terminal_internal_nonzero_only": ratio_nonzero,
        "left_right_terminal_log2_imbalance": (
            float(np.log2((left.mean() + 0.5) / (right.mean() + 0.5)))
        ),
    }


def _trimmed_mean(values: np.ndarray, fraction: float = 0.1) -> float:
    ordered = np.sort(values)
    trim = int(np.floor(len(ordered) * fraction))
    kept = ordered[trim : len(ordered) - trim] if trim else ordered
    return float(kept.mean())


def _bootstrap_interval(
    values: np.ndarray, statistic: str, seed_text: str
) -> tuple:
    seed = int(hashlib.sha256(seed_text.encode("utf-8")).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    sampled = rng.choice(values, size=(BOOTSTRAP_REPLICATES, len(values)), replace=True)
    if statistic == "median":
        estimates = np.median(sampled, axis=1)
    elif statistic == "geometric_mean":
        estimates = np.exp(np.mean(np.log(sampled), axis=1))
    else:
        raise ValueError(f"Unknown bootstrap statistic: {statistic}")
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def summarize_species(elements: pd.DataFrame) -> pd.DataFrame:
    branches = {
        "all_5plus_domain_elements": pd.Series(True, index=elements.index),
        "coverage_ge_80pct": (
            elements["left_ltr_positive_coverage_fraction"].ge(0.8)
            & elements["right_ltr_positive_coverage_fraction"].ge(0.8)
            & elements["internal_positive_coverage_fraction"].ge(0.8)
        ),
        "tesorter_complete_yes": elements["Complete"].astype(str).str.lower().eq("yes"),
        "coverage_ge_80pct_and_complete_yes": (
            elements["left_ltr_positive_coverage_fraction"].ge(0.8)
            & elements["right_ltr_positive_coverage_fraction"].ge(0.8)
            & elements["internal_positive_coverage_fraction"].ge(0.8)
            & elements["Complete"].astype(str).str.lower().eq("yes")
        ),
    }
    rows: List[Dict[str, object]] = []
    for branch, mask in branches.items():
        selected = elements.loc[mask].copy()
        for species, group in selected.groupby("species", sort=True):
            values = group["ratio_terminal_internal_all_positions"].dropna().to_numpy()
            values = values[np.isfinite(values) & (values > 0)]
            if len(values) == 0:
                continue
            median_ci = _bootstrap_interval(values, "median", f"{species}:{branch}:median")
            geometric_ci = _bootstrap_interval(
                values, "geometric_mean", f"{species}:{branch}:geometric"
            )
            mean_value = float(values.mean())
            median_value = float(np.median(values))
            geometric_mean = float(np.exp(np.mean(np.log(values))))
            winsorized = np.clip(values, *np.quantile(values, [0.05, 0.95]))
            loo_means = [np.delete(values, index).mean() for index in range(len(values))]
            loo_medians = [np.median(np.delete(values, index)) for index in range(len(values))]
            rows.append(
                {
                    "species": species,
                    "te_sra_accession": group["te_sra_accession"].iloc[0],
                    "te_assembly_accession": group["te_assembly_accession"].iloc[0],
                    "analysis_branch": branch,
                    "n_elements": len(values),
                    "ratio_arithmetic_mean": mean_value,
                    "ratio_median": median_value,
                    "ratio_geometric_mean": geometric_mean,
                    "ratio_trimmed_mean_10pct": _trimmed_mean(values),
                    "ratio_winsorized_mean_5pct": float(winsorized.mean()),
                    "ratio_q05": float(np.quantile(values, 0.05)),
                    "ratio_q95": float(np.quantile(values, 0.95)),
                    "ratio_max": float(values.max()),
                    "fraction_ratio_gt_1": float(np.mean(values > 1)),
                    "median_bootstrap_ci_low": median_ci[0],
                    "median_bootstrap_ci_high": median_ci[1],
                    "geometric_mean_bootstrap_ci_low": geometric_ci[0],
                    "geometric_mean_bootstrap_ci_high": geometric_ci[1],
                    "max_abs_leave_one_out_mean_shift": float(
                        np.max(np.abs(np.asarray(loo_means) - mean_value))
                    ),
                    "max_abs_leave_one_out_median_shift": float(
                        np.max(np.abs(np.asarray(loo_medians) - median_value))
                    ),
                    "largest_element_fraction_of_ratio_sum": float(
                        values.max() / values.sum()
                    ),
                    "median_terminal_positive_coverage_fraction": float(
                        group["terminal_positive_coverage_fraction"].median()
                    ),
                    "median_internal_positive_coverage_fraction": float(
                        group["internal_positive_coverage_fraction"].median()
                    ),
                }
            )
    return pd.DataFrame(rows).sort_values(["analysis_branch", "species"]).reset_index(
        drop=True
    )


def _load_analysis_elements() -> tuple:
    panel_species = pd.read_csv(PANEL)["species"].map(_canonical_species)
    if len(panel_species) != 18 or panel_species.duplicated().any():
        raise ValueError("Final TE/genome panel must contain 18 unique species")
    lookup = pd.read_csv(LOOKUP, sep="\t")
    lookup["species"] = lookup["Species"].map(_canonical_species)
    lookup = lookup.loc[lookup["species"].isin(set(panel_species))].copy()
    if len(lookup) != 18 or lookup["species"].duplicated().any():
        raise ValueError("Active TE lookup does not map the final 18 species one-to-one")

    frames = []
    raw_paths = []
    missing_resources = []
    for row in lookup.sort_values("species").itertuples(index=False):
        path = INPUT_DIR / f"{row.Genome_Accension}_tabout.csv"
        if not path.exists():
            missing_resources.append(
                {
                    "species": row.species,
                    "te_sra_accession": row.SRA_Accension,
                    "te_assembly_accession": row.Genome_Accension,
                    "expected_path": _portable(path),
                }
            )
            continue
        frame = pd.read_csv(path, sep="\t", index_col=False)
        numeric_columns = [
            "element start",
            "element end",
            "element length",
            "lLTR length",
            "rLTR length",
        ]
        converted = frame[numeric_columns].apply(pd.to_numeric, errors="coerce")
        if converted.isna().any().any():
            raise ValueError(f"{path} has non-numeric required LTR coordinates")
        frame[numeric_columns] = converted.astype(int)
        frame = frame.loc[frame["element length"].ge(MIN_ELEMENT_LENGTH)].copy()
        frame["species"] = row.species
        frame["te_sra_accession"] = row.SRA_Accension
        frame["te_assembly_accession"] = row.Genome_Accension
        frames.append(frame)
        raw_paths.append(path)
    elements = pd.concat(frames, ignore_index=True)
    tesorter = pd.read_csv(TESORTER, sep="\t")
    merged = merge_elements_with_tesorter(elements, tesorter)
    merged["domain_count"] = merged["Domains"].apply(_count_domains)
    selected = merged.loc[merged["domain_count"].isin(DOMAIN_COUNTS)].copy()
    if selected.empty:
        raise RuntimeError("No exact-joined 5/6-domain LTR elements remain")
    return selected, merged, raw_paths, missing_resources


def _add_depth_metrics(elements: pd.DataFrame) -> pd.DataFrame:
    metrics = []
    for _, row in elements.iterrows():
        path = INPUT_DIR / f"{row['element_id']}.fa.depth.txt"
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
        depth = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["depth_sequence_id", "position", "depth"],
        )
        metric = calculate_zero_aware_depth_metrics(
            depth,
            row["element length"],
            row["lLTR length"],
            row["rLTR length"],
        )
        metric["depth_file"] = _portable(path)
        metric["depth_file_sha256"] = _sha256(path)
        metrics.append(metric)
    result = pd.concat(
        [elements.reset_index(drop=True), pd.DataFrame(metrics)], axis=1
    )
    result["coverage_ge_80pct"] = (
        result["left_ltr_positive_coverage_fraction"].ge(0.8)
        & result["right_ltr_positive_coverage_fraction"].ge(0.8)
        & result["internal_positive_coverage_fraction"].ge(0.8)
    )
    return result


def _add_historical_comparison(elements: pd.DataFrame) -> pd.DataFrame:
    historical = pd.read_csv(HISTORICAL_MASTER, sep="\t")
    historical["species"] = historical["species"].map(_canonical_species)
    historical["element_id"] = _element_ids(historical)
    historical = historical.loc[
        historical["domain_count"].isin(DOMAIN_COUNTS),
        ["element_id", "ratio_terminal_internal"],
    ].copy()
    if historical["element_id"].duplicated().any():
        raise ValueError("Historical filtered branch has duplicate exact element IDs")
    historical = historical.rename(
        columns={"ratio_terminal_internal": "historical_nonzero_only_ratio"}
    )
    compared = elements.merge(historical, on="element_id", how="left", validate="one_to_one")
    compared["nonzero_ratio_reproduction_abs_error"] = (
        compared["ratio_terminal_internal_nonzero_only"]
        - compared["historical_nonzero_only_ratio"]
    ).abs()
    return compared


def _historical_join_stats() -> Dict[str, int]:
    historical = pd.read_csv(HISTORICAL_MASTER, sep="\t")
    historical["element_id"] = _element_ids(historical)
    return {
        "historical_rows": len(historical),
        "historical_unique_element_ids": int(historical["element_id"].nunique()),
        "historical_extra_rows_from_contig_join": int(
            len(historical) - historical["element_id"].nunique()
        ),
        "historical_element_annotation_mismatches": int(
            historical["element_id"].ne(historical["#TE"]).sum()
        ),
    }


def _write_figures(elements: pd.DataFrame, species: pd.DataFrame) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    paths = []

    order = (
        elements.groupby("species")["log2_ratio_terminal_internal_all_positions"]
        .median()
        .sort_values()
        .index
    )
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.boxplot(
        data=elements,
        y="species",
        x="log2_ratio_terminal_internal_all_positions",
        order=order,
        color="#9ecae1",
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=elements,
        y="species",
        x="log2_ratio_terminal_internal_all_positions",
        order=order,
        color="#1f2937",
        alpha=0.35,
        size=2.5,
        ax=ax,
    )
    ax.axvline(0, color="#b91c1c", linestyle="--", linewidth=1)
    ax.set(
        xlabel="Element log2(terminal depth / internal depth), zeros retained",
        ylabel="Species",
        title="LTR terminal:internal depth distributions (screening proxy)",
    )
    fig.tight_layout()
    path = FIGURE_DIR / "ectopic_element_log2_ratio_by_species_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    base = species.loc[
        species["analysis_branch"].eq("all_5plus_domain_elements")
    ].copy()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(base["ratio_median"], base["ratio_arithmetic_mean"], s=45, color="#2563eb")
    low = min(base["ratio_median"].min(), base["ratio_arithmetic_mean"].min())
    high = max(base["ratio_median"].max(), base["ratio_arithmetic_mean"].max())
    ax.plot([low, high], [low, high], linestyle="--", color="#6b7280")
    for row in base.nlargest(4, "max_abs_leave_one_out_mean_shift").itertuples():
        ax.annotate(row.species, (row.ratio_median, row.ratio_arithmetic_mean), fontsize=8)
    ax.set(
        xlabel="Species median ratio",
        ylabel="Species arithmetic mean ratio",
        title="Arithmetic means are sensitive to rare high-ratio elements",
    )
    fig.tight_layout()
    path = FIGURE_DIR / "ectopic_mean_vs_median_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(
        data=elements,
        x="internal_positive_coverage_fraction",
        y="terminal_positive_coverage_fraction",
        hue="species",
        legend=False,
        alpha=0.55,
        s=25,
        ax=ax,
    )
    ax.axhline(0.8, color="#b91c1c", linestyle="--", linewidth=1)
    ax.axvline(0.8, color="#b91c1c", linestyle="--", linewidth=1)
    ax.set(
        xlabel="Internal-region positive-position fraction",
        ylabel="Combined terminal positive-position fraction",
        title="Per-element depth coverage support",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
    )
    fig.tight_layout()
    path = FIGURE_DIR / "ectopic_coverage_qc_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        base["n_elements"],
        base["max_abs_leave_one_out_mean_shift"],
        s=45,
        color="#7c3aed",
    )
    for row in base.nlargest(5, "max_abs_leave_one_out_mean_shift").itertuples():
        ax.annotate(row.species, (row.n_elements, row.max_abs_leave_one_out_mean_shift), fontsize=8)
    ax.set(
        xlabel="Elements per species",
        ylabel="Maximum leave-one-element-out shift in arithmetic mean",
        title="Influence of individual LTR candidates on species means",
    )
    fig.tight_layout()
    path = FIGURE_DIR / "ectopic_mean_influence_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)
    return paths


def _render_report(
    elements: pd.DataFrame,
    species: pd.DataFrame,
    join_all: pd.DataFrame,
    missing_resources: List[Dict[str, object]],
    historical_join_stats: Dict[str, int],
) -> str:
    base = species.loc[species["analysis_branch"].eq("all_5plus_domain_elements")]
    worst = base.loc[base["max_abs_leave_one_out_mean_shift"].idxmax()]
    coverage_pass = int(elements["coverage_ge_80pct"].sum())
    corrected_delta = (
        elements["ratio_terminal_internal_all_positions"]
        - elements["ratio_terminal_internal_nonzero_only"]
    ).abs()
    return f"""# Corrected ectopic-recombination depth audit: final panel

This branch rebuilds element identity and depth summaries for only the declared 18-species TE/genome panel. It preserves historical outputs unchanged.

## What was corrected

- TEsorter annotations are joined by exact `sequence_start_end`, not contig name alone.
- Zero-depth positions remain in terminal and internal regional means.
- Per-element left-LTR, right-LTR, terminal, and internal positive-position coverage is explicit.
- Species summaries include arithmetic mean, median, geometric mean, trimmed/winsorized means, bootstrap intervals, and leave-one-element-out influence.
- 5/6-domain elements are reported under all-element, >=80% coverage, TEsorter-complete, and combined sensitivity branches.

## Audit result

- Exact 5/6-domain elements: {len(elements)} across {elements.species.nunique()} of the 18 panel species.
- Species with no retained element-level estimate: {18 - elements.species.nunique()}.
- Missing assembly-level `tabout` resources: {', '.join(item['species'] for item in missing_resources) if missing_resources else 'none'}.
- Elements passing >=80% positive-position coverage in both LTRs and the internal region: {coverage_pass}/{len(elements)} ({coverage_pass / len(elements):.1%}).
- Maximum absolute change caused by retaining zeros instead of deleting them: {corrected_delta.max():.6f} in the element ratio.
- Maximum nonzero-only reproduction error against the historical stored ratio: {elements.nonzero_ratio_reproduction_abs_error.max():.3e}.
- The historical contig-only TEsorter join created {historical_join_stats['historical_extra_rows_from_contig_join']} extra row and {historical_join_stats['historical_element_annotation_mismatches']} coordinate-mismatched annotation; that row has fewer than five domains and did not enter the current 5/6-domain predictor.
- Strongest arithmetic-mean influence remains *D. {worst.species}*: removing one element can shift its species mean by {worst.max_abs_leave_one_out_mean_shift:.3f}.

## Interpretation boundary

Terminal:internal read depth is a screening proxy for excess LTR-like sequence relative to internal sequence. It is not yet a measured ectopic-recombination rate or a direct solo-LTR count. Mapping command, reference construction, multimapper handling, MAPQ, secondary/supplementary alignment policy, duplicate policy, and independent validation against structurally called solo/intact LTRs are absent. No corrected ectopic variable is approved for confirmatory path analysis until those are resolved.

## Validation figures

- `results/figures/corrected/ectopic/ectopic_element_log2_ratio_by_species_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_mean_vs_median_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_coverage_qc_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_mean_influence_analysis18_v1.png`
"""


def main() -> None:
    occupied = [path for path in (ELEMENT_OUTPUT, SPECIES_OUTPUT, MANIFEST_OUTPUT) if path.exists()]
    if occupied:
        raise FileExistsError(
            "Non-destructive ectopic rebuild requires unused outputs: "
            + ", ".join(str(path) for path in occupied)
        )
    selected, all_joined, raw_paths, missing_resources = _load_analysis_elements()
    elements = _add_depth_metrics(selected)
    elements = _add_historical_comparison(elements)
    species = summarize_species(elements)
    historical_join_stats = _historical_join_stats()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    elements.to_csv(ELEMENT_OUTPUT, index=False)
    species.to_csv(SPECIES_OUTPUT, index=False)
    figure_paths = _write_figures(elements, species)
    REPORT_OUTPUT.write_text(
        _render_report(
            elements,
            species,
            all_joined,
            missing_resources,
            historical_join_stats,
        ),
        encoding="utf-8",
    )
    manifest = {
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "analysis_status": "sensitivity_only_not_validated_ectopic_rate",
        "eligible_for_confirmatory_path_analysis": False,
        "n_panel_species": 18,
        "n_species_with_5plus_domain_elements": int(elements["species"].nunique()),
        "n_elements": len(elements),
        "n_elements_coverage_ge_80pct": int(elements["coverage_ge_80pct"].sum()),
        "missing_panel_resources": missing_resources,
        "exact_tesorter_join": True,
        "historical_join_audit": historical_join_stats,
        "zero_depth_positions_retained": True,
        "out_of_panel_species_processed": False,
        "panel": _portable(PANEL),
        "panel_sha256": _sha256(PANEL),
        "lookup": _portable(LOOKUP),
        "lookup_sha256": _sha256(LOOKUP),
        "tesorter": _portable(TESORTER),
        "tesorter_sha256": _sha256(TESORTER),
        "historical_master": _portable(HISTORICAL_MASTER),
        "historical_master_sha256": _sha256(HISTORICAL_MASTER),
        "assembly_element_inputs": [
            {"path": _portable(path), "sha256": _sha256(path)} for path in raw_paths
        ],
        "element_output": _portable(ELEMENT_OUTPUT),
        "element_output_sha256": _sha256(ELEMENT_OUTPUT),
        "species_output": _portable(SPECIES_OUTPUT),
        "species_output_sha256": _sha256(SPECIES_OUTPUT),
        "report": _portable(REPORT_OUTPUT),
        "report_sha256": _sha256(REPORT_OUTPUT),
        "figures": [
            {"path": _portable(path), "sha256": _sha256(path)} for path in figure_paths
        ],
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "coverage_sensitivity_threshold": 0.8,
        "blocking_provenance": [
            "read_mapping_command_and_reference",
            "multimapper_policy",
            "mapq_threshold",
            "secondary_and_supplementary_alignment_policy",
            "duplicate_policy",
            "structural_solo_ltr_validation",
        ],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {ELEMENT_OUTPUT}")
    print(f"Wrote {SPECIES_OUTPUT}")
    print(f"Wrote {MANIFEST_OUTPUT}")
    print(f"Wrote {REPORT_OUTPUT}")
    for path in figure_paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()

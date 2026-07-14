#!/usr/bin/env python3
"""Build final-panel TE diversity and compositional-PCA audit products."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
PANEL = PROJECT_ROOT / "path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv"
ORDER = PROJECT_ROOT / "results/data/dnaPipeTE_order_breakdown.csv"
SUPERFAMILY = PROJECT_ROOT / "results/data/dnaPipeTE_superfamily_breakdown.csv"
HISTORICAL_ORDER_DIVERSITY = PROJECT_ROOT / "results/data/diversity_order_stats.csv"
HISTORICAL_SUPERFAMILY_DIVERSITY = (
    PROJECT_ROOT / "results/data/diversity_superfamily_stats.csv"
)
MASS = (
    PROJECT_ROOT
    / "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv"
)
OUTPUT_DIR = PROJECT_ROOT / "results/data/corrected/diversity_pca"
FIGURE_DIR = PROJECT_ROOT / "results/figures/corrected/diversity_pca"
AUDIT_DIR = PROJECT_ROOT / "plans/publication-readiness-deep-audit"
DIVERSITY_OUTPUT = OUTPUT_DIR / "te_diversity_mass_sensitivity_analysis18_v1.csv"
COMPOSITION_OUTPUT = OUTPUT_DIR / "te_composition_matrices_analysis18_v1.csv"
CLR_OUTPUT = OUTPUT_DIR / "te_pca_clr_matrices_analysis18_v1.csv"
SCORES_OUTPUT = OUTPUT_DIR / "te_pca_scores_analysis18_v1.csv"
LOADINGS_OUTPUT = OUTPUT_DIR / "te_pca_loadings_analysis18_v1.csv"
VARIANCE_OUTPUT = OUTPUT_DIR / "te_pca_variance_analysis18_v1.csv"
STABILITY_OUTPUT = OUTPUT_DIR / "te_pca_stability_analysis18_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "te_diversity_pca_analysis18_v1.manifest.json"
REPORT_OUTPUT = AUDIT_DIR / "te_diversity_pca_corrected_analysis18_v1.md"


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


def calculate_diversity_indices(values: np.ndarray) -> Dict[str, float]:
    """Return standard entropy and Hill-number definitions for a composition."""

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) == 0 or values.sum() <= 0:
        return {
            "observed_richness": 0,
            "shannon_entropy": 0.0,
            "simpson_dominance": 1.0,
            "gini_simpson": 0.0,
            "hill_q1": 0.0,
            "hill_q2": 0.0,
            "pielou_evenness": 0.0,
        }
    p = values / values.sum()
    richness = len(p)
    shannon = float(-np.sum(p * np.log(p)))
    dominance = float(np.sum(p**2))
    return {
        "observed_richness": richness,
        "shannon_entropy": shannon,
        "simpson_dominance": dominance,
        "gini_simpson": 1.0 - dominance,
        "hill_q1": float(np.exp(shannon)),
        "hill_q2": float(1.0 / dominance),
        "pielou_evenness": float(shannon / np.log(richness)) if richness > 1 else 0.0,
    }


def build_mass_aware_composition(
    breakdown: pd.DataFrame,
    retained_fraction: pd.Series,
    unresolved_fraction: pd.Series,
) -> pd.DataFrame:
    """Scale a closed classified composition and append unresolved mass."""

    classified = breakdown.apply(pd.to_numeric, errors="coerce")
    if classified.isna().any().any() or classified.lt(0).any().any():
        raise ValueError("Breakdown contains invalid composition values")
    classified = classified.div(classified.sum(axis=1), axis=0)
    retained = retained_fraction.reindex(classified.index)
    unresolved = unresolved_fraction.reindex(classified.index)
    if retained.isna().any() or unresolved.isna().any():
        raise ValueError("Mass ledger is missing composition species")
    if not (retained + unresolved).sub(1).abs().le(1e-10).all():
        raise ValueError("Retained plus unresolved fractions do not sum to one")
    result = classified.mul(retained, axis=0)
    result["Unresolved"] = unresolved
    if not result.sum(axis=1).sub(1).abs().le(1e-10).all():
        raise RuntimeError("Mass-aware composition does not conserve unit mass")
    return result


def replace_zeros_and_clr(
    composition: pd.DataFrame, method: str
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """Apply an explicit zero replacement, reclose, and CLR-transform."""

    table = composition.astype(float).copy()
    if table.lt(0).any().any() or table.sum(axis=1).le(0).any():
        raise ValueError("Composition has negative values or empty rows")
    table = table.div(table.sum(axis=1), axis=0)
    zero_mask = table.eq(0)
    n_zeros = int(zero_mask.sum().sum())
    replacements: List[float] = []
    if n_zeros:
        if method == "half_global_min_positive":
            delta = float(table.where(table.gt(0)).min().min() / 2)
            table = table.mask(zero_mask, delta)
            replacements = [delta]
        elif method == "half_feature_min_positive":
            for column in table:
                mask = table[column].eq(0)
                if not mask.any():
                    continue
                positives = table.loc[~mask, column]
                if positives.empty:
                    raise ValueError(f"Feature {column} is entirely zero")
                delta = float(positives.min() / 2)
                table.loc[mask, column] = delta
                replacements.append(delta)
        else:
            raise ValueError(f"Unknown zero replacement method: {method}")
        table = table.div(table.sum(axis=1), axis=0)
    elif method not in {"half_global_min_positive", "half_feature_min_positive", "not_needed"}:
        raise ValueError(f"Unknown zero replacement method: {method}")
    clr = np.log(table)
    clr = clr.sub(clr.mean(axis=1), axis=0)
    return table, clr, {
        "method": "not_needed" if n_zeros == 0 else method,
        "n_zeros_replaced": n_zeros,
        "replacement_min": min(replacements) if replacements else 0.0,
        "replacement_max": max(replacements) if replacements else 0.0,
    }


def _run_pca(clr: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    centered = clr.to_numpy() - clr.to_numpy().mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    n_components = min(len(clr) - 1, clr.shape[1] - 1)
    singular = singular[:n_components]
    scores = u[:, :n_components] * singular
    loadings = vt[:n_components].T
    for index in range(n_components):
        anchor = np.argmax(np.abs(loadings[:, index]))
        if loadings[anchor, index] < 0:
            loadings[:, index] *= -1
            scores[:, index] *= -1
    eigenvalues = singular**2 / (len(clr) - 1)
    variance = eigenvalues / eigenvalues.sum()
    pc_names = [f"PC{index + 1}" for index in range(n_components)]
    return {
        "scores": pd.DataFrame(scores, index=clr.index, columns=pc_names),
        "loadings": pd.DataFrame(loadings, index=clr.columns, columns=pc_names),
        "variance": pd.DataFrame(
            {
                "pc": pc_names,
                "eigenvalue": eigenvalues,
                "proportion_variance": variance,
                "cumulative_variance": np.cumsum(variance),
            }
        ),
    }


def _pca_stability(clr: pd.DataFrame, full: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for species in clr.index:
        loo = _run_pca(clr.drop(index=species))
        for pc in [name for name in ("PC1", "PC2") if name in loo["loadings"]]:
            full_loading = full["loadings"][pc]
            loo_loading = loo["loadings"][pc].reindex(full_loading.index)
            loading_cor = float(np.corrcoef(full_loading, loo_loading)[0, 1])
            sign = 1 if loading_cor >= 0 else -1
            loading_cor = abs(loading_cor)
            common = clr.index.difference([species])
            score_cor = float(
                np.corrcoef(
                    full["scores"].loc[common, pc],
                    loo["scores"].loc[common, pc] * sign,
                )[0, 1]
            )
            rows.append(
                {
                    "species_left_out": species,
                    "pc": pc,
                    "loading_correlation": loading_cor,
                    "score_correlation": abs(score_cor),
                }
            )
    return pd.DataFrame(rows)


def _load_breakdown(path: Path, species: List[str]) -> pd.DataFrame:
    table = pd.read_csv(path, index_col=0)
    table.index = pd.Index([_canonical_species(value) for value in table.index], name="species")
    selected = table.reindex(species)
    if selected.isna().any().any():
        raise ValueError(f"{path} is missing final-panel values")
    return selected.apply(pd.to_numeric, errors="raise")


def _diversity_table(
    compositions: Dict[Tuple[str, str], pd.DataFrame], mass: pd.DataFrame
) -> pd.DataFrame:
    resources = mass.set_index("species")[["te_sra_accession", "te_assembly_accession"]]
    rows = []
    for (level, mode), composition in compositions.items():
        for species, values in composition.iterrows():
            row = {
                "species": species,
                "te_level": level,
                "composition_mode": mode,
                **resources.loc[species].to_dict(),
                **calculate_diversity_indices(values.to_numpy()),
            }
            rows.append(row)
    return pd.DataFrame(rows)


def _historical_reproduction(diversity: pd.DataFrame) -> Dict[str, float]:
    maximum = {}
    for level, path in [
        ("order", HISTORICAL_ORDER_DIVERSITY),
        ("superfamily", HISTORICAL_SUPERFAMILY_DIVERSITY),
    ]:
        historical = pd.read_csv(path, index_col=0)
        historical.index = [_canonical_species(value) for value in historical.index]
        current = diversity.loc[
            diversity["te_level"].eq(level)
            & diversity["composition_mode"].eq("classified_conditional")
        ].set_index("species")
        pairs = {
            "gini_simpson": "Simpson_Diversity",
            "shannon_entropy": "Shannon_Diversity",
            "pielou_evenness": "Pielou_Evenness",
        }
        maximum[level] = max(
            float((current[column] - historical.loc[current.index, old]).abs().max())
            for column, old in pairs.items()
        )
    return maximum


def _run_all_pca(compositions: Dict[Tuple[str, str], pd.DataFrame]):
    specs = [
        ("order_classified_primary", ("order", "classified_conditional"), 1, "primary_descriptive"),
        ("order_mass_aware_sensitivity", ("order", "mass_aware_unresolved_bin"), 1, "sensitivity"),
        ("superfamily_classified_presence3", ("superfamily", "classified_conditional"), 3, "supplementary"),
        ("superfamily_mass_aware_presence3", ("superfamily", "mass_aware_unresolved_bin"), 3, "sensitivity"),
        ("superfamily_classified_presence5", ("superfamily", "classified_conditional"), 5, "sensitivity"),
    ]
    score_rows = []
    loading_rows = []
    variance_rows = []
    stability_rows = []
    clr_rows = []
    manifests = []
    pca_cache = {}
    for analysis_id, key, min_presence, role in specs:
        composition = compositions[key]
        presence = composition.gt(0).sum(axis=0)
        selected = composition.loc[:, presence.ge(min_presence)].copy()
        retained_mass = selected.sum(axis=1)
        selected = selected.div(retained_mass, axis=0)
        methods = ["half_global_min_positive"]
        if selected.eq(0).any().any():
            methods.append("half_feature_min_positive")
        for method in methods:
            _, clr, replacement = replace_zeros_and_clr(selected, method)
            result = _run_pca(clr)
            stability = _pca_stability(clr, result)
            pca_cache[(analysis_id, replacement["method"])] = result

            scores = result["scores"].reset_index().rename(columns={"index": "species"})
            scores.insert(0, "zero_replacement", replacement["method"])
            scores.insert(0, "analysis_id", analysis_id)
            score_rows.append(scores)
            loadings = result["loadings"].reset_index().rename(columns={"index": "feature"})
            loadings.insert(0, "zero_replacement", replacement["method"])
            loadings.insert(0, "analysis_id", analysis_id)
            loading_rows.append(loadings)
            variance = result["variance"].copy()
            variance.insert(0, "zero_replacement", replacement["method"])
            variance.insert(0, "analysis_id", analysis_id)
            variance_rows.append(variance)
            stability.insert(0, "zero_replacement", replacement["method"])
            stability.insert(0, "analysis_id", analysis_id)
            stability_rows.append(stability)
            clr_long = clr.rename_axis("species").reset_index().melt(
                id_vars="species", var_name="feature", value_name="clr_value"
            )
            clr_long.insert(0, "zero_replacement", replacement["method"])
            clr_long.insert(0, "analysis_id", analysis_id)
            clr_rows.append(clr_long)
            manifests.append(
                {
                    "analysis_id": analysis_id,
                    "role": role,
                    "te_level": key[0],
                    "composition_mode": key[1],
                    "min_presence_species": min_presence,
                    "n_species": len(selected),
                    "n_features": selected.shape[1],
                    "median_retained_feature_mass": float(retained_mass.median()),
                    "zero_replacement": replacement["method"],
                    "n_zeros_replaced": replacement["n_zeros_replaced"],
                    "replacement_min": replacement["replacement_min"],
                    "replacement_max": replacement["replacement_max"],
                    "pc1_variance": float(result["variance"].iloc[0]["proportion_variance"]),
                    "pc2_variance": float(result["variance"].iloc[1]["proportion_variance"]),
                }
            )
    return {
        "scores": pd.concat(score_rows, ignore_index=True),
        "loadings": pd.concat(loading_rows, ignore_index=True),
        "variance": pd.concat(variance_rows, ignore_index=True),
        "stability": pd.concat(stability_rows, ignore_index=True),
        "clr": pd.concat(clr_rows, ignore_index=True),
        "manifest": pd.DataFrame(manifests),
        "cache": pca_cache,
    }


def _write_figures(diversity: pd.DataFrame, pca: Dict[str, object]) -> List[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    paths = []

    metrics = ["gini_simpson", "shannon_entropy", "pielou_evenness"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for row_index, level in enumerate(["order", "superfamily"]):
        base = diversity.loc[diversity["te_level"].eq(level)]
        wide = base.pivot(index="species", columns="composition_mode", values=metrics)
        for column_index, metric in enumerate(metrics):
            ax = axes[row_index, column_index]
            x = wide[(metric, "classified_conditional")]
            y = wide[(metric, "mass_aware_unresolved_bin")]
            ax.scatter(x, y, color="#2563eb")
            low, high = min(x.min(), y.min()), max(x.max(), y.max())
            ax.plot([low, high], [low, high], "--", color="#6b7280")
            ax.set_title(f"{level}: {metric.replace('_', ' ')}")
            ax.set_xlabel("Conditional on classified mass")
            ax.set_ylabel("Unresolved-bin sensitivity")
    fig.suptitle("TE diversity sensitivity to unresolved classification mass", fontweight="bold")
    fig.tight_layout()
    path = FIGURE_DIR / "te_diversity_unresolved_mass_sensitivity_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    def plot_scores(analysis_id: str, path_name: str, title: str):
        manifest = pca["manifest"].loc[pca["manifest"]["analysis_id"].eq(analysis_id)].iloc[0]
        method = manifest["zero_replacement"]
        scores = pca["scores"].loc[
            pca["scores"]["analysis_id"].eq(analysis_id)
            & pca["scores"]["zero_replacement"].eq(method)
        ]
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(scores["PC1"], scores["PC2"], s=45, color="#0f766e")
        for row in scores.itertuples():
            ax.annotate(row.species, (row.PC1, row.PC2), fontsize=8, xytext=(3, 3), textcoords="offset points")
        ax.set(
            xlabel=f"PC1 ({manifest.pc1_variance:.1%})",
            ylabel=f"PC2 ({manifest.pc2_variance:.1%})",
            title=title,
        )
        ax.set_aspect("equal", adjustable="datalim")
        fig.tight_layout()
        output = FIGURE_DIR / path_name
        fig.savefig(output, dpi=300)
        plt.close(fig)
        paths.append(output)

    plot_scores(
        "order_classified_primary",
        "te_order_clr_pca_scores_analysis18_v1.png",
        "Order-level CLR PCA (classified composition; descriptive)",
    )
    plot_scores(
        "superfamily_classified_presence3",
        "te_superfamily_clr_pca_scores_analysis18_v1.png",
        "Superfamily CLR PCA (presence >=3; supplementary)",
    )

    loadings = pca["loadings"].loc[
        pca["loadings"]["analysis_id"].eq("order_classified_primary")
    ].drop_duplicates("feature")
    plot_loadings = loadings.set_index("feature")[["PC1", "PC2"]].stack().reset_index()
    plot_loadings.columns = ["feature", "pc", "loading"]
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.barplot(data=plot_loadings, y="feature", x="loading", hue="pc", ax=ax)
    ax.set_title("Order-level CLR PCA loadings")
    fig.tight_layout()
    path = FIGURE_DIR / "te_order_clr_pca_loadings_analysis18_v1.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    paths.append(path)

    sensitivity = pca["scores"].loc[
        pca["scores"]["analysis_id"].eq("superfamily_classified_presence3")
    ]
    methods = sensitivity["zero_replacement"].unique()
    if len(methods) == 2:
        left = sensitivity.loc[sensitivity["zero_replacement"].eq(methods[0])].set_index("species")
        right = sensitivity.loc[sensitivity["zero_replacement"].eq(methods[1])].set_index("species")
        fig, ax = plt.subplots(figsize=(8, 7))
        for species in left.index:
            ax.plot([left.loc[species, "PC1"], right.loc[species, "PC1"]], [left.loc[species, "PC2"], right.loc[species, "PC2"]], color="#9ca3af", linewidth=0.8)
        ax.scatter(left["PC1"], left["PC2"], label=methods[0], color="#2563eb")
        ax.scatter(right["PC1"], right["PC2"], label=methods[1], color="#dc2626", marker="x")
        ax.legend()
        ax.set(title="Superfamily PCA sensitivity to zero replacement", xlabel="PC1", ylabel="PC2")
        fig.tight_layout()
        path = FIGURE_DIR / "te_superfamily_pca_zero_replacement_sensitivity_analysis18_v1.png"
        fig.savefig(path, dpi=300)
        plt.close(fig)
        paths.append(path)
    return paths


def _render_report(diversity: pd.DataFrame, pca: Dict[str, object], reproduction: Dict[str, float]) -> str:
    correlations = []
    for level in ["order", "superfamily"]:
        table = diversity.loc[diversity["te_level"].eq(level)]
        for metric in ["gini_simpson", "shannon_entropy", "pielou_evenness"]:
            wide = table.pivot(index="species", columns="composition_mode", values=metric)
            rho = spearmanr(wide["classified_conditional"], wide["mass_aware_unresolved_bin"]).statistic
            correlations.append((level, metric, rho))
    stability = pca["stability"].groupby(["analysis_id", "zero_replacement", "pc"])[
        ["loading_correlation", "score_correlation"]
    ].min()
    primary = pca["manifest"].loc[pca["manifest"]["analysis_id"].eq("order_classified_primary")].iloc[0]
    corr_lines = "\n".join(
        f"- {level} {metric}: Spearman rho = {rho:.3f}." for level, metric, rho in correlations
    )
    return f"""# Corrected TE diversity and compositional PCA audit

This branch is restricted to the declared 18-species TE/genome panel and preserves all historical diversity/PCA outputs.

## Diversity definitions

- `Simpson_Diversity` in the historical tables is reproduced as Gini-Simpson `1 - sum(p_i^2)`; it is not a finite-count-corrected Simpson estimator.
- Shannon entropy uses natural logarithms.
- Pielou evenness is `H / log(S_observed)` and is conditional on positive classified categories.
- Hill numbers `exp(H)` (q=1) and `1 / sum(p_i^2)` (q=2) are now included.
- Biological indices conditional on classified mass are kept separate from a technical sensitivity that appends unresolved mass as one non-biological bin.

Maximum historical reproduction error is {reproduction['order']:.3e} for order and {reproduction['superfamily']:.3e} for superfamily.

Rank sensitivity to unresolved mass:
{corr_lines}

## PCA

The primary descriptive ordination is order-level CLR PCA conditional on classified order mass ({int(primary.n_features)} features; PC1 {primary.pc1_variance:.1%}, PC2 {primary.pc2_variance:.1%}). Superfamily PCA is supplementary because p approaches/exceeds n, sparse zeros require replacement, and rare-feature retention changes geometry. Both global-half-minimum and feature-half-minimum zero replacements are exported as sensitivities, along with exact CLR matrices, loadings, scores, variance, and leave-one-species-out stability.

Minimum leave-one-out loading/score correlations by analysis are available in `te_pca_stability_analysis18_v1.csv`; no PCA axis is approved as a causal variable without the phylogenetic and zero-replacement checks.

## PERMANOVA verdict

The current clade PERMANOVA is **not approved for inference**. It permutes species labels freely even though groups are phylogenetic clades, violating exchangeability under phylogenetic covariance. Its Bray-Curtis and CLR sensitivity calculations and beta-dispersion checks are useful descriptively, but p-values must not support a clade claim. A phylogenetically valid simulation/permutation or comparative multivariate model is required.

## Validation figures

- `results/figures/corrected/diversity_pca/te_diversity_unresolved_mass_sensitivity_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_clr_pca_scores_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_superfamily_clr_pca_scores_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_clr_pca_loadings_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_superfamily_pca_zero_replacement_sensitivity_analysis18_v1.png`
"""


def main() -> None:
    occupied = [
        path
        for path in (
            DIVERSITY_OUTPUT,
            COMPOSITION_OUTPUT,
            CLR_OUTPUT,
            SCORES_OUTPUT,
            LOADINGS_OUTPUT,
            VARIANCE_OUTPUT,
            STABILITY_OUTPUT,
            MANIFEST_OUTPUT,
        )
        if path.exists()
    ]
    if occupied:
        raise FileExistsError("Non-destructive diversity/PCA audit requires unused outputs")
    species = pd.read_csv(PANEL)["species"].map(_canonical_species).tolist()
    if len(species) != 18 or len(set(species)) != 18:
        raise ValueError("Final panel must contain 18 unique species")
    mass = pd.read_csv(MASS).sort_values("species").reset_index(drop=True)
    if set(mass["species"]) != set(species):
        raise ValueError("Mass ledger differs from final panel")
    mass_index = mass.set_index("species")
    order = _load_breakdown(ORDER, species)
    superfamily = _load_breakdown(SUPERFAMILY, species)
    compositions = {
        ("order", "classified_conditional"): order.div(order.sum(axis=1), axis=0),
        ("order", "mass_aware_unresolved_bin"): build_mass_aware_composition(
            order,
            mass_index["order_retained_fraction"],
            mass_index["order_unresolved_fraction"],
        ),
        ("superfamily", "classified_conditional"): superfamily.div(
            superfamily.sum(axis=1), axis=0
        ),
        ("superfamily", "mass_aware_unresolved_bin"): build_mass_aware_composition(
            superfamily,
            mass_index["superfamily_retained_fraction"],
            mass_index["superfamily_unresolved_fraction"],
        ),
    }
    diversity = _diversity_table(compositions, mass)
    reproduction = _historical_reproduction(diversity)
    pca = _run_all_pca(compositions)

    composition_long = []
    for (level, mode), table in compositions.items():
        long = table.rename_axis("species").reset_index().melt(
            id_vars="species", var_name="feature", value_name="proportion"
        )
        long.insert(0, "composition_mode", mode)
        long.insert(0, "te_level", level)
        composition_long.append(long)
    composition_long = pd.concat(composition_long, ignore_index=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    diversity.to_csv(DIVERSITY_OUTPUT, index=False)
    composition_long.to_csv(COMPOSITION_OUTPUT, index=False)
    pca["clr"].to_csv(CLR_OUTPUT, index=False)
    pca["scores"].to_csv(SCORES_OUTPUT, index=False)
    pca["loadings"].to_csv(LOADINGS_OUTPUT, index=False)
    pca["variance"].to_csv(VARIANCE_OUTPUT, index=False)
    pca["stability"].to_csv(STABILITY_OUTPUT, index=False)
    figure_paths = _write_figures(diversity, pca)
    REPORT_OUTPUT.write_text(_render_report(diversity, pca, reproduction), encoding="utf-8")
    outputs = [
        DIVERSITY_OUTPUT,
        COMPOSITION_OUTPUT,
        CLR_OUTPUT,
        SCORES_OUTPUT,
        LOADINGS_OUTPUT,
        VARIANCE_OUTPUT,
        STABILITY_OUTPUT,
    ]
    manifest = {
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "n_species": 18,
        "out_of_panel_species_processed": False,
        "diversity_status": "approved_descriptive_conditional_on_classified_mass",
        "pca_status": "approved_descriptive_with_sensitivity_not_causal",
        "permanova_status": "not_approved_phylogenetic_exchangeability_violation",
        "historical_diversity_reproduction_max_abs_error": reproduction,
        "panel": _portable(PANEL),
        "panel_sha256": _sha256(PANEL),
        "mass_input": _portable(MASS),
        "mass_input_sha256": _sha256(MASS),
        "order_input": _portable(ORDER),
        "order_input_sha256": _sha256(ORDER),
        "superfamily_input": _portable(SUPERFAMILY),
        "superfamily_input_sha256": _sha256(SUPERFAMILY),
        "pca_analyses": pca["manifest"].to_dict("records"),
        "outputs": [{"path": _portable(path), "sha256": _sha256(path)} for path in outputs],
        "report": _portable(REPORT_OUTPUT),
        "report_sha256": _sha256(REPORT_OUTPUT),
        "figures": [{"path": _portable(path), "sha256": _sha256(path)} for path in figure_paths],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    for path in outputs + [MANIFEST_OUTPUT, REPORT_OUTPUT] + figure_paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()

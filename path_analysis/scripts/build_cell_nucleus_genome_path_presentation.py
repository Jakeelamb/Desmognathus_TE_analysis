#!/usr/bin/env python3
"""Render publication figures for the all-finalized-species path notebook."""

from __future__ import annotations

import hashlib
import json
import textwrap
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import nbformat as nbf
import numpy as np
import pandas as pd

import build_cell_nucleus_genome_path_notebook as design


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = design.OUTPUT_DIR
FIGURE_DIR = (
    PROJECT_ROOT
    / "results"
    / "figures"
    / "research_review"
    / "cell_nucleus_genome_path"
)
NOTEBOOK_PATH = (
    PROJECT_ROOT
    / "notebooks"
    / "research_review"
    / "07_cell_nucleus_genome_path_analysis.ipynb"
)
R_SCRIPT_PATH = (
    PROJECT_ROOT
    / "path_analysis"
    / "scripts"
    / "run_cell_nucleus_genome_phylogenetic_path_analysis.R"
)
REFERENCE_AUDIT_PATH = (
    PROJECT_ROOT / "path_analysis" / "CELL_NUCLEUS_GENOME_PATH_REFERENCE_METHOD_AUDIT.md"
)

PRIMARY_RANKING_PATH = INPUT_DIR / "primary_model_ranking.csv"
PRIMARY_DSEP_PATH = INPUT_DIR / "primary_dsep_tests.csv"
PAIRWISE_PATH = INPUT_DIR / "pairwise_phylogenetic_regressions.csv"
CORRELATION_PATH = INPUT_DIR / "raw_trait_correlations.csv"
MECHANISM_RANKING_PATH = INPUT_DIR / "user_mechanism_equivalence_ranking.csv"
MECHANISM_EDGES_PATH = INPUT_DIR / "user_mechanism_edges.csv"
EVOLUTIONARY_PATH = INPUT_DIR / "evolutionary_model_sensitivity_rankings.csv"
LOO_PATH = INPUT_DIR / "leave_one_species_out_rankings.csv.gz"
MEASUREMENT_PATH = INPUT_DIR / "measurement_bootstrap_rankings.csv.gz"
TREE_PATH = INPUT_DIR / "tree_sensitivity_rankings.csv.gz"
SIMULATION_RANKING_PATH = INPUT_DIR / "simulation_rankings.csv.gz"
SIMULATION_SUMMARY_PATH = INPUT_DIR / "simulation_summary.csv"
FAILURE_PATH = INPUT_DIR / "analysis_failures.csv"
ANALYSIS_MANIFEST_PATH = INPUT_DIR / "analysis_manifest.json"
STABILITY_SUMMARY_PATH = INPUT_DIR / "model_stability_summary.csv"
CONCLUSION_PATH = INPUT_DIR / "conclusion_summary.json"
FIGURE_MANIFEST_PATH = INPUT_DIR / "cell_nucleus_genome_path_figure_manifest.json"

FIGURE_PATHS = {
    "raw_relationships": FIGURE_DIR / "01_raw_trait_relationships.png",
    "model_ranking": FIGURE_DIR / "02_primary_model_ranking.png",
    "equivalent_mechanisms": FIGURE_DIR / "03_markov_equivalent_user_mechanisms.png",
    "candidate_gallery": FIGURE_DIR / "04_candidate_equivalence_class_gallery.png",
    "model_stability": FIGURE_DIR / "05_model_stability.png",
    "bridge_contrast": FIGURE_DIR / "06_nucleus_vs_cell_bridge_contrast.png",
    "sensitivity_heatmap": FIGURE_DIR / "07_structural_sensitivity_heatmap.png",
    "simulation": FIGURE_DIR / "08_simulation_calibration.png",
}

MODEL_ORDER = [
    "nucleus_bridge",
    "cell_bridge",
    "nucleus_cell_only",
    "nucleus_collider",
    "cell_collider",
    "genome_bridge",
    "genome_nucleus_only",
    "genome_cell_only",
    "genome_collider",
    "independent",
]
MODEL_COLORS = {
    "nucleus_bridge": "#2E7D6E",
    "cell_bridge": "#D17C33",
    "genome_bridge": "#5576A3",
    "nucleus_collider": "#8A6FA8",
    "cell_collider": "#B67887",
    "genome_collider": "#708599",
    "nucleus_cell_only": "#699B72",
    "genome_nucleus_only": "#7194B8",
    "genome_cell_only": "#B29A64",
    "independent": "#8A8F94",
    "saturated": "#C24A45",
}
MODEL_LABELS = {
    "nucleus_bridge": "Nucleus bridge",
    "cell_bridge": "Cell bridge",
    "genome_bridge": "Genome bridge",
    "nucleus_collider": "Nucleus collider",
    "cell_collider": "Cell collider",
    "genome_collider": "Genome collider",
    "nucleus_cell_only": "Nucleus–cell only",
    "genome_nucleus_only": "Genome–nucleus only",
    "genome_cell_only": "Genome–cell only",
    "independent": "Independent",
    "saturated": "Saturated (untestable)",
}
CANDIDATE_GALLERY_EXCLUDED_CLASSES = {"genome_cell_only", "saturated"}
TRAIT_LABELS = {
    "genome_size": "Genome size",
    "nucleus_size": "Nucleus size",
    "cell_size": "Cell size",
}
ABBREVIATIONS = {
    "genome_size": "GS",
    "nucleus_size": "NS",
    "cell_size": "CS",
}
PNG_METADATA = {"Software": "Desmognathus_TE audited path analysis"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
        metadata=PNG_METADATA,
    )
    fig.savefig(
        path.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
        metadata={"Creator": "Desmognathus_TE audited path analysis"},
    )
    plt.close(fig)


def load_analysis_tables() -> dict[str, Any]:
    paths = {
        "traits": INPUT_DIR / "primary_traits_all_finalized.csv",
        "dags": INPUT_DIR / "all_three_trait_dags.csv",
        "primary": PRIMARY_RANKING_PATH,
        "dsep": PRIMARY_DSEP_PATH,
        "pairwise": PAIRWISE_PATH,
        "correlations": CORRELATION_PATH,
        "mechanisms": MECHANISM_RANKING_PATH,
        "mechanism_edges": MECHANISM_EDGES_PATH,
        "evolutionary": EVOLUTIONARY_PATH,
        "loo": LOO_PATH,
        "measurement": MEASUREMENT_PATH,
        "trees": TREE_PATH,
        "simulation_rankings": SIMULATION_RANKING_PATH,
        "simulation": SIMULATION_SUMMARY_PATH,
        "failures": FAILURE_PATH,
    }
    missing = [path for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing fitted path-analysis output(s): " + ", ".join(map(str, missing))
        )
    tables = {name: pd.read_csv(path, low_memory=False) for name, path in paths.items()}
    tables["manifest"] = json.loads(ANALYSIS_MANIFEST_PATH.read_text())
    return tables


def _winner_summary(rankings: pd.DataFrame, analysis: str) -> pd.DataFrame:
    if rankings.empty:
        return pd.DataFrame()
    winners = rankings.loc[rankings["rank"].eq(1), ["replicate_id", "model"]]
    competitive = (
        rankings.loc[rankings["admissible_competitive"].astype(bool)]
        .groupby("model")["replicate_id"]
        .nunique()
    )
    total = winners["replicate_id"].nunique()
    counts = winners["model"].value_counts()
    models = sorted(set(rankings["model"]), key=lambda value: MODEL_ORDER.index(value))
    return pd.DataFrame(
        {
            "analysis": analysis,
            "model": models,
            "n_fits": total,
            "top_count": [int(counts.get(model, 0)) for model in models],
            "top_rate": [float(counts.get(model, 0) / total) for model in models],
            "admissible_competitive_count": [
                int(competitive.get(model, 0)) for model in models
            ],
            "admissible_competitive_rate": [
                float(competitive.get(model, 0) / total) for model in models
            ],
        }
    )


def build_stability_summary(tables: dict[str, Any]) -> pd.DataFrame:
    tree_bootstrap = tables["trees"].loc[
        tables["trees"]["analysis_type"].eq("published_bootstrap_tree_sensitivity")
    ]
    groups = [
        _winner_summary(tables["measurement"], "measurement_bootstrap_primary"),
        _winner_summary(tables["loo"], "leave_one_species_out_24"),
        _winner_summary(tree_bootstrap, "published_bootstrap_trees_18"),
        _winner_summary(tables["evolutionary"], "evolutionary_model_choice_primary"),
    ]
    summary = pd.concat(groups, ignore_index=True)
    summary.to_csv(STABILITY_SUMMARY_PATH, index=False)
    return summary


def _bridge_contrast(rankings: pd.DataFrame) -> pd.DataFrame:
    pivot = rankings.pivot_table(
        index="replicate_id", columns="model", values="CICc", aggfunc="first"
    )
    required = {"nucleus_bridge", "cell_bridge"}
    if not required.issubset(pivot.columns):
        raise ValueError("Ranking table is missing a bridge equivalence class")
    contrast = pivot["nucleus_bridge"] - pivot["cell_bridge"]
    return contrast.rename("nucleus_minus_cell_CICc").reset_index()


def build_conclusion_summary(
    tables: dict[str, Any], stability: pd.DataFrame
) -> dict[str, Any]:
    primary = tables["primary"].sort_values("rank")
    competitive = primary.loc[primary["admissible_competitive"].astype(bool), "model"]
    measurement_contrast = _bridge_contrast(tables["measurement"])
    loo_contrast = _bridge_contrast(tables["loo"])
    tree_bootstrap = tables["trees"].loc[
        tables["trees"]["analysis_type"].eq("published_bootstrap_tree_sensitivity")
    ]
    top_tree = tree_bootstrap.loc[tree_bootstrap["rank"].eq(1), "model"].value_counts()
    evolutionary_winners = (
        tables["evolutionary"]
        .loc[tables["evolutionary"]["rank"].eq(1), ["evolutionary_model", "model"]]
        .set_index("evolutionary_model")["model"]
        .to_dict()
    )
    measurement_winners = tables["measurement"].loc[
        tables["measurement"]["rank"].eq(1), "model"
    ]
    loo_winners = tables["loo"].loc[tables["loo"]["rank"].eq(1), "model"]
    result = {
        "primary_top_class": str(primary.iloc[0]["model"]),
        "primary_top_global_p": float(primary.iloc[0]["p"]),
        "primary_top_cicc": float(primary.iloc[0]["CICc"]),
        "primary_top_weight": float(primary.iloc[0]["w"]),
        "primary_second_class": str(primary.iloc[1]["model"]),
        "primary_second_delta_cicc": float(primary.iloc[1]["delta_CICc"]),
        "primary_second_weight": float(primary.iloc[1]["w"]),
        "primary_admissible_competitive_classes": competitive.tolist(),
        "measurement_bootstrap_replicates": int(measurement_contrast.shape[0]),
        "measurement_nucleus_bridge_top_rate": float(
            measurement_winners.eq("nucleus_bridge").mean()
        ),
        "measurement_cell_bridge_top_rate": float(
            measurement_winners.eq("cell_bridge").mean()
        ),
        "measurement_median_nucleus_minus_cell_cicc": float(
            measurement_contrast["nucleus_minus_cell_CICc"].median()
        ),
        "leave_one_out_fits": int(loo_contrast.shape[0]),
        "leave_one_out_nucleus_bridge_top_rate": float(
            loo_winners.eq("nucleus_bridge").mean()
        ),
        "leave_one_out_cell_bridge_top_rate": float(
            loo_winners.eq("cell_bridge").mean()
        ),
        "published_bootstrap_tree_fits": int(top_tree.sum()),
        "published_bootstrap_tree_top_counts": {
            str(key): int(value) for key, value in top_tree.items()
        },
        "evolutionary_model_top_classes": evolutionary_winners,
        "causal_direction_identified": False,
        "user_mechanisms_markov_equivalent": True,
        "shared_testable_claim": (
            "genome_size is independent of cell_size conditional on nucleus_size"
        ),
        "iod_nucleus_coupling_caveat": (
            "The genome estimate is derived from nuclear integrated optical density, "
            "which contains a nuclear-area component; genome–nucleus evidence is not "
            "an independent causal validation."
        ),
        "reportable_verdict": (
            f"The data support a connected three-trait system; the point estimate "
            f"favors {_friendly_model(str(primary.iloc[0]['model']))} "
            f"(second-class delta CICc = {float(primary.iloc[1]['delta_CICc']):.2f}). "
            "The cross-sectional data still cannot orient Markov-equivalent causal arrows."
        ),
    }
    CONCLUSION_PATH.write_text(json.dumps(result, indent=2) + "\n")
    return result


def _friendly_model(model: str) -> str:
    return MODEL_LABELS.get(model, model.replace("_", " ").title())


def _natural_join(labels: list[str]) -> str:
    if not labels:
        return ""
    if len(labels) == 1:
        return labels[0]
    if len(labels) == 2:
        return f"{labels[0]} and {labels[1]}"
    return f"{', '.join(labels[:-1])}, and {labels[-1]}"


def primary_ranking_footer(ranking: pd.DataFrame) -> str:
    """Describe support from the current fit instead of a historical result."""
    supported = ranking["globally_supported"].fillna(False).astype(bool)
    competitive = supported & ranking["delta_CICc"].le(2)
    supported_noncompetitive = supported & ~competitive

    competitive_labels = [
        _friendly_model(model)
        for model in ranking.loc[competitive, "model"].astype(str)
    ]
    noncompetitive_labels = [
        _friendly_model(model)
        for model in ranking.loc[supported_noncompetitive, "model"].astype(str)
    ]

    if len(competitive_labels) == 1:
        support_sentence = (
            f"{competitive_labels[0]} is the sole globally supported class "
            "within ΔCICc ≤ 2."
        )
    elif competitive_labels:
        support_sentence = (
            f"{_natural_join(competitive_labels)} are globally supported classes "
            "within ΔCICc ≤ 2."
        )
    else:
        support_sentence = "No class is globally supported within ΔCICc ≤ 2."

    if noncompetitive_labels:
        remainder_sentence = (
            f"{_natural_join(noncompetitive_labels)} pass global fit but are not competitive."
        )
    else:
        remainder_sentence = "All other testable classes fail global fit."

    return (
        f"{support_sentence} {remainder_sentence} "
        "The saturated class has no rankable conditional-independence test."
    )


def render_raw_relationships(tables: dict[str, Any]) -> None:
    traits = tables["traits"]
    correlations = tables["correlations"].set_index("comparison")
    pairwise = tables["pairwise"].set_index("relation_id")
    specs = [
        (
            "genome_nucleus",
            "genome_to_nucleus",
            "genome_size_pg",
            "nucleus_area_um2",
            "Genome size (pg; D. fuscus = 16.36)",
            "Nucleus area (µm²)",
            "#3C78A8",
        ),
        (
            "genome_cell",
            "genome_to_cell",
            "genome_size_pg",
            "cell_area_um2",
            "Genome size (pg; D. fuscus = 16.36)",
            "Cell area (µm²)",
            "#D16D45",
        ),
        (
            "nucleus_cell",
            "nucleus_to_cell",
            "nucleus_area_um2",
            "cell_area_um2",
            "Nucleus area (µm²)",
            "Cell area (µm²)",
            "#765B9E",
        ),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.4))
    for panel_index, (comparison, relation, x_col, y_col, x_label, y_label, color) in enumerate(specs):
        ax = axes[panel_index]
        x = traits[x_col].to_numpy(float)
        y = traits[y_col].to_numpy(float)
        ax.scatter(x, y, s=55, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        for index, row in traits.reset_index(drop=True).iterrows():
            dx = 3 if index % 2 == 0 else -3
            align = "left" if dx > 0 else "right"
            ax.annotate(
                row["tree_tip"],
                (row[x_col], row[y_col]),
                xytext=(dx, 4 + (index % 3) * 2),
                textcoords="offset points",
                ha=align,
                va="bottom",
                fontsize=6.5,
                color="#43484C",
                fontstyle="italic",
            )
        regression = pairwise.loc[relation]
        if isinstance(regression, pd.DataFrame):
            regression = regression.iloc[0]
        slope = float(regression["standardized_coefficient"])
        slope_low = float(regression["ci_low"])
        slope_high = float(regression["ci_high"])
        log_x = np.log10(x)
        log_y = np.log10(y)
        x_grid = np.linspace(x.min() * 0.985, x.max() * 1.015, 250)
        z_x = (np.log10(x_grid) - log_x.mean()) / log_x.std(ddof=1)

        def backtransform(beta: float) -> np.ndarray:
            z_y = beta * z_x
            return 10 ** (log_y.mean() + log_y.std(ddof=1) * z_y)

        estimate = backtransform(slope)
        low_line = backtransform(slope_low)
        high_line = backtransform(slope_high)
        ax.fill_between(
            x_grid,
            np.minimum(low_line, high_line),
            np.maximum(low_line, high_line),
            color=color,
            alpha=0.12,
            linewidth=0,
        )
        ax.plot(x_grid, estimate, color=color, linewidth=2.1)
        corr = correlations.loc[comparison]
        ax.text(
            0.03,
            0.97,
            (
                f"raw Spearman ρ = {corr['spearman_rho']:.2f} "
                f"(P = {corr['spearman_p']:.3f})\n"
                f"PGLS βstd = {slope:.2f} "
                f"[{slope_low:.2f}, {slope_high:.2f}]\n"
                f"P = {regression['p_value']:.3f}; λ = "
                f"{regression['fitted_phylogenetic_parameter']:.1g}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9.2,
            bbox={"facecolor": "white", "edgecolor": "#D8DADD", "alpha": 0.92},
        )
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(color="#D8DBDE", alpha=0.45, linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "All raw species estimates: genome size, nucleus size, and cell size",
        fontsize=16,
        fontweight="bold",
        y=1.01,
    )
    fig.text(
        0.5,
        -0.015,
        "Points are all finalized species. Lines are standardized log10 PGLS fits transformed back to raw units; shading varies the slope over its approximate 95% interval.",
        ha="center",
        fontsize=9,
        color="#4B5055",
    )
    fig.tight_layout()
    _save_figure(fig, FIGURE_PATHS["raw_relationships"])


def render_primary_ranking(tables: dict[str, Any]) -> None:
    ranking = tables["primary"].sort_values("rank").reset_index(drop=True)
    labels = [_friendly_model(model) for model in ranking["model"]]
    colors = [MODEL_COLORS[model] for model in ranking["model"]]
    y = np.arange(len(ranking))
    fig, (delta_ax, weight_ax) = plt.subplots(
        1, 2, figsize=(14.5, 7.3), gridspec_kw={"width_ratios": [1.2, 1]}
    )
    delta_ax.barh(y, ranking["delta_CICc"], color=colors, alpha=0.82, height=0.62)
    delta_ax.axvline(2, color="#B03A37", linestyle="--", linewidth=1.4)
    for index, row in ranking.iterrows():
        marker = "pass" if bool(row["globally_supported"]) else "reject"
        delta_ax.text(
            max(float(row["delta_CICc"]), 0) + 0.18,
            index,
            f"{row['delta_CICc']:.2f}   global P={row['p']:.3f} {marker}",
            va="center",
            fontsize=8.4,
            color="#33383C" if bool(row["globally_supported"]) else "#A13B38",
        )
    delta_ax.set_yticks(y, labels)
    delta_ax.invert_yaxis()
    delta_ax.set_xlabel("ΔCICc (lower is better; dashed line = 2)")
    delta_ax.set_title("Small-sample model comparison", fontweight="bold")
    delta_ax.grid(axis="x", color="#D8DBDE", alpha=0.5)
    delta_ax.spines[["top", "right", "left"]].set_visible(False)

    weight_ax.barh(y, ranking["w"], color=colors, alpha=0.82, height=0.62)
    for index, value in enumerate(ranking["w"]):
        weight_ax.text(value + 0.008, index, f"{value:.3f}", va="center", fontsize=8.6)
    weight_ax.set_yticks(y, [""] * len(y))
    weight_ax.invert_yaxis()
    weight_ax.set_xlabel("CICc weight")
    weight_ax.set_title("Relative support within the tested set", fontweight="bold")
    weight_ax.grid(axis="x", color="#D8DBDE", alpha=0.5)
    weight_ax.spines[["top", "right", "left"]].set_visible(False)
    fig.suptitle(
        "All-finalized-species phylogenetic path ranking",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.015,
        primary_ranking_footer(ranking),
        ha="center",
        fontsize=9.2,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    _save_figure(fig, FIGURE_PATHS["model_ranking"])


def _draw_box_node(
    ax: plt.Axes,
    position: tuple[float, float],
    label: str,
    *,
    width: float = 0.18,
    height: float = 0.15,
) -> None:
    x, y = position
    rectangle = Rectangle(
        (x - width / 2, y - height / 2),
        width,
        height,
        facecolor="white",
        edgecolor="#111111",
        linewidth=1.5,
        zorder=4,
    )
    ax.add_patch(rectangle)
    ax.text(x, y, label, ha="center", va="center", fontsize=15, zorder=5)


def _draw_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    label: str | None = None,
    label_offset: float = 0.08,
    linewidth: float = 2.7,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=18,
        linewidth=linewidth,
        color="#111111",
        shrinkA=28,
        shrinkB=28,
        zorder=3,
    )
    ax.add_patch(arrow)
    if label:
        midpoint = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        ax.text(
            midpoint[0],
            midpoint[1] + label_offset,
            label,
            ha="center",
            va="center",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
            zorder=6,
        )


def render_equivalent_mechanisms(tables: dict[str, Any]) -> None:
    ranking = tables["mechanisms"].set_index("model")
    edges = tables["mechanism_edges"]
    mechanisms = [
        ("genome_to_nucleus_to_cell", "Genome → nucleus → cell"),
        ("cell_to_nucleus_to_genome", "Cell → nucleus → genome"),
        ("nucleus_to_genome_and_cell", "Nucleus → genome and cell"),
    ]
    positions = {
        "genome_size": (0.16, 0.54),
        "nucleus_size": (0.50, 0.54),
        "cell_size": (0.84, 0.54),
    }
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8))
    for ax, (mechanism, title) in zip(axes, mechanisms):
        for trait, position in positions.items():
            _draw_box_node(ax, position, ABBREVIATIONS[trait])
        for _, edge in edges.loc[edges["mechanism"].eq(mechanism)].iterrows():
            label = (
                f"β={edge['standardized_coefficient']:.2f}\n"
                f"SE={edge['standard_error']:.2f}"
            )
            _draw_arrow(
                ax,
                positions[str(edge["parent"])],
                positions[str(edge["child"])],
                label=label,
                label_offset=0.10 if edge["parent"] != "cell_size" else -0.12,
            )
        row = ranking.loc[mechanism]
        ax.text(
            0.5,
            0.15,
            (
                f"global P={row['p']:.3f}   CICc={row['CICc']:.2f}\n"
                f"Δ={row['delta_CICc']:.2f}   duplicated-set weight={row['w']:.3f}"
            ),
            ha="center",
            va="center",
            fontsize=9.2,
            color="#353A3E",
        )
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    fig.suptitle(
        "Three biological stories, one observational equivalence class",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.91,
        "All imply GS independent of CS given NS and therefore receive identical Fisher C, global P, and CICc. Coefficients change orientation only because the arrows were assumed.",
        ha="center",
        fontsize=9.5,
        color="#4B5055",
    )
    fig.text(
        0.5,
        0.025,
        "GS, calibrated genome-size estimate; NS, median nucleus area of vetted top-50 cells; CS, median cell area of vetted top-50 cells. Equal 1/3 weights are an artifact of listing the same equivalence class three times.",
        ha="center",
        fontsize=9,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.88))
    _save_figure(fig, FIGURE_PATHS["equivalent_mechanisms"])


def _parse_dag_id(dag_id: str) -> list[tuple[str, str]]:
    if dag_id == "no_edges":
        return []
    mapping = {"GS": "genome_size", "NS": "nucleus_size", "CS": "cell_size"}
    edges = []
    for token in dag_id.split(";"):
        source, target = token.split("->")
        edges.append((mapping[source], mapping[target]))
    return edges


def candidate_gallery_rows(dags: pd.DataFrame) -> pd.DataFrame:
    """Return the biologically retained DAG panels shown in the gallery.

    The fitted analysis remains exhaustive.  This is only the manuscript-facing
    display contract: the genome-cell-only skeleton is not biologically plausible
    for this study, the saturated class is untestable, and both focal chain
    orientations are shown even though they belong to one equivalence class.
    """

    representatives = (
        dags.sort_values(["n_edges", "dag_id"])
        .groupby("equivalence_class", sort=False)
        .first()
        .reset_index()
    )
    representatives = representatives.loc[
        ~representatives["equivalence_class"].isin(
            CANDIDATE_GALLERY_EXCLUDED_CLASSES
        )
    ].copy()
    representatives["gallery_id"] = representatives["equivalence_class"]
    representatives["gallery_title"] = representatives["equivalence_class"].map(
        _friendly_model
    )

    reverse_chain_id = design.USER_MECHANISMS["cell_to_nucleus_to_genome"]
    reverse_chain = dags.loc[dags["dag_id"].eq(reverse_chain_id)].copy()
    if len(reverse_chain) != 1:
        raise ValueError("Expected one Cell -> nucleus -> genome DAG")
    reverse_chain["gallery_id"] = "cell_to_nucleus_to_genome"
    reverse_chain["gallery_title"] = "Cell → nucleus → genome"
    representatives = pd.concat([representatives, reverse_chain], ignore_index=True)

    class_order = [
        "independent",
        "genome_nucleus_only",
        "nucleus_cell_only",
        "nucleus_bridge",
        "cell_to_nucleus_to_genome",
        "cell_bridge",
        "genome_bridge",
        "nucleus_collider",
        "cell_collider",
        "genome_collider",
    ]
    representatives["order"] = representatives["gallery_id"].map(
        {value: index for index, value in enumerate(class_order)}
    )
    return representatives.sort_values("order").reset_index(drop=True)


def render_candidate_gallery(tables: dict[str, Any]) -> None:
    representatives = candidate_gallery_rows(tables["dags"].copy())
    positions = {
        "genome_size": (0.50, 0.78),
        "nucleus_size": (0.22, 0.37),
        "cell_size": (0.78, 0.37),
    }
    fig, axes = plt.subplots(2, 5, figsize=(19, 8.5))
    flat_axes = axes.ravel()
    for ax, (_, row) in zip(flat_axes, representatives.iterrows()):
        for trait, position in positions.items():
            _draw_box_node(ax, position, ABBREVIATIONS[trait], width=0.19, height=0.13)
        for source, target in _parse_dag_id(str(row["dag_id"])):
            _draw_arrow(
                ax,
                positions[source],
                positions[target],
                linewidth=2.0,
            )
        ax.set_title(str(row["gallery_title"]), fontsize=11, fontweight="bold")
        basis = str(row["basis_claim"])
        basis = (
            basis.replace("genome_size", "GS")
            .replace("nucleus_size", "NS")
            .replace("cell_size", "CS")
            .replace("_||_", "independent of")
        )
        ax.text(
            0.5,
            0.08,
            textwrap.fill(basis, 34),
            ha="center",
            va="center",
            fontsize=7.8,
            color="#4B5055",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
    for ax in flat_axes[len(representatives) :]:
        ax.axis("off")
    fig.suptitle(
        "Biologically retained three-trait DAG gallery",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.955,
        "Genome-cell-only is excluded as biologically implausible here; the saturated class is excluded because it is untestable. Both focal nucleus-bridge chain orientations are shown.",
        ha="center",
        fontsize=9.5,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save_figure(fig, FIGURE_PATHS["candidate_gallery"])


def _plot_winner_rates(
    ax: plt.Axes,
    ranking: pd.DataFrame,
    title: str,
    denominator_label: str,
) -> None:
    winners = ranking.loc[ranking["rank"].eq(1), "model"]
    rates = winners.value_counts(normalize=True)
    ordered = [model for model in MODEL_ORDER if model in rates.index]
    values = [float(rates[model]) for model in ordered]
    y = np.arange(len(ordered))
    ax.barh(y, values, color=[MODEL_COLORS[model] for model in ordered], height=0.62)
    for index, value in enumerate(values):
        ax.text(value + 0.015, index, f"{100 * value:.1f}%", va="center", fontsize=8.6)
    ax.set_yticks(y, [_friendly_model(model) for model in ordered])
    ax.invert_yaxis()
    ax.set_xlim(0, 1.08)
    ax.set_xlabel(f"Top-ranked frequency ({denominator_label}; n={len(winners)})")
    ax.set_title(title, fontweight="bold")
    ax.grid(axis="x", color="#D8DBDE", alpha=0.5)
    ax.spines[["top", "right", "left"]].set_visible(False)


def render_model_stability(tables: dict[str, Any]) -> None:
    tree_bootstrap = tables["trees"].loc[
        tables["trees"]["analysis_type"].eq("published_bootstrap_tree_sensitivity")
    ]
    fig, axes = plt.subplots(2, 2, figsize=(15.5, 10.5))
    _plot_winner_rates(
        axes[0, 0],
        tables["measurement"],
        "Measurement resampling (all finalized species)",
        "bootstrap draws",
    )
    _plot_winner_rates(
        axes[0, 1],
        tables["loo"],
        "Leave one species out (19 species)",
        "omissions",
    )
    _plot_winner_rates(
        axes[1, 0],
        tree_bootstrap,
        "Published trees (fixed 18-species set)",
        "trees",
    )
    evolutionary = tables["evolutionary"].loc[
        tables["evolutionary"]["rank"].eq(1), ["evolutionary_model", "model", "w"]
    ]
    evolutionary = evolutionary.reset_index(drop=True)
    y = np.arange(len(evolutionary))
    axes[1, 1].barh(
        y,
        evolutionary["w"],
        color=[MODEL_COLORS[model] for model in evolutionary["model"]],
        height=0.62,
    )
    for index, row in evolutionary.iterrows():
        axes[1, 1].text(
            row["w"] + 0.012,
            index,
            f"{_friendly_model(row['model'])}  w={row['w']:.2f}",
            va="center",
            fontsize=8.8,
        )
    axes[1, 1].set_yticks(y, evolutionary["evolutionary_model"])
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlim(0, max(0.55, evolutionary["w"].max() + 0.18))
    axes[1, 1].set_xlabel("Weight of top class")
    axes[1, 1].set_title("Evolutionary residual model (all finalized species)", fontweight="bold")
    axes[1, 1].grid(axis="x", color="#D8DBDE", alpha=0.5)
    axes[1, 1].spines[["top", "right", "left"]].set_visible(False)
    fig.suptitle(
        "The structural winner is not stable across reasonable perturbations",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.015,
        "Focal-tree and published-tree fits are identical at n=18 because fitted lambda is approximately zero; the 20-to-18 species-set change, not bootstrap topology, drives that flip.",
        ha="center",
        fontsize=9,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    _save_figure(fig, FIGURE_PATHS["model_stability"])


def render_bridge_contrast(tables: dict[str, Any]) -> None:
    measurement = _bridge_contrast(tables["measurement"])
    loo = _bridge_contrast(tables["loo"])
    primary = _bridge_contrast(tables["primary"])["nucleus_minus_cell_CICc"].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw={"width_ratios": [1, 1.25]})
    values = measurement["nucleus_minus_cell_CICc"]
    axes[0].hist(values, bins=24, color="#6C8FA8", alpha=0.86, edgecolor="white")
    axes[0].axvline(0, color="#222222", linewidth=1.2)
    axes[0].axvline(primary, color="#B03A37", linestyle="--", linewidth=2)
    axes[0].text(
        0.03,
        0.97,
        (
            f"Nucleus bridge favored: {(values < 0).mean():.1%}\n"
            f"Cell bridge favored: {(values > 0).mean():.1%}\n"
            f"Point estimate: {primary:+.2f}"
        ),
        transform=axes[0].transAxes,
        va="top",
        fontsize=9.5,
        bbox={"facecolor": "white", "edgecolor": "#D8DADD", "alpha": 0.9},
    )
    axes[0].set_xlabel("CICc(nucleus bridge) − CICc(cell bridge)")
    axes[0].set_ylabel("Measurement-bootstrap replicates")
    axes[0].set_title("Measurement uncertainty", fontweight="bold")
    axes[0].grid(axis="y", color="#D8DBDE", alpha=0.45)
    axes[0].spines[["top", "right"]].set_visible(False)

    loo = loo.sort_values("nucleus_minus_cell_CICc")
    colors = np.where(
        loo["nucleus_minus_cell_CICc"] < 0,
        MODEL_COLORS["nucleus_bridge"],
        MODEL_COLORS["cell_bridge"],
    )
    y = np.arange(len(loo))
    axes[1].barh(y, loo["nucleus_minus_cell_CICc"], color=colors, height=0.68)
    axes[1].axvline(0, color="#222222", linewidth=1.2)
    axes[1].axvline(primary, color="#B03A37", linestyle="--", linewidth=1.6)
    axes[1].set_yticks(y, [f"omit {value}" for value in loo["replicate_id"]])
    axes[1].set_xlabel("CICc(nucleus bridge) − CICc(cell bridge)")
    axes[1].set_title("Species influence", fontweight="bold")
    axes[1].grid(axis="x", color="#D8DBDE", alpha=0.45)
    axes[1].spines[["top", "right", "left"]].set_visible(False)
    fig.suptitle(
        "Nucleus-bridge versus cell-bridge support",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.015,
        "Negative values favor nucleus as the bridge; positive values favor cell as the bridge. The dashed line is the full 24-species point estimate.",
        ha="center",
        fontsize=9.2,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    _save_figure(fig, FIGURE_PATHS["bridge_contrast"])


def render_sensitivity_heatmap(tables: dict[str, Any]) -> None:
    rows: list[pd.DataFrame] = []
    for model, label in [
        ("lambda", "Lambda · focal 20"),
        ("BM", "BM · focal 20"),
        ("OUfixedRoot", "OU fixed root · focal 20"),
        ("OUrandomRoot", "OU random root · focal 20"),
    ]:
        subset = tables["evolutionary"].loc[
            tables["evolutionary"]["evolutionary_model"].eq(model)
        ].copy()
        subset["sensitivity_label"] = label
        rows.append(subset)
    for tree_id, label in [
        ("focal_tree_18", "Focal tree · 18 species"),
        ("published_main_18", "Published tree · 18 species"),
    ]:
        subset = tables["trees"].loc[tables["trees"]["tree_id"].eq(tree_id)].copy()
        subset["sensitivity_label"] = label
        rows.append(subset)
    data = pd.concat(rows, ignore_index=True)
    row_order = [frame["sensitivity_label"].iloc[0] for frame in rows]
    matrix = data.pivot(index="sensitivity_label", columns="model", values="delta_CICc")
    support = data.pivot(index="sensitivity_label", columns="model", values="globally_supported")
    matrix = matrix.reindex(index=row_order, columns=MODEL_ORDER)
    support = support.reindex(index=row_order, columns=MODEL_ORDER)
    fig, ax = plt.subplots(figsize=(15.5, 6.8))
    image = ax.imshow(np.minimum(matrix.to_numpy(float), 10), cmap="YlGnBu", aspect="auto", vmin=0, vmax=10)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix.iloc[row_index, column_index]
            if not np.isfinite(value):
                continue
            rejected = not bool(support.iloc[row_index, column_index])
            label = f"{value:.1f}" + (" ×" if rejected else "")
            ax.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=7.7,
                color="white" if min(value, 10) > 6 else "#222222",
            )
    ax.set_xticks(np.arange(len(MODEL_ORDER)), [_friendly_model(model) for model in MODEL_ORDER], rotation=42, ha="right")
    ax.set_yticks(np.arange(len(row_order)), row_order)
    ax.set_title(
        "ΔCICc across evolutionary models, trees, and species sets",
        fontsize=15,
        fontweight="bold",
        pad=14,
    )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.02)
    colorbar.set_label("ΔCICc (values ≥10 share the darkest color)")
    fig.text(
        0.5,
        0.015,
        "× marks a globally rejected model (Fisher-C P < 0.05). The 18-species rows exclude D. brimleyorum and D. folkertsi.",
        ha="center",
        fontsize=9.2,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _save_figure(fig, FIGURE_PATHS["sensitivity_heatmap"])


def render_simulation(tables: dict[str, Any]) -> None:
    summary = tables["simulation"].copy()
    class_order = ["independent", "nucleus_bridge", "cell_bridge", "genome_bridge"]
    regime_order = [regime for regime in ["iid", "brownian"] if regime in set(summary["residual_regime"])]
    metrics = [
        ("true_class_top_rate", "True class ranked first"),
        ("true_class_delta_le_2_rate", "True class retained at ΔCICc ≤ 2"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.6), sharey=True)
    x = np.arange(len(class_order))
    width = 0.34
    for ax, (metric, title) in zip(axes, metrics):
        for regime_index, regime in enumerate(regime_order):
            values = []
            for class_name in class_order:
                match = summary.loc[
                    summary["generating_class"].eq(class_name)
                    & summary["residual_regime"].eq(regime),
                    metric,
                ]
                values.append(float(match.iloc[0]))
            positions = x + (regime_index - (len(regime_order) - 1) / 2) * width
            bars = ax.bar(
                positions,
                values,
                width=width,
                label="IID residuals" if regime == "iid" else "Brownian residuals",
                color="#4D83A6" if regime == "iid" else "#B56D46",
                alpha=0.86,
            )
            for bar, value in zip(bars, values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + 0.025,
                    f"{100 * value:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=8.3,
                )
        ax.set_xticks(x, [_friendly_model(value) for value in class_order], rotation=24, ha="right")
        ax.set_ylim(0, 1.08)
        ax.set_ylabel("Proportion of simulations")
        ax.set_title(title, fontweight="bold")
        ax.grid(axis="y", color="#D8DBDE", alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[1].legend(frameon=False, loc="upper right")
    n_per = int(summary["n_simulations"].iloc[0])
    fig.suptitle(
        f"Actual-tree simulation: how often can n = {len(tables['traits'])} recover the generating class?",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.015,
        f"{n_per} simulations per generating class × residual regime. Bridge effect sizes match the observed standardized bivariate paths; every dataset is analyzed with the primary lambda model set.",
        ha="center",
        fontsize=9.2,
        color="#4B5055",
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    _save_figure(fig, FIGURE_PATHS["simulation"])


def render_all_figures(tables: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    render_raw_relationships(tables)
    render_primary_ranking(tables)
    render_equivalent_mechanisms(tables)
    render_candidate_gallery(tables)
    render_model_stability(tables)
    render_bridge_contrast(tables)
    render_sensitivity_heatmap(tables)
    render_simulation(tables)


def validate_full_release(tables: dict[str, Any]) -> None:
    """Reject quick-run or incomplete inputs before producing share-facing figures."""

    manifest = tables["manifest"]
    expected = {
        "run_scope": "full_release",
        "primary_species": 24,
        "measurement_bootstrap_replicates": 250,
        "published_bootstrap_trees_analyzed": 200,
        "simulation_replicates_per_class_and_regime": 100,
        "failure_count": 0,
        "causal_direction_identified": False,
        "user_mechanisms_markov_equivalent": True,
    }
    mismatches = {
        key: {"expected": value, "observed": manifest.get(key)}
        for key, value in expected.items()
        if manifest.get(key) != value
    }
    if mismatches:
        raise RuntimeError(f"Common-support path release gate failed: {mismatches}")
    if not tables["failures"].empty:
        raise RuntimeError("Common-support path analysis contains failed fits")
    primary = tables["primary"].sort_values("rank")
    if len(primary) != 10 or not bool(primary.iloc[0]["admissible_competitive"]):
        raise RuntimeError("Primary ranking must contain ten classes and an admissible winner")


def build() -> dict[str, Any]:
    tables = load_analysis_tables()
    validate_full_release(tables)
    stability = build_stability_summary(tables)
    conclusion = build_conclusion_summary(tables, stability)
    render_all_figures(tables)

    figure_records = []
    for role, png_path in FIGURE_PATHS.items():
        pdf_path = png_path.with_suffix(".pdf")
        for path in (png_path, pdf_path):
            if not path.exists():
                raise FileNotFoundError(path)
        figure_records.append(
            {
                "role": role,
                "png": str(png_path.relative_to(PROJECT_ROOT)),
                "png_sha256": sha256_file(png_path),
                "pdf": str(pdf_path.relative_to(PROJECT_ROOT)),
                "pdf_sha256": sha256_file(pdf_path),
            }
        )

    manifest = {
        "analysis_id": "cell_nucleus_genome_path_figures_v1",
        "source_analysis": str(ANALYSIS_MANIFEST_PATH.relative_to(PROJECT_ROOT)),
        "source_analysis_sha256": sha256_file(ANALYSIS_MANIFEST_PATH),
        "statistical_script": str(R_SCRIPT_PATH.relative_to(PROJECT_ROOT)),
        "statistical_script_sha256": sha256_file(R_SCRIPT_PATH),
        "reference_method_audit": str(REFERENCE_AUDIT_PATH.relative_to(PROJECT_ROOT)),
        "reference_method_audit_sha256": sha256_file(REFERENCE_AUDIT_PATH),
        "source_script": str(Path(__file__).resolve().relative_to(PROJECT_ROOT)),
        "release_gates_passed": True,
        "primary_top_class": conclusion["primary_top_class"],
        "primary_second_class": conclusion["primary_second_class"],
        "primary_second_delta_cicc": conclusion["primary_second_delta_cicc"],
        "causal_direction_identified": False,
        "figures": figure_records,
        "derived_tables": {
            "stability_summary": str(STABILITY_SUMMARY_PATH.relative_to(PROJECT_ROOT)),
            "stability_summary_sha256": sha256_file(STABILITY_SUMMARY_PATH),
            "conclusion_summary": str(CONCLUSION_PATH.relative_to(PROJECT_ROOT)),
            "conclusion_summary_sha256": sha256_file(CONCLUSION_PATH),
        },
    }
    FIGURE_MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def _markdown(text: str, cell_id: str) -> nbf.NotebookNode:
    cell = nbf.v4.new_markdown_cell(textwrap.dedent(text).strip() + "\n")
    cell["id"] = cell_id
    return cell


def _code(text: str, cell_id: str) -> nbf.NotebookNode:
    cell = nbf.v4.new_code_cell(textwrap.dedent(text).strip() + "\n")
    cell["id"] = cell_id
    return cell


def build_notebook(conclusion: dict[str, Any]) -> nbf.NotebookNode:
    analysis_manifest = json.loads(ANALYSIS_MANIFEST_PATH.read_text())
    primary_species = int(analysis_manifest["primary_species"])
    primary_top = _friendly_model(str(conclusion["primary_top_class"]))
    second = _friendly_model(str(conclusion["primary_second_class"]))
    second_status = (
        "remains competitive"
        if conclusion["primary_second_delta_cicc"] <= 2
        else "is not competitive at the prespecified delta CICc <= 2 threshold"
    )
    measurement_nucleus = 100 * conclusion["measurement_nucleus_bridge_top_rate"]
    measurement_cell = 100 * conclusion["measurement_cell_bridge_top_rate"]
    loo_nucleus = 100 * conclusion["leave_one_out_nucleus_bridge_top_rate"]
    loo_cell = 100 * conclusion["leave_one_out_cell_bridge_top_rate"]
    published_counts = ", ".join(
        f"{_friendly_model(model)} {count}/{conclusion['published_bootstrap_tree_fits']}"
        for model, count in conclusion["published_bootstrap_tree_top_counts"].items()
    )
    dsep_table = pd.read_csv(PRIMARY_DSEP_PATH)
    bridge_p = {
        str(row.candidate_model): float(row.component_p)
        for row in dsep_table.loc[
            dsep_table["candidate_model"].isin(["nucleus_bridge", "cell_bridge"])
        ].itertuples(index=False)
    }
    evolutionary_counts = ", ".join(
        f"{model}: {_friendly_model(winner)}"
        for model, winner in conclusion["evolutionary_model_top_classes"].items()
    )
    notebook = nbf.v4.new_notebook()
    notebook.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        "analysis": {
            "id": "cell_nucleus_genome_phylogenetic_path_analysis_v1",
            "causal_direction_identified": False,
            "primary_species": primary_species,
        },
    }
    notebook.cells = [
        _markdown(
            f"""
            <style>
            .jp-Notebook {{ max-width: 1500px; margin: auto; }}
            table.dataframe {{ font-size: 12px; }}
            .verdict {{ border-left: 7px solid #2E7D6E; padding: 0.8em 1.1em; background: #F3F8F6; }}
            .warning {{ border-left: 7px solid #B44A46; padding: 0.8em 1.1em; background: #FBF4F3; }}
            </style>

            # Cell–nucleus–genome phylogenetic path analysis

            **Finalized {primary_species}-species *Desmognathus* panel · calibrated genome size in pg ·
            vetted top-50 cell/nucleus morphology · exhaustive three-node DAG audit**

            <div class="verdict">
            <b>Bottom line.</b> The traits form a connected positive system, but these
            data do not identify a unique causal pathway. {primary_top} ranks first
            (global P = {conclusion['primary_top_global_p']:.3f}, weight =
            {conclusion['primary_top_weight']:.3f}), while {second} {second_status}
            (ΔCICc = {conclusion['primary_second_delta_cicc']:.2f}, weight =
            {conclusion['primary_second_weight']:.3f}). The three proposed arrow
            directions are Markov-equivalent and therefore cannot be separated by
            cross-sectional phylogenetic covariance.
            </div>

            This notebook is the canonical review surface. It exposes the raw
            species data, every distinguishable three-node structure, primary path
            coefficients, global-fit tests, CICc ranking, measurement uncertainty,
            leave-one-species-out influence, evolutionary-model and tree sensitivity,
            and actual-tree simulation calibration.
            """,
            "path-analysis-01",
        ),
        _markdown("## 1. Exact analysis inputs and finalized estimands", "path-analysis-02"),
        _code(
            """
            from pathlib import Path
            import json
            import pandas as pd
            from IPython.display import Image, display, FileLink, Markdown

            pd.set_option("display.max_columns", 80)
            pd.set_option("display.max_rows", 300)
            pd.set_option("display.width", 180)

            def find_project_root(start=None):
                current = (start or Path.cwd()).resolve()
                for candidate in [current, *current.parents]:
                    if (candidate / "paths.yaml").exists():
                        return candidate
                raise FileNotFoundError("Could not locate the Desmognathus_TE project root")

            ROOT = find_project_root()
            REVIEW_DATA_DIR = ROOT / "results/data/research_review"
            FROZEN_REGISTRY_PATH = REVIEW_DATA_DIR / "frozen_input_registry.csv"
            frozen_registry = pd.read_csv(FROZEN_REGISTRY_PATH)
            required_registry_columns = {"source_path", "frozen_path", "source_sha256", "frozen_sha256"}
            assert required_registry_columns.issubset(frozen_registry.columns)
            assert frozen_registry.source_sha256.eq(frozen_registry.frozen_sha256).all()
            FROZEN_ARTIFACTS = dict(zip(frozen_registry.source_path, frozen_registry.frozen_path))

            def resolve_artifact(relative):
                relative = str(relative)
                source = ROOT / relative
                if source.exists():
                    return source
                frozen_relative = FROZEN_ARTIFACTS.get(relative)
                if frozen_relative is None:
                    raise FileNotFoundError(relative)
                frozen = ROOT / frozen_relative
                if not frozen.exists():
                    raise FileNotFoundError(frozen)
                return frozen

            INPUT_DIR = ROOT / "results/data/research_review/cell_nucleus_genome_path"
            REPORT_DIR = INPUT_DIR
            FIG_DIR = ROOT / "results/figures/research_review/cell_nucleus_genome_path"

            traits = pd.read_csv(INPUT_DIR / "primary_traits_all_finalized.csv")
            dags = pd.read_csv(INPUT_DIR / "all_three_trait_dags.csv")
            primary = pd.read_csv(INPUT_DIR / "primary_model_ranking.csv")
            dsep = pd.read_csv(INPUT_DIR / "primary_dsep_tests.csv")
            pairwise = pd.read_csv(INPUT_DIR / "pairwise_phylogenetic_regressions.csv")
            mechanisms = pd.read_csv(INPUT_DIR / "user_mechanism_equivalence_ranking.csv")
            mechanism_edges = pd.read_csv(INPUT_DIR / "user_mechanism_edges.csv")
            evolutionary = pd.read_csv(INPUT_DIR / "evolutionary_model_sensitivity_rankings.csv")
            loo = pd.read_csv(INPUT_DIR / "leave_one_species_out_rankings.csv.gz")
            measurement = pd.read_csv(INPUT_DIR / "measurement_bootstrap_rankings.csv.gz")
            trees = pd.read_csv(INPUT_DIR / "tree_sensitivity_rankings.csv.gz")
            simulation = pd.read_csv(INPUT_DIR / "simulation_summary.csv")
            failures = pd.read_csv(INPUT_DIR / "analysis_failures.csv")
            fit_warnings = pd.read_csv(INPUT_DIR / "primary_fit_warnings.csv")
            manifest = json.loads((INPUT_DIR / "analysis_manifest.json").read_text())
            conclusion = json.loads((INPUT_DIR / "conclusion_summary.json").read_text())

            assert len(traits) == manifest["primary_species"] and traits.species.is_unique
            assert len(dags) == 25 and dags.equivalence_class.nunique() == 11
            assert primary.model.nunique() == 10
            assert failures.empty

            display(pd.DataFrame({
                "quantity": [
                    "Primary species", "Labeled DAGs", "Markov-equivalence classes",
                    "Testable classes", "Measurement-bootstrap fits",
                    "Published bootstrap-tree fits", "Fit failures"
                ],
                "value": [
                    len(traits), len(dags), dags.equivalence_class.nunique(),
                    primary.model.nunique(), measurement.replicate_id.nunique(),
                    trees.loc[trees.analysis_type.eq("published_bootstrap_tree_sensitivity"), "replicate_id"].nunique(),
                    len(failures),
                ],
            }))
            """,
            "path-analysis-03",
        ),
        _markdown(
            """
            Genome size is the image-quality-matched nuclear-IOD estimate calibrated
            as `species IOD / D. fuscus IOD × 16.36 pg`. Thus every plotted and modeled
            genome value is in **picograms**, not a relative-IOD axis. A shared
            multiplicative calibration does not change standardized log-scale path
            coefficients or model ranking, although image-level measurement uncertainty
            does and is propagated below.

            Nucleus and cell sizes are the species medians of the manually vetted 50
            largest valid cell masks and their corresponding nucleus masks. These are
            paired **upper-tail morphology estimands**, not mean erythrocyte size.
            `D. aeneus`, `D. ochrophaeus`, and `D. wrighti` retain limited-overlap
            image-quality flags. All three are included in the 24-species path panel,
            with observed nucleus and cell counts reported for interpretation; no
            species is phylogenetically filled.
            """,
            "path-analysis-04",
        ),
        _code(
            """
            display(traits.style.format({
                "genome_size_pg": "{:.3f}",
                "nucleus_area_um2": "{:.3f}",
                "cell_area_um2": "{:.3f}",
            }))
            """,
            "path-analysis-05",
        ),
        _markdown("## 2. What the raw data say", "path-analysis-06"),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "01_raw_trait_relationships.png"), width=1450))
            """,
            "path-analysis-07",
        ),
        _markdown(
            """
            All three pairwise associations are positive. The path question is harder:
            which association remains after conditioning on the third trait, and is that
            conclusion stable after phylogenetic correction and measurement uncertainty?
            Simple pairwise significance cannot orient a causal arrow.
            """,
            "path-analysis-08",
        ),
        _code(
            """
            display(pairwise[[
                "relation_id", "outcome", "predictor", "standardized_coefficient",
                "standard_error", "ci_low", "ci_high", "p_value",
                "fitted_phylogenetic_parameter", "r_squared", "n_species"
            ]].style.format({
                "standardized_coefficient": "{:.3f}", "standard_error": "{:.3f}",
                "ci_low": "{:.3f}", "ci_high": "{:.3f}", "p_value": "{:.4f}",
                "fitted_phylogenetic_parameter": "{:.2g}", "r_squared": "{:.3f}"
            }))
            """,
            "path-analysis-09",
        ),
        _markdown("## 3. Method mirrored from the reference study", "path-analysis-10"),
        _markdown(
            """
            The workflow mirrors Yu et al. (2020): state candidate DAGs, derive their
            d-separation claims, test those claims by phylogenetic generalized least
            squares, combine component probabilities with Fisher's C, reject models with
            global P < 0.05, and rank the surviving structures with small-sample CICc.
            Models at ΔCICc ≤ 2 are treated as competitive.

            Primary sources: [Yu et al. reference paper](https://doi.org/10.1111/jzo.12755),
            [phylopath method paper](https://doi.org/10.7717/peerj.4718), and
            [original phylogenetic confirmatory path framework](https://doi.org/10.1111/j.1558-5646.2012.01790.x).

            With only three variables there are 25 labeled DAGs but 11
            Markov-equivalence classes. Ten have testable independence claims. Six fully
            connected DAGs collapse into one saturated class with no d-separation claim,
            so `phylopath` cannot score that class. Because this analysis exhaustively
            screens every distinguishable class rather than testing a small preregistered
            model set, it is an **exploratory structural screen**, not confirmatory proof
            of causality.

            The displayed gallery is intentionally narrower than that exhaustive fitted
            universe: it omits the biologically implausible genome-cell-only skeleton and
            the untestable saturated class, and it explicitly shows both genome ->
            nucleus -> cell and cell -> nucleus -> genome.
            """,
            "path-analysis-11",
        ),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "04_candidate_equivalence_class_gallery.png"), width=1300))
            class_summary = (
                dags.groupby(["equivalence_class", "basis_claim", "testable_by_dsep"], dropna=False)
                .agg(n_labeled_dags=("dag_id", "size"), min_edges=("n_edges", "min"))
                .reset_index()
                .sort_values(["min_edges", "equivalence_class"])
            )
            display(class_summary)
            """,
            "path-analysis-12",
        ),
        _markdown("## 4. Primary model comparison", "path-analysis-13"),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "02_primary_model_ranking.png"), width=1250))
            display(primary[[
                "rank", "model", "k", "q", "C", "p", "CICc", "delta_CICc", "w",
                "globally_supported", "admissible_competitive"
            ]].style.format({
                "C": "{:.3f}", "p": "{:.4f}", "CICc": "{:.3f}",
                "delta_CICc": "{:.3f}", "w": "{:.3f}"
            }))
            """,
            "path-analysis-14",
        ),
        _markdown(
            f"""
            The point estimate favors {primary_top}. {second} is ΔCICc =
            {conclusion['primary_second_delta_cicc']:.2f} behind and {second_status}.
            The nucleus–cell-only
            model is globally non-rejected but falls outside ΔCICc ≤ 2. The independence,
            collider, and genome-bridge structures are rejected or weakly supported.

            The decisive d-separation tests are symmetric in their warning:

            - **Nucleus bridge:** does genome still explain cell after nucleus is in the
              model? P = {bridge_p['nucleus_bridge']:.3f}.
            - **Cell bridge:** does genome still explain nucleus after cell is in the
              model? P = {bridge_p['cell_bridge']:.3f}.

            Non-significance means those structures are *not rejected*; it does not prove
            the corresponding independence relation. All six primary lambda component
            fits reached the lambda boundary (approximately zero), so the BM and OU
            sensitivity rows are part of the result rather than an optional appendix.
            """,
            "path-analysis-15",
        ),
        _code(
            """
            display(dsep[[
                "candidate_model", "independence_claim", "component_p",
                "fitted_phylogenetic_parameter"
            ]].style.format({"component_p": "{:.4f}", "fitted_phylogenetic_parameter": "{:.2g}"}))
            display(fit_warnings)
            """,
            "path-analysis-16",
        ),
        _markdown("## 5. The three proposed causal stories", "path-analysis-17"),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "03_markov_equivalent_user_mechanisms.png"), width=1450))
            display(mechanisms[["model", "C", "p", "CICc", "delta_CICc", "w"]])
            display(mechanism_edges)
            """,
            "path-analysis-18",
        ),
        _markdown(
            """
            <div class="warning">
            <b>Hard identifiability limit.</b> Genome → nucleus → cell, cell → nucleus →
            genome, and nucleus → both genome and cell have the same skeleton, no
            collider, and the same testable statement: genome and cell are independent
            conditional on nucleus. They therefore receive exactly the same Fisher C,
            global P, and CICc. No amount of refitting these same cross-sectional traits
            can choose among the three directions.
            </div>

            The displayed coefficients are valid standardized associations under each
            assumed arrow orientation. They are not evidence that one orientation caused
            the others.
            """,
            "path-analysis-19",
        ),
        _markdown("## 6. Ruthless stability tests", "path-analysis-20"),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "05_model_stability.png"), width=1300))
            display(Image(filename=str(FIG_DIR / "06_nucleus_vs_cell_bridge_contrast.png"), width=1300))
            """,
            "path-analysis-21",
        ),
        _markdown(
            f"""
            The structural label is fragile:

            - Across 250 paired measurement-bootstrap draws, nucleus bridge ranks first
              {measurement_nucleus:.1f}% of the time and cell bridge {measurement_cell:.1f}%.
            - Across the {primary_species} leave-one-species-out fits, the rates are
              {loo_nucleus:.1f}% and {loo_cell:.1f}%, respectively.
            - In the 18-species published-tree bootstrap panel: {published_counts}.
            - Evolutionary-model winners are: {evolutionary_counts}.

            The focal and published trees give identical rankings on the same 18 species,
            and all 200 published bootstrap trees agree because fitted lambda is near
            zero. The 18-species flip is therefore driven by dropping `D. brimleyorum`
            and `D. folkertsi`, not by alternative bootstrap topologies.

            This is not a technical footnote: the answer to “which trait is the bridge?”
            changes under reasonable versions of the same analysis.
            """,
            "path-analysis-22",
        ),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "07_structural_sensitivity_heatmap.png"), width=1350))
            point_tree_rows = trees.loc[
                trees.analysis_type.eq("tree_species_set_sensitivity")
                & trees["rank"].le(3),
                ["tree_id", "n_species", "rank", "model", "p", "delta_CICc", "w"]
            ]
            display(point_tree_rows)
            """,
            "path-analysis-23",
        ),
        _markdown("## 7. Small-sample simulation calibration", "path-analysis-24"),
        _code(
            """
            display(Image(filename=str(FIG_DIR / "08_simulation_calibration.png"), width=1300))
            display(simulation.style.format({
                "true_class_top_rate": "{:.1%}",
                "true_class_delta_le_2_rate": "{:.1%}",
                "true_class_global_fit_rate": "{:.1%}",
                "unique_true_class_rate": "{:.1%}",
                "any_bridge_top_rate": "{:.1%}",
            }))
            """,
            "path-analysis-25",
        ),
        _markdown(
            f"""
            The simulation uses the actual focal tree, n = {primary_species}, observed bridge-sized
            coefficients, IID and Brownian residual regimes, and the same ten-class
            lambda analysis. It answers whether this study design can reliably recover a
            known generating class—not whether any simulated class is biologically true.
            There are 100 simulations per generating class and residual regime.
            Recovery is imperfect, and retaining several competitive models is common.
            That validates a model-uncertainty conclusion rather than a forced winner.
            """,
            "path-analysis-26",
        ),
        _markdown("## 8. Biological interpretation and final decision", "path-analysis-27"),
        _markdown(
            f"""
            | Question | Evidence-based answer |
            |---|---|
            | Are genome, nucleus, and cell size associated? | **Yes.** All three raw and simple lambda-PGLS associations are positive in the finalized 24-species panel. |
            | Is a connected bridge structure better than independence? | **Yes.** Independence is globally rejected and has ΔCICc > 15. |
            | Which bridge is best supported? | **{primary_top}.** The point estimate ranks it first; {second} has ΔCICc = {conclusion['primary_second_delta_cicc']:.2f}. |
            | Can these data choose genome → nucleus → cell, the reverse chain, or nucleus → both? | **No.** They are Markov-equivalent. |
            | Do the effect sizes make biological sense? | **As positive allometric associations, yes.** Larger-genome species tend to have larger nuclei and cells, and nucleus–cell coupling is strongest. Causal direction is not established. |
            | What is defensible to report? | A positive three-trait association favoring {primary_top} in this panel, with causal arrow direction explicitly unidentified. |

            ### Measurement caveat that matters for causality

            The genome-size estimate is derived from nuclear integrated optical density;
            integrated density includes a nuclear-area component. Therefore the observed
            genome–nucleus link is partly coupled by the measurement construction and is
            not an independent validation that genome size caused nucleus size. The paired
            measurement bootstrap propagates sampling variation but cannot remove that
            structural coupling.

            ### Final verdict

            **The analysis resolves association and bridge structure more strongly than
            causal direction.** The strongest honest statement is that genome size,
            nucleus size, and cell size covary and that {primary_top} is favored in this
            finalized panel. None of the Markov-equivalent arrow directions can be selected
            from these cross-sectional data alone.
            """,
            "path-analysis-28",
        ),
        _markdown("## 9. Reproducibility, failures, and downloadable outputs", "path-analysis-29"),
        _code(
            """
            display(pd.DataFrame([manifest]).T.rename(columns={0: "value"}))
            display(failures)
            assert failures.empty

            output_files = sorted(
                [path for path in REPORT_DIR.iterdir() if path.is_file()],
                key=lambda path: path.name,
            )
            display(Markdown("### Report files"))
            for path in output_files:
                display(FileLink(path))
            """,
            "path-analysis-30",
        ),
        _markdown(
            """
            Statistical engine: `phylopath` 1.3.1 and `phylolm` 2.6.5 in the repository's
            Dusky environment. Primary transformations are log10 followed by within-fit
            z-scoring. The exact command, package versions, replicate counts, and failure
            gate are recorded in `analysis_manifest.json`; hashes of every PNG/PDF path
            figure are recorded in `cell_nucleus_genome_path_figure_manifest.json`, and
            the executed notebook hash is recorded in `research_review_manifest.json`.
            """,
            "path-analysis-31",
        ),
    ]
    return notebook


def write_notebook(conclusion: dict[str, Any]) -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(build_notebook(conclusion), NOTEBOOK_PATH)


def main() -> None:
    figure_manifest = build()
    conclusion = json.loads(CONCLUSION_PATH.read_text())
    write_notebook(conclusion)
    print(json.dumps(figure_manifest, indent=2))
    print(f"Notebook source: {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()

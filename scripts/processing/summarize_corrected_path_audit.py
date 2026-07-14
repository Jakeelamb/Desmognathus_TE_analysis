#!/usr/bin/env python3
"""Freeze corrected path-model evidence, figures, and publication gates.

The analysis is intentionally fail-closed. Relative nuclear IOD is never
renamed genome size, terminal:internal mapping is never renamed an ectopic-
recombination rate, and stable model selection is not treated as causation.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "results/data/corrected/path_analysis"
FIGURES = ROOT / "results/figures/corrected/path_analysis"
AUDIT = ROOT / "plans/publication-readiness-deep-audit"

RANKINGS = {
    "data": DATA / "corrected_path_data_sensitivity_rankings_analysis18_v1.csv",
    "trees": DATA / "corrected_path_tree_sensitivity_rankings_analysis18_v1.csv",
    "loo": DATA / "corrected_path_leave_one_out_rankings_analysis18_v1.csv",
}
EDGES = {
    "data": DATA / "corrected_path_data_sensitivity_best_edges_analysis18_v1.csv",
    "trees": DATA / "corrected_path_tree_sensitivity_best_edges_analysis18_v1.csv",
    "loo": DATA / "corrected_path_leave_one_out_best_edges_analysis18_v1.csv",
}
FAILURES = {
    "data": DATA / "corrected_path_data_sensitivity_failures_analysis18_v1.csv",
    "trees": DATA / "corrected_path_tree_sensitivity_failures_analysis18_v1.csv",
    "loo": DATA / "corrected_path_leave_one_out_failures_analysis18_v1.csv",
}
SIMULATION = DATA / "corrected_path_simulation_calibration_analysis18_v1.csv"
SIMULATION_EDGES = DATA / "corrected_path_simulation_edge_calibration_analysis18_v1.csv"

STABILITY_OUTPUT = DATA / "corrected_path_model_stability_analysis18_v1.csv"
ANCHOR_OUTPUT = DATA / "corrected_path_anchor_model_comparison_analysis18_v1.csv"
ANCHOR_EDGE_OUTPUT = DATA / "corrected_path_anchor_edges_analysis18_v1.csv"
TREE_EDGE_OUTPUT = DATA / "corrected_path_tree_edge_uncertainty_analysis18_v1.csv"
GATE_OUTPUT = DATA / "publication_release_gate_matrix_analysis18_v1.csv"
MANIFEST_OUTPUT = DATA / "corrected_path_release_audit_analysis18_v1.manifest.json"
REPORT_OUTPUT = AUDIT / "corrected_path_analysis_audit_analysis18_v1.md"

WEIGHT_FIGURE = FIGURES / "corrected_path_measurement_model_weights_analysis18_v1.png"
MODEL_FIGURE = FIGURES / "corrected_path_anchor_model_comparison_analysis18_v1.png"
DAG_FIGURE = FIGURES / "corrected_path_anchor_dag_analysis18_v1.png"
TREE_FIGURE = FIGURES / "corrected_path_tree_uncertainty_analysis18_v1.png"
LOO_FIGURE = FIGURES / "corrected_path_leave_one_out_influence_analysis18_v1.png"
SIMULATION_FIGURE = FIGURES / "corrected_path_simulation_calibration_analysis18_v1.png"
GATE_FIGURE = FIGURES / "publication_release_gate_matrix_analysis18_v1.png"

ANCHOR_MORPHOLOGY = "image_balanced_selected50"
ANCHOR_IOD = "image_qc_pass"

FAMILY_LABELS = {
    "te_iod": "TE composition → relative IOD",
    "iod_morphology": "Relative IOD → morphology",
    "integrated": "Integrated TE–IOD–morphology",
    "terminal_internal_iod": "Terminal:internal proxy → relative IOD",
}

MODEL_LABELS = {
    "proxy_null": "Null",
    "morphology_null": "Null",
    "integrated_null": "Null",
    "ltr_only": "LTR only",
    "evenness_only": "Evenness only",
    "additive_composition": "LTR + evenness",
    "mediated_evenness": "LTR → evenness → IOD",
    "iod_to_nucleus": "IOD → nucleus",
    "iod_to_cell": "IOD → cell",
    "iod_nucleus_cell_chain": "IOD → nucleus → cell",
    "iod_morphology_only": "IOD–morphology only",
    "te_evenness_chain": "LTR → evenness → IOD → nucleus → cell",
    "te_additive_chain": "LTR + evenness → IOD chain",
    "te_nucleus_bypass": "Evenness nucleus bypass",
    "te_cell_bypass": "Evenness cell bypass",
    "terminal_internal_only": "Terminal:internal only",
    "ltr_terminal_internal_chain": "LTR → terminal:internal → IOD",
    "terminal_internal_evenness_additive": "Terminal:internal + evenness",
    "te_terminal_internal_state": "TE state → terminal:internal → IOD",
}

NODE_LABELS = {
    "ltr": "LTR:LINE\nlog-ratio",
    "even": "TE Pielou\nevenness",
    "iod": "Relative nuclear\nIOD proxy",
    "ns": "Nucleus area",
    "cs": "Cell area",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def fmt(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    if value != 0 and abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def require_inputs() -> None:
    paths = list(RANKINGS.values()) + list(EDGES.values()) + list(FAILURES.values())
    paths += [SIMULATION, SIMULATION_EDGES]
    missing = [portable(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing corrected path audit inputs: " + ", ".join(missing))


def count_failure_rows(path: Path) -> int:
    text = path.read_text().strip()
    if not text:
        return 0
    return len(pd.read_csv(path))


def load_tables() -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    rankings = {phase: pd.read_csv(path) for phase, path in RANKINGS.items()}
    edges = {phase: pd.read_csv(path) for phase, path in EDGES.items()}
    for frame in rankings.values():
        frame["model_label"] = frame["model"].map(MODEL_LABELS).fillna(frame["model"])
    return rankings, edges


def build_stability(rankings: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for phase, frame in rankings.items():
        winners = frame[frame["rank"].eq(1)].copy()
        for family, group in winners.groupby("family", sort=True):
            counts = group["model"].value_counts()
            modal = str(counts.index[0])
            rows.append(
                {
                    "phase": phase,
                    "family": family,
                    "family_label": FAMILY_LABELS[family],
                    "n_fits": len(group),
                    "modal_top_model": modal,
                    "modal_top_model_label": MODEL_LABELS.get(modal, modal),
                    "modal_top_model_rate": float(counts.iloc[0] / len(group)),
                    "global_fit_pass_rate": float(group["global_fit_pass"].mean()),
                    "unique_top_gate_rate": float(
                        group["release_gate_status"].eq("supported_unique_top_model").mean()
                    ),
                    "minimum_top_weight": float(group["w"].min()),
                    "median_top_weight": float(group["w"].median()),
                    "maximum_top_weight": float(group["w"].max()),
                    "minimum_global_p": float(group["p"].min()),
                    "maximum_global_p": float(group["p"].max()),
                    "absolute_genome_size_used": False,
                    "publication_claim_allowed": False,
                }
            )
    return pd.DataFrame(rows)


def anchor_rows(rankings: pd.DataFrame) -> pd.DataFrame:
    out = rankings[
        rankings["tree_id"].eq("published_main")
        & rankings["morphology_estimator"].eq(ANCHOR_MORPHOLOGY)
        & rankings["iod_subset"].eq(ANCHOR_IOD)
    ].copy()
    if out.groupby("family")["fit_id"].nunique().ne(1).any():
        raise RuntimeError("Anchor selection did not resolve one fit per family")
    return out.sort_values(["family", "rank"])


def anchor_edges(edges: pd.DataFrame) -> pd.DataFrame:
    return edges[
        edges["tree_id"].eq("published_main")
        & edges["morphology_estimator"].eq(ANCHOR_MORPHOLOGY)
        & edges["iod_subset"].eq(ANCHOR_IOD)
    ].copy().sort_values(["family", "parent", "child"])


def summarize_tree_edges(edges: pd.DataFrame) -> pd.DataFrame:
    return (
        edges.groupby(["family", "best_model", "parent", "child"], as_index=False)
        .agg(
            n_trees=("tree_id", "nunique"),
            coefficient_median=("coefficient", "median"),
            coefficient_q025=("coefficient", lambda x: x.quantile(0.025)),
            coefficient_q975=("coefficient", lambda x: x.quantile(0.975)),
            minimum_coefficient=("coefficient", "min"),
            maximum_coefficient=("coefficient", "max"),
            sign_stability_rate=("coefficient", lambda x: max((x > 0).mean(), (x < 0).mean())),
        )
        .sort_values(["family", "parent", "child"])
    )


def release_gates() -> pd.DataFrame:
    rows = [
        (
            "LTR deletion / ectopic recombination",
            "BLOCKED_AS_RATE",
            "Terminal:internal mapping is a deletion-footprint proxy; two species lack assembly mapping and every proxy path family is globally rejected.",
            "Exploratory terminal:internal deletion-footprint proxy only.",
        ),
        (
            "TE diversity indices",
            "APPROVED_DESCRIPTIVE",
            "Definitions, denominators, unresolved mass, order/superfamily levels, and sensitivity tables are explicit.",
            "Descriptive Shannon, Gini-Simpson, Hill, richness, and evenness patterns.",
        ),
        (
            "TE compositional PCA",
            "APPROVED_DESCRIPTIVE",
            "Order-level CLR PCA includes zero handling, scores/loadings, leave-one-out stability, phylogenetic PCA, and 200-tree sensitivity.",
            "Descriptive ordination; no unrestricted PERMANOVA group claim.",
        ),
        (
            "Time-tree sensitivity",
            "APPROVED_FOR_SENSITIVITY",
            "Published main tree plus 200 published time trees are propagated; focal-tree source/calibration metadata remain incomplete.",
            "Tree-set sensitivity, with source limitation stated.",
        ),
        (
            "Cell/nucleus linkage",
            "APPROVED_PROVENANCE",
            "All 900 frozen final-panel pairs retain one-to-one cell/nucleus linkage.",
            "Linked upper-tail morphometry provenance.",
        ),
        (
            "Segmentation generalization",
            "BLOCKED",
            "Held-out annotations cover no final-panel species; production cell masks overcount objects and nucleus validation is training-source-overlapping.",
            "No final-panel segmentation-accuracy claim.",
        ),
        (
            "Morphology trait",
            "SENSITIVITY_ONLY",
            "Top-50 estimates an upper tail; ranks are robust to largest-100, but specimen support and manual-review depth are uneven.",
            "Upper-tail cell/nucleus area sensitivity, not typical erythrocyte size.",
        ),
        (
            "Absolute genome size",
            "BLOCKED",
            "No documented same-batch DNA-stoichiometric reference calibration; IOD equals area times mean OD and is image-dependent.",
            "Relative nuclear-IOD proxy only; never pg or C-value.",
        ),
        (
            "Phylogenetic path analysis",
            "EXPLORATORY_ONLY",
            "Model selection is stable, but actual-tree null simulations yield non-negligible unique false selections; candidate graphs are post-audit, IOD is not genome size, algebraic area coupling remains, and measurement error is not propagated.",
            "Exploratory phylogenetically corrected association-model sensitivity; no causal claim.",
        ),
    ]
    return pd.DataFrame(rows, columns=["component", "status", "evidence", "permitted_language"])


def plot_model_weights(data_rankings: pd.DataFrame) -> None:
    estimators = [
        "composite_selected50",
        "image_balanced_selected50",
        "literal_largest50_eligible",
        "literal_largest100_eligible",
        "eligible_candidate_pool",
        "manual_keep_pool",
    ]
    estimator_labels = {
        "composite_selected50": "Frozen composite 50",
        "image_balanced_selected50": "Equal-image frozen 50",
        "literal_largest50_eligible": "Literal largest 50",
        "literal_largest100_eligible": "Literal largest 100",
        "eligible_candidate_pool": "All eligible",
        "manual_keep_pool": "Manual keeps",
    }
    iod_order = ["all_selected", "image_qc_pass", "high_qc"]
    iod_labels = ["All selected nuclei", "Image QC pass", "High-quality nuclei"]
    winners = data_rankings[
        data_rankings["rank"].eq(1) & data_rankings["tree_id"].eq("published_main")
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8), constrained_layout=True)
    for axis, family in zip(axes, ["iod_morphology", "integrated"]):
        selected = winners[winners["family"].eq(family)]
        matrix = selected.pivot(index="morphology_estimator", columns="iod_subset", values="w")
        matrix = matrix.reindex(index=estimators, columns=iod_order)
        sns.heatmap(
            matrix,
            ax=axis,
            cmap="viridis",
            vmin=0,
            vmax=1,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar=axis is axes[-1],
            cbar_kws={"label": "Top-model weight"} if axis is axes[-1] else None,
        )
        axis.set_title(FAMILY_LABELS[family])
        axis.set_xlabel("Relative-IOD image subset")
        axis.set_ylabel("")
        axis.set_xticklabels(iod_labels, rotation=25, ha="right")
        axis.set_yticklabels([estimator_labels[x] for x in estimators], rotation=0)
    fig.suptitle("Top-model support across microscopy measurement specifications", fontsize=14)
    fig.savefig(WEIGHT_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_anchor_models(anchor: pd.DataFrame) -> None:
    families = ["te_iod", "iod_morphology", "integrated", "terminal_internal_iod"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for axis, family in zip(axes.flat, families):
        group = anchor[anchor["family"].eq(family)].sort_values("w")
        colors = np.where(group["global_fit_pass"], "#367C6B", "#B64E4E")
        axis.barh(group["model_label"], group["w"], color=colors, alpha=0.9)
        for y, (_, row) in enumerate(group.iterrows()):
            axis.text(
                min(row["w"] + 0.015, 0.92),
                y,
                f"ΔCICc={row['delta_CICc']:.1f}; p={fmt(row['p'])}",
                va="center",
                fontsize=8,
            )
        axis.set_xlim(0, 1.02)
        axis.set_xlabel("CICc model weight")
        axis.set_title(FAMILY_LABELS[family])
        axis.grid(axis="x", alpha=0.2)
    fig.suptitle("Anchor candidate-model comparison (published main tree)", fontsize=14)
    fig.savefig(MODEL_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_anchor_dag(anchor_edge_frame: pd.DataFrame) -> None:
    group = anchor_edge_frame[anchor_edge_frame["family"].eq("integrated")]
    expected = {("ltr", "even"), ("even", "iod"), ("iod", "ns"), ("ns", "cs")}
    found = set(zip(group["parent"], group["child"]))
    if found != expected:
        raise RuntimeError(f"Integrated anchor edges changed: {found}")
    positions = {node: (index, 0) for index, node in enumerate(["ltr", "even", "iod", "ns", "cs"])}
    fig, axis = plt.subplots(figsize=(14, 4.2), constrained_layout=True)
    axis.set_xlim(-0.5, 4.5)
    axis.set_ylim(-1.0, 1.0)
    axis.axis("off")
    for node, (x, y) in positions.items():
        color = "#F1B56B" if node == "iod" else "#E9EFF4"
        axis.text(
            x,
            y,
            NODE_LABELS[node],
            ha="center",
            va="center",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.5", "fc": color, "ec": "#334155", "lw": 1.4},
        )
    for _, row in group.iterrows():
        x1, y1 = positions[row["parent"]]
        x2, y2 = positions[row["child"]]
        color = "#8B4C8C" if row["coefficient"] < 0 else "#246B8E"
        axis.annotate(
            "",
            xy=(x2 - 0.34, y2),
            xytext=(x1 + 0.34, y1),
            arrowprops={
                "arrowstyle": "-|>",
                "lw": 1.5 + 3 * abs(row["coefficient"]),
                "color": color,
                "shrinkA": 4,
                "shrinkB": 4,
            },
        )
        axis.text(
            (x1 + x2) / 2,
            0.32,
            f"β={row['coefficient']:.2f}\n95% approx. [{row['approx_ci_low']:.2f}, {row['approx_ci_high']:.2f}]",
            ha="center",
            va="center",
            fontsize=9,
            color=color,
        )
    axis.text(
        2,
        -0.72,
        "Exploratory association graph. Relative IOD is not genome size, and IOD contains nuclear area by construction.",
        ha="center",
        fontsize=10,
        color="#8A3B3B",
    )
    axis.set_title("Anchor integrated path: standardized direct estimates", fontsize=14)
    fig.savefig(DAG_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_tree_uncertainty(tree_edges: pd.DataFrame) -> None:
    group = tree_edges[tree_edges["family"].eq("integrated")].copy()
    group["edge"] = group["parent"].map(NODE_LABELS).str.replace("\n", " ") + " → " + group[
        "child"
    ].map(NODE_LABELS).str.replace("\n", " ")
    order = list(dict.fromkeys(group["edge"]))
    fig, axis = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    sns.violinplot(data=group, y="edge", x="coefficient", order=order, inner=None, color="#8CB9C8", ax=axis)
    sns.boxplot(
        data=group,
        y="edge",
        x="coefficient",
        order=order,
        width=0.18,
        showfliers=False,
        boxprops={"facecolor": "white", "zorder": 3},
        ax=axis,
    )
    axis.axvline(0, color="#555", lw=1, ls="--")
    axis.set_xlabel("Standardized coefficient across published time trees")
    axis.set_ylabel("")
    axis.set_title("Integrated-path branch-length uncertainty (main + 200 time trees)")
    fig.savefig(TREE_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_loo(loo_edges: pd.DataFrame, anchor_edge_frame: pd.DataFrame) -> None:
    group = loo_edges[loo_edges["family"].eq("integrated")].copy()
    group["edge"] = group["parent"] + "→" + group["child"]
    anchor = anchor_edge_frame[anchor_edge_frame["family"].eq("integrated")].copy()
    anchor["edge"] = anchor["parent"] + "→" + anchor["child"]
    anchor_map = anchor.set_index("edge")["coefficient"]
    group["coefficient_change"] = group["coefficient"] - group["edge"].map(anchor_map)
    matrix = group.pivot(index="omitted_species", columns="edge", values="coefficient_change")
    matrix = matrix.reindex(sorted(matrix.index))
    limit = max(abs(matrix.min().min()), abs(matrix.max().max()), 0.05)
    fig, axis = plt.subplots(figsize=(8.5, 8), constrained_layout=True)
    sns.heatmap(
        matrix,
        cmap="vlag",
        center=0,
        vmin=-limit,
        vmax=limit,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        cbar_kws={"label": "Coefficient change from full panel"},
        ax=axis,
    )
    axis.set_xlabel("Integrated-path edge")
    axis.set_ylabel("Omitted species")
    axis.set_title("Leave-one-species-out influence on direct paths")
    fig.savefig(LOO_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_simulation(simulation: pd.DataFrame, edge_calibration: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    families = ["te_iod", "iod_morphology", "integrated"]
    for column, family in enumerate(families):
        axis = axes[0, column]
        group = simulation[simulation["family"].eq(family)].copy()
        null = group[group["scenario"].eq("independent_null")].iloc[0]
        signal = group[group["scenario"].eq("observed_chain")].sort_values("effect_scale")
        x = np.arange(4)
        values = np.r_[null["false_unique_nonnull_rate"], signal["expected_unique_supported_rate"]]
        lows = np.r_[null["primary_rate_ci_low"], signal["primary_rate_ci_low"]]
        highs = np.r_[null["primary_rate_ci_high"], signal["primary_rate_ci_high"]]
        axis.errorbar(
            x,
            values,
            yerr=np.vstack([values - lows, highs - values]),
            fmt="o-",
            color="#326A81",
            capsize=4,
            lw=2,
        )
        axis.set_xticks(x, ["Null\nfalse +", "0.5×\nrecovery", "1×\nrecovery", "1.5×\nrecovery"])
        axis.set_ylim(-0.03, 1.03)
        axis.set_ylabel("Rate")
        axis.set_title(FAMILY_LABELS[family])
        axis.grid(axis="y", alpha=0.25)

        edge_axis = axes[1, column]
        edge_group = edge_calibration[edge_calibration["family"].eq(family)].copy()
        edge_group["edge"] = edge_group["parent"] + "→" + edge_group["child"]
        for edge, subset in edge_group.groupby("edge"):
            edge_axis.plot(
                subset["effect_scale"],
                subset["interval_coverage_rate"],
                marker="o",
                label=edge,
            )
        edge_axis.axhline(0.95, color="#555", ls="--", lw=1, label="Nominal 0.95")
        edge_axis.set_xticks([0.5, 1.0, 1.5])
        edge_axis.set_ylim(-0.03, 1.03)
        edge_axis.set_xlabel("Observed coefficient scale")
        edge_axis.set_ylabel("Approximate 95% interval coverage")
        edge_axis.grid(axis="y", alpha=0.25)
        edge_axis.legend(fontsize=7, loc="lower right")
    fig.suptitle("Actual-tree path-model calibration at n=18", fontsize=14)
    fig.savefig(SIMULATION_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_gates(gates: pd.DataFrame) -> None:
    color_map = {
        "APPROVED_DESCRIPTIVE": "#3F8D6A",
        "APPROVED_FOR_SENSITIVITY": "#4D879C",
        "APPROVED_PROVENANCE": "#4D879C",
        "SENSITIVITY_ONLY": "#D39A3A",
        "EXPLORATORY_ONLY": "#D0703B",
        "BLOCKED_AS_RATE": "#B85252",
        "BLOCKED": "#B85252",
    }
    fig, axis = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    axis.set_xlim(0, 1)
    axis.set_ylim(-0.5, len(gates) - 0.5)
    axis.axis("off")
    for index, row in gates.reset_index(drop=True).iterrows():
        y = len(gates) - index - 1
        axis.add_patch(plt.Rectangle((0.01, y - 0.36), 0.98, 0.72, color="#F7F8FA", ec="#D9DEE5"))
        axis.add_patch(plt.Rectangle((0.01, y - 0.36), 0.025, 0.72, color=color_map[row["status"]]))
        axis.text(0.055, y + 0.12, row["component"], va="center", fontsize=10.5, weight="bold")
        axis.text(0.055, y - 0.16, row["permitted_language"], va="center", fontsize=8.5, color="#4B5563")
        axis.text(
            0.97,
            y,
            row["status"].replace("_", " "),
            va="center",
            ha="right",
            fontsize=9,
            color=color_map[row["status"]],
            weight="bold",
        )
    axis.set_title("Publication release gates: what the current evidence permits", fontsize=14, pad=14)
    fig.savefig(GATE_FIGURE, dpi=220, bbox_inches="tight")
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    display = frame[columns].copy()
    def render(value: object) -> str:
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return fmt(float(value), 3)
        return str(value).replace("|", "\\|").replace("\n", " ")

    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(render(row[column]) for column in columns) + " |"
        for _, row in display.iterrows()
    ]
    return "\n".join([header, separator, *rows])


def write_report(
    rankings: dict[str, pd.DataFrame],
    stability: pd.DataFrame,
    anchor: pd.DataFrame,
    anchor_edge_frame: pd.DataFrame,
    simulation: pd.DataFrame,
    edge_simulation: pd.DataFrame,
    gates: pd.DataFrame,
) -> None:
    data_winners = rankings["data"][rankings["data"]["rank"].eq(1)]
    n_fits = sum(frame["fit_id"].nunique() for frame in rankings.values())
    n_failures = sum(count_failure_rows(path) for path in FAILURES.values())
    anchor_top = anchor[anchor["rank"].eq(1)].copy()
    anchor_top["family"] = anchor_top["family"].map(FAMILY_LABELS)
    anchor_top["model"] = anchor_top["model"].map(MODEL_LABELS)
    anchor_top["global_p"] = anchor_top["p"].map(fmt)
    anchor_top["weight"] = anchor_top["w"].map(lambda x: fmt(x, 3))
    anchor_top["gate"] = anchor_top["release_gate_status"].str.replace("_", " ")

    integrated_edges = anchor_edge_frame[anchor_edge_frame["family"].eq("integrated")].copy()
    integrated_edges["edge"] = integrated_edges["parent"] + " → " + integrated_edges["child"]
    integrated_edges["estimate"] = integrated_edges.apply(
        lambda row: f"{row['coefficient']:.3f} [{row['approx_ci_low']:.3f}, {row['approx_ci_high']:.3f}]",
        axis=1,
    )

    calibration_display = simulation.copy()
    calibration_display["family"] = calibration_display["family"].map(FAMILY_LABELS)
    calibration_display["scenario"] = calibration_display.apply(
        lambda row: "independent null" if row["scenario"] == "independent_null" else f"observed chain ×{row['effect_scale']:g}",
        axis=1,
    )
    calibration_display["recovery_or_false_rate"] = calibration_display.apply(
        lambda row: row["false_unique_nonnull_rate"]
        if row["effect_scale"] == 0
        else row["expected_unique_supported_rate"],
        axis=1,
    )
    calibration_display["rate_95_ci"] = calibration_display.apply(
        lambda row: f"{row['recovery_or_false_rate']:.3f} [{row['primary_rate_ci_low']:.3f}, {row['primary_rate_ci_high']:.3f}]",
        axis=1,
    )
    null_calibration = simulation[simulation["scenario"].eq("independent_null")].set_index("family")
    observed_calibration = simulation[
        simulation["scenario"].eq("observed_chain") & simulation["effect_scale"].eq(1.0)
    ].set_index("family")

    report = f"""# Corrected phylogenetic path-analysis audit — final 18 species (v1)

## Verdict

The corrected path implementation is technically reproducible and substantially stronger than the historical analysis: **{n_fits} fits completed with {n_failures} failures** across 18 measurement specifications, the published main time tree, 200 published time-tree replicates, and leave-one-species-out analyses. Every candidate set now has an explicit null, every basis set is exported, rejected global models are fail-closed, and all outputs identify the image node as a **relative nuclear-IOD proxy**, never absolute genome size.

The path analysis is nevertheless **not approved for a publication claim about genome size or causal mechanisms**. The decisive reasons are upstream measurement validity, not software execution: relative IOD is uncalibrated and algebraically includes nuclear area; segmentation has no held-out final-panel species validation; the largest-50 trait is an upper-tail estimand; within-species measurement error is not propagated; and the candidate DAGs were formalized during this audit rather than genuinely preregistered before viewing the data. The terminal:internal family is additionally rejected by its global-fit test in every full-panel specification.

## Anchor result

Anchor specification: equal-image frozen selected-50 morphology, image-QC-passing relative IOD, and the published main Stewart–Wiens time tree.

{markdown_table(anchor_top, ['family', 'model', 'global_p', 'weight', 'gate'])}

The selected integrated graph has these standardized direct estimates (approximate Wald intervals from `phylopath`):

{markdown_table(integrated_edges, ['edge', 'estimate', 'interval_excludes_zero'])}

The `even → IOD` interval includes zero even though the mediated graph has overwhelming CICc weight. Model weight compares whole conditional-independence structures; it is not a substitute for uncertainty on a particular arrow. The positive `IOD → nucleus area` path is not an independent biological validation because the stored IOD is exactly nuclear pixel area multiplied by mean optical density.

## Robustness coverage

- Data specification: {len(data_winners)} family/specification/tree fits. The integrated TE-evenness chain ranked first in {int((data_winners.query("family == 'integrated'")['model'] == 'te_evenness_chain').sum())}/{len(data_winners.query("family == 'integrated'"))}; the IOD–nucleus–cell chain ranked first in {int((data_winners.query("family == 'iod_morphology'")['model'] == 'iod_nucleus_cell_chain').sum())}/{len(data_winners.query("family == 'iod_morphology'"))}.
- Tree uncertainty: each family was refit on the published main tree plus 200 published time trees. The same top model was selected on every tree for all four families; the terminal:internal winner still failed global fit on every tree.
- Species influence: the same TE-IOD, morphology, and integrated winners remained top after omitting each species. Terminal:internal models failed global fit in all 16 analyzable omissions.
- Measurement uncertainty remains incompletely propagated. Switching estimands and IOD QC subsets is a sensitivity analysis, not an error-aware hierarchical SEM.

{markdown_table(stability, ['phase', 'family_label', 'n_fits', 'modal_top_model_label', 'modal_top_model_rate', 'global_fit_pass_rate', 'minimum_top_weight', 'maximum_top_weight'])}

## Actual-tree simulation calibration

Independent Brownian residual traits and the observed chain were simulated on the published final-18 tree. Null scenarios use 200 replicates per family; signal scenarios use 100 replicates per family at 0.5×, 1×, and 1.5× the observed standardized coefficients. “Recovery” requires the generating model to rank first and pass the unique-support gate. “False +” uses the matching strict definition: a non-null model ranks first and receives the unique-support gate under independent null traits. Wilson 95% intervals quantify Monte Carlo/binomial uncertainty. Coefficient bias and approximate interval coverage are exported for every generating edge.

{markdown_table(calibration_display, ['family', 'scenario', 'n_successful', 'modal_top_model', 'rate_95_ci', 'mean_top_weight'])}

Simulation is a diagnostic under the stated Brownian data-generating process, not proof that the real DAG is causal. Poor recovery would block model discrimination; good recovery cannot rescue invalid trait meaning or unmodeled measurement error.

The strict unique false-selection rates were {null_calibration.loc['te_iod', 'false_unique_nonnull_rate']:.1%} for TE–IOD, {null_calibration.loc['iod_morphology', 'false_unique_nonnull_rate']:.1%} for IOD–morphology, and {null_calibration.loc['integrated', 'false_unique_nonnull_rate']:.1%} for the integrated family. A non-null model ranked first under the null even more often ({null_calibration.loc['te_iod', 'false_nonnull_top_rate']:.1%}, {null_calibration.loc['iod_morphology', 'false_nonnull_top_rate']:.1%}, and {null_calibration.loc['integrated', 'false_nonnull_top_rate']:.1%}, respectively), usually with model-set competition. At the observed coefficient scale, strict generating-model recovery was {observed_calibration.loc['te_iod', 'expected_unique_supported_rate']:.1%}, {observed_calibration.loc['iod_morphology', 'expected_unique_supported_rate']:.1%}, and only {observed_calibration.loc['integrated', 'expected_unique_supported_rate']:.1%}. Thus the integrated pattern is empirically stable in the observed dataset but only moderately identifiable at `n=18`; the false-selection surface is too large for confirmatory language.

Across generating edges and signal scales, approximate 95% interval coverage ranged from {edge_simulation['interval_coverage_rate'].min():.3f} to {edge_simulation['interval_coverage_rate'].max():.3f}; mean coefficient bias ranged from {edge_simulation['mean_bias'].min():.3f} to {edge_simulation['mean_bias'].max():.3f}. See the calibration CSV for each edge and effect scale.

## Publication release matrix

{markdown_table(gates, ['component', 'status', 'permitted_language'])}

## Comparison with leading methods

- `phylopath` requires a biologically justified model set, basis-set/global-fit evidence, CICc comparison, and explicit model uncertainty (van der Bijl 2018, [doi:10.7717/peerj.4718](https://doi.org/10.7717/peerj.4718)). Those computational reporting surfaces are now present.
- The complete model comparison is propagated over all 200 time trees supplied for comparative analyses by Stewart and Wiens (2025, [doi:10.1016/j.ympev.2024.108272](https://doi.org/10.1016/j.ympev.2024.108272)), matching current tree-uncertainty expectations.
- Actual-tree Monte Carlo addresses the small-panel model-selection problem highlighted by Boettiger, Coop, and Ralph (2012, [doi:10.1111/j.1558-5646.2011.01574.x](https://doi.org/10.1111/j.1558-5646.2011.01574.x)), but a hierarchical error-aware model remains the stronger next analysis.
- The image assay does not meet Feulgen densitometry standards in Hardie, Gregory, and Hebert (2002, [doi:10.1177/002215540205000601](https://doi.org/10.1177/002215540205000601)) or the same-batch reference practice used in salamanders by Mueller et al. (2008, [doi:10.1016/j.zool.2007.07.010](https://doi.org/10.1016/j.zool.2007.07.010).)
- Recent salamander TE-diversity studies report clearly defined Shannon/Gini-Simpson metrics and phylogenetically controlled trait tests; matching their quality means transparent definitions and uncertainty, not forcing the same direction of result (Decena-Segarra and Rovito 2024, [doi:10.1093/molbev/msae225](https://doi.org/10.1093/molbev/msae225)).

## Figures and exact data

- `{portable(MODEL_FIGURE)}` — all anchor candidate models, weights, global fit.
- `{portable(DAG_FIGURE)}` — standardized anchor direct paths and approximate intervals.
- `{portable(WEIGHT_FIGURE)}` — morphology/IOD specification sensitivity.
- `{portable(TREE_FIGURE)}` — direct paths across published time trees.
- `{portable(LOO_FIGURE)}` — species influence on the integrated path.
- `{portable(SIMULATION_FIGURE)}` — recovery, false selection, and interval coverage.
- `{portable(GATE_FIGURE)}` — permitted-language release matrix.
- `{portable(ANCHOR_OUTPUT)}`, `{portable(ANCHOR_EDGE_OUTPUT)}`, `{portable(STABILITY_OUTPUT)}`, `{portable(TREE_EDGE_OUTPUT)}` — review-sized exact tables.
- All raw model rankings, basis sets, fitted best edges, simulation replicates, and manifests remain in `{portable(DATA)}/`.

## Remaining work required for confirmatory promotion

1. Establish an absolute genome-size trait with a documented DNA-stoichiometric assay, same-batch reference standard, OD equation, slide/batch/specimen replication, and external validation; otherwise reframe the paper around relative IOD and morphology without genome-size claims.
2. Produce specimen/slide-held-out and final-panel-species segmentation truth with object-level split/merge/miss metrics, then rerun the morphometry from the validated production model.
3. Make specimen-balanced central tendency primary and keep largest-50 as a named upper-tail sensitivity; acquire biological replicates where only one image/specimen exists.
4. Propagate specimen-level uncertainty and TE uncertainty through an error-aware phylogenetic SEM; repeat tree and species influence analyses.
5. Treat the current DAG family as exploratory. Freeze a biologically justified confirmatory set before collecting/reprocessing the validation data.

No historical output was deleted or overwritten by this audit.
"""
    REPORT_OUTPUT.write_text(report)


def main() -> None:
    require_inputs()
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    rankings, edges = load_tables()
    stability = build_stability(rankings)
    anchor = anchor_rows(rankings["data"])
    anchor_edge_frame = anchor_edges(edges["data"])
    tree_edge_summary = summarize_tree_edges(edges["trees"])
    gates = release_gates()
    simulation = pd.read_csv(SIMULATION)
    edge_simulation = pd.read_csv(SIMULATION_EDGES)

    stability.to_csv(STABILITY_OUTPUT, index=False)
    anchor.to_csv(ANCHOR_OUTPUT, index=False)
    anchor_edge_frame.to_csv(ANCHOR_EDGE_OUTPUT, index=False)
    tree_edge_summary.to_csv(TREE_EDGE_OUTPUT, index=False)
    gates.to_csv(GATE_OUTPUT, index=False)

    plot_model_weights(rankings["data"])
    plot_anchor_models(anchor)
    plot_anchor_dag(anchor_edge_frame)
    plot_tree_uncertainty(edges["trees"])
    plot_loo(edges["loo"], anchor_edge_frame)
    plot_simulation(simulation, edge_simulation)
    plot_gates(gates)
    write_report(rankings, stability, anchor, anchor_edge_frame, simulation, edge_simulation, gates)

    inputs = list(RANKINGS.values()) + list(EDGES.values()) + [SIMULATION, SIMULATION_EDGES]
    outputs = [
        STABILITY_OUTPUT,
        ANCHOR_OUTPUT,
        ANCHOR_EDGE_OUTPUT,
        TREE_EDGE_OUTPUT,
        GATE_OUTPUT,
        REPORT_OUTPUT,
        WEIGHT_FIGURE,
        MODEL_FIGURE,
        DAG_FIGURE,
        TREE_FIGURE,
        LOO_FIGURE,
        SIMULATION_FIGURE,
        GATE_FIGURE,
    ]
    manifest = {
        "analysis_id": "corrected_path_release_audit_analysis18_v1",
        "inputs": [{"path": portable(path), "sha256": sha256(path)} for path in inputs],
        "outputs": [portable(path) for path in outputs],
        "n_path_fits": int(sum(frame["fit_id"].nunique() for frame in rankings.values())),
        "n_fit_failures": int(sum(count_failure_rows(path) for path in FAILURES.values())),
        "absolute_genome_size_used": False,
        "causal_claim_allowed": False,
        "historical_outputs_overwritten": False,
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote corrected path release audit with {manifest['n_path_fits']} fits.")


if __name__ == "__main__":
    main()

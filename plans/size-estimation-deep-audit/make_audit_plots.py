#!/usr/bin/env python3
"""Generate read-only audit figures for the size-estimation deep audit."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
FIG_DIR = SCRIPT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def species_key(value: object) -> str:
    text = str(value).strip()
    for prefix in ["D. ", "D.", "Desmognathus "]:
        if text.startswith(prefix):
            text = text[len(prefix) :]
            break
    return text.strip().lower()


def species_label(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        return text
    if text.startswith("Desmognathus "):
        return "D. " + text.split(" ", 1)[1]
    return "D. " + text


def read_inputs() -> pd.DataFrame:
    derived = PROJECT_ROOT / "path_analysis" / "data" / "derived"
    external = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"
    cellprofiler_dir = (
        Path.home()
        / "Projects"
        / "cellprofiler_test"
        / "output"
        / "runs"
        / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
        / "verified_species_dataset_top50_latest"
    )

    master = pd.read_csv(derived / "path_input_master.csv", low_memory=False)
    sensitivity = pd.read_csv(external / "cellprofiler_genome_sensitivity.csv", low_memory=False)
    image_qc = pd.read_csv(external / "cellprofiler_image_iod_quality_summary.csv", low_memory=False)
    verified_species = pd.read_csv(cellprofiler_dir / "species_estimates_verified.csv", low_memory=False)
    selected_pairs = pd.read_csv(
        cellprofiler_dir / "selected_high_quality_linked_pairs_with_iod_qc.csv.gz",
        low_memory=False,
    )

    master["species_key"] = master["species"].map(species_key)
    sensitivity["species_key"] = sensitivity["species"].map(species_key)
    image_qc["species_key"] = image_qc["species"].map(species_key)
    verified_species["species_key"] = verified_species["species"].map(species_key)
    selected_pairs["species_key"] = selected_pairs["species"].map(species_key)

    shift = (
        sensitivity[sensitivity["subset_name"].ne("all_selected")]
        .pivot(index="species_key", columns="subset_name", values="pct_shift_vs_all_selected")
        .rename(
            columns={
                "image_qc_pass": "image_qc_shift_pct",
                "high_iod_qc": "high_iod_qc_shift_pct",
            }
        )
        .reset_index()
    )
    qc_counts = (
        image_qc.pivot_table(
            index="species_key",
            columns="image_iod_qc_status",
            values="filename",
            aggfunc="count",
            fill_value=0,
        )
        .add_prefix("n_image_qc_")
        .reset_index()
    )
    pair_stats = (
        selected_pairs.groupby("species_key", sort=True)
        .agg(
            pair_rows=("species_key", "size"),
            pair_images=("filename", "nunique"),
            median_nuc_iod=("nuc_iod", "median"),
            median_nuc_mean_od=("nuc_mean_od", "median"),
            median_iod_quality=("iod_final_quality_score", "median"),
            median_pair_quality=("quality_score", "median"),
            qc_pass_rows=("image_iod_qc_pass", "sum"),
        )
        .reset_index()
    )

    keep_cols = [
        "species_key",
        "species",
        "has_tree_tip",
        "has_te",
        "has_genome",
        "has_morphology",
        "genome_size_pg",
        "morph_cell_area_um2",
        "morph_nucleus_area_um2",
        "morph_nc_ratio",
        "genome_result_status",
        "genome_support_tier",
        "genome_support_warnings",
    ]
    data = master[[col for col in keep_cols if col in master.columns]].copy()
    data = data[data["has_genome"].fillna(False) | data["has_morphology"].fillna(False)].copy()
    data = data.merge(shift, on="species_key", how="left")
    data = data.merge(qc_counts, on="species_key", how="left")
    data = data.merge(pair_stats, on="species_key", how="left")
    data = data.merge(
        verified_species[
            [
                "species_key",
                "linked_n_selected_manual_total",
                "linked_n_selected_auto",
                "estimated_genome_pg_source",
                "genome_primary_missing",
                "genome_primary_n_selected_pairs",
                "genome_primary_n_selected_images",
                "genome_primary_support_label",
            ]
        ],
        on="species_key",
        how="left",
    )
    data["display_species"] = data["species"].map(species_label)
    data.to_csv(SCRIPT_DIR / "current_size_qc_audit_table.csv", index=False)
    return data


def terminal_order_from_tree(data: pd.DataFrame) -> list[str]:
    try:
        from Bio import Phylo
    except Exception:
        return sorted(data["species_key"].tolist())
    tree_path = PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"
    tree = Phylo.read(tree_path, "newick")
    tree.ladderize()
    order = [species_key(tip.name) for tip in tree.get_terminals()]
    present = set(data["species_key"])
    return [item for item in order if item in present]


def draw_tree_axis(ax, ordered_keys: list[str]) -> None:
    try:
        from Bio import Phylo
    except Exception:
        ax.set_axis_off()
        return

    tree_path = PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"
    tree = Phylo.read(tree_path, "newick")
    tree.ladderize()
    wanted = set(ordered_keys)
    for tip in list(tree.get_terminals()):
        if species_key(tip.name) not in wanted:
            tree.prune(tip)

    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    y_by_key = {key: idx for idx, key in enumerate(ordered_keys)}

    def clade_y(clade):
        if clade.is_terminal():
            return y_by_key[species_key(clade.name)]
        ys = [clade_y(child) for child in clade.clades]
        return float(np.mean(ys))

    def draw_clade(clade):
        x = depths[clade]
        y = clade_y(clade)
        if clade.clades:
            child_ys = [clade_y(child) for child in clade.clades]
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="#333333", linewidth=0.8)
            for child in clade.clades:
                cx = depths[child]
                cy = clade_y(child)
                ax.plot([x, cx], [cy, cy], color="#333333", linewidth=0.8)
                draw_clade(child)
        else:
            ax.plot([x], [y], marker="o", markersize=2.5, color="#333333")

    draw_clade(tree.root)
    ax.set_ylim(-0.5, len(ordered_keys) - 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("branch length")
    ax.set_title("Pruned tree")
    ax.set_yticks([])
    ax.grid(axis="x", alpha=0.18)


def normalized(values: pd.Series) -> pd.Series:
    vals = pd.to_numeric(values, errors="coerce")
    finite = vals[np.isfinite(vals)]
    if finite.empty or finite.max() == finite.min():
        return pd.Series(0.5, index=values.index)
    return (vals - finite.min()) / (finite.max() - finite.min())


def plot_phylogeny_traits(data: pd.DataFrame) -> None:
    order = terminal_order_from_tree(data)
    ordered = data.set_index("species_key").loc[order].reset_index()
    trait_specs = [
        ("genome_size_pg", "Genome pg", "viridis"),
        ("morph_cell_area_um2", "Cell area", "Blues"),
        ("morph_nucleus_area_um2", "Nucleus area", "Purples"),
        ("morph_nc_ratio", "N:C", "YlGn"),
        ("image_qc_shift_pct", "QC shift %", "RdBu_r"),
        ("linked_n_selected_auto", "Auto rows", "Oranges"),
    ]

    fig = plt.figure(figsize=(15, max(7, 0.32 * len(ordered))))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.4, 3.0], wspace=0.24)
    tree_ax = fig.add_subplot(grid[0, 0])
    heat_ax = fig.add_subplot(grid[0, 1])
    draw_tree_axis(tree_ax, order)

    mat = []
    for col, _label, _cmap in trait_specs:
        if col == "image_qc_shift_pct":
            vals = pd.to_numeric(ordered[col], errors="coerce")
            finite = vals[np.isfinite(vals)]
            if finite.empty:
                normed = pd.Series(0.5, index=ordered.index)
            else:
                max_abs = max(abs(float(finite.min())), abs(float(finite.max())), 1.0)
                normed = (vals + max_abs) / (2 * max_abs)
        else:
            normed = normalized(ordered[col])
        mat.append(normed.to_numpy(float))
    matrix = np.vstack(mat).T
    heat_ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    heat_ax.set_yticks(np.arange(len(ordered)))
    heat_ax.set_yticklabels(ordered["display_species"], fontsize=8)
    heat_ax.set_xticks(np.arange(len(trait_specs)))
    heat_ax.set_xticklabels([label for _col, label, _cmap in trait_specs], rotation=35, ha="right")
    heat_ax.set_title("Size and QC traits ordered by phylogeny")
    for i in range(len(ordered)):
        for j, (col, _label, _cmap) in enumerate(trait_specs):
            value = ordered.loc[i, col]
            if pd.notna(value):
                precision = 2 if col == "morph_nc_ratio" else 1
                heat_ax.text(
                    j,
                    i,
                    f"{float(value):.{precision}f}",
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white",
                )
    fig.suptitle("Current Cell/Nucleus/Genome Evidence Against Phylogeny", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "phylogeny_size_qc_audit.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_variable_relationships(data: pd.DataFrame) -> None:
    vars_ = [
        ("genome_size_pg", "Genome pg"),
        ("morph_cell_area_um2", "Cell area"),
        ("morph_nucleus_area_um2", "Nucleus area"),
        ("morph_nc_ratio", "N:C"),
        ("image_qc_shift_pct", "QC shift %"),
        ("linked_n_selected_auto", "Auto rows"),
    ]
    clean = data[["display_species"] + [col for col, _label in vars_]].copy()
    for col, _label in vars_:
        clean[col] = pd.to_numeric(clean[col], errors="coerce")
    n = len(vars_)
    fig, axes = plt.subplots(n, n, figsize=(16, 16))
    for i, (y_col, y_label) in enumerate(vars_):
        for j, (x_col, x_label) in enumerate(vars_):
            ax = axes[i, j]
            if i == j:
                vals = clean[x_col].dropna()
                ax.hist(vals, bins=min(8, max(3, len(vals) // 2)), color="#6C757D", alpha=0.8)
                ax.set_ylabel("")
            else:
                sub = clean[[x_col, y_col, "display_species"]].dropna()
                ax.scatter(sub[x_col], sub[y_col], s=26, color="#315C7A", alpha=0.72)
                if len(sub) >= 5:
                    rho = sub[[x_col, y_col]].corr(method="spearman").iloc[0, 1]
                    ax.text(0.05, 0.9, f"rho={rho:.2f}", transform=ax.transAxes, fontsize=8)
            if i == n - 1:
                ax.set_xlabel(x_label, fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel(y_label, fontsize=8)
            else:
                ax.set_yticklabels([])
            ax.grid(alpha=0.15)
    fig.suptitle("Current Size, Genome, QC, and Selection Relationships", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "size_qc_variable_relationships.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_qc_sensitivity(data: pd.DataFrame) -> None:
    sub = data.dropna(subset=["image_qc_shift_pct"]).copy()
    sub = sub.sort_values("image_qc_shift_pct")
    fig, ax = plt.subplots(figsize=(11, max(5, 0.28 * len(sub))))
    colors = np.where(sub["estimated_genome_pg_source"].eq("all_selected_fallback"), "#B94747", "#315C7A")
    ax.barh(sub["display_species"], sub["image_qc_shift_pct"], color=colors, alpha=0.86)
    ax.axvline(0, color="#222222", linewidth=1)
    ax.set_xlabel("image-QC-pass genome estimate shift vs all-selected (%)")
    ax.set_title("Direction and magnitude of genome-estimate QC sensitivity")
    ax.grid(axis="x", alpha=0.22)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "genome_qc_sensitivity_ranked.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = read_inputs()
    plot_phylogeny_traits(data)
    plot_variable_relationships(data)
    plot_qc_sensitivity(data)
    print(f"Wrote audit table and figures under {SCRIPT_DIR}")


if __name__ == "__main__":
    main()

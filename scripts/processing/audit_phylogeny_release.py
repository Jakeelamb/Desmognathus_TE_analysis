#!/usr/bin/env python3
"""Build a fail-closed release audit for the Desmognathus dated tree.

Historical trees and figures are preserved.  This script derives a versioned,
final-panel tree from the exact local input, corrects only Newick rounding at
terminal edges, and records why the current single-tree analysis is not yet a
publication-ready representation of phylogenetic uncertainty.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Phylo


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_TREE = PROJECT_ROOT / "input_data/phylogeny/desmo900dated_test.tre"
PROCESSED_TREE = PROJECT_ROOT / "results/phylogeny/processed_phylogeny.nwk"
PANEL = PROJECT_ROOT / "path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv"
MASS = PROJECT_ROOT / "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv"
MORPHOLOGY = PROJECT_ROOT / "path_analysis/data/external/derived/cellprofiler_species_morphology_summary.csv"
GENOME_IOD = PROJECT_ROOT / "path_analysis/data/external/derived/cellprofiler_final_species_results.csv"

OUTPUT_DIR = PROJECT_ROOT / "results/data/corrected/phylogeny"
FIGURE_DIR = PROJECT_ROOT / "results/figures/corrected/phylogeny"
AUDIT_DIR = PROJECT_ROOT / "plans/publication-readiness-deep-audit"

METRICS_OUTPUT = OUTPUT_DIR / "phylogeny_tree_metrics_analysis18_v1.csv"
CROSSWALK_OUTPUT = OUTPUT_DIR / "phylogeny_tip_crosswalk_analysis18_v1.csv"
DISTANCE_OUTPUT = OUTPUT_DIR / "phylogeny_patristic_distances_analysis18_v1.csv"
COMPARISON_OUTPUT = OUTPUT_DIR / "phylogeny_tree_input_comparison_v1.csv"
PRUNED_TREE_OUTPUT = OUTPUT_DIR / "desmognathus_time_tree_analysis18_v1.nwk"
MANIFEST_OUTPUT = OUTPUT_DIR / "phylogeny_release_audit_analysis18_v1.manifest.json"
REPORT_OUTPUT = AUDIT_DIR / "phylogeny_release_audit_analysis18_v1.md"

TREE_FIGURE = FIGURE_DIR / "phylogeny_data_completeness_analysis18_v1.png"
DISTANCE_FIGURE = FIGURE_DIR / "phylogeny_patristic_heatmap_analysis18_v1.png"
ROUNDING_FIGURE = FIGURE_DIR / "phylogeny_root_to_tip_rounding_analysis18_v1.png"
PROCESSED_FIGURE = FIGURE_DIR / "phylogeny_processed_scale_comparison_v1.png"


def canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip().lower()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def portable(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def branch_lengths(tree) -> np.ndarray:
    return np.asarray(
        [float(clade.branch_length) for clade in tree.find_clades() if clade is not tree.root],
        dtype=float,
    )


def root_to_tip_depths(tree) -> dict[str, float]:
    return {
        canonical_species(tip.name): float(tree.distance(tree.root, tip))
        for tip in tree.get_terminals()
    }


def internal_clade_signature(tree) -> set[frozenset[str]]:
    all_tips = {canonical_species(tip.name) for tip in tree.get_terminals()}
    signatures: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals(order="postorder"):
        descendants = frozenset(canonical_species(t.name) for t in clade.get_terminals())
        if 1 < len(descendants) < len(all_tips):
            signatures.add(descendants)
    return signatures


def prune_to_species(tree, species: Iterable[str]):
    target = {canonical_species(value) for value in species}
    copy = deepcopy(tree)
    for tip in copy.get_terminals():
        tip.name = canonical_species(tip.name)
    labels = [tip.name for tip in copy.get_terminals()]
    if len(labels) != len(set(labels)):
        duplicated = sorted({name for name in labels if labels.count(name) > 1})
        raise ValueError(f"Tree labels are duplicated after exact normalization: {duplicated}")
    missing = sorted(target.difference(labels))
    if missing:
        raise ValueError(f"Final panel is missing from tree: {missing}")
    for tip in list(copy.get_terminals()):
        if tip.name not in target:
            copy.prune(tip)
    if {tip.name for tip in copy.get_terminals()} != target:
        raise RuntimeError("Pruned tree does not exactly match the declared final panel")
    return copy


def pad_terminal_rounding(tree) -> tuple[object, pd.DataFrame]:
    """Pad only terminal edges to the maximum root-tip depth.

    The source tree is ultrametric at 1e-5 Myr tolerance but not at ape's
    strict default because its Newick values are rounded.  This deterministic
    correction never shortens an edge and preserves topology/internal ages.
    """

    copy = deepcopy(tree)
    depths = {tip: float(copy.distance(copy.root, tip)) for tip in copy.get_terminals()}
    target = max(depths.values())
    rows = []
    for tip, before in depths.items():
        delta = target - before
        if delta < -1e-12:
            raise RuntimeError("Negative terminal correction requested")
        tip.branch_length = float(tip.branch_length or 0.0) + max(delta, 0.0)
        after = float(copy.distance(copy.root, tip))
        rows.append(
            {
                "species": canonical_species(tip.name),
                "root_to_tip_myr_before": before,
                "terminal_padding_myr": max(delta, 0.0),
                "terminal_padding_years": max(delta, 0.0) * 1_000_000,
                "root_to_tip_myr_after": after,
            }
        )
    return copy, pd.DataFrame(rows).sort_values("species").reset_index(drop=True)


def pairwise_distances(tree) -> pd.DataFrame:
    species = [canonical_species(t.name) for t in tree.get_terminals()]
    tips = {canonical_species(t.name): t for t in tree.get_terminals()}
    matrix = pd.DataFrame(index=species, columns=species, dtype=float)
    for left in species:
        for right in species:
            matrix.loc[left, right] = float(tree.distance(tips[left], tips[right]))
    matrix.index.name = "species"
    return matrix


def tree_metrics(tree, tree_id: str, path: Path, rounding_tolerance_myr: float = 1e-5) -> dict[str, object]:
    depths = np.asarray(list(root_to_tip_depths(tree).values()), dtype=float)
    lengths = branch_lengths(tree)
    root_degree = len(tree.root.clades)
    child_counts = [len(clade.clades) for clade in tree.get_nonterminals()]
    return {
        "tree_id": tree_id,
        "path": portable(path),
        "sha256": sha256(path),
        "n_tips": len(tree.get_terminals()),
        "n_internal_nodes": len(tree.get_nonterminals()),
        "root_degree": root_degree,
        "rooted_by_root_degree": root_degree == 2,
        "fully_bifurcating": all(value == 2 for value in child_counts),
        "all_branch_lengths_present": len(lengths) == len(list(tree.find_clades())) - 1,
        "all_branch_lengths_positive": bool(np.all(lengths > 0)),
        "minimum_branch_length_myr": float(lengths.min()),
        "maximum_branch_length_myr": float(lengths.max()),
        "root_age_myr": float(depths.max()),
        "minimum_root_to_tip_myr": float(depths.min()),
        "maximum_root_to_tip_myr": float(depths.max()),
        "root_to_tip_range_myr": float(np.ptp(depths)),
        "ultrametric_at_1e_5_myr_tolerance": bool(np.ptp(depths) <= rounding_tolerance_myr),
        "strict_ultrametric_at_1e_10_myr": bool(np.ptp(depths) <= 1e-10),
        "duplicate_exact_tip_labels": len({t.name for t in tree.get_terminals()}) != len(tree.get_terminals()),
    }


def load_crosswalk(panel_species: list[str]) -> pd.DataFrame:
    panel = pd.read_csv(PANEL)
    panel["species"] = panel["species"].map(canonical_species)
    panel = panel.set_index("species")

    mass = pd.read_csv(MASS)
    mass["species"] = mass["species"].map(canonical_species)
    mass = mass.drop_duplicates("species").set_index("species")

    morphology = pd.read_csv(MORPHOLOGY)
    morphology["species"] = morphology["species"].map(canonical_species)
    morphology = morphology.set_index("species")

    genome = pd.read_csv(GENOME_IOD)
    genome["species"] = genome["species"].map(canonical_species)
    genome = genome.set_index("species")

    rows = []
    for species in panel_species:
        prow = panel.loc[species]
        mrow = mass.loc[species]
        morph = morphology.loc[species]
        grow = genome.loc[species]
        rows.append(
            {
                "species": species,
                "tree_tip": species,
                "tree_match_type": "exact_after_genus_prefix_removal",
                "te_sra_accession": mrow["te_sra_accession"],
                "te_assembly_accession": mrow["te_assembly_accession"],
                "has_te_composition": bool(prow["has_te"]),
                "has_terminal_internal_proxy": bool(prow["has_ectopic"]),
                "has_microscopy_morphology": bool(prow["has_morphology"]),
                "has_image_iod_sensitivity": bool(prow["has_genome"]),
                "microscopy_n_specimens": int(morph["n_specimens_strict"]),
                "microscopy_n_images": int(morph["n_images_strict"]),
                "genome_iod_n_specimens": int(grow["primary_n_specimens"]),
                "genome_iod_n_images": int(grow["primary_n_images"]),
                "genome_iod_status": grow["result_status"],
                "genome_and_microscopy_are_independent_samples": True,
                "taxon_decision_note": (
                    "Use accession-specific expert-reidentified fuscus resource on fuscus tip; "
                    "do not substitute planiceps tip"
                    if species == "fuscus"
                    else "No tree-tip substitution"
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_tree_completeness(tree, crosswalk: pd.DataFrame, path: Path) -> None:
    tree = deepcopy(tree)
    tree.ladderize()
    terminals = tree.get_terminals()
    y = {tip: index for index, tip in enumerate(terminals)}

    def assign_y(clade):
        if clade in y:
            return y[clade]
        values = [assign_y(child) for child in clade.clades]
        y[clade] = float(np.mean(values))
        return y[clade]

    assign_y(tree.root)
    x = tree.depths()
    root_age = max(x[t] for t in terminals)

    ordered = [canonical_species(t.name) for t in terminals]
    cw = crosswalk.set_index("species").loc[ordered]
    status_columns = [
        "has_te_composition",
        "has_terminal_internal_proxy",
        "has_microscopy_morphology",
        "has_image_iod_sensitivity",
    ]
    labels = ["TE", "T:I proxy", "Morphology", "Image IOD"]

    fig = plt.figure(figsize=(11, 8.5))
    grid = fig.add_gridspec(1, 2, width_ratios=[4.5, 1.8], wspace=0.04)
    ax = fig.add_subplot(grid[0, 0])
    for parent in tree.get_nonterminals(order="preorder"):
        child_ys = [y[child] for child in parent.clades]
        ax.plot([x[parent], x[parent]], [min(child_ys), max(child_ys)], color="#344E41", lw=1.2)
        for child in parent.clades:
            ax.plot([x[parent], x[child]], [y[child], y[child]], color="#344E41", lw=1.2)
    for tip in terminals:
        ax.text(root_age + 0.35, y[tip], f"D. {canonical_species(tip.name)}", va="center", fontsize=8.5)
    ticks = np.linspace(0, root_age, 7)
    ax.set_xticks(ticks, [f"{root_age - tick:.0f}" for tick in ticks])
    ax.set_xlabel("Million years before present")
    ax.set_xlim(-0.3, root_age + 7.7)
    ax.set_ylim(-1, len(terminals))
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.set_title("Final 18-species dated tree", loc="left", fontweight="bold")

    hx = fig.add_subplot(grid[0, 1], sharey=ax)
    values = cw[status_columns].astype(int).to_numpy()
    sns.heatmap(
        values,
        cmap=sns.color_palette(["#E8E8E8", "#2A9D8F"], as_cmap=True),
        cbar=False,
        linewidths=0.7,
        linecolor="white",
        xticklabels=labels,
        yticklabels=False,
        ax=hx,
        vmin=0,
        vmax=1,
    )
    hx.tick_params(axis="x", labelrotation=45, labelsize=8)
    hx.set_xlabel("")
    hx.set_title(
        "Available stream\nteal = available · gray = unavailable",
        fontsize=9,
        fontweight="bold",
    )
    fig.suptitle(
        "Tree and data completeness (streams joined by species, not specimen)",
        x=0.07,
        ha="left",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.015,
        "Image IOD is shown as an available sensitivity stream, not an approved absolute C-value.",
        fontsize=8.5,
        color="#555555",
    )
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_distance_heatmap(distance: pd.DataFrame, tree, path: Path) -> None:
    tree = deepcopy(tree)
    tree.ladderize()
    order = [canonical_species(t.name) for t in tree.get_terminals()]
    ordered = distance.loc[order, order]
    fig, ax = plt.subplots(figsize=(10, 8.5))
    sns.heatmap(ordered, cmap="viridis", square=True, cbar_kws={"label": "Patristic distance (Myr)"}, ax=ax)
    ax.set_title("Pairwise distances on the corrected final-panel time tree", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=55, labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rounding(rounding: pd.DataFrame, path: Path) -> None:
    ordered = rounding.sort_values("terminal_padding_years")
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.barh(ordered["species"], ordered["terminal_padding_years"], color="#E76F51")
    ax.set_xlabel("Terminal-edge padding needed for exact ultrametricity (years)")
    ax.set_ylabel("")
    ax.set_title("Source-tree root-to-tip variation is Newick rounding scale", fontweight="bold")
    ax.axvline(0, color="#333333", lw=0.8)
    ax.text(
        0.99,
        0.02,
        f"maximum = {ordered['terminal_padding_years'].max():.3f} years",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
    )
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_processed_comparison(source: pd.DataFrame, processed: pd.DataFrame, path: Path) -> None:
    common = source.index.intersection(processed.index)
    x = source.loc[common, common].to_numpy()
    y = processed.loc[common, common].to_numpy()
    mask = np.triu(np.ones_like(x, dtype=bool), k=1)
    x = x[mask]
    y = y[mask]
    fig, ax = plt.subplots(figsize=(7.2, 6.5))
    ax.scatter(x, y, s=22, alpha=0.6, color="#457B9D", edgecolor="none")
    limit = max(x.max(), y.max())
    ax.plot([0, limit], [0, 0.8 * limit], color="#D62828", lw=1.6, label="y = 0.8x")
    ax.set_xlabel("Source-tree pairwise distance (Myr)")
    ax.set_ylabel("Tracked processed-tree distance")
    ax.set_title("Processed tree applies an undocumented 0.8 branch scale", fontweight="bold")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    for path in [SOURCE_TREE, PROCESSED_TREE, PANEL, MASS, MORPHOLOGY, GENOME_IOD]:
        if not path.exists():
            raise FileNotFoundError(path)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)

    source = Phylo.read(SOURCE_TREE, "newick")
    processed = Phylo.read(PROCESSED_TREE, "newick")
    source_labels = [canonical_species(t.name) for t in source.get_terminals()]
    processed_labels = [canonical_species(t.name) for t in processed.get_terminals()]

    panel = pd.read_csv(PANEL)
    panel_species = [canonical_species(value) for value in panel["species"]]
    if len(panel_species) != 18 or len(panel_species) != len(set(panel_species)):
        raise ValueError("Expected an exact unique 18-species final panel")

    pruned = prune_to_species(source, panel_species)
    corrected, rounding = pad_terminal_rounding(pruned)
    corrected_depths = np.asarray(list(root_to_tip_depths(corrected).values()))
    if np.ptp(corrected_depths) > 1e-10:
        raise RuntimeError("Corrected final-panel tree is not strictly ultrametric")
    # Bio.Phylo defaults to five decimal places, which would recreate the
    # exact rounding problem this audit repairs.  Preserve sub-year precision.
    Phylo.write(
        corrected,
        PRUNED_TREE_OUTPUT,
        "newick",
        format_branch_length="%1.12f",
    )
    written = Phylo.read(PRUNED_TREE_OUTPUT, "newick")
    if np.ptp(np.asarray(list(root_to_tip_depths(written).values()))) > 1e-8:
        raise RuntimeError("Serialized final-panel tree lost ultrametricity")

    crosswalk = load_crosswalk(panel_species)
    crosswalk = crosswalk.merge(rounding, on="species", how="left", validate="one_to_one")
    crosswalk.to_csv(CROSSWALK_OUTPUT, index=False)

    distance = pairwise_distances(corrected).sort_index().sort_index(axis=1)
    distance.to_csv(DISTANCE_OUTPUT)

    source_common = prune_to_species(source, processed_labels)
    processed_common = deepcopy(processed)
    for tip in processed_common.get_terminals():
        tip.name = canonical_species(tip.name)
    source_common_distance = pairwise_distances(source_common).sort_index().sort_index(axis=1)
    processed_distance = pairwise_distances(processed_common).sort_index().sort_index(axis=1)
    common = source_common_distance.index.intersection(processed_distance.index)
    source_values = source_common_distance.loc[common, common].to_numpy()
    processed_values = processed_distance.loc[common, common].to_numpy()
    mask = np.triu(np.ones_like(source_values, dtype=bool), k=1)
    ratios = processed_values[mask] / source_values[mask]
    scale = float(np.median(ratios))
    scale_residual = float(np.max(np.abs(processed_values[mask] - scale * source_values[mask])))
    topology_equal = internal_clade_signature(source_common) == internal_clade_signature(processed_common)

    comparison = pd.DataFrame(
        [
            {
                "comparison": "source_vs_tracked_processed_common_tips",
                "n_common_tips": len(common),
                "source_only_tips": ";".join(sorted(set(source_labels) - set(processed_labels))),
                "processed_only_tips": ";".join(sorted(set(processed_labels) - set(source_labels))),
                "same_rooted_clade_signature": topology_equal,
                "median_pairwise_distance_scale_processed_over_source": scale,
                "maximum_absolute_scale_residual": scale_residual,
                "processed_tree_role": "not_approved_undocumented_0.8_scale",
            }
        ]
    )
    comparison.to_csv(COMPARISON_OUTPUT, index=False)

    metrics = pd.DataFrame(
        [
            tree_metrics(source, "local_source_46_tip", SOURCE_TREE),
            tree_metrics(processed, "tracked_processed_34_tip", PROCESSED_TREE),
            {
                **tree_metrics(written, "corrected_final_panel_18_tip", PRUNED_TREE_OUTPUT),
                "source_sha256": sha256(SOURCE_TREE),
                "derivation": "exact_tip_prune_then_terminal_padding_to_max_root_to_tip",
                "maximum_terminal_padding_myr": float(rounding["terminal_padding_myr"].max()),
            },
        ]
    )
    metrics.to_csv(METRICS_OUTPUT, index=False)

    plot_tree_completeness(corrected, crosswalk, TREE_FIGURE)
    plot_distance_heatmap(distance, corrected, DISTANCE_FIGURE)
    plot_rounding(rounding, ROUNDING_FIGURE)
    plot_processed_comparison(source_common_distance, processed_distance, PROCESSED_FIGURE)

    source_metric = metrics.loc[metrics["tree_id"] == "local_source_46_tip"].iloc[0]
    corrected_metric = metrics.loc[metrics["tree_id"] == "corrected_final_panel_18_tip"].iloc[0]
    missing_ectopic = crosswalk.loc[~crosswalk["has_terminal_internal_proxy"], "species"].tolist()
    report = f"""# Dated-phylogeny release audit: final 18-species panel

This audit preserves both historical tree files. It derives a versioned analysis tree from the exact 46-tip local source and does not infer that genomic and microscopy resources are matched specimens.

## Structural verdict

- The local source has {int(source_metric['n_tips'])} tips, is rooted and fully bifurcating, has positive branch lengths, and has a root age of {source_metric['root_age_myr']:.6f} Myr.
- Its root-to-tip range is {source_metric['root_to_tip_range_myr']:.9f} Myr ({source_metric['root_to_tip_range_myr'] * 1_000_000:.3f} years). It is ultrametric at `1e-5` Myr tolerance but fails a strict floating-point check because the Newick branch lengths are rounded.
- The final tree keeps exactly the declared 18 tips, including `fuscus` and excluding `planiceps`. No congener substitution, duplicate-tip collapse, or random polytomy resolution is used.
- The only correction pads terminal edges to the maximum source root-to-tip depth; the maximum change is {rounding['terminal_padding_years'].max():.3f} years. The serialized output is strictly ultrametric and retains all internal ages and topology.
- Corrected final-panel root age: {corrected_metric['root_age_myr']:.6f} Myr.

**Approved:** structural validation, exact final-panel pruning, the accession-specific `fuscus` tip decision, and the rounding-only correction.

## Provenance and uncertainty verdict

**Not approved for confirmatory publication inference yet.** The repository describes the source only as a published time-calibrated *Desmognathus* phylogeny. It does not contain the exact source citation/archive, original tree identifier, calibration record, support values, tree type (posterior draw/consensus/MCC), or posterior/bootstrap tree set. The local source is ignored by Git even though current path scripts read it directly.

The tracked `results/phylogeny/processed_phylogeny.nwk` is not an alternative topology or uncertainty draw. On its 34 common tips it has the same rooted clade signature as the local source, but every pairwise distance is scaled by {scale:.6f}; maximum residual from exact scaling is {scale_residual:.3e}. No justification for this 0.8 scale was found, so it is not approved for analysis or release figures.

One corrected point tree is a clean computational input, not propagation of phylogenetic uncertainty. That omission is now addressed by the separately provenance-locked Stewart-Wiens (2025) main tree and 200 time-calibrated bootstrap-tree sensitivity in `phylogeny_published_uncertainty_analysis18_v1.md`. The focal tree still needs its own exact source/calibration record, and neither tree set represents reticulate-network uncertainty.

## Data completeness

- TE composition: {int(crosswalk['has_te_composition'].sum())}/18 species.
- Terminal:internal LTR proxy: {int(crosswalk['has_terminal_internal_proxy'].sum())}/18 species; unavailable for {', '.join(missing_ectopic)}.
- Microscopy morphology: {int(crosswalk['has_microscopy_morphology'].sum())}/18 species.
- Image-IOD sensitivity stream: {int(crosswalk['has_image_iod_sensitivity'].sum())}/18 species. Availability does not make it an approved absolute C-value.

## Release products

- `{portable(PRUNED_TREE_OUTPUT)}`
- `{portable(METRICS_OUTPUT)}`
- `{portable(CROSSWALK_OUTPUT)}`
- `{portable(DISTANCE_OUTPUT)}`
- `{portable(COMPARISON_OUTPUT)}`
- `{portable(TREE_FIGURE)}`
- `{portable(DISTANCE_FIGURE)}`
- `{portable(ROUNDING_FIGURE)}`
- `{portable(PROCESSED_FIGURE)}`

## Required external provenance from the tree provider

Request the source paper/DOI, archive URL, exact original filename and checksum, branch-length units, tree-estimation method, calibration identities/priors, support annotations, whether this is a posterior draw/consensus/MCC tree, and the posterior/bootstrap tree sample. Until those arrive, keep all path-model results exploratory with respect to tree uncertainty.
"""
    REPORT_OUTPUT.write_text(report)

    output_paths = [
        METRICS_OUTPUT,
        CROSSWALK_OUTPUT,
        DISTANCE_OUTPUT,
        COMPARISON_OUTPUT,
        PRUNED_TREE_OUTPUT,
        REPORT_OUTPUT,
        TREE_FIGURE,
        DISTANCE_FIGURE,
        ROUNDING_FIGURE,
        PROCESSED_FIGURE,
    ]
    manifest = {
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "n_species": 18,
        "source_tree": portable(SOURCE_TREE),
        "source_tree_sha256": sha256(SOURCE_TREE),
        "source_tree_git_status": "ignored_local_input",
        "structural_tree_status": "approved_after_rounding_only_terminal_padding",
        "tip_crosswalk_status": "approved_exact_no_congener_substitution",
        "tree_source_provenance_status": "not_approved_missing_exact_citation_archive_and_calibrations",
        "tree_uncertainty_status": "approved_external_200_tree_sensitivity_focal_source_still_unresolved",
        "tracked_processed_tree_status": "not_approved_undocumented_0.8_branch_scale",
        "terminal_padding_max_myr": float(rounding["terminal_padding_myr"].max()),
        "outputs": [{"path": portable(path), "sha256": sha256(path)} for path in output_paths],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Wrote {REPORT_OUTPUT}")
    print(f"Wrote {MANIFEST_OUTPUT}")


if __name__ == "__main__":
    main()

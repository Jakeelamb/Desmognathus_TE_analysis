#!/usr/bin/env python3
"""Legacy builder for the superseded analysis18 presentation notebook tree.

The canonical collaborator bundle is built by
``scripts/processing/build_research_review_notebooks.py``. This file preserves
the historical notebook-generation code for provenance only.

The notebooks consume frozen results and never launch dnaPipeTE, RepeatMasker,
RepeatModeler, TEsorter, segmentation inference, model training, or path-model
fitting.  They are review surfaces for collaborators, not hidden pipelines.
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = ROOT / "notebooks/legacy_analysis18"

FINAL18 = [
    "amphileucus",
    "anicetus",
    "apalachicolae",
    "auriculatus",
    "bairdi",
    "campi",
    "fuscus",
    "gvnigeusgwotli",
    "intermedius",
    "kanawha",
    "marmoratus",
    "mavrokoilius",
    "monticola",
    "ocoee",
    "perlapsus",
    "tilleyi",
    "valtos",
    "welteri",
]


def md(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def setup_cell(domain: str):
    return code(
        f"""
from pathlib import Path
import json
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import HTML, IFrame, Image, Markdown, display

ROOT = Path.cwd().resolve()
for candidate in [ROOT, ROOT.parent, ROOT.parent.parent]:
    if (candidate / "results/data/corrected").exists():
        ROOT = candidate
        break
else:
    raise RuntimeError("Open this notebook from the Desmognathus_TE repository")

CELL_ROOT = ROOT.parent / "cellprofiler_test"
FIGURE_DIR = ROOT / "results/figures/legacy_analysis18/{domain}"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
FINAL18 = {FINAL18!r}
TREE_ORDER = [
    "gvnigeusgwotli", "amphileucus", "kanawha", "mavrokoilius", "marmoratus", "intermedius",
    "welteri", "fuscus", "auriculatus", "valtos", "bairdi", "campi", "tilleyi", "anicetus",
    "perlapsus", "monticola", "ocoee", "apalachicolae",
]
assert set(TREE_ORDER) == set(FINAL18)

pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 120)
sns.set_theme(style="ticks", context="notebook", font_scale=1.0)
plt.rcParams.update({{
    "figure.dpi": 120,
    "savefig.dpi": 240,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.family": "DejaVu Sans",
}})

def read_csv(relative, **kwargs):
    path = ROOT / relative
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)

def show_figure(relative, width=1100):
    path = ROOT / relative
    if not path.exists():
        raise FileNotFoundError(path)
    display(Image(filename=str(path), width=width))

def artifact_table(paths):
    rows = []
    for relative in paths:
        path = ROOT / relative
        rows.append({{
            "artifact": relative,
            "exists": path.exists(),
            "size_mb": round(path.stat().st_size / 1024**2, 3) if path.exists() else np.nan,
        }})
    result = pd.DataFrame(rows)
    if not result.exists.all():
        raise FileNotFoundError(result.loc[~result.exists, "artifact"].tolist())
    return result

def save_show(fig, filename):
    path = FIGURE_DIR / filename
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    display(fig)
    plt.close(fig)
    print(f"Saved: {{path.relative_to(ROOT)}}")
    return path

print(f"Repository: {{ROOT}}")
print(f"Presentation figure directory: {{FIGURE_DIR.relative_to(ROOT)}}")
"""
    )


def notebook(title: str, domain: str, cells: List[object]):
    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.8"},
        "desmognathus_analysis_id": "analysis18_v1",
        "presentation_domain": domain,
        "expensive_upstream_tools_executed": False,
    }
    nb["cells"] = [
        md(
            f"""
# {title}

**Final 18-species collaborator audit · `analysis18_v1`**

This notebook is a lightweight presentation and inspection surface. It loads frozen tables, manifests, trees, HTML reports, and figures. It does **not** rerun genome annotation, repeat discovery, read mapping, segmentation, model training, or phylogenetic path fitting.

The aim is to make every displayed claim traceable to a file while keeping historical results visibly separate from corrected release-facing results.
"""
        ),
        setup_cell(domain),
        *cells,
    ]
    return nb


def phylogeny_notebook():
    cells = [
        md(
            """
## 1. Frozen tree sources and release status

The focal collaborator tree, published optimal tree, and 200 published time trees are separate hypotheses. The final-18 tree is an exact tip prune with only a ≤2-year terminal correction for Newick rounding. No congener substitution is used.
"""
        ),
        code(
            """
artifacts = [
    "input_data/phylogeny/desmo900dated_test.tre",
    "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk",
    "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk",
    "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex",
    "results/data/corrected/phylogeny/phylogeny_source_registry_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv",
]
display(artifact_table(artifacts))
source_registry = read_csv("results/data/corrected/phylogeny/phylogeny_source_registry_v1.csv")
tree_metrics = read_csv("results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv")
display(source_registry)
display(tree_metrics)
"""
        ),
        md(
            """
## 2. Reproduce the trimming decision in memory

This cell demonstrates the biological unit and pruning rule without estimating a tree. It prunes a copy of the 46-tip focal source to exact final-panel names, then compares rooted clade signatures with the frozen corrected tree. Branch lengths are not rewritten here; the frozen writer and manifest own the rounding correction.
"""
        ),
        code(
            """
from copy import deepcopy
from Bio import Phylo

def canonical_tip(label):
    text = str(label).strip().replace(" ", "_")
    if text.startswith("D._"):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    elif text.lower().startswith("desmognathus_"):
        text = text.split("_", 1)[1]
    return text.lower()

def rooted_clade_signature(tree):
    all_tips = frozenset(canonical_tip(t.name) for t in tree.get_terminals())
    signature = set()
    for clade in tree.get_nonterminals():
        descendants = frozenset(canonical_tip(t.name) for t in clade.get_terminals())
        if 1 < len(descendants) < len(all_tips):
            signature.add(descendants)
    return signature

source_tree = Phylo.read(ROOT / "input_data/phylogeny/desmo900dated_test.tre", "newick")
pruned = deepcopy(source_tree)
dropped = sorted(canonical_tip(t.name) for t in pruned.get_terminals() if canonical_tip(t.name) not in FINAL18)
trim_audit = pd.DataFrame([
    {
        "source_tip": terminal.name,
        "canonical_species": canonical_tip(terminal.name),
        "action": "retain" if canonical_tip(terminal.name) in FINAL18 else "prune",
        "reason": (
            "exact member of frozen final-18 species panel"
            if canonical_tip(terminal.name) in FINAL18
            else "outside frozen final-18 species panel"
        ),
        "evidence": "exact canonical-tip set comparison; no congener substitution",
        "final_tip": canonical_tip(terminal.name) if canonical_tip(terminal.name) in FINAL18 else "",
    }
    for terminal in source_tree.get_terminals()
]).sort_values(["action", "canonical_species"])
assert trim_audit.source_tip.is_unique
assert (trim_audit.action == "retain").sum() == 18
for terminal in list(pruned.get_terminals()):
    if canonical_tip(terminal.name) not in FINAL18:
        pruned.prune(terminal)
for terminal in pruned.get_terminals():
    terminal.name = canonical_tip(terminal.name)

frozen = Phylo.read(ROOT / "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk", "newick")
assert sorted(t.name for t in pruned.get_terminals()) == sorted(FINAL18)
assert sorted(canonical_tip(t.name) for t in frozen.get_terminals()) == sorted(FINAL18)
assert rooted_clade_signature(pruned) == rooted_clade_signature(frozen)

print(f"Source tips: {len(source_tree.get_terminals())}")
print(f"Dropped tips: {len(dropped)}")
print(f"Final tips: {len(pruned.get_terminals())}")
print("Rooted topology matches frozen final-18 tree: PASS")
display(trim_audit)

fig, axes = plt.subplots(1, 2, figsize=(17, 12), gridspec_kw={"width_ratios": [1.35, 1]})
Phylo.draw(
    source_tree, axes=axes[0], do_show=False, show_confidence=False,
    label_func=lambda clade: canonical_tip(clade.name) if clade.name else None,
)
Phylo.draw(
    pruned, axes=axes[1], do_show=False, show_confidence=False,
    label_func=lambda clade: clade.name if clade.name else None,
)
axes[0].set_title("Source focal tree (46 tips; before exact panel prune)")
axes[1].set_title("Frozen analysis tree (18 tips; after prune)")
for axis in axes:
    axis.set_ylabel("")
    axis.set_xlabel("Branch time (Myr)")
fig.suptitle("Exact phylogenetic tip trimming: source versus final analysis panel", fontsize=15)
fig.tight_layout()
save_show(fig, "phylogeny_source46_vs_final18_trim_analysis18_v1.png")
"""
        ),
        md("## 3. Taxon/accession crosswalk and independent evidence streams"),
        code(
            """
crosswalk = read_csv("results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv")
assert set(crosswalk.species) == set(FINAL18)
columns = [
    "species", "tree_tip", "tree_match_type", "te_sra_accession", "te_assembly_accession",
    "has_terminal_internal_proxy", "microscopy_n_specimens", "microscopy_n_images",
    "genome_and_microscopy_are_independent_samples", "taxon_decision_note",
]
display(crosswalk[columns].sort_values("species"))
"""
        ),
        md("## 4. Tree structure, topology, and data completeness figures"),
        code(
            """
for figure in [
    "results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_data_completeness_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_root_to_tip_rounding_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_processed_scale_comparison_v1.png",
]:
    show_figure(figure)
"""
        ),
        md("## 5. Published time-tree uncertainty"),
        code(
            """
uncertainty = read_csv("results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv")
topology = read_csv("results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv")
display(uncertainty.describe().T)
display(topology)
for figure in [
    "results/figures/corrected/phylogeny/phylogeny_root_age_uncertainty_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_published_bootstrap_density_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md(
            """
## 6. Interpretation and paper-facing boundary

- **Approved:** exact final-18 matching, positive bifurcating branches, published main-tree sensitivity, and propagation over all 200 published time trees.
- **Open:** the collaborator focal tree's original publication/archive identifier, tree type, and calibration record.
- **Do not use:** the tracked 34-tip processed tree as an independent hypothesis; it is the same common-tip topology with an undocumented exact 0.8 branch-length scale.

Figure design follows the primary-paper expectation that trees show explicit tips, branch-time scale, support/uncertainty, and trait completeness. See Stewart & Wiens (2025, DOI `10.1016/j.ympev.2024.108272`) and the [primary-literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).

Exact audit: [phylogeny_release_audit_analysis18_v1.md](../../plans/publication-readiness-deep-audit/phylogeny_release_audit_analysis18_v1.md).
"""
        ),
    ]
    return notebook("Phylogenetic tree trimming, provenance, and uncertainty", "phylogeny", cells)


def data_tables_notebook():
    cells = [
        md(
            """
## 1. What counts as a final-panel row

The genomic, microscopy, and tree streams are joined at the species level but do not represent the same physical specimen. dnaPipeTE retry suffixes are computational run aliases. The genomic *fuscus* resource is accession-explicit and placed on the *fuscus* tip; there is no *planiceps* analysis row.
"""
        ),
        code(
            """
panel = read_csv("path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv")
readiness = read_csv("path_analysis/data/derived/analysis_species_readiness.csv")
crosswalk = read_csv("results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv")
support = read_csv("results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv")
assert set(panel.species.str.replace("D.", "", regex=False).str.lower()) == set(FINAL18)
assert set(crosswalk.species) == set(FINAL18)
assert set(support.species) == set(FINAL18)
assert panel.species.is_unique and crosswalk.species.is_unique and support.species.is_unique
assert crosswalk.tree_tip.is_unique
display(crosswalk[[
    "species", "te_sra_accession", "te_assembly_accession", "tree_tip",
    "has_te_composition", "has_terminal_internal_proxy", "has_microscopy_morphology",
    "genome_and_microscopy_are_independent_samples", "taxon_decision_note",
]].sort_values("species"))
"""
        ),
        md("## 2. Data completeness and biological replication"),
        code(
            """
completeness = crosswalk.set_index("species")[[
    "has_te_composition", "has_terminal_internal_proxy", "has_microscopy_morphology",
    "has_image_iod_sensitivity",
]].astype(float)
support_indexed = support.set_index("species")
completeness["≥2 microscopy specimens"] = (support_indexed.n_specimens >= 2).astype(float)
completeness = completeness.loc[TREE_ORDER]
completeness.columns = [
    "TE\\ncomposition", "LTR terminal:\\ninternal", "Microscopy\\nmorphology",
    "Relative\\nnuclear IOD", "≥2 microscopy\\nspecimens",
]
annotations = pd.DataFrame(index=completeness.index, columns=completeness.columns, dtype=object)
for column in completeness.columns[:-1]:
    annotations[column] = np.where(completeness[column].eq(1), "available", "missing")
annotations["≥2 microscopy\\nspecimens"] = np.where(
    completeness["≥2 microscopy\\nspecimens"].eq(1), "≥2 specimens", "1 specimen"
)

fig, ax = plt.subplots(figsize=(10.5, 7.4))
sns.heatmap(
    completeness,
    cmap=sns.color_palette(["#d9dde3", "#167d73"], as_cmap=True),
    vmin=0, vmax=1, linewidths=0.7, linecolor="white", cbar=False,
    annot=annotations, fmt="", ax=ax,
)
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("Final-18 evidence completeness and biological replication")
ax.tick_params(axis="x", rotation=0)
save_show(fig, "final18_evidence_completeness_analysis18_v1.png")
"""
        ),
        md("## 3. Source manifests and immutable identifiers"),
        code(
            """
source_manifest = read_csv("path_analysis/data/templates/source_manifest.csv")
source_registry = read_csv("path_analysis/data/derived/source_file_registry.csv")
taxonomy = read_csv("path_analysis/data/templates/species_taxonomy_crosswalk.csv")
genomic = source_manifest[source_manifest.evidence_stream.fillna("").str.contains("genom", case=False)]
display(genomic[[c for c in [
    "source_id", "short_citation", "sra_accessions", "assembly_accession",
    "public_taxon_name", "analysis_taxon_name", "tree_tip", "analysis_role",
    "inclusion_decision", "decision_basis_source_id",
] if c in genomic.columns]])
display(source_registry[["source_id", "local_path", "sha256", "exists", "origin"]].head(30))
display(taxonomy[taxonomy.current_species.isin(FINAL18)])
"""
        ),
        md("## 4. Corrected artifact catalog—load, do not recompute"),
        code(
            """
selected_tables = [
    "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv",
    "results/data/corrected/divergence/divergence_summary_statistics_by_species_analysis18_v1.csv",
    "results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv",
    "results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv",
]
catalog = artifact_table(selected_tables)
shapes = []
for relative in selected_tables:
    frame = read_csv(relative)
    shapes.append({"artifact": relative, "rows": len(frame), "columns": len(frame.columns)})
display(catalog.merge(pd.DataFrame(shapes), on="artifact"))
"""
        ),
        md("## 5. Pipeline map and analysis boundaries"),
        code(
            """
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(16, 6.8))
ax.axis("off")
boxes = [
    (0.02, 0.66, 0.14, 0.15, "WGS reads\\nSRX accession", "#d9edf0"),
    (0.22, 0.66, 0.19, 0.15, "dnaPipeTE + custom library\\n(frozen upstream)", "#d9edf0"),
    (0.47, 0.66, 0.20, 0.15, "TE composition, diversity,\\ndivergence, LTR proxies", "#d9edf0"),
    (0.02, 0.17, 0.17, 0.15, "Blood-smear images\\nindependent specimens", "#f5e4ce"),
    (0.24, 0.17, 0.18, 0.15, "Cell/nucleus masks + QC\\n(frozen upstream)", "#f5e4ce"),
    (0.48, 0.17, 0.19, 0.15, "Upper-tail morphology +\\nrelative nuclear IOD", "#f5e4ce"),
    (0.73, 0.39, 0.12, 0.20, "Exact taxon\\ncrosswalk +\\nfinal-18 trees", "#e5dff2"),
    (0.875, 0.39, 0.105, 0.20, "Exploratory\\nphylogenetic path\\nsensitivity", "#e5dff2"),
]
for x, y, width, height, label, color in boxes:
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, transform=ax.transAxes,
        boxstyle="round,pad=0.012", facecolor=color, edgecolor="#425466", linewidth=1.4,
    ))
    ax.text(x + width / 2, y + height / 2, label, transform=ax.transAxes,
            ha="center", va="center", fontsize=9.5)
arrows = [
    ((0.16,0.735),(0.22,0.735)), ((0.41,0.735),(0.47,0.735)),
    ((0.19,0.245),(0.24,0.245)), ((0.42,0.245),(0.48,0.245)),
    ((0.67,0.735),(0.73,0.54)), ((0.67,0.245),(0.73,0.44)),
    ((0.85,0.49),(0.875,0.49)),
]
for start, end in arrows:
    ax.annotate("", xy=end, xytext=start, xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "-|>", "color": "#425466", "lw": 1.8})
ax.set_title("Analysis18 evidence streams: accession-explicit, non-interchangeable samples", fontsize=15)
save_show(fig, "analysis18_pipeline_evidence_streams_v1.png")
"""
        ),
        md("## 6. Release gates"),
        code(
            """
gates = read_csv("results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv")
display(gates[["component", "status", "evidence", "permitted_language"]])
show_figure("results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png")
"""
        ),
        md(
            """
## 7. Presentation checklist

- Start with the accession/tip crosswalk and state that genome and microscopy specimens are independent.
- Treat a row as one species, not one read file, cell, TE family, or retry directory.
- Load corrected artifacts by exact versioned path. Historical similarly named tables are provenance, not release inputs.
- Preserve `SRX20497025 / GCA_032353935.1` as the current expert-reidentified *fuscus* genomic resource and do not create a *planiceps* row.
- Show replication counts whenever morphology or IOD appears.

The tabular style follows Nature-family expectations for transparent sample counts, machine-readable supplements, and explicit missingness. See the [figure literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).
"""
        ),
    ]
    return notebook("Analysis data tables, identity, and provenance", "data_tables", cells)


def repeat_notebook():
    cells = [
        md(
            """
## 1. Pipeline overview—what is frozen

Genomic WGS reads were subsampled for dnaPipeTE, assembled into repeat components, annotated with the custom TE library, and parsed into relative composition tables. RepeatMasker hit-native classification and divergence were subsequently corrected non-destructively. This notebook reads those products; it never invokes any upstream bioinformatics executable.

Retry suffixes such as `2`, `3`, or `R2` identify alternative computational runs of the same SRX—not additional samples and never additional species.
"""
        ),
        code(
            """
artifacts = [
    "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv",
    "results/data/corrected/dnapipete/dnapipete_absolute_load_sensitivity_analysis18_v1.csv",
    "results/data/corrected/dnapipete/dnaPipeTE_order_breakdown_mass_accounted_analysis18_v2.csv",
    "results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv.manifest.json",
    "results/data/corrected/repeat_landscape/repeatmasker_divergence_landscape_analysis18_v1.csv",
    "results/data/corrected/repeat_landscape/repeatmasker_divergence_landscape_analysis18_v1.manifest.json",
    "results/data/corrected/divergence/divergence_summary_statistics_by_species_analysis18_v1.csv",
    "results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv",
]
display(artifact_table(artifacts))
"""
        ),
        md("## 2. dnaPipeTE relative composition and configured-load sensitivity"),
        code(
            """
mass = read_csv("results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv")
load = read_csv("results/data/corrected/dnapipete/dnapipete_absolute_load_sensitivity_analysis18_v1.csv")
composition = read_csv("results/data/corrected/dnapipete/dnaPipeTE_order_breakdown_mass_accounted_analysis18_v2.csv")
assert set(mass.species) == set(FINAL18)
assert set(load.species) == set(FINAL18)

orders = ["DIRS", "Helitron", "LINE", "LTR", "Maverick", "PLE", "SINE", "TIR", "Transposon Derivatives", "YR"]
stacked = composition.set_index("species")[orders].div(100)
stacked = stacked.mul(composition.set_index("species").order_retained_fraction, axis=0)
stacked["Unresolved"] = composition.set_index("species").order_unresolved_fraction
stacked = stacked.loc[TREE_ORDER]
assert np.allclose(stacked.sum(axis=1), 1)

fig, axes = plt.subplots(1, 2, figsize=(16, 7.2), gridspec_kw={"width_ratios": [1.55, 1]})
stacked.plot(kind="barh", stacked=True, ax=axes[0], width=0.82, colormap="tab20")
axes[0].set_xlabel("Fraction of repeat-aligned mass")
axes[0].set_ylabel("")
axes[0].set_title("Mass-accounted TE order composition")
axes[0].legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)

load_plot = load.set_index("species").loc[TREE_ORDER]
y = np.arange(len(load_plot))
axes[1].scatter(load_plot.repeat_aligned_fraction, y, label="Repeat aligned", color="#244f73", s=45)
axes[1].scatter(load_plot.te_classified_fraction, y, label="TE classified", color="#d07a32", s=45)
for yi, (a, b) in enumerate(zip(load_plot.repeat_aligned_fraction, load_plot.te_classified_fraction)):
    axes[1].plot([b, a], [yi, yi], color="#aeb8c2", lw=1)
axes[1].set_yticks(y, load_plot.index)
axes[1].set_xlabel("Configured 1.5-Gb quantification-sample fraction")
axes[1].set_title("Absolute-load sensitivity (not confirmatory)")
axes[1].legend(frameon=False)
fig.suptitle("Final-18 dnaPipeTE composition and mass accounting", fontsize=15)
fig.tight_layout()
save_show(fig, "repeat_mass_and_composition_analysis18_v1.png")
display(mass.sort_values("superfamily_unresolved_fraction", ascending=False))
"""
        ),
        md(
            """
## 3. RepeatMasker classification correction

The historical merge inherited one dnaPipeTE contig classification for every RepeatMasker hit. The corrected branch uses each hit's native RepeatMasker label and preserves hit and aligned-base totals. The full 5.13-million-row corrected table is not loaded here; its manifest is the review surface, followed by a 20-row preview.
"""
        ),
        code(
            """
manifest_path = ROOT / "results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv.manifest.json"
manifest = json.loads(manifest_path.read_text())
audit = read_csv("plans/publication-readiness-deep-audit/repeatmasker_hit_classification_audit.csv")
sensitivity = read_csv("plans/publication-readiness-deep-audit/repeatmasker_corrected_feature_sensitivity_summary.csv")
display(pd.DataFrame({
    "metric": ["corrected hits", "inclusive aligned bp", "species", "classified hits", "unmapped hits"],
    "value": [manifest["n_hits"], manifest["inclusive_aligned_bp"], manifest["analysis_species_count"],
              manifest["classification_status_counts"]["classified"], manifest["classification_status_counts"]["unmapped"]],
}))
display(audit.T.rename(columns={0: "historical merge audit"}))
display(sensitivity)
preview = pd.read_csv(ROOT / manifest["output_path"], nrows=20)
display(preview)
"""
        ),
        md(
            """
## 4. Corrected RepeatMasker divergence landscape and fingerprint

The landscape denominator is **inclusive corrected RepeatMasker query-coordinate aligned bp within each species**. It is not percent of the genome, assembly, dnaPipeTE sample, or inferred genome size. RepeatMasker hits that overlap in query coordinates are not interval-deduplicated, matching the conserved hit-level manifest. The `50+` bin is an overflow bin.
"""
        ),
        code(
            """
landscape = read_csv("results/data/corrected/repeat_landscape/repeatmasker_divergence_landscape_analysis18_v1.csv")
landscape_manifest = json.loads((ROOT / "results/data/corrected/repeat_landscape/repeatmasker_divergence_landscape_analysis18_v1.manifest.json").read_text())
assert landscape_manifest["validation"]["n_species"] == 18
assert landscape_manifest["validation"]["hit_count"] == manifest["n_hits"]
assert landscape_manifest["validation"]["hit_bp"] == manifest["inclusive_aligned_bp"]

major_orders = ["LTR", "LINE", "TIR", "DIRS", "SINE", "Unclassified"]
landscape["display_order"] = landscape.order.where(landscape.order.isin(major_orders), "Other")
display_orders = ["LTR", "LINE", "TIR", "DIRS", "SINE", "Other", "Unclassified"]
colors = {
    "LTR": "#3b6fb6", "LINE": "#d4673f", "TIR": "#56a270", "DIRS": "#9a67ad",
    "SINE": "#d6a84b", "Other": "#8c9aa6", "Unclassified": "#d9dde2",
}

fig, axes = plt.subplots(6, 3, figsize=(14.5, 17), sharex=True, sharey=True)
for ax, species_name in zip(axes.flat, TREE_ORDER):
    panel = (
        landscape[landscape.species == species_name]
        .groupby(["divergence_bin_start_pct", "display_order"], as_index=False)["percent_species_hit_bp"]
        .sum()
        .pivot(index="divergence_bin_start_pct", columns="display_order", values="percent_species_hit_bp")
        .reindex(index=range(51), columns=display_orders, fill_value=0)
        .fillna(0)
    )
    ax.stackplot(
        panel.index,
        *[panel[column].values for column in display_orders],
        colors=[colors[column] for column in display_orders],
        labels=display_orders,
        linewidth=0,
    )
    ax.set_title(f"D. {species_name}", fontsize=10, fontstyle="italic", loc="left")
    ax.axvline(50, color="#5a6570", lw=0.6, ls=":")
    ax.grid(axis="y", color="#e7eaee", linewidth=0.6)
    ax.set_xlim(0, 50)
for ax in axes[-1, :]:
    ax.set_xlabel("RepeatMasker divergence from consensus (%)\\n(50 = 50+ overflow)")
for ax in axes[:, 0]:
    ax.set_ylabel("% of species aligned hit bp")
handles = [plt.Rectangle((0, 0), 1, 1, color=colors[name]) for name in display_orders]
fig.legend(handles, display_orders, loc="upper center", ncol=7, frameon=False, bbox_to_anchor=(0.5, 0.996))
fig.suptitle("Final-18 hit-native RepeatMasker divergence landscapes", fontsize=16, y=1.018)
fig.tight_layout(rect=(0, 0, 1, 0.965))
save_show(fig, "repeatmasker_divergence_landscapes_analysis18_v1.png")

divergence = read_csv("results/data/corrected/divergence/divergence_summary_statistics_by_species_analysis18_v1.csv")
order_all = divergence[(divergence.group_level == "order") & (divergence.threshold == -1.0)].copy()
main_orders = ["DIRS", "Helitron", "LINE", "LTR", "Maverick", "PLE", "SINE", "TIR", "YR"]
order_all["species"] = order_all.Desmognathus_Species.str.replace("D.", "", regex=False)
matrix = order_all[order_all.group_name.isin(main_orders)].pivot(
    index="species", columns="group_name", values="percent_divergence_median"
).reindex(index=TREE_ORDER, columns=main_orders)

fig, ax = plt.subplots(figsize=(10.5, 7.8))
sns.heatmap(matrix, cmap="mako", linewidths=0.45, linecolor="white", ax=ax,
            cbar_kws={"label": "Median RepeatMasker divergence (%)"})
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("Hit-native RepeatMasker divergence fingerprint")
save_show(fig, "repeatmasker_divergence_fingerprint_analysis18_v1.png")
coverage = read_csv("results/data/corrected/divergence/divergence_threshold_context_coverage_by_species_analysis18_v1.csv")
display(coverage.sort_values("threshold_context_hit_fraction"))
"""
        ),
        md("## 5. Diversity indices and compositional PCA"),
        code(
            """
diversity = read_csv("results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv")
primary = diversity[(diversity.te_level == "order") & (diversity.composition_mode == "classified_conditional")]
display(primary[["species", "observed_richness", "shannon_entropy", "gini_simpson", "hill_q1", "hill_q2", "pielou_evenness"]])
for figure in [
    "results/figures/corrected/diversity_pca/te_order_clr_pca_scores_analysis18_v1.png",
    "results/figures/corrected/diversity_pca/te_order_clr_pca_loadings_analysis18_v1.png",
    "results/figures/corrected/diversity_pca/te_diversity_unresolved_mass_sensitivity_analysis18_v1.png",
    "results/figures/corrected/diversity_pca/te_order_ordinary_vs_phylogenetic_pca_analysis18_v1.png",
    "results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_loading_intervals_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md("## 6. Paired-LTR divergence branch (separate ascertainment question)"),
        code(
            """
ltr = read_csv("results/data/ltr_age/ltr_age_species_summary.csv")
ltr["species_clean"] = ltr.species.str.replace("D.", "", regex=False)
focal_ltr = ltr[ltr.species_clean.isin(FINAL18)].copy()
display(focal_ltr[["species", "gca_accession", "n_pairs_attempted", "n_pairs_estimated", "n_high_confidence_pairs_estimated", "median_k2p_distance", "median_comparable_sites"]].sort_values("n_high_confidence_pairs_estimated"))
"""
        ),
        md(
            """
## 7. Interpretation and paper-facing boundary

- **Approved descriptively:** final-18 relative TE composition, explicit unresolved mass, Shannon/Gini-Simpson/Hill indices, order-level CLR PCA, phylogenetic PCA, and 200-tree PCA sensitivity.
- **Sensitivity only:** configured-denominator absolute load; the historical container, custom-library checksum, runtime log, and repeated sampling are incomplete.
- **Corrected:** hit-native RepeatMasker classification. Alignment `%del/%ins` describe gaps relative to a repeat consensus—not genomic deletion rates.
- **Landscape denominator:** percent of inclusive corrected RepeatMasker query-coordinate aligned bp within species; it is not genome occupancy and overlapping hit intervals are not deduplicated.
- **Supplementary:** superfamily PCA because zeros and rare-category choices affect its geometry.
- **Not promoted:** unrestricted clade PERMANOVA.
- **Paired-LTR branch:** describes detectable, overwhelmingly Gypsy paired elements; it is not a genome-wide insertion-age census.

Nature-style TE figures generally combine a mass-accounted composition panel, an explicit repeat-landscape/divergence panel, raw species points, phylogenetic control, and sensitivity to classification depth. See Wang et al. (2025, *Communications Biology*, DOI `10.1038/s42003-025-08127-3`) and the [literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).

Exact audit: [te_diversity_pca_corrected_analysis18_v1.md](../../plans/publication-readiness-deep-audit/te_diversity_pca_corrected_analysis18_v1.md).
"""
        ),
    ]
    return notebook("Repeat analysis, diversity indices, and compositional PCA", "repeat_analysis", cells)


def ectopic_notebook():
    cells = [
        md(
            """
## 1. Estimand first

The stored quantity is terminal LTR read depth divided by internal-region read depth for structurally annotated LTR elements. It can be discussed as a **terminal:internal deletion-footprint/mapping proxy**. It is not a direct solo-LTR count, an orthologous deletion event, or an ectopic-recombination rate.
"""
        ),
        code(
            """
artifacts = [
    "results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv",
    "results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv",
    "results/data/corrected/ectopic/ectopic_analysis18_v1.manifest.json",
    "plans/publication-readiness-deep-audit/ectopic_recombination_corrected_analysis18_v1.md",
]
display(artifact_table(artifacts))
elements = read_csv("results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv")
species = read_csv("results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv")
print(f"Elements: {len(elements)}; species: {elements.species.nunique()}")
print(f"Coverage ≥80%: {elements.coverage_ge_80pct.sum()}/{len(elements)}")
display(elements[[
    "species", "te_sra_accession", "te_assembly_accession", "element_id", "Complete",
    "domain_count", "terminal_positive_coverage_fraction", "internal_positive_coverage_fraction",
    "ratio_terminal_internal_all_positions", "ratio_terminal_internal_nonzero_only",
]].head(20))
"""
        ),
        md("## 2. Reproduce the per-element ratio and historical comparison"),
        code(
            """
rederived = elements.mean_depth_terminal_all_positions / elements.mean_depth_internal_all_positions
error = np.nanmax(np.abs(rederived - elements.ratio_terminal_internal_all_positions))
assert error < 1e-12
assert elements.nonzero_ratio_reproduction_abs_error.max() < 1e-12
print(f"Maximum zero-aware ratio reconstruction error: {error:.3e}")
print(f"Explicit zero-depth positions retained: {int(elements.depth_positions_explicit_zero.sum()):,}")
print("Per-element arithmetic identity and historical nonzero-only reproduction: PASS")
"""
        ),
        md(
            """
## 3. Robust species summaries and conditional uncertainty

Bootstrap intervals resample detected elements within a species. They quantify conditional variation among the ascertained candidates; they do not estimate biological between-specimen or between-population uncertainty, because each genomic branch is accession-specific and mapping provenance is incomplete.
"""
        ),
        code(
            """
primary = species[species.analysis_branch == "all_5plus_domain_elements"].copy()
assert primary.species.nunique() == 16
display(primary[[
    "species", "n_elements", "ratio_arithmetic_mean", "ratio_median", "ratio_geometric_mean",
    "median_bootstrap_ci_low", "median_bootstrap_ci_high", "ratio_max",
    "max_abs_leave_one_out_mean_shift", "largest_element_fraction_of_ratio_sum",
]].sort_values("ratio_median", ascending=False))
for figure in [
    "results/figures/corrected/ectopic/ectopic_element_log2_ratio_by_species_analysis18_v1.png",
    "results/figures/corrected/ectopic/ectopic_mean_vs_median_analysis18_v1.png",
    "results/figures/corrected/ectopic/ectopic_mean_influence_analysis18_v1.png",
    "results/figures/corrected/ectopic/ectopic_coverage_qc_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md("## 4. Branch and missing-species sensitivity"),
        code(
            """
branch_table = species.pivot(index="species", columns="analysis_branch", values="ratio_median")
display(branch_table)
missing = sorted(set(FINAL18) - set(primary.species))
print("Final-panel species lacking the assembly-dependent proxy:", missing)
assert missing == ["kanawha", "valtos"]
"""
        ),
        md("## 5. Does the proxy support the proposed path family?"),
        code(
            """
rankings = read_csv("results/data/corrected/path_analysis/corrected_path_data_sensitivity_rankings_analysis18_v1.csv")
terminal = rankings[(rankings.family == "terminal_internal_iod") & (rankings["rank"] == 1)]
display(terminal[["tree_id", "iod_subset", "model", "p", "CICc", "w", "release_gate_status"]])
assert terminal.release_gate_status.eq("blocked_no_globally_supported_model").all()
print("All full-panel terminal:internal candidate families fail global fit: PASS")
"""
        ),
        md(
            """
## 6. Interpretation and paper-facing boundary

- Exact ≥5/6-domain elements, zero retention, coverage metrics, robust location estimates, bootstrap intervals, and influence diagnostics are reproducible.
- The arithmetic mean is outlier-sensitive; the median/geometric mean and full element distribution must accompany it.
- `kanawha` and `valtos` have no assembly-dependent proxy and are never silently imputed.
- Mapping reference, MAPQ, multimapper, duplicate, and secondary-alignment provenance still need recovery.
- A structural intact-versus-solo LTR analysis or manually validated subset is required before ectopic-recombination language.
- The corrected phylogenetic path family is globally rejected, so no mechanism claim is supported by current path-model fit.

Gold-standard figures pair the full element distribution with structural cartoons/calls, uncertainty, mappability/QC, and sensitivity to high-confidence elements. See Frahry et al. (2015, DOI `10.1007/s00239-014-9663-7`) and the [literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).

Exact audit: [ectopic_recombination_corrected_analysis18_v1.md](../../plans/publication-readiness-deep-audit/ectopic_recombination_corrected_analysis18_v1.md).
"""
        ),
    ]
    return notebook("LTR terminal:internal deletion-footprint proxy", "ectopic_proxy", cells)


def cell_notebook():
    cells = [
        md(
            """
## 1. Measurement layers must remain separate

1. **Segmentation and linkage:** Are cell/nucleus masks and one-to-one pairs credible?
2. **Morphometry:** Which upper-tail or central-tendency cell/nucleus estimand is being reported?
3. **Image density:** Is nuclear IOD comparable across images and calibrated to DNA content?

The current release supports linked upper-tail morphometry and relative nuclear-IOD sensitivity. It does **not** support absolute genome size in picograms.
"""
        ),
        code(
            """
artifacts = [
    "results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv",
    "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk",
]
display(artifact_table(artifacts))
if not CELL_ROOT.exists():
    raise FileNotFoundError(f"Microscopy workspace missing: {CELL_ROOT}")
"""
        ),
        md("## 2. Live HTML mask and analysis viewers"),
        code(
            """
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import threading

class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

if "_viewer_server" not in globals():
    handler = partial(QuietHandler, directory=str(ROOT.parent))
    _viewer_server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    _viewer_thread = threading.Thread(target=_viewer_server.serve_forever, daemon=True)
    _viewer_thread.start()

port = _viewer_server.server_address[1]
base = f"http://127.0.0.1:{port}/cellprofiler_test"
viewers = pd.DataFrame([
    {"viewer": "Cell–nucleus linkage", "path": "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linkage/index.html", "status": "diagnostic linkage/QC"},
    {"viewer": "Tile annotation gallery", "path": "output/tile_annotation_bundle_v1/index.html", "status": "training-domain tiles; inspect scope"},
    {"viewer": "Raw / overlay / mask comparison", "path": "output/tile_bootstrap_review_v1/index.html", "status": "model-QC gallery"},
    {"viewer": "Reviewed species statistics", "path": "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linked_species_stats_reviewed/index.html", "status": "historical pg labels—not approved"},
    {"viewer": "Historical cell phylogenetic analysis", "path": "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/phylogenetic_analysis/index.html", "status": "historical—not release inference"},
    {"viewer": "Historical integrated analysis", "path": "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/integrated_multinode_phylo/index.html", "status": "historical—not release inference"},
])
viewers["url"] = viewers.path.map(lambda path: f"{base}/{path}")
for path in viewers.path:
    if not (CELL_ROOT / path).exists():
        raise FileNotFoundError(CELL_ROOT / path)

cards = "".join(
    f'<p><a href="{row.url}" target="_blank"><strong>{row.viewer}</strong></a> — {row.status}</p>'
    for row in viewers.itertuples()
)
display(HTML(f'<div style="padding:14px;border:1px solid #ccd6df;border-radius:12px">{cards}</div>'))
display(IFrame(src=viewers.iloc[0].url, width="100%", height=650))
print("If an embedded viewer is blank after reopening the notebook, rerun only this lightweight cell.")
"""
        ),
        md("## 3. Static mask examples and validation domain"),
        code(
            """
prediction_manifest = pd.read_csv(CELL_ROOT / "output/tile_bootstrap_review_v1/predictions_manifest.csv")
display(prediction_manifest[["species", "filename", "tile_name", "n_masks", "coverage", "median_area_px", "is_annotated"]])
examples = prediction_manifest.drop_duplicates("species").head(3)
fig, axes = plt.subplots(len(examples), 2, figsize=(12, 4 * len(examples)))
if len(examples) == 1:
    axes = np.array([axes])
for row_index, row in enumerate(examples.itertuples()):
    raw = plt.imread(CELL_ROOT / "output/tile_bootstrap_review_v1" / row.preview_raw_path)
    overlay = plt.imread(CELL_ROOT / "output/tile_bootstrap_review_v1" / row.preview_overlay_path)
    axes[row_index, 0].imshow(raw, cmap="gray")
    axes[row_index, 1].imshow(overlay)
    axes[row_index, 0].set_title(f"{row.species}: raw tile")
    axes[row_index, 1].set_title(f"{row.species}: predicted cell masks")
    for axis in axes[row_index]:
        axis.axis("off")
fig.suptitle("Static segmentation examples (viewer remains the full audit surface)", fontsize=14)
fig.tight_layout()
save_show(fig, "cell_mask_examples_analysis18_v1.png")
"""
        ),
        md(
            """
## 4. Segmentation benchmark and linkage support

> **Publication gate:** the available segmentation benchmark contains two tiles, represents 0 of 18 focal species, and the nucleus train/test split overlaps at the source-image level. The plots below are diagnostic. They do not establish final-panel generalization.
"""
        ),
        code(
            """
validation = read_csv("results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv")
support = read_csv("results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv")
assert set(support.species) == set(FINAL18)
assert support.n_frozen_pairs.sum() == 900
display(validation)
display(support.sort_values(["n_specimens", "n_images", "manual_keep_fraction"]))
show_figure("results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png")
"""
        ),
        md("## 5. What does the selected-50 trait estimate?"),
        code(
            """
estimands = read_csv("results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv")
rank_stability = read_csv("results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv")
bootstrap = read_csv("results/data/corrected/microscopy/microscopy_selected50_hierarchical_bootstrap_analysis18_v1.csv")
display(estimands)
display(rank_stability.sort_values(["metric", "spearman_rho"], ascending=[True, False]))
display(bootstrap)
for figure in [
    "results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md("## 6. Relative nuclear-IOD quality and algebraic dependence"),
        code(
            """
iod = read_csv("results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv")
image_qc = read_csv("results/data/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.csv")
display(iod.sort_values("relative_iod_to_fuscus_image_qc_pass", ascending=False))
display(image_qc)
for figure in [
    "results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md(
            """
## 7. Pairwise morphology and relative nuclear-IOD relationships

The next plot uses the frozen equal-image selected-50 morphology and image-QC relative IOD. It reports raw Spearman correlation and a Brownian-covariance PGLS sensitivity on the published final-18 tree. This is lightweight refitting of three pairwise summaries—not a rerun of segmentation or the path analysis. IOD remains a candidate image phenotype, not absolute genome size.
"""
        ),
        code(
            """
from Bio import Phylo
from scipy.stats import spearmanr, t as student_t

path_input = read_csv("results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv")
traits = path_input[
    (path_input.morphology_estimator == "image_balanced_selected50") &
    (path_input.iod_subset == "image_qc_pass")
].copy().set_index("species").loc[FINAL18]

tree = Phylo.read(ROOT / "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk", "newick")
tip_map = {terminal.name.replace("D.", ""): terminal for terminal in tree.get_terminals()}
if set(tip_map) != set(FINAL18):
    raise RuntimeError("Published tree and trait panel differ")
root = tree.root
covariance = np.zeros((len(FINAL18), len(FINAL18)))
for i, left in enumerate(FINAL18):
    for j, right in enumerate(FINAL18):
        ancestor = tree.common_ancestor(tip_map[left], tip_map[right])
        covariance[i, j] = tree.distance(root, ancestor)
covariance += np.eye(len(FINAL18)) * 1e-8

def standardize(values):
    values = np.asarray(values, dtype=float)
    return (values - values.mean()) / values.std(ddof=1)

def bm_pgls(x, y):
    xz, yz = standardize(np.log(x)), standardize(np.log(y))
    design = np.column_stack([np.ones(len(xz)), xz])
    inverse = np.linalg.pinv(covariance)
    bread = np.linalg.pinv(design.T @ inverse @ design)
    beta = bread @ design.T @ inverse @ yz
    residual = yz - design @ beta
    dof = len(yz) - design.shape[1]
    sigma2 = float(residual.T @ inverse @ residual / dof)
    se = np.sqrt(np.diag(bread * sigma2))
    critical = student_t.ppf(0.975, dof)
    slope_p = 2 * student_t.sf(abs(beta[1] / se[1]), dof)
    return xz, yz, beta, se, (beta[1] - critical * se[1], beta[1] + critical * se[1]), slope_p

pairs = [
    ("relative_nuclear_iod_proxy", "nuc_area_um2", "Relative nuclear IOD", "Nucleus area"),
    ("relative_nuclear_iod_proxy", "cell_area_um2", "Relative nuclear IOD", "Cell area"),
    ("nuc_area_um2", "cell_area_um2", "Nucleus area", "Cell area"),
]
fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
results = []
for axis, (xcol, ycol, xlabel, ylabel) in zip(axes, pairs):
    xz, yz, beta, se, ci, p_pgls = bm_pgls(traits[xcol].values, traits[ycol].values)
    rho, p_raw = spearmanr(traits[xcol], traits[ycol])
    grid = np.linspace(xz.min() - 0.1, xz.max() + 0.1, 100)
    axis.scatter(xz, yz, s=48, color="#236a73", edgecolor="white", linewidth=0.7, zorder=3)
    axis.plot(grid, beta[0] + beta[1] * grid, color="#d47a32", lw=2.2)
    for species_name, xvalue, yvalue in zip(FINAL18, xz, yz):
        axis.annotate(species_name[:4], (xvalue, yvalue), xytext=(3, 3), textcoords="offset points", fontsize=6.5, alpha=0.8)
    axis.axhline(0, color="#d9dee3", lw=0.8)
    axis.axvline(0, color="#d9dee3", lw=0.8)
    axis.set_xlabel(f"ln {xlabel} (standardized)")
    axis.set_ylabel(f"ln {ylabel} (standardized)")
    axis.set_title(
        f"n={len(traits)}; Spearman ρ={rho:.2f}\\n"
        f"BM (λ=1) PGLS β={beta[1]:.2f} [{ci[0]:.2f}, {ci[1]:.2f}], p={p_pgls:.3g}"
    )
    results.append({"x": xlabel, "y": ylabel, "spearman_rho": rho, "spearman_p": p_raw,
                    "n_species": len(traits), "bm_pgls_beta": beta[1], "bm_pgls_se": se[1],
                    "bm_pgls_p": p_pgls, "ci_low": ci[0], "ci_high": ci[1]})
fig.suptitle("Species-level morphology and relative nuclear-IOD relationships", fontsize=15)
fig.tight_layout()
save_show(fig, "cell_trait_pairwise_bm_pgls_analysis18_v1.png")
display(pd.DataFrame(results))
"""
        ),
        md("## 8. Trait tracks across the phylogeny"),
        code(
            """
from copy import deepcopy
from matplotlib.colors import TwoSlopeNorm

track_tree = deepcopy(tree)
track_tree.ladderize()
depths = track_tree.depths()
terminals = track_tree.get_terminals()
species_order = [tip.name.replace("D.", "") for tip in terminals]
y_position = {tip: len(terminals) - index - 1 for index, tip in enumerate(terminals)}

def node_y(clade):
    if clade in y_position:
        return y_position[clade]
    values = [node_y(child) for child in clade.clades]
    y_position[clade] = float(np.mean(values))
    return y_position[clade]

node_y(track_tree.root)
track_values = traits.loc[species_order, ["relative_nuclear_iod_proxy", "nuc_area_um2", "cell_area_um2"]].copy()
track_values = np.log(track_values)
track_values = (track_values - track_values.mean()) / track_values.std(ddof=1)

fig = plt.figure(figsize=(14, 9))
grid = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=0.04)
tree_ax = fig.add_subplot(grid[0, 0])
heat_ax = fig.add_subplot(grid[0, 1])
for clade in track_tree.find_clades(order="preorder"):
    x = depths[clade]
    y = y_position[clade]
    if clade.clades:
        child_ys = [y_position[child] for child in clade.clades]
        tree_ax.plot([x, x], [min(child_ys), max(child_ys)], color="#44515c", lw=1.2)
        for child in clade.clades:
            tree_ax.plot([x, depths[child]], [y_position[child], y_position[child]], color="#44515c", lw=1.2)
for tip in terminals:
    species_name = tip.name.replace("D.", "")
    tree_ax.text(depths[tip] + 0.12, y_position[tip], species_name, va="center", fontsize=9, style="italic")
tree_ax.set_ylim(-1, len(terminals))
tree_ax.set_xlim(0, max(depths.values()) * 1.25)
tree_ax.set_xlabel("Time from root (Myr)")
tree_ax.set_yticks([])
tree_ax.set_title("Published main time tree")

image = heat_ax.imshow(track_values.values, aspect="auto", cmap="vlag", norm=TwoSlopeNorm(vcenter=0, vmin=-2.2, vmax=2.2))
heat_ax.set_yticks(np.arange(len(species_order)), [""] * len(species_order))
heat_ax.set_xticks([0, 1, 2], ["Relative\\nIOD", "Nucleus\\narea", "Cell\\narea"])
heat_ax.tick_params(length=0)
for i in range(len(species_order)):
    for j in range(3):
        heat_ax.text(j, i, f"{track_values.iloc[i, j]:.1f}", ha="center", va="center", fontsize=7,
                     color="white" if abs(track_values.iloc[i, j]) > 1.1 else "#24313a")
heat_ax.set_title("Standardized ln traits")
fig.colorbar(image, ax=heat_ax, fraction=0.055, pad=0.04, label="Within-trait z score")
fig.suptitle("Relative nuclear-IOD, nucleus area, and cell area across the final-18 phylogeny", fontsize=15)
save_show(fig, "cell_trait_phylogeny_tracks_analysis18_v1.png")
raw_track_table = traits.loc[species_order, ["relative_nuclear_iod_proxy", "nuc_area_um2", "cell_area_um2"]].copy()
raw_track_table = raw_track_table.join(
    support.set_index("species")[["n_specimens", "n_images", "n_frozen_pairs"]]
)
display(raw_track_table)
"""
        ),
        md(
            """
## 9. Interpretation and paper-facing boundary

- **Approved provenance:** 900 frozen linked cell/nucleus pairs; one-to-one and physical-linkage checks pass.
- **Morphology:** call the selected trait an upper-tail or equal-image selected-50 sensitivity, not typical erythrocyte size.
- **Segmentation:** no held-out final-panel species benchmark exists. The archived cell benchmark and current nucleus benchmark are diagnostic, not publication-level generalization evidence.
- **IOD:** exactly equals nuclear pixel area × mean optical density in the stored table. Without documented DNA-stoichiometric staining, same-batch reference standards, and external validation, it is relative nuclear IOD—not pg or C-value.
- **Replication:** three final species have one microscopy image/specimen; those conditional intervals do not estimate intraspecific biological variance.
- **Historical HTML reports:** useful for inspecting masks and pipeline history, but any hard-coded 16.36-pg results are not approved release evidence.

Nature Methods conventions support raw/truth/prediction overlays, object-level metrics, leakage-proof held-out splits, and stratified error analysis (Stringer et al. 2021, DOI `10.1038/s41592-020-01018-x`; Reinke et al. 2024, DOI `10.1038/s41592-023-02150-0`). Comparable salamander cytometry requires same-batch DNA standards and optical-density calibration. See the [literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).

Exact audit: [microscopy_release_audit_analysis18_v1.md](../../plans/publication-readiness-deep-audit/microscopy_release_audit_analysis18_v1.md).
"""
        ),
    ]
    return notebook("Cell segmentation, linked morphometry, and relative nuclear IOD", "cell_models", cells)


def path_notebook():
    cells = [
        md(
            """
## 1. Release boundary

The corrected path branch uses a relative nuclear-IOD proxy, never the historical `genome_size_pg` column. Candidate DAGs were formalized during the audit after historical results existed, so this is exploratory model comparison rather than a prospectively confirmatory causal analysis.
"""
        ),
        code(
            """
gates = read_csv("results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv")
display(gates[["component", "status", "permitted_language"]])
show_figure("results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png")
"""
        ),
        md("## 2. Frozen input cube and transformations"),
        code(
            """
inputs = read_csv("results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv")
specifications = read_csv("results/data/corrected/path_analysis/corrected_path_measurement_specifications_analysis18_v1.csv")
assert inputs.shape[0] == 18 * 6 * 3
assert set(inputs.species) == set(FINAL18)
assert "genome_size_pg" not in inputs.columns
display(specifications)
display(inputs.head(18))
"""
        ),
        md(
            """
## 3. Candidate-model comparison and basis sets

> **Shared-measurement gate:** relative nuclear IOD is nuclear pixel area × mean optical density, so an IOD–nucleus-area edge is mechanically coupled. Likewise, LTR:LINE log-ratio and Pielou evenness are both functions of the same closed TE composition. CICc support for those edges describes coherence among derived variables; it is not independent evidence for either causal direction.
"""
        ),
        code(
            """
anchor = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv")
anchor["fisher_C_df"] = 2 * anchor.q
edges = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv")
basis = read_csv("results/data/corrected/path_analysis/corrected_path_data_sensitivity_basis_sets_analysis18_v1.csv")
anchor_basis = basis[(basis.tree_id == "published_main") &
                     (basis.morphology_estimator == "image_balanced_selected50") &
                     (basis.iod_subset == "image_qc_pass")]
display(anchor[["family", "n_species", "rank", "model", "q", "fisher_C_df", "C", "p", "CICc", "delta_CICc", "w", "release_gate_status"]])
display(anchor_basis)
display(edges)
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.png")
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_dag_analysis18_v1.png")
"""
        ),
        md("## 4. Measurement, tree, and species-influence sensitivity"),
        code(
            """
stability = read_csv("results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv")
tree_edges = read_csv("results/data/corrected/path_analysis/corrected_path_tree_edge_uncertainty_analysis18_v1.csv")
loo = read_csv("results/data/corrected/path_analysis/corrected_path_leave_one_out_rankings_analysis18_v1.csv")
display(stability)
display(tree_edges)
display(loo[loo["rank"] == 1].groupby(["family", "model"]).agg(
    n_omissions=("omitted_species", "nunique"), min_global_p=("p", "min"),
    min_weight=("w", "min"), max_weight=("w", "max")))
for figure in [
    "results/figures/corrected/path_analysis/corrected_path_measurement_model_weights_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_tree_uncertainty_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_leave_one_out_influence_analysis18_v1.png",
]:
    show_figure(figure)
"""
        ),
        md("## 5. Actual-tree selection calibration"),
        code(
            """
simulation = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv")
edge_calibration = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_edge_calibration_analysis18_v1.csv")
assert simulation.n_failures.sum() == 0
assert simulation.n_attempted.sum() == 1500
display(simulation)
display(edge_calibration)
show_figure("results/figures/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.png")
"""
        ),
        md("## 6. Compact evidence summary"),
        code(
            """
anchor_top = anchor[anchor["rank"] == 1].set_index("family")
null_sim = simulation[simulation.scenario == "independent_null"].set_index("family")
observed_sim = simulation[(simulation.scenario == "observed_chain") & (simulation.effect_scale == 1.0)].set_index("family")
summary = pd.DataFrame({
    "anchor_top_model": anchor_top.model,
    "anchor_global_p": anchor_top.p,
    "anchor_weight": anchor_top.w,
    "strict_unique_false_selection": null_sim.false_unique_nonnull_rate,
    "observed_scale_exact_recovery": observed_sim.expected_unique_supported_rate,
})
display(summary)
"""
        ),
        md(
            """
## 7. Interpretation and paper-facing boundary

- 958 corrected fits completed with no fitting failure: 84 data/tree specifications, 804 published-tree fits, and 70 species omissions.
- The TE-evenness chain and IOD–nucleus–cell chain are empirically stable in these data.
- Stability is not causal identification. The integrated generating model is recovered uniquely in only 50% of observed-scale simulations at `n=18`.
- Strict unique false-selection rates under independent Brownian traits are 19% for TE–IOD, 20% for IOD–morphology, and 9.5% for the integrated family.
- The `evenness → IOD` approximate interval crosses zero.
- The terminal:internal family is globally rejected in every full-panel specification.
- IOD algebraically contains nuclear area, so the IOD→nucleus path is mechanically coupled and cannot validate a nucleotypic mechanism.
- Measurement error is explored through estimator/subset sensitivity but not jointly propagated in an error-aware hierarchical SEM.

Leading phylogenetic path papers expect a small a priori DAG set, printed basis sets, global fit, CICc/model weights, edge uncertainty, tree uncertainty, influence analysis, and actual-tree calibration. Those computational surfaces are present; the upstream trait-validity and prospective-model gates are not. See van der Bijl (2018, DOI `10.7717/peerj.4718`) and the [literature benchmark](../../plans/publication-readiness-deep-audit/notebook_figure_primary_literature_benchmark.md).

Exact audit: [corrected_path_analysis_audit_analysis18_v1.md](../../plans/publication-readiness-deep-audit/corrected_path_analysis_audit_analysis18_v1.md).
"""
        ),
    ]
    return notebook("Corrected phylogenetic path-analysis audit", "path_analysis", cells)


def descriptive_panel_notebook(title: str, analysis_id: str, domain: str, panel: str, artifacts: list[str], figure: str, narrative: str):
    """Build a lightweight, executed-review notebook for a non-path panel."""

    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.8"},
        "desmognathus_analysis_id": analysis_id,
        "presentation_domain": domain,
        "expensive_upstream_tools_executed": False,
    }
    nb["cells"] = [
        md(f"""
# {title}

**Independent descriptive panel · `{analysis_id}`**

This notebook is deliberately outside the integrated path18 scope. It loads
frozen downstream products only; it does not run repeat annotation, dnaPipeTE,
segmentation, image inference, or phylogenetic models.
"""),
        setup_cell(domain),
        md("## 1. Panel contract and frozen inputs\n\n" + narrative),
        code(f"""
artifacts = {artifacts!r}
display(artifact_table(artifacts))
panel = read_csv({panel!r})
print(f"Declared species: {{panel.species.nunique()}}")
display(panel)
"""),
        md("## 2. Derived data and source manifest"),
        code(f"""
manifest_path = ROOT / {artifacts[-1]!r}
manifest = json.loads(manifest_path.read_text())
display(pd.DataFrame([manifest]))
tables = [path for path in artifacts if path.endswith('.csv') and 'panel' not in path]
for relative in tables:
    table = read_csv(relative)
    print(f"{{relative}}: {{table.shape[0]}} rows, {{table.shape[1]}} columns")
    display(table.head(8))
"""),
        md("## 3. Presentation figure"),
        code(f"""
show_figure({figure!r})
"""),
        md("## 4. Scope boundary"),
        code("""
study_manifest = json.loads((ROOT / "path_analysis/data/derived/panels/study_species_panels_v1.manifest.json").read_text())
display(pd.DataFrame([{
    "TE resource panel": study_manifest["te_resource_panel"]["n_species"],
    "linked-cell panel": study_manifest["cell_linked_panel"]["n_species"],
    "integrated path panel": study_manifest["integrated_path_panel"]["n_species"],
    "microscopy traits imputed": study_manifest["microscopy_traits_imputed_for_te_only_species"],
}]))
"""),
    ]
    return nb


def te34_notebook():
    return descriptive_panel_notebook(
        "TE34 repeat landscapes, diversity, and compositional PCA",
        "te34_descriptive_bundle_v1",
        "repeat_analysis_te34",
        "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
        [
            "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
            "results/data/corrected/te34/dnapipete_mass_accounting_te34_v1.csv",
            "results/data/corrected/te34/te_diversity_mass_sensitivity_te34_v1.csv",
            "results/data/corrected/te34/te_pca_scores_te34_v1.csv",
            "results/data/corrected/te34/repeatmasker_hit_inventory_te34_v1.csv",
            "results/data/corrected/te34/repeatmasker_divergence_landscape_te34_v1.csv",
            "results/data/corrected/te34/te34_descriptive_bundle_v1.manifest.json",
        ],
        "results/figures/corrected/te34/te34_order_clr_pca_v1.png",
        "TE34 is the complete active vetted genomic-resource panel. Its repeat abundance is normalized within the stated data stream; no Cell21 trait is required or imputed.",
    )


def te34_replicate_averaged_notebook():
    return descriptive_panel_notebook(
        "TE34 repeat landscapes, diversity, and PCA: same-species runs averaged",
        "te34_replicate_averaged_v1",
        "repeat_analysis_te34_replicate_averaged",
        "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
        [
            "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
            "results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv",
            "results/data/corrected/te34_replicate_averaged/te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv",
            "results/data/corrected/te34_replicate_averaged/te_pca_scores_te34_replicate_averaged_v1.csv",
            "results/data/corrected/te34_replicate_averaged/repeatmasker_hit_inventory_te34_replicate_averaged_v1.csv",
            "results/data/corrected/te34_replicate_averaged/repeatmasker_divergence_landscape_te34_replicate_averaged_v1.csv",
            "results/data/corrected/te34_replicate_averaged/te34_replicate_averaged_v1.manifest.json",
        ],
        "results/figures/corrected/te34_replicate_averaged/te34_replicate_averaged_order_clr_pca_v1.png",
        "Same-species dnaPipeTE retries are equal-weight averaged at the run level before species-level composition, diversity, and PCA. They never become extra species or doubled sequence mass. This is the preferred TE34 presentation branch.",
    )


def cell21_notebook():
    return descriptive_panel_notebook(
        "Cell21 linked morphometry and relative nuclear IOD",
        "cell21_descriptive_bundle_v1",
        "cell_models_cell21",
        "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
        [
            "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
            "results/data/corrected/cell21/cell_linked_traits_cell21_v1.csv",
            "results/data/corrected/cell21/cell21_descriptive_bundle_v1.manifest.json",
        ],
        "results/figures/corrected/cell21/cell21_cell_nucleus_relative_iod_v1.png",
        "Cell21 is the current linked cell/nucleus and image-IOD panel. Relative nuclear IOD is an image phenotype, not absolute genome size; the 18-species path analysis uses only the evidence-complete overlap.",
    )


def legacy_descriptive_te_notebook():
    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.8"},
        "desmognathus_analysis_id": "legacy_descriptive_figure_bundle_v1",
        "presentation_domain": "legacy_descriptive_te",
        "expensive_upstream_tools_executed": False,
    }
    figures = [
        "results/figures/legacy_descriptive/dnapipete_input_quality_te34_legacy_descriptive_v1.png",
        "results/figures/talk_iv_restored_composition/Order_average.png",
        "results/figures/talk_iv_restored_composition/Superfamily_average.png",
        "results/figures/talk_iv_restored_pca/pca_scree_plot.png",
        "results/figures/talk_iv_restored_pca/pca_by_species.png",
        "results/figures/talk_iv_restored_pca/cluster_elbow_plot.png",
        "results/figures/talk_iv_restored_pca/cluster_silhouette_plot.png",
        "results/figures/talk_iv_restored_pca/pca_cluster_plot.png",
        "results/figures/legacy_descriptive/te_diversity_indices_te34_legacy_descriptive_v1.png",
        "results/figures/talk_iv_restored_ectopic/elements_5plus_domains.png",
    ]
    artifacts = [
        "results/data/legacy_descriptive/legacy_descriptive_figure_bundle_v1.manifest.json",
        "results/data/legacy_descriptive/dnapipete_run_quality_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/dnapipete_species_quality_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/te_superfamily_pca_variance_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/te_superfamily_pca_scores_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/te_superfamily_cluster_metrics_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/te_diversity_indices_te34_legacy_descriptive_v1.csv",
        "results/data/legacy_descriptive/ectopic_anova_te34_legacy_descriptive_v1.txt",
    ]
    nb["cells"] = [
        md("""
# Restored TE34 figures, clustering, and LTR statistics

**Lab-meeting figure suite · `legacy_descriptive_figure_bundle_v1`**

This notebook restores the high-information descriptive plots from the earlier
Research Talk IV workflow: dnaPipeTE input quality, TE composition, PCA scree
and clustering diagnostics, diversity, and LTR terminal:internal figures with
their statistical outputs. It reads frozen local products only.
"""),
        setup_cell("legacy_descriptive_te"),
        md("## 1. Frozen figure bundle and source tables"),
        code(f"""
artifacts = {artifacts!r}
display(artifact_table(artifacts))
manifest = json.loads((ROOT / artifacts[0]).read_text())
display(pd.DataFrame([manifest]))
"""),
        md("## 2. Sequencing-input and component quality"),
        code(f"""
run_quality = read_csv("results/data/legacy_descriptive/dnapipete_run_quality_te34_legacy_descriptive_v1.csv")
species_quality = read_csv("results/data/legacy_descriptive/dnapipete_species_quality_te34_legacy_descriptive_v1.csv")
display(run_quality)
display(species_quality)
show_figure({figures[0]!r})
"""),
        md("## 3. Original Talk IV PCA, elbow selection, and clustering"),
        code(f"""
variance = read_csv("results/figures/talk_iv_restored_pca/pca_variance_summary.csv")
scores = read_csv("results/figures/talk_iv_restored_pca/pca_scores.csv")
cluster_metrics = read_csv("results/figures/talk_iv_restored_pca/clustering_metrics.csv")
display(variance)
display(cluster_metrics)
for relative in {figures[1:8]!r}:
    show_figure(relative)
"""),
        md("## 4. TE diversity and LTR terminal:internal statistics"),
        code(f"""
diversity = read_csv("results/data/legacy_descriptive/te_diversity_indices_te34_legacy_descriptive_v1.csv")
display(diversity)
show_figure({figures[8]!r})
anova_text = (ROOT / "results/data/legacy_descriptive/talk_iv_ectopic_stats/anova_results.txt").read_text()
display(Markdown("### One-way ANOVA and Kruskal-Wallis outputs\\n```text\\n" + anova_text + "\\n```"))
show_figure({figures[9]!r})
"""),
        md("""
## 5. How this connects to the corrected notebooks

This restores the original descriptive communication layer. The TE34
replicate-averaged notebook records its active-resource and retry policy; the
corrected LTR notebook provides element-level coverage and influence checks;
the path18 notebook remains the only combined TE/cell comparative analysis.
"""),
    ]
    return nb


NOTEBOOKS = {
    "01_phylogeny_tree_trimming_analysis18_v1.ipynb": phylogeny_notebook,
    "02_data_tables_and_provenance_analysis18_v1.ipynb": data_tables_notebook,
    "03_repeat_landscapes_diversity_pca_analysis18_v1.ipynb": repeat_notebook,
    "04_ltr_deletion_proxy_analysis18_v1.ipynb": ectopic_notebook,
    "05_cell_models_morphometry_iod_analysis18_v1.ipynb": cell_notebook,
    "06_phylogenetic_path_analysis18_v1.ipynb": path_notebook,
    "03_repeat_landscapes_diversity_pca_te34_v1.ipynb": te34_notebook,
    "03_repeat_landscapes_diversity_pca_te34_replicate_averaged_v1.ipynb": te34_replicate_averaged_notebook,
    "05_cell_models_morphometry_iod_cell21_v1.ipynb": cell21_notebook,
    "03_te34_legacy_descriptive_figures_v1.ipynb": legacy_descriptive_te_notebook,
}


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for filename, builder in NOTEBOOKS.items():
        output = OUTPUT_DIR / filename
        built = builder()
        nbf.write(built, output)
        print(f"Wrote {output.relative_to(ROOT)} with {len(built['cells'])} cells")


if __name__ == "__main__":
    main()

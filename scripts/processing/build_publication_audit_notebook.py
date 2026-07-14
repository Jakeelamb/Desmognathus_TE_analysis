#!/usr/bin/env python3
"""Build the portable final-18 publication-audit review notebook."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def build_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.8"},
        "desmognathus_analysis_id": "publication_audit_analysis18_v1",
    }
    nb["cells"] = [
        markdown(
            """
# *Desmognathus* TE–microscopy–phylogeny publication audit

Final 18-species review workbench, analysis version `analysis18_v1`.

This notebook compiles the corrected, non-destructive audit products. It is designed for collaborators to inspect species-level anomalies, change filters, and trace every displayed result to an exact CSV and source script. Historical outputs are preserved outside `results/data/corrected/`.

> **Current release boundary:** TE diversity/CLR ordination and time-tree sensitivity are suitable for descriptive use. The terminal:internal metric is a deletion-footprint proxy, not an ectopic-recombination rate. Morphometry is upper-tail sensitivity with incomplete final-panel segmentation validation. Relative nuclear IOD is **not absolute genome size**. Therefore the path models remain exploratory.
"""
        ),
        markdown("## 0. Setup and immutable final panel"),
        code(
            """
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import Image, Markdown, display

ROOT = Path.cwd().resolve()
if not (ROOT / "results/data/corrected").exists():
    ROOT = ROOT.parent
assert (ROOT / "results/data/corrected").exists(), "Run from the repository root or notebooks/"

CORRECTED = ROOT / "results/data/corrected"
FIGURES = ROOT / "results/figures/corrected"
AUDIT = ROOT / "plans/publication-readiness-deep-audit"
pd.set_option("display.max_columns", 60)
pd.set_option("display.max_rows", 100)
sns.set_theme(style="whitegrid")

FINAL18 = [
    "amphileucus", "anicetus", "apalachicolae", "auriculatus", "bairdi", "campi",
    "fuscus", "gvnigeusgwotli", "intermedius", "kanawha", "marmoratus", "mavrokoilius",
    "monticola", "ocoee", "perlapsus", "tilleyi", "valtos", "welteri",
]

def read_csv(relative):
    path = ROOT / relative
    assert path.exists(), path
    return pd.read_csv(path)

def show_figure(relative, width=1050):
    path = ROOT / relative
    assert path.exists(), path
    display(Image(filename=str(path), width=width))

path_input = read_csv("results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv")
assert sorted(path_input.species.unique()) == sorted(FINAL18)
assert path_input.shape[0] == 18 * 6 * 3
print(f"Repository root: {ROOT}")
print(f"Final panel: {path_input.species.nunique()} species; {path_input.shape[0]} measurement-specification rows")
"""
        ),
        markdown(
            """
Species identity is trait-specific. The genomic `fuscus` row uses the current short-read comparative resource (`SRX20497025 / GCA_032353935.1`, publicly labeled *planiceps* but expert-reidentified as *fuscus*). The separately published chromosome-level `GCA_050004315.1` resource is validation-only and was not used for the current analyses. Microscopy specimens are independent biological samples and are not accession substitutions. dnaPipeTE directory suffixes such as `2` or `3` are computational retries, not additional species or replicates.
"""
        ),
        markdown("## 1. Claim and release gates"),
        code(
            """
gates = read_csv("results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv")
display(gates[["component", "status", "permitted_language"]])
show_figure("results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png")
"""
        ),
        markdown(
            "Full rationale: [`corrected_path_analysis_audit_analysis18_v1.md`](../plans/publication-readiness-deep-audit/corrected_path_analysis_audit_analysis18_v1.md) and [`upstream_methods_primary_literature_benchmark.md`](../plans/publication-readiness-deep-audit/upstream_methods_primary_literature_benchmark.md)."
        ),
        markdown("## 2. dnaPipeTE mass accounting and TE composition"),
        code(
            """
mass = read_csv("results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv")
display(mass.sort_values("species"))

composition = read_csv("results/data/corrected/dnapipete/dnaPipeTE_order_breakdown_mass_accounted_analysis18_v2.csv")
te_orders = ["DIRS", "Helitron", "LINE", "LTR", "Maverick", "PLE", "SINE", "TIR", "Transposon Derivatives", "YR"]
plot_frame = composition.set_index("species")[te_orders].div(100)
plot_frame = plot_frame.mul(composition.set_index("species")["order_retained_fraction"], axis=0)
plot_frame["Unresolved"] = composition.set_index("species")["order_unresolved_fraction"]
assert np.allclose(plot_frame.sum(axis=1), 1)
plot_frame = plot_frame.loc[[s for s in FINAL18 if s in plot_frame.index]]
ax = plot_frame.plot(kind="bar", stacked=True, figsize=(15, 6), width=0.85, colormap="tab20")
ax.set_ylabel("Mass-accounted fraction")
ax.set_xlabel("")
ax.set_title("dnaPipeTE order composition with unresolved mass retained")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
plt.tight_layout()
plt.show()
"""
        ),
        markdown(
            "Audit sources: `scripts/processing/audit_dnapipete_mass_accounting.py`, `scripts/processing/build_dnapipete_absolute_load_sensitivity.py`, and the retry audit report. Retry folders resolve to one canonical species resource before analysis."
        ),
        markdown("## 3. LTR terminal:internal deletion-footprint proxy"),
        code(
            """
ectopic = read_csv("results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv")
display(ectopic.sort_values(["analysis_branch", "ratio_median"], ascending=[True, False]))
show_figure("results/figures/corrected/ectopic/ectopic_element_log2_ratio_by_species_analysis18_v1.png")
show_figure("results/figures/corrected/ectopic/ectopic_mean_influence_analysis18_v1.png")
"""
        ),
        markdown(
            """
Interpretation is deliberately narrow: exact ≥5/6-domain LTR elements, coverage QC, zero retention, robust summaries, and mean-influence diagnostics are reproducible. The ratio does not count orthologous solo-LTR/full-length events and lacks two assemblies, so it is not an ectopic-recombination rate. The corrected path family using this proxy fails global fit in every full-panel specification.

Source: `scripts/processing/build_corrected_ectopic_recombination.py`.
"""
        ),
        markdown("## 4. TE diversity and compositional ordination"),
        code(
            """
diversity = read_csv("results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv")
order_primary = diversity[(diversity.te_level == "order") & (diversity.composition_mode == "classified_conditional")]
display(order_primary[["species", "observed_richness", "shannon_entropy", "gini_simpson", "hill_q1", "hill_q2", "pielou_evenness"]].sort_values("pielou_evenness"))

variance = read_csv("results/data/corrected/diversity_pca/te_pca_variance_analysis18_v1.csv")
display(variance.head(10))
show_figure("results/figures/corrected/diversity_pca/te_order_clr_pca_scores_analysis18_v1.png")
show_figure("results/figures/corrected/diversity_pca/te_order_clr_pca_loadings_analysis18_v1.png")
show_figure("results/figures/corrected/diversity_pca/te_diversity_unresolved_mass_sensitivity_analysis18_v1.png")
"""
        ),
        code(
            """
ppca = read_csv("results/data/corrected/diversity_pca/te_order_phylogenetic_pca_comparison_analysis18_v1.csv")
tree_pca = read_csv("results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_metrics_analysis18_v1.csv")
display(ppca)
display(tree_pca.describe(include="all").T)
show_figure("results/figures/corrected/diversity_pca/te_order_ordinary_vs_phylogenetic_pca_analysis18_v1.png")
show_figure("results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_loading_intervals_analysis18_v1.png")
"""
        ),
        markdown(
            """
Approved scope: descriptive Shannon, Gini-Simpson, Hill, richness/evenness, and CLR ordination at declared taxonomic levels. The historical “Simpson” ambiguity is corrected. The order PCA has explicit closure and CLR handling, full scores/loadings/scree, leave-one-out stability, phylogenetic PCA, and 200-tree sensitivity. An unrestricted PERMANOVA is not promoted because no defensible replicated grouping factor was prespecified.

Sources: `scripts/processing/build_corrected_te_diversity_pca.py`, `audit_order_phylogenetic_pca.R`, and `audit_order_pca_tree_uncertainty.R`.
"""
        ),
        markdown("## 5. Phylogeny provenance and uncertainty"),
        code(
            """
tree_metrics = read_csv("results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv")
topology = read_csv("results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv")
crosswalk = read_csv("results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv")
display(tree_metrics)
display(topology)
display(crosswalk)
show_figure("results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png")
show_figure("results/figures/corrected/phylogeny/phylogeny_data_completeness_analysis18_v1.png")
show_figure("results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png")
"""
        ),
        markdown(
            """
The exact final-18 focal tree is structurally valid after correcting terminal rounding. Comparative sensitivity uses the published Stewart–Wiens main tree plus all 200 supplied time trees. Published-tree provenance is strong; the focal tree's exact archive/citation and calibration history remain unresolved and must be supplied before treating it as a primary publication tree.

Sources: `scripts/processing/audit_phylogeny_release.py` and `audit_published_tree_uncertainty.R`.
"""
        ),
        markdown("## 6. Microscopy estimand, segmentation, and relative IOD"),
        code(
            """
support = read_csv("results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv")
estimands = read_csv("results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv")
rank_stability = read_csv("results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv")
display(support.sort_values(["n_specimens", "n_images", "species"]))
display(rank_stability.sort_values(["metric", "spearman_rho"], ascending=[True, False]))
show_figure("results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png")
"""
        ),
        code(
            """
validation = read_csv("results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv")
iod = read_csv("results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv")
image_qc = read_csv("results/data/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.csv")
display(validation)
display(iod.sort_values("relative_iod_to_fuscus_image_qc_pass", ascending=False))
show_figure("results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png")
"""
        ),
        markdown(
            """
The frozen data preserve 900 one-cell/one-nucleus pairs. They are not 900 independent biological replicates: cells are nested in images/specimens/species. The largest-50 workflow estimates an upper tail; it is not “least biased” typical size. Only 324/900 frozen pairs are explicit manual keeps, and some species have one image/specimen.

Absolute genome size is blocked. The stored IOD is exactly `area_px × mean_OD`; no same-batch Feulgen/reference calibration is documented, and only one *fuscus* image passes the image-QC subset. Use `relative_nuclear_iod_proxy` only.

Source: `scripts/processing/audit_microscopy_release.py`.
"""
        ),
        markdown("## 7. Corrected phylogenetic path models"),
        code(
            """
anchor = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv")
anchor_edges = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv")
stability = read_csv("results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv")
basis = read_csv("results/data/corrected/path_analysis/corrected_path_data_sensitivity_basis_sets_analysis18_v1.csv")

display(anchor[["family", "rank", "model", "q", "C", "p", "CICc", "delta_CICc", "w", "release_gate_status"]])
display(anchor_edges)
display(stability)
display(basis[(basis.tree_id == "published_main") & (basis.morphology_estimator == "image_balanced_selected50") & (basis.iod_subset == "image_qc_pass")])
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.png")
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_dag_analysis18_v1.png")
"""
        ),
        code(
            """
tree_rankings = read_csv("results/data/corrected/path_analysis/corrected_path_tree_sensitivity_rankings_analysis18_v1.csv")
loo_rankings = read_csv("results/data/corrected/path_analysis/corrected_path_leave_one_out_rankings_analysis18_v1.csv")

tree_top = tree_rankings[tree_rankings["rank"] == 1]
loo_top = loo_rankings[loo_rankings["rank"] == 1]
display(tree_top.groupby(["family", "model"]).agg(n_trees=("tree_id", "nunique"), min_weight=("w", "min"), max_weight=("w", "max"), min_global_p=("p", "min")))
display(loo_top.groupby(["family", "model"]).agg(n_omissions=("omitted_species", "nunique"), min_weight=("w", "min"), max_weight=("w", "max"), min_global_p=("p", "min")))
show_figure("results/figures/corrected/path_analysis/corrected_path_measurement_model_weights_analysis18_v1.png")
show_figure("results/figures/corrected/path_analysis/corrected_path_tree_uncertainty_analysis18_v1.png")
show_figure("results/figures/corrected/path_analysis/corrected_path_leave_one_out_influence_analysis18_v1.png")
"""
        ),
        markdown("## 8. Actual-tree model-selection and interval calibration"),
        code(
            """
simulation = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv")
edge_calibration = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_edge_calibration_analysis18_v1.csv")
display(simulation)
display(edge_calibration)
show_figure("results/figures/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.png")
"""
        ),
        markdown(
            """
Simulation uses the actual published final-18 tree, 200 independent-Brownian null replicates per family, and 100 signal replicates per family at 0.5×/1×/1.5× observed coefficients. It reports false supported non-null selection, exact unique model recovery, coefficient bias, RMSE, and approximate interval coverage. This diagnoses information content under the specified generator; it cannot validate the real causal direction.

Sources: `scripts/processing/audit_corrected_path_models.R`, `simulate_corrected_path_calibration.R`, and `summarize_corrected_path_audit.py`.
"""
        ),
        markdown("## 9. Automatically generated weird-result queue"),
        code(
            """
issues = []

low_support = support[(support.n_specimens < 2) | (support.n_images < 2)].copy()
if len(low_support):
    low_support.insert(0, "issue", "<2 morphology specimens or images")
    issues.append(low_support[["issue", "species", "n_images", "n_specimens"]])

missing_ectopic = path_input[path_input.terminal_internal_ratio_median.isna()][["species"]].drop_duplicates()
if len(missing_ectopic):
    missing_ectopic.insert(0, "issue", "no terminal:internal assembly proxy")
    issues.append(missing_ectopic)

low_manual = support.nsmallest(8, "manual_keep_fraction").copy()
low_manual.insert(0, "issue", "lowest explicit manual-keep fractions")
issues.append(low_manual[["issue", "species", "manual_keep_fraction", "n_manual_keep", "n_model_ranked_unlabeled"]])

uncertain_edges = anchor_edges[~anchor_edges.interval_excludes_zero].copy()
if len(uncertain_edges):
    uncertain_edges.insert(0, "issue", "anchor edge interval includes zero")
    issues.append(uncertain_edges[["issue", "family", "parent", "child", "coefficient", "approx_ci_low", "approx_ci_high"]])

terminal_rejection = anchor[(anchor.family == "terminal_internal_iod") & (anchor["rank"] == 1)].copy()
terminal_rejection.insert(0, "issue", "terminal:internal family top model rejected globally")
issues.append(terminal_rejection[["issue", "model", "p", "w", "release_gate_status"]])

for issue_table in issues:
    display(issue_table)
"""
        ),
        markdown(
            """
This queue is intentionally redundant with the figures: it makes the most consequential rows easy to discuss together. Add collaborator notes in a copy of this notebook; keep corrected source CSVs immutable and regenerate them through the listed scripts.
"""
        ),
        markdown("## 10. Artifact and reproducibility index"),
        code(
            """
manifests = sorted((ROOT / "results/data/corrected").rglob("*.manifest.json"))
manifest_index = []
for path in manifests:
    payload = json.loads(path.read_text())
    manifest_index.append({
        "manifest": str(path.relative_to(ROOT)),
        "analysis_id": payload.get("analysis_id", ""),
        "absolute_genome_size_used": payload.get("absolute_genome_size_used", ""),
        "publication_or_causal_claim": payload.get("publication_claim_allowed", payload.get("causal_claim_allowed", "")),
    })
display(pd.DataFrame(manifest_index))

reports = sorted(AUDIT.glob("*.md"))
display(pd.DataFrame({"audit_report": [str(path.relative_to(ROOT)) for path in reports]}))
"""
        ),
        markdown(
            """
### Rebuild order

1. `scripts/run_in_dusky.sh python scripts/processing/build_corrected_ectopic_recombination.py`
2. `scripts/run_in_dusky.sh python scripts/processing/build_corrected_te_diversity_pca.py`
3. Run the phylogenetic PCA/tree-sensitivity R scripts and `audit_phylogeny_release.py`.
4. `scripts/run_in_dusky.sh python scripts/processing/audit_microscopy_release.py`
5. `scripts/run_in_dusky.sh python scripts/processing/build_corrected_path_inputs.py`
6. `scripts/run_in_dusky.sh Rscript scripts/processing/audit_corrected_path_models.R --phase all`
7. `scripts/run_in_dusky.sh Rscript scripts/processing/simulate_corrected_path_calibration.R`
8. `scripts/run_in_dusky.sh python scripts/processing/summarize_corrected_path_audit.py`
9. Rebuild and execute this notebook, then run `scripts/processing/audit_publication_notebook.py` and `scripts/run_tests.sh`.

The publication-facing analysis should cite exact corrected artifacts and report versions, not historical similarly named tables.
"""
        ),
    ]
    return nb


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    nb = build_notebook()
    nbf.write(nb, OUTPUT)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} with {len(nb['cells'])} cells")


if __name__ == "__main__":
    main()

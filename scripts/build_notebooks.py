#!/usr/bin/env python3
"""Generate the three small, editable data-review notebooks."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]


def notebook(title: str, boundary: str, cells: list[str]) -> nbf.NotebookNode:
    book = nbf.v4.new_notebook()
    book["metadata"] = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3"},
    }
    book["cells"] = [
        nbf.v4.new_markdown_cell(
            f"# {title}\n\n{boundary}\n\n"
            "This notebook reads only the compact canonical tables in this analysis "
            "directory. Change or add cells freely; the paper figures are built "
            "separately in R/ggplot."
        ),
        nbf.v4.new_code_cell(
            "from pathlib import Path\n"
            "import pandas as pd\n\n"
            "ROOT = next(\n"
            "    candidate for candidate in (Path.cwd(), *Path.cwd().parents)\n"
            "    if (candidate / 'pyproject.toml').exists()\n"
            ")\n"
            "pd.set_option('display.max_columns', 100)"
        ),
        *[nbf.v4.new_code_cell(source) for source in cells],
    ]
    return book


NOTEBOOKS = {
    "analyses/01_transposable_elements/explore.ipynb": notebook(
        "Transposable-element data",
        "TE34 is the accession-linked genomic-resource panel. The LTR "
        "terminal:internal statistic is a mapping/deletion-footprint proxy, not a "
        "direct ectopic-recombination or DNA-loss rate.",
        [
            (
                "import json\n\n"
                "DATA = ROOT / 'analyses/01_transposable_elements/data'\n"
                "panel = pd.read_csv(ROOT / 'data/identity/te34_panel.csv')\n"
                "diversity = pd.read_csv(DATA / 'te_diversity.csv')\n"
                "prevalence = pd.read_csv(DATA / 'te_feature_prevalence.csv')\n"
                "clr = pd.read_csv(DATA / 'te_pca_clr_matrix.csv')\n"
                "pca = pd.read_csv(DATA / 'te_pca_scores.csv')\n"
                "variance = pd.read_csv(DATA / 'te_pca_variance.csv')\n"
                "loadings = pd.read_csv(DATA / 'te_pca_loadings.csv')\n"
                "diagnostics = pd.read_csv(DATA / 'te_pca_clustering_diagnostics.csv')\n"
                "stability = pd.read_csv(DATA / 'te_pca_clustering_stability.csv')\n"
                "candidates = pd.read_csv(DATA / 'te_pca_candidate_cluster_assignments.csv')\n"
                "trait_tests = pd.read_csv(DATA / 'te_pca_trait_association_tests.csv')\n"
                "phylogeny_tests = pd.read_csv(DATA / 'te_pca_phylogenetic_signal_tests.csv')\n"
                "manifest = json.loads((DATA / 'te_pca_clustering_analysis_manifest.json').read_text())\n"
                "landscape = pd.read_csv(DATA / 'repeatmasker_divergence_landscape.csv')\n"
                "ltr = pd.read_csv(DATA / 'ltr_element_metrics.csv')\n\n"
                "assert (panel['species'].nunique(), len(pca), len(variance), len(loadings)) "
                "== (34, 34, 23, 24)\n"
                "assert manifest['selected_k'] is None\n"
                "assert set(candidates['candidate_status']) == {'rejected'}\n"
                "{'TE34 species': panel['species'].nunique(), "
                "'PCA features': len(loadings), "
                "'PCA axes': len(variance), "
                "'landscape rows': len(landscape), 'usable LTR elements': len(ltr)}"
            ),
            (
                "diversity.sort_values(['te_level', 'shannon_entropy'], "
                "ascending=[True, False]).head(20)"
            ),
            (
                "display(\n"
                "    pca[['species', 'PC1', 'PC2']].sort_values('PC1'),\n"
                "    variance[['PC', 'variance_explained', "
                "'cumulative_variance_explained']].head(10),\n"
                "    loadings[['superfamily', 'PC1', 'PC2']].sort_values('PC1'),\n"
                ")"
            ),
            (
                "pd.Series({\n"
                "    'analysis_id': manifest['analysis_id'],\n"
                "    'analysis_status': manifest['analysis_status'],\n"
                "    'species': manifest['pca_geometry']['species'],\n"
                "    'ubiquitous_superfamilies': "
                "manifest['pca_geometry']['ubiquitous_superfamilies'],\n"
                "    'PC1_PC2_variance': manifest['pca_geometry']['pc1_pc2_variance'],\n"
                "    'gap_selected_k': "
                "manifest['clustering']['gap_primary_rule']['selected_k'],\n"
                "    'accepted_k': manifest['clustering']['accepted_k'],\n"
                "    'conclusion': manifest['clustering']['conclusion'],\n"
                "    'tree_provenance': "
                "manifest['phylogenetic_signal']['provenance_status'],\n"
                "}, name='final PCA audit')"
            ),
            (
                "display(\n"
                "    diagnostics.loc[diagnostics['representation'].eq('all_23_pcs'), [\n"
                "        'k', 'mean_silhouette', 'minimum_cluster_size', "
                "'singleton_count',\n"
                "        'cluster_sizes_ascending', 'gap_supported', "
                "'passes_feature_stability',\n"
                "        'passes_species_stability', 'passes_six_pc_sensitivity',\n"
                "        'accepted_cluster_solution',\n"
                "    ]],\n"
                "    stability[['k', 'stability_type', 'perturbation', "
                "'n_replicates',\n"
                "        'ari_q10', 'ari_median', 'ari_q90', "
                "'fraction_ari_at_least_0_80']],\n"
                ")"
            ),
            (
                "trait_tests[['test_id', 'analysis_role', 'predictor', 'n_species',\n"
                "    'statistic_name', 'statistic', 'variance_explained_r2',\n"
                "    'p_value', 'adjusted_p_value', 'interpretation_boundary']]"
            ),
            (
                "phylogeny_tests[['test_id', 'analysis_role', 'representation',\n"
                "    'statistic_name', 'statistic', 'p_value', 'branch_length_use',\n"
                "    'interpretation_boundary']]"
            ),
        ],
    ),
    "analyses/02_morphology/explore.ipynb": notebook(
        "Cell, nucleus, and nuclear-IOD data",
        "Path24 contains finalized reviewed image objects. Relative nuclear IOD is "
        "an image-derived phenotype, not an independently validated absolute genome "
        "size. Genomic and microscopy records are joined only at species level.",
        [
            (
                "DATA = ROOT / 'analyses/02_morphology/data'\n"
                "panel = pd.read_csv(ROOT / 'data/identity/path24_panel.csv')\n"
                "objects = pd.read_csv(DATA / 'cell_nucleus_objects.csv', low_memory=False)\n"
                "species = pd.read_csv(DATA / 'cell_nucleus_species_estimates.csv')\n"
                "iod_objects = pd.read_csv(DATA / 'nuclear_iod_objects.csv', low_memory=False)\n"
                "iod_species = pd.read_csv(DATA / 'relative_nuclear_iod_species.csv')\n\n"
                "assert (panel['species'].nunique(), objects.shape[0], "
                "iod_objects.shape[0]) == (24, 1152, 805)\n"
                "{'Path24 species': panel['species'].nunique(), "
                "'reviewed cell-nucleus objects': len(objects), "
                "'reviewed IOD nuclei': len(iod_objects)}"
            ),
            "species.merge(iod_species, on='species', suffixes=('_morphology', '_iod'))",
            (
                "panel.loc[panel['manual_review_inclusion'], "
                "['species', 'decision_basis', 'has_finalized_genome_iod', "
                "'has_finalized_cell_morphology']]"
            ),
        ],
    ),
    "analyses/03_phylogenetic_path/explore.ipynb": notebook(
        "Phylogenetic path data",
        "The final fitted model is the three-trait Path24 analysis: relative nuclear "
        "IOD, nucleus area, and cell area. It is exploratory and does not uniquely "
        "identify causal arrow direction. SVL and life-history tables are retained "
        "for inspection but are not predictors in this fitted Path24 model.",
        [
            (
                "DATA = ROOT / 'analyses/03_phylogenetic_path/data'\n"
                "traits = pd.read_csv(DATA / 'path24_traits.csv')\n"
                "organismal = pd.read_csv(DATA / 'organismal_traits.csv')\n"
                "models = pd.read_csv(DATA / 'path_model_comparison.csv')\n"
                "edges = pd.read_csv(DATA / 'path_edge_estimates.csv')\n"
                "overlap = pd.read_csv(DATA / 'te_relative_iod_overlap.csv')\n\n"
                "assert (traits['tree_tip'].nunique(), len(overlap)) == (24, 21)\n"
                "{'Path24 species': traits['tree_tip'].nunique(), "
                "'TE-IOD overlap': len(overlap), 'candidate classes': len(models)}"
            ),
            (
                "organismal_for_join = organismal.rename("
                "columns={'species': 'organismal_species'})\n"
                "traits.merge(organismal_for_join, left_on='tree_tip', "
                "right_on='organismal_species', how='left', validate='one_to_one')[[\n"
                "'species', 'relative_iod_index', 'nucleus_area_um2', "
                "'cell_area_um2', 'max_svl_mm', 'reproductive_strategy']].sort_values('species')"
            ),
            (
                "organismal[['species', 'max_svl_mm', 'development_mode', "
                "'reproductive_strategy', 'aquaticity_index', "
                "'microhabitat_class']].sort_values('species')"
            ),
            "models.sort_values('delta_CICc').head(10)",
            "edges",
        ],
    ),
}


def main() -> None:
    for relative, book in NOTEBOOKS.items():
        path = ROOT / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        nbf.write(book, path)
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()

# Manuscript Figure Plan

This file translates the frozen TE analysis state into a practical main-text versus supplement asset plan.

## Main Text Figures

- `Figure 1`: `path_analysis/results/te_genome_primary_mediumplus_best_model.pdf`
  Purpose: Baseline path model for the observed-only TE plus genome family.
  Why it belongs in the main text: Shows the cleanest mechanistic TE result before adding ectopic or organismal competition.
- `Figure 2`: `path_analysis/results/te_genome_ectopic_organismal_primary_mediumplus_best_model.pdf`
  Purpose: Combined-family best model when ectopic support and body size are allowed to compete.
  Why it belongs in the main text: Captures the most important updated inference: body size remains while ectopic_index drops.
- `Figure 3`: `results/figures/pca/superfamily_primary_scores_pc1_pc2.png`
  Purpose: Primary TE ordination at the superfamily level.
  Why it belongs in the main text: Uses the canonical compositional PCA and provides one visually interpretable summary of TE composition structure.

## Supplementary Figures

- `Figure S1`: `path_analysis/results/te_genome_primary_strict_body_best_model.pdf`
  Purpose: Strict-body robustness check for the baseline TE family.
  Why it belongs in the supplement: Confirms winner stability under adult-oriented body-size restriction.
- `Figure S2`: `path_analysis/results/te_genome_ectopic_primary_mediumplus_best_model.pdf`
  Purpose: Best ectopic-only family before organismal competition.
  Why it belongs in the supplement: Useful to show the intermediate result that ectopic support is strongest only before body-size competition.
- `Figure S3`: `path_analysis/results/te_genome_organismal_primary_mediumplus_best_model.pdf`
  Purpose: Observed-only TE plus organismal best model.
  Why it belongs in the supplement: Documents the body-size result directly in the organismal family.
- `Figure S4`: `results/figures/pca/superfamily_primary_ppca_phylomorphospace.png`
  Purpose: Phylogenetic PCA view of the main superfamily ordination.
  Why it belongs in the supplement: Keeps the phylogenetically structured ordination available without making it the main visual.
- `Figure S5`: `results/figures/pca/superfamily_primary_scree_plot.png`
  Purpose: Variance explained for the primary superfamily PCA.
  Why it belongs in the supplement: Supports the ordination methods and PC1/PC2 interpretation.
- `Figure S6`: `results/figures/phylo_signal/blombergs_k_signal.png`
  Purpose: Trait-wise Blomberg's K summary from the focused phylogenetic-signal workflow.
  Why it belongs in the supplement: Documents broader TE signal patterns without crowding the main text.
- `Figure S7`: `results/figures/phylo_signal/phylogenetic_correlogram_moran_per_trait.png`
  Purpose: Distance-binned Moran's I correlogram summary.
  Why it belongs in the supplement: Extends the signal results beyond global K and lambda statistics.
- `Figure S8`: `results/figures/landscape/te_order_distribution.png`
  Purpose: Across-species TE order landscape summary.
  Why it belongs in the supplement: Provides a compact entry point to the full landscape rebuild without showing dozens of per-species curves.
- `Figure S9`: `results/figures/phylo_landscape/phylogeny_with_landscape.png`
  Purpose: Phylogeny-linked landscape heatmap summary.
  Why it belongs in the supplement: Connects the landscape rebuild to the comparative tree in a single supplemental graphic.

## Supplementary Tables

- `Table S1`: `paper_freeze/manuscript_primary_results_summary.csv`
  Purpose: Primary and robustness path-model winner summary.
- `Table S2`: `results/tables/pca/te_pca_analysis_manifest.csv`
  Purpose: Canonical PCA analysis manifest and variance explained.
- `Table S3`: `results/tables/pca/te_pca_phylogenetic_signal.csv`
  Purpose: Phylogenetic signal statistics for PCA axes.
- `Table S4`: `paper_freeze/key_file_manifest.csv`
  Purpose: Checksummed freeze manifest for critical inputs and derived tables.

## Recommended Narrative Order

- Start with the baseline TE family and the observed-only medium-plus subset as the clean primary comparative result.
- Then show the combined ectopic-plus-organismal family to make the key updated point: body size persists while ectopic_index drops out.
- Use the superfamily compositional PCA as the single TE ordination figure, but keep it framed as supplementary structure rather than the primary causal result.

## Primary Result Snapshot

- `te_genome_primary_mediumplus`: winner `mediated_evenness` at `n = 27`; Baseline TE model favors mediated_evenness and supports a stable LTR-balance to TE-evenness relationship.
- `te_genome_ectopic_primary_mediumplus`: winner `ectopic_only` at `n = 24`; Adding ectopic support favors ectopic_only before organismal competition is introduced.
- `te_genome_organismal_primary_mediumplus`: winner `body_size_additive` at `n = 27`; Body size is the only organismal covariate that meaningfully improves the TE-genome family.
- `te_genome_ectopic_organismal_primary_mediumplus`: winner `te_body_size_baseline` at `n = 24`; When ectopic support and body size compete directly, the best model drops ectopic_index and keeps TE plus body size.

## Exclusions From Main Text

- Morphology-linked families remain sensitivity-only because the genome-size side is still provisional.
- Strict-body reruns should stay supplemental because they are robustness checks on species inclusion, not the primary observed-only dataset.
- Order-level all-feature ordinations and their pPCA variants should stay supplemental because they are less interpretable than the superfamily primary ordination.


# Order-level PCA sensitivity across 200 published time trees

The corrected order composition was projected from CLR to a 9-dimensional orthonormal ILR basis and refit with Brownian-covariance phylogenetic PCA on the Stewart-Wiens (2025) optimal time tree plus all 200 time-calibrated bootstrap trees. Axes were matched to the collaborator-focal pPCA using loading correlations; signs were aligned for reporting only.

## Stability results

- Axis swaps required: 0/200 bootstrap trees.
- Minimum ordinary-versus-phylogenetic species-score correlation: 0.988 for axis 1 and 0.910 for axis 2.
- Minimum ordinary-versus-phylogenetic loading correlation: 0.955 for axis 1 and 0.723 for axis 2.
- Minimum published-tree versus collaborator-focal species-score correlation: 0.995 for pPC1 and 0.971 for pPC2.
- Median variance explained across trees: 51.4% for pPC1 and 24.8% for pPC2.
- 1 of 20 feature-by-axis 95% loading intervals cross zero; these loadings must not be given direction-stable biological interpretations.

## Verdict

**Approved as a descriptive ordination across fixed-tree and published tree-time uncertainty.** The dominant order-composition geometry, especially pPC1, is robust to shared ancestry and the 200 published time trees. This still does not promote PCA axes to causal path variables: direct declared log-ratios or a multivariate phylogenetic model remain preferable for mechanistic inference, and the bootstrap set does not represent reticulation uncertainty.

## Outputs

- `results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_metrics_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_scores_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_score_summary_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_loadings_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_pca_tree_uncertainty_loading_summary_analysis18_v1.csv`
- `results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_score_intervals_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_loading_intervals_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_correlation_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_pca_tree_uncertainty_variance_analysis18_v1.png`

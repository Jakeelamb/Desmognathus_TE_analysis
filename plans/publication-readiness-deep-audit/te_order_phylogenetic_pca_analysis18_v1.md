# Order-level phylogenetic PCA sensitivity

The approved descriptive order CLR-PCA was refit on the exact corrected 18-tip time tree. Because CLR coordinates are singular by construction, the CLR matrix was projected into a 9-dimensional orthonormal ILR basis before `phytools::phyl.pca(method = "BM", mode = "cov")`; loadings were then projected back into CLR feature space.

## Results

- pPC1 explains 53.5% and pPC2 explains 23.6% of phylogenetic ILR variance (77.1% cumulative).
- Ordinary versus phylogenetic species-score correlations are 0.983 for axis 1 and 0.913 for axis 2 after sign alignment; axis signs have no biological meaning.
- Ordinary versus back-projected loading correlations are 0.954 for axis 1 and 0.742 for axis 2. The dominant axis is strongly preserved; axis 2 is recognizably similar but more sensitive to phylogenetic covariance.
- Minimum leave-one-species-out loading/score correlations are 0.994/0.998 for pPC1 and 0.987/0.985 for pPC2.

## Verdict

**Approved as a fixed-tree descriptive sensitivity.** The principal order-composition gradient is not an artifact of ignoring shared ancestry. This does not promote either PCA axis to a causal path variable: the tree source/posterior uncertainty remains unresolved, and direct prespecified log-ratios or multivariate phylogenetic models are preferable for inference.

## Outputs

- `results/data/corrected/diversity_pca/te_order_phylogenetic_pca_scores_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_phylogenetic_pca_loadings_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_phylogenetic_pca_variance_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_phylogenetic_pca_comparison_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_order_phylogenetic_pca_stability_analysis18_v1.csv`
- `results/figures/corrected/diversity_pca/te_order_phylogenetic_pca_scores_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_ordinary_vs_phylogenetic_pca_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_phylogenetic_pca_loadings_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_phylogenetic_pca_stability_analysis18_v1.png`

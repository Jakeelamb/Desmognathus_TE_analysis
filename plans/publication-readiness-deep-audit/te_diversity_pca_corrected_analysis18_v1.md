# Corrected TE diversity and compositional PCA audit

This branch is restricted to the declared 18-species TE/genome panel and preserves all historical diversity/PCA outputs.

## Diversity definitions

- `Simpson_Diversity` in the historical tables is reproduced as Gini-Simpson `1 - sum(p_i^2)`; it is not a finite-count-corrected Simpson estimator.
- Shannon entropy uses natural logarithms.
- Pielou evenness is `H / log(S_observed)` and is conditional on positive classified categories.
- Hill numbers `exp(H)` (q=1) and `1 / sum(p_i^2)` (q=2) are now included.
- Biological indices conditional on classified mass are kept separate from a technical sensitivity that appends unresolved mass as one non-biological bin.

Maximum historical reproduction error is 3.331e-16 for order and 4.441e-16 for superfamily.

Rank sensitivity to unresolved mass:
- order gini_simpson: Spearman rho = 0.975.
- order shannon_entropy: Spearman rho = 0.965.
- order pielou_evenness: Spearman rho = 0.965.
- superfamily gini_simpson: Spearman rho = 0.986.
- superfamily shannon_entropy: Spearman rho = 0.979.
- superfamily pielou_evenness: Spearman rho = 0.930.

## PCA

The primary descriptive ordination is order-level CLR PCA conditional on classified order mass (10 features; PC1 56.7%, PC2 28.9%). Superfamily PCA is supplementary because p approaches/exceeds n, sparse zeros require replacement, and rare-feature retention changes geometry. Both global-half-minimum and feature-half-minimum zero replacements are exported as sensitivities, along with exact CLR matrices, loadings, scores, variance, and leave-one-species-out stability.

Minimum leave-one-out loading/score correlations by analysis are available in `te_pca_stability_analysis18_v1.csv`; no PCA axis is approved as a causal variable without the phylogenetic and zero-replacement checks.

## PERMANOVA verdict

The current clade PERMANOVA is **not approved for inference**. It permutes species labels freely even though groups are phylogenetic clades, violating exchangeability under phylogenetic covariance. Its Bray-Curtis and CLR sensitivity calculations and beta-dispersion checks are useful descriptively, but p-values must not support a clade claim. A phylogenetically valid simulation/permutation or comparative multivariate model is required.

## Validation figures

- `results/figures/corrected/diversity_pca/te_diversity_unresolved_mass_sensitivity_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_clr_pca_scores_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_superfamily_clr_pca_scores_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_order_clr_pca_loadings_analysis18_v1.png`
- `results/figures/corrected/diversity_pca/te_superfamily_pca_zero_replacement_sensitivity_analysis18_v1.png`

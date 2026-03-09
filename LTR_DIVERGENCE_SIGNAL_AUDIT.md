# LTR Divergence Signal Audit

This note summarizes what the sequence-derived paired-LTR divergence layer is
saying biologically after the assembly-backed extraction workflow completed.

## Bottom Line

- Successful sequence-derived LTR divergence estimates: `1291` pairs across `30` species.
- The estimated set is overwhelmingly `Gypsy` (`1283 / 1291` successful pairs).
- The recommended high-confidence subset (`Complete = yes` and `5+` domains) has a slightly lower median K2P distance (`0.066201`) than the full estimated set (`0.070564`).
- Species-level median LTR divergence shows weak phylogenetic structure (`K = 0.3591`, `p = 0.464`; `lambda = 0.0001`, `p = 1`).
- Species-level median LTR divergence does not show a strong monotonic relationship with the existing TE/path summary features tested here.

## Species Extremes

Highest median K2P species:
- `D.aureatus`: median K2P `0.089795` from `27` estimated pairs
- `D.ocoee`: median K2P `0.081763` from `38` estimated pairs
- `D.monticola`: median K2P `0.079247` from `29` estimated pairs
- `D.campi`: median K2P `0.076760` from `44` estimated pairs
- `D.marmoratus`: median K2P `0.076287` from `15` estimated pairs

Lowest median K2P species:
- `D.aeneus`: median K2P `0.053658` from `33` estimated pairs
- `D.conanti`: median K2P `0.057500` from `29` estimated pairs
- `D.balsameus`: median K2P `0.064368` from `38` estimated pairs
- `D.orestes`: median K2P `0.064740` from `40` estimated pairs
- `D.mavrokoilius`: median K2P `0.066286` from `34` estimated pairs

## Feature Correlations

- `ectopic_log10_mean_ratio`: Spearman rho `0.291`, `p = 0.1190`, `n = 30`
- `ltr_line_logratio`: Spearman rho `0.229`, `p = 0.2228`, `n = 30`
- `order_pielou`: Spearman rho `-0.100`, `p = 0.5978`, `n = 30`
- `ltr_divergence_p90`: Spearman rho `0.094`, `p = 0.6215`, `n = 30`
- `weighted_te_divergence_p90`: Spearman rho `0.001`, `p = 0.9953`, `n = 30`

## Interpretation

- The new LTR divergence layer appears to add information that is not trivially redundant with the current TE composition and ectopic summary features.
- Because the signal is weakly phylogenetically structured and not strongly tied to the existing TE/path features, it looks more species-specific than deeply clade-conserved at the current species coverage.
- This branch is now suitable for supplementary comparative use as a sequence-derived LTR recency/divergence layer, but absolute insertion ages still require an externally justified substitution-rate calibration.

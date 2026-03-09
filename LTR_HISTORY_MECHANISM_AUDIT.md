# LTR History Mechanism Audit

## Purpose

This note tests whether the paired-LTR history layer helps explain current TE
state, or whether it is mostly orthogonal to the current TE composition and
landscape summaries.

## Upstream inputs

- `path_analysis/data/derived/ltr_history_features.csv`
- `path_analysis/data/derived/te_path_features.csv`
- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`

## Outputs

- `results/data/ltr_history_mechanism/ltr_history_mechanism_tests.csv`

## Coverage

- Species in the merged mechanism audit: `30`

## Predictor interpretation

- `ltr_history_median_k2p_distance`: oldness of retained paired LTRs
- `ltr_history_n_pairs_estimated`: amount of recoverable paired-LTR substrate
- `ltr_history_pairs_per_million_ltr_bp`: recoverable paired-LTR density after
  normalizing by total LTR landscape mass
- `ltr_history_n_pairs_high_confidence`: stricter support subset

## Main result

The paired-LTR history branch splits into two different signals:

- retained pair age is mostly independent of the current TE state
- retained pair abundance/density tracks landscape recency more than it tracks
  TE balance or ectopic ratio

Key checks:

- `ltr_history_median_k2p_distance` vs `weighted_te_divergence_p90`: rho = 0.001, BH p = 0.995, n = 30
- `ltr_history_median_k2p_distance` vs `all_recent_mass_frac_0_5`: rho = -0.008, BH p = 0.992, n = 30
- `ltr_history_n_pairs_estimated` vs `all_recent_mass_frac_0_5`: rho = 0.599, BH p = 0.0186, n = 30
- `ltr_history_n_pairs_estimated` vs `all_old_tail_frac_20plus`: rho = -0.493, BH p = 0.0281, n = 30
- `ltr_history_pairs_per_million_ltr_bp` vs `all_recent_mass_frac_0_5`: rho = 0.561, BH p = 0.022, n = 30
- `ltr_history_pairs_per_million_ltr_bp` vs `ltr_recent_mass_frac_0_5`: rho = 0.505, BH p = 0.0281, n = 30
- `ltr_history_n_pairs_estimated` vs `ectopic_log10_mean_ratio`: rho = 0.043, BH p = 0.967, n = 30

## Interpretation

- The age/oldness axis (`ltr_history_median_k2p_distance`) does not explain the
  current broad TE divergence or landscape recency summaries.
- The abundance axis (`ltr_history_n_pairs_estimated`) does connect to current
  TE tempo: species with more recoverable paired LTRs tend to show younger,
  more recent landscapes and lower old-tail burden.
- That signal survives simple normalization by total LTR landscape mass, so it
  is not only a trivial "more LTR sequence means more pairs" artifact.
- The ectopic proxy stays largely separate from the paired-LTR mechanism layer.
- The strict high-confidence pair count is weaker than the full recoverable
  count, so this abundance link is real enough to note but should remain a
  secondary mechanism result rather than a headline claim.

## Bottom line

The paired-LTR history layer is not one thing. Retained pair age behaves like a
mostly orthogonal historical axis, while retained pair abundance behaves more
like a recency/recoverability signal tied to the present-day TE landscape. That
means the current repo supports a nuanced mechanism story: old retained LTR
history and current TE tempo are related only weakly, but the amount of intact
paired-LTR substrate does carry information about how recent the present TE
landscape looks.

## Audit summary

- BH-significant tests for paired-LTR age: `0`
- BH-significant tests for paired-LTR count: `4`
- BH-significant tests for paired-LTR density: `4`

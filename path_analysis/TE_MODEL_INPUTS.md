# TE Model Inputs

This note records the recommended TE predictors for the comparative path-analysis stage.

## Recommended TE blocks

- `ltr_line_logratio`
  Primary composition axis. Higher values indicate relatively greater LTR representation versus LINEs.
- `order_pielou`
  Order-level evenness axis, conditional on order-classified mass. Higher
  values indicate a more even TE order composition. Use with the unresolved-bin
  rank sensitivity; do not present it as absolute-genome TE diversity.
- order-level CLR composition scores
  Descriptive ordination only. The corrected order CLR PCA is stable to
  leave-one-species-out analysis, fixed-tree ILR-pPCA, and the 200 published
  Stewart-Wiens time trees. Minimum ordinary-versus-pPCA score correlations are
  0.988/0.910 for axes 1/2. Despite that stability, do not promote axes to causal
  predictors; use declared log-ratios or a direct multivariate phylogenetic model.
- superfamily CLR composition scores
  Supplementary only because zero replacement and p >= n materially affect
  the geometry; do not promote these scores to the confirmatory path set.
- `weighted_te_divergence_p90`
  Recommended turnover or age-structure proxy for the main model set.
- `weighted_te_deletions_p90`
  RepeatMasker alignment-gap sensitivity statistic, not a direct DNA-loss rate; do not combine with `weighted_te_divergence_p90` unless collinearity is checked explicitly.
- corrected terminal:internal LTR depth summaries
  Exploratory mapping/structure proxy only. The historical `ectopic_log10_mean_ratio` arithmetic mean is not recommended. If this layer is shown, use the zero-aware species median or mean element-level log-ratio with bootstrap intervals and influence diagnostics from `results/data/corrected/ectopic/`.

## Why this panel

- It avoids overloading the path model with many correlated TE fractions.
- It aims for one observed variable per TE concept block, but the alignment-gap and terminal:internal layers remain sensitivity-only rather than validated loss or ectopic-rate constructs.
- Diversity variables are explicitly conditional on classified mass; unresolved
  mass is carried in a separate sensitivity rather than silently discarded.
- It maps directly onto the candidate DAG families already scaffolded in `path_analysis`.

## Data products

- `data/derived/te_path_features.csv`
  Full reusable TE feature table with provenance columns.
- `data/derived/te_model_feature_panel.csv`
  Compact analysis-facing TE panel for model assembly and overlap checks.

## Presentation guidance

- Use the compact panel for path models and main-text comparative tables.
- Reserve the full order and superfamily breakdowns for supplements, QC, and feature-engineering transparency.
- Treat `weighted_te_divergence_p90` and `weighted_te_deletions_p90` as alternative mechanistic summaries unless the final model set justifies including both.

# TE Model Inputs

This note records the recommended TE predictors for the comparative path-analysis stage.

## Recommended TE blocks

- `ltr_line_logratio`
  Primary composition axis. Higher values indicate relatively greater LTR representation versus LINEs.
- `order_pielou`
  Order-level evenness axis. Higher values indicate a more even TE order composition.
- `weighted_te_divergence_p90`
  Recommended turnover or age-structure proxy for the main model set.
- `weighted_te_deletions_p90`
  Recommended DNA-loss proxy for sensitivity analyses, not for the same small model as `weighted_te_divergence_p90` unless collinearity is checked explicitly.
- `ectopic_log10_mean_ratio`
  Ectopic recombination proxy for the extended `te_genome_ectopic` family.

## Why this panel

- It avoids overloading the path model with many correlated TE fractions.
- It keeps one observed variable per TE concept block: composition, evenness, turnover or loss, and ectopic constraint.
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

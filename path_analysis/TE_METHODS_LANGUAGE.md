# TE Methods Language

This file is the analysis-facing wording guide for the TE layer used in the
current `path_analysis/` workflow. It translates the audit conclusions in
`TE_DATA_AUDIT.md` into short methods-ready text that does not accidentally
overstate what the current TE tables mean.

## Recommended Short Methods Text

Use this wording if you need a compact paragraph for the main Methods section:

> We used the repo-local species-level TE summary tables stored in
> `results/data/` as the canonical TE input layer for comparative analyses.
> Order-level TE composition was summarized with Shannon entropy,
> Gini-Simpson diversity (`1 - sum(p_i^2)`), and Pielou evenness. For the
> path-analysis families, we
> used order-level Pielou evenness as the primary evenness variable because it
> is the cleanest interpretable evenness metric in the current stored tables.
> Terminal:internal LTR depth was retained only as an exploratory mapping and
> structure proxy. Corrected element summaries used exact-coordinate TEsorter
> joins and included zero-depth positions. Species-level robustness was
> evaluated with medians, geometric means, bootstrap intervals, coverage
> thresholds, and leave-one-element-out influence. The ratio was not
> interpreted as a measured ectopic-recombination rate.

## Recommended Detailed Methods Text

Use this if the analysis notes or supplement needs more explicit semantics:

> The TE predictors used in the comparative path-analysis layer were taken from
> the repo-local TE summary tables rather than regenerated from raw HPC
> outputs. Order-level TE composition was summarized from the stored
> `dnaPipeTE_order_breakdown.csv` and diversity summary tables. The current
> `Simpson_Diversity` values in the stored order-level diversity table are
> standard Gini-Simpson values, `1 - sum(p_i^2)`, calculated after normalizing
> positive relative abundances. They are not finite-count-corrected estimates.
> Because the stored evenness tables are
> internally consistent and directly traceable, we treated them as the
> analysis TE input layer and used order-level Pielou evenness as the main
> evenness predictor in the path models.
>
> The historical ectopic table used nonzero-only regional means and an
> arithmetic species mean. The corrected sensitivity branch joins each LTR to
> TEsorter by full `sequence_start_end`, includes explicit zero-depth positions,
> reports left-LTR/right-LTR/internal positive-position coverage, and provides
> all-element, >=80% coverage, TEsorter-complete, and combined branches. The
> historical arithmetic mean remains available for provenance but is not the
> recommended estimator.

## Explicit Variable Semantics To Preserve

If these variables are named in the results, describe them this way:

- `order_pielou`
  Order-level Pielou evenness derived from the stored order diversity table.
- `order_simpson`
  Gini-Simpson diversity, `1 - sum(p_i^2)`, from the stored order diversity
  table. State the convention because other Simpson indices reverse or
  rescale the interpretation.
- `ectopic_n_rows_total`
  Total filtered ectopic rows staged for a species.
- `ectopic_n_elements`
  Count of ectopic rows with usable non-missing ratio support.
- `ectopic_mean_ratio`
  Historical arithmetic mean terminal:internal ratio; pre-audit and
  outlier-sensitive, not recommended for confirmatory inference.
- `ectopic_complete_fraction`
  Fraction of known `Complete` calls that are `yes`, excluding `unknown`
  records from the denominator.
- `ectopic_n_complete_known`
  Number of ectopic rows contributing to `ectopic_complete_fraction`.

## Language To Avoid

Avoid these shortcuts in the results:

- "Simpson diversity" without stating that the stored value is Gini-Simpson
  `1 - sum(p_i^2)`.
- "All ectopic elements" if the statistic excludes missing-ratio rows.
- "Complete fraction across all rows" if `Complete = unknown` rows are excluded.
- Any wording that implies the original upstream HPC TE analyses were rerun for
  this audit pass.
- "Ectopic-recombination rate", "solo-LTR formation rate", or "deletion rate"
  for terminal:internal depth without independent structural validation.

## Canonical Files Behind This Language

The wording above is anchored to these current files:

1. `results/data/dnaPipeTE_order_breakdown.csv`
2. `results/data/diversity_order_stats.csv`
3. `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
4. `results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv`
5. `results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv`
6. `path_analysis/data/derived/te_path_features.csv`
7. `path_analysis/data/derived/te_model_feature_panel.csv`

## Practical Recommendation

For analysis writeups:

- use `order_pielou` as the narrative evenness metric
- report Gini-Simpson alongside Shannon, richness, and dominance/evenness when
  making a diversity claim
- report ectopic denominator, coverage, estimator, bootstrap, and influence
  semantics whenever this exploratory layer appears

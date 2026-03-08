# TE Methods Language

This file is the manuscript-facing wording guide for the TE layer used in the
current `path_analysis/` workflow. It translates the audit conclusions in
`TE_DATA_AUDIT.md` into short methods-ready text so the paper does not
accidentally overstate what the current TE tables mean.

## Recommended Short Methods Text

Use this wording if you need a compact paragraph for the main Methods section:

> We used the frozen species-level TE summary tables stored in
> `results/data/` as the canonical TE input layer for comparative analyses.
> Order-level TE composition was summarized with Shannon diversity, corrected
> Simpson diversity, and Pielou evenness. For the path-analysis families, we
> used order-level Pielou evenness as the primary evenness variable because it
> is the cleanest interpretable evenness metric in the current stored tables.
> Ectopic-recombination summaries were derived from the filtered species table
> already present in the repository; rows with missing
> `ratio_terminal_internal` were excluded from ratio-based summaries, and
> `ectopic_complete_fraction` was calculated only from rows with known
> `Complete = yes` or `no`.

## Recommended Detailed Methods Text

Use this if the chapter or supplement needs more explicit semantics:

> The TE predictors used in the comparative path-analysis layer were taken from
> the frozen repo-local TE summary tables rather than regenerated from raw HPC
> outputs. Order-level TE composition was summarized from the stored
> `dnaPipeTE_order_breakdown.csv` and diversity summary tables. The current
> `Simpson_Diversity` values in the stored order-level diversity table are the
> corrected Simpson values already present in the canonical file, not a newly
> recomputed raw Gini-Simpson index. Because the stored evenness tables are
> internally consistent and directly traceable, we treated them as the
> manuscript TE input layer and used order-level Pielou evenness as the main
> evenness predictor in the path models.
>
> Ectopic-recombination support was summarized from the filtered ectopic table
> already present in `results/data/`. Species-level ectopic support metrics were
> defined explicitly as follows: `ectopic_n_rows_total` counts all staged
> filtered rows for a species, `ectopic_n_elements` counts only rows with usable
> non-missing `ratio_terminal_internal` values, `ectopic_mean_ratio` is
> calculated only from those usable ratio-bearing rows, and
> `ectopic_complete_fraction` is calculated only from rows with known
> `Complete = yes` or `Complete = no`. Rows with `Complete = unknown` are
> excluded from the denominator, and `ectopic_n_complete_known` records that
> denominator explicitly.

## Explicit Variable Semantics To Preserve

If these variables are named in the manuscript, describe them this way:

- `order_pielou`
  Order-level Pielou evenness derived from the stored order diversity table.
- `order_simpson`
  Corrected Simpson diversity from the stored order diversity table; do not
  call this the raw Gini-Simpson value without clarification.
- `ectopic_n_rows_total`
  Total filtered ectopic rows staged for a species.
- `ectopic_n_elements`
  Count of ectopic rows with usable non-missing ratio support.
- `ectopic_mean_ratio`
  Mean terminal:internal ratio across usable ratio-bearing ectopic rows only.
- `ectopic_complete_fraction`
  Fraction of known `Complete` calls that are `yes`, excluding `unknown`
  records from the denominator.
- `ectopic_n_complete_known`
  Number of ectopic rows contributing to `ectopic_complete_fraction`.

## Language To Avoid

Avoid these shortcuts in the manuscript:

- "Simpson diversity" with no qualifier if you mean the stored corrected value.
- "All ectopic elements" if the statistic excludes missing-ratio rows.
- "Complete fraction across all rows" if `Complete = unknown` rows are excluded.
- Any wording that implies the original upstream HPC TE analyses were rerun for
  this audit pass.

## Canonical Files Behind This Language

The wording above is anchored to these current files:

1. `results/data/dnaPipeTE_order_breakdown.csv`
2. `results/data/diversity_order_stats.csv`
3. `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
4. `path_analysis/data/derived/te_path_features.csv`
5. `path_analysis/data/derived/te_model_feature_panel.csv`

## Practical Recommendation

For the paper:

- use `order_pielou` as the narrative evenness metric
- mention corrected Simpson only if needed for completeness
- report ectopic denominator semantics whenever `ectopic_n_elements` or
  `ectopic_complete_fraction` appear in text, tables, or supplements

# LTR History Path Integration

This note records how the sequence-derived paired-LTR history layer was
integrated into the path-analysis workspace and what it currently supports.

## Source Basis

- Sequence-divergence workflow:
  - [`LTR_AGE_AUDIT.md`](/home/jake/Projects/Desmognathus_TE/LTR_AGE_AUDIT.md)
  - [`LTR_DIVERGENCE_SIGNAL_AUDIT.md`](/home/jake/Projects/Desmognathus_TE/LTR_DIVERGENCE_SIGNAL_AUDIT.md)
- Primary-literature rate calibration:
  - [`LTR_SUBSTITUTION_RATE_CALIBRATION.md`](/home/jake/Projects/Desmognathus_TE/LTR_SUBSTITUTION_RATE_CALIBRATION.md)
  - Herrick and Sclavi 2014: https://doi.org/10.1016/j.crpv.2014.06.002
  - Crawford 2003: https://doi.org/10.1007/s00239-003-2513-7
- Derived path-analysis inputs:
  - [`ltr_history_features.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/ltr_history_features.csv)
  - [`path_input_master.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/path_input_master.csv)
  - [`analysis_panel_summary.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_panel_summary.csv)

## Modeling Choice

- The path-analysis variable is the species-level median paired-LTR `K2P`
  divergence, not a single hard-dated age estimate.
- This is deliberate: for any fixed positive substitution rate, calibrated age is
  a constant multiple of `K2P`, so the standardized path-model covariate is
  mathematically identical whether expressed as `K2P` or calibrated age.
- The literature calibration is therefore used for biological interpretation of
  timescale, not to define a separate model predictor.

## Added Data Layer

- [`prepare_ltr_history_features.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/prepare_ltr_history_features.py)
  builds a canonical species-level LTR-history feature table from the audited
  LTR divergence and calibration outputs.
- The current feature table carries:
  - species-level `median_k2p_distance`
  - estimated-pair counts and high-confidence pair counts
  - low/central/high calibrated median age summaries
  - explicit source IDs, paths, and checksums

## Panel Overlap

- `te_genome_ltr_history_primary_mediumplus`: `24` species
- `te_genome_ltr_history_primary_strict_body`: `15` species
- The TE+genome species currently missing LTR-history coverage are `kanawha`,
  `lycos`, and `valtos`.

## Sensitivity Family

Files:

- [`te_genome_ltr_history_primary_mediumplus_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ltr_history_primary_mediumplus_model_ranking.csv)
- [`te_genome_ltr_history_primary_strict_body_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ltr_history_primary_strict_body_model_ranking.csv)
- [`te_genome_ltr_history_primary_mediumplus_best_model_edges.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ltr_history_primary_mediumplus_best_model_edges.csv)
- [`te_genome_ltr_history_primary_strict_body_best_model_edges.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ltr_history_primary_strict_body_best_model_edges.csv)

Candidate models:

- `te_baseline`
- `history_direct`
- `history_additive`
- `history_to_balance`
- `history_to_evenness`

## Current Result

- In the `n = 24` medium-plus panel, `history_additive` wins clearly
  (`delta_CICc = 4.16` over `te_baseline`, model weight `0.81`).
- The best medium-plus model keeps the existing `ltr_balance -> te_evenness`
  path and adds a positive direct `ltr_history -> gs` edge (`0.357 +/- 0.171`).
- In the stricter `n = 15` body-filtered panel, `te_baseline` ranks first, but
  `history_additive` remains close (`delta_CICc = 0.65`) and `history_direct`
  is also competitive (`delta_CICc = 1.49`).
- The `te_evenness -> gs` edge remains weak in the history family, just as it
  does in the baseline TE family.

## Interpretation Boundary

- The historical LTR axis is informative enough to keep as a sensitivity family.
- It is not yet strong enough to replace the current primary TE-family results.
- The clean interpretation is that older retained paired-LTR histories may
  explain some additional genome-size variation in the broader medium-plus set,
  but that signal is not fully stable once the strictest body-size filter is
  applied.
- Because the underlying LTR layer is overwhelmingly `Gypsy`-dominated, this is
  best understood as a Gypsy-heavy historical axis rather than a global TE-age
  clock for all TE classes.

## Practical Recommendation

- Keep `te_genome_ltr_history` as a sensitivity family.
- Use calibrated age ranges from
  [`LTR_SUBSTITUTION_RATE_CALIBRATION.md`](/home/jake/Projects/Desmognathus_TE/LTR_SUBSTITUTION_RATE_CALIBRATION.md)
  in the manuscript text, but keep the path-model predictor itself as
  divergence-based `ltr_history`.

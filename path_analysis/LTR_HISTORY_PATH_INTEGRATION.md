# LTR History Path Integration

This note records how the sequence-derived paired-LTR history layer was
integrated into the path-analysis workspace and what it currently supports.

## Source Basis

- Sequence-divergence workflow:
  - `results/reports/LTR_AGE_AUDIT.md`
  - `results/reports/LTR_DIVERGENCE_SIGNAL_AUDIT.md`
- Primary-literature rate calibration:
  - `LTR_SUBSTITUTION_RATE_CALIBRATION.md`
  - Herrick and Sclavi 2014: https://doi.org/10.1016/j.crpv.2014.06.002
  - Crawford 2003: https://doi.org/10.1007/s00239-003-2513-7
- Derived path-analysis inputs:
  - `path_analysis/data/derived/ltr_history_features.csv`
  - `path_analysis/data/derived/path_input_master.csv`
  - `path_analysis/data/derived/analysis_panel_summary.csv`

## Modeling Choice

- The path-analysis variable is the species-level median paired-LTR `K2P`
  divergence, not a single hard-dated age estimate.
- This is deliberate: for any fixed positive substitution rate, calibrated age is
  a constant multiple of `K2P`, so the standardized path-model covariate is
  mathematically identical whether expressed as `K2P` or calibrated age.
- The literature calibration is therefore used for biological interpretation of
  timescale, not to define a separate model predictor.

## Added Data Layer

- `path_analysis/scripts/prepare_ltr_history_features.py`
  builds a canonical species-level LTR-history feature table from the audited
  LTR divergence and calibration outputs.
- The current feature table carries:
  - species-level `median_k2p_distance`
  - estimated-pair counts and high-confidence pair counts
  - low/central/high calibrated median age summaries
  - explicit source IDs, paths, and checksums

## Panel Overlap

- `te_genome_ltr_history_primary_mediumplus`: `16` species
- `te_genome_ltr_history_primary_strict_body`: `8` species
- Relative to the current TE+genome medium-plus panel, the species missing
  LTR-history coverage are `kanawha` and `valtos`.
- Relative to the current strict-body panel, the species missing LTR-history
  coverage is `kanawha`.

## Sensitivity Family

Files:

- `path_analysis/results/te_genome_ltr_history_primary_mediumplus_model_ranking.csv`
- `path_analysis/results/te_genome_ltr_history_primary_strict_body_model_ranking.csv`
- `path_analysis/results/te_genome_ltr_history_primary_mediumplus_best_model_edges.csv`
- `path_analysis/results/te_genome_ltr_history_primary_strict_body_best_model_edges.csv`

Candidate models:

- `te_baseline`
- `history_direct`
- `history_additive`
- `history_to_balance`
- `history_to_evenness`

## Current Result Boundary

- The source layer has been rebuilt from local paired-LTR sequence extraction
  and substitution-rate calibration outputs.
- Current panel sizes are `n = 16` for the medium-plus LTR-history panel and
  `n = 8` for the strict-body LTR-history panel.
- Model rankings for this family need to be rerun before carrying forward older
  winner claims.

## Interpretation Boundary

- The historical LTR axis is source-traceable enough to keep as a sensitivity
  family.
- It should not replace the current primary TE-family results unless refreshed
  model rankings support that promotion.
- Because the underlying LTR layer is overwhelmingly `Gypsy`-dominated, this is
  best understood as a Gypsy-heavy historical axis rather than a global TE-age
  clock for all TE classes.

## Practical Recommendation

- Keep `te_genome_ltr_history` as a sensitivity family.
- Use calibrated age ranges from
  `LTR_SUBSTITUTION_RATE_CALIBRATION.md`
  in the results text, but keep the path-model predictor itself as
  divergence-based `ltr_history`.

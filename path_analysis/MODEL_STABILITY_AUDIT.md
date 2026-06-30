# Model Stability Audit

## Purpose

This note audits whether the current primary medium-plus path-model winners
are robust to single-species removal and named-clade jackknife exclusion.

## Upstream inputs

- `path_analysis/data/derived/panels/*.csv`
- `path_analysis/results/*primary_mediumplus_model_ranking.csv`
- `path_analysis/results/*primary_mediumplus_best_model_edges.csv`
- `results/data/permanova/species_clade_assignments.csv`
- `input_data/phylogeny/desmo900dated_test.tre`

## Outputs

- `path_analysis/data/derived/model_stability_detail.csv`
- `path_analysis/data/derived/model_stability_summary.csv`

## te_genome_primary_mediumplus

- Baseline winner: `mediated_evenness`
- Winner retained across all runs: `31 / 31` (1.000)
- Leave-one-species-out retention: `27 / 27` (1.000)
- Clade-jackknife retention: `4 / 4` (1.000)
- Alternative winners observed: `0`
- Median winner margin (delta CICc to runner-up): `28.966`

Tracked edge stability:
- `ltr_balance -> te_evenness`: baseline `negative`, present in 1.000 of runs, same sign in 1.000 of runs
- `te_evenness -> gs`: baseline `negative`, present in 1.000 of runs, same sign in 0.903 of runs

## te_genome_organismal_primary_mediumplus

- Baseline winner: `body_size_additive`
- Winner retained across all runs: `16 / 31` (0.516)
- Leave-one-species-out retention: `14 / 27` (0.519)
- Clade-jackknife retention: `2 / 4` (0.500)
- Alternative winners observed: `2`
- Median winner margin (delta CICc to runner-up): `0.824`

Tracked edge stability:
- `ltr_balance -> te_evenness`: baseline `negative`, present in 1.000 of runs, same sign in 1.000 of runs
- `te_evenness -> gs`: baseline `negative`, present in 1.000 of runs, same sign in 0.935 of runs
- `body_size -> gs`: baseline `positive`, present in 0.516 of runs, same sign in 0.516 of runs

## te_genome_ectopic_organismal_primary_mediumplus

- Baseline winner: `te_body_size_baseline`
- Winner retained across all runs: `28 / 28` (1.000)
- Leave-one-species-out retention: `24 / 24` (1.000)
- Clade-jackknife retention: `4 / 4` (1.000)
- Alternative winners observed: `1`
- Median winner margin (delta CICc to runner-up): `5.542`

Tracked edge stability:
- `ltr_balance -> te_evenness`: baseline `negative`, present in 1.000 of runs, same sign in 1.000 of runs
- `te_evenness -> gs`: baseline `negative`, present in 1.000 of runs, same sign in 0.929 of runs
- `body_size -> gs`: baseline `positive`, present in 1.000 of runs, same sign in 1.000 of runs

## te_genome_ltr_history_primary_mediumplus

- Baseline winner: `history_additive`
- Winner retained across all runs: `26 / 28` (0.929)
- Leave-one-species-out retention: `23 / 24` (0.958)
- Clade-jackknife retention: `3 / 4` (0.750)
- Alternative winners observed: `2`
- Median winner margin (delta CICc to runner-up): `3.699`

Tracked edge stability:
- `ltr_balance -> te_evenness`: baseline `negative`, present in 1.000 of runs, same sign in 1.000 of runs
- `te_evenness -> gs`: baseline `positive`, present in 1.000 of runs, same sign in 0.750 of runs
- `ltr_history -> gs`: baseline `positive`, present in 0.929 of runs, same sign in 0.929 of runs

## Bottom line

Treat families as primary only if winner identity is broadly retained and
their baseline edge signs remain stable under both leave-one-out and named
clade jackknife perturbations.

# Phylogenetic Outlier Audit

## Purpose

This note identifies species that sit unusually high or low for key TE and
LTR-history variables relative to nearby phylogenetic expectation.

## Upstream inputs

- `input_data/phylogeny/desmo900dated_test.tre`
- `path_analysis/data/derived/path_input_master.csv`
- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`
- `results/data/permanova/species_clade_assignments.csv`

## Outputs

- `results/data/phylo_residuals/phylo_residual_summary.csv`
- `results/data/phylo_residuals/phylo_outlier_species.csv`

## Method

For each trait, ancestral states were reconstructed on the time-calibrated
tree with `phytools::fastAnc()`. Species-level outlier scores are terminal
branch deviations from the reconstructed parent state, scaled by branch
length and then z-scored across species within trait.

This is a branch-local outlier screen, not a claim about full causal
direction or adaptive explanation.

## Strongest branch-level deviations

- `all_recent_mass_frac_0_5`: positive `fuscus` (z = 2.975), negative `aeneus` (z = -2.495), lambda = 1.000
- `ectopic_log10_mean_ratio`: positive `intermedius` (z = 3.566), negative `mavrokoilius` (z = -2.351), lambda = 0.000
- `ltr_line_logratio`: positive `catahoula` (z = 2.546), negative `aeneus` (z = -1.533), lambda = 1.000
- `weighted_te_divergence_p90`: positive `aeneus` (z = 2.170), negative `fuscus` (z = -2.289), lambda = 1.000
- `all_old_tail_frac_20plus`: positive `aeneus` (z = 3.938), negative `santeetlah` (z = -1.719), lambda = 1.000
- `ltr_history_n_pairs_estimated`: positive `tilleyi` (z = 2.331), negative `auriculatus` (z = -1.780), lambda = 0.000
- `ltr_history_median_k2p_distance`: positive `aureatus` (z = 2.262), negative `conanti` (z = -1.745), lambda = 0.000
- `order_pielou`: positive `lycos` (z = 1.470), negative `catahoula` (z = -2.666), lambda = 1.000

## Recurrent outlier species

- `aeneus`: `3` traits
- `catahoula`: `2` traits
- `fuscus`: `2` traits
- `aureatus`: `1` traits
- `auriculatus`: `1` traits
- `cheaha`: `1` traits
- `intermedius`: `1` traits
- `mavrokoilius`: `1` traits
- `orestes`: `1` traits
- `tilleyi`: `1` traits

## Bottom line

The TE landscape and LTR-history layers are not only clade-level patterns.
They also contain species-level terminal deviations that can anchor concrete
biological discussion in the paper.

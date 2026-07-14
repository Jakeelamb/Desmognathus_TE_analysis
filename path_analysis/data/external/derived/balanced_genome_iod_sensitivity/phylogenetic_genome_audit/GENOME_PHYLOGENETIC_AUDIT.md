# Genome Phylogenetic Deviation Audit

## Purpose

This audit ranks species that most deserve image-level review because their
estimated genome sizes are far from `D. fuscus`, far from local
phylogenetic expectation, or sharply different from their nearest measured
relative on the time-calibrated tree.

## Inputs

- Tree: `/home/jake/Projects/Desmognathus_TE/input_data/phylogeny/desmo900dated_test.tre`
- Genome estimates: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/balanced_genome_iod_species_estimates.csv`
- Current-vs-balanced comparison: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/balanced_genome_iod_comparison.csv`

## Outputs

- Species audit table: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/phylogenetic_genome_audit/genome_phylo_species_audit.csv`
- Phylogenetic signal summary: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/phylogenetic_genome_audit/genome_phylo_signal_summary.csv`
- Branch jump table: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/phylogenetic_genome_audit/genome_phylo_branch_jumps.csv`
- Browser dashboard: `/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/balanced_genome_iod_sensitivity/phylogenetic_genome_audit/index.html`

## Method

Genome size was analyzed on the log scale. For each estimate panel,
`phytools::fastAnc()` reconstructed internal states on the dated tree.
Terminal outlier scores are observed tip values minus reconstructed parent
states, scaled by terminal branch length and then z-scored within panel.
Nearest-relative jumps compare each species to the closest measured species
by patristic distance. This is an audit-prioritization screen, not a final
claim that any species is biologically wrong.

## balanced_qc_curated signal snapshot

- Species: `21`
- Blomberg K: `0.500`; p = `0.477`
- Pagel lambda: `0.000`; p = `1.000`
- Spearman genome vs distance from D. fuscus: `0.109`
- Terminal |z| >= 2 species: `2`
- Species >= 50% from D. fuscus: `1`
- Species >= 25% from nearest measured relative: `3`

## Largest deviations from D. fuscus

species | genome_pg | vs_fuscus | terminal_z | nearest | nearest_jump | triggers
--- | --- | --- | --- | --- | --- | ---
D. valtos | 28.95 | 77.0% | 2.38 | D. bairdi | 18.6% | fuscus_delta_ge_50pct; terminal_phylo_z_ge_2; method_shift_ge_40pct
D. bairdi | 24.42 | 49.3% | 0.24 | D. valtos | -15.7% | method_shift_ge_40pct; low_independent_support
D. brimleyorum | 23.66 | 44.6% | 1.25 | D. ocoee | 9.2% | method_shift_ge_40pct
D. intermedius | 23.18 | 41.7% | 1.94 | D. marmoratus | 21.5% | method_shift_ge_40pct
D. kanawha | 22.46 | 37.3% | 1.67 | D. mavrokoilius | 31.7% | nearest_relative_jump_ge_25pct; method_shift_ge_40pct
D. ocoee | 21.65 | 32.4% | 0.87 | D. apalachicolae | 7.6% | method_shift_ge_40pct
D. tilleyi | 21.63 | 32.2% | 0.56 | D. anicetus | 8.9% | method_shift_ge_40pct
D. folkertsi | 21.60 | 32.0% | 0.20 | D. amphileucus | 0.4% | method_shift_ge_40pct
D. amphileucus | 21.51 | 31.5% | 0.20 | D. gvnigeusgwotli | 5.3% | method_shift_ge_40pct
D. gvnigeusgwotli | 20.43 | 24.9% | -0.41 | D. amphileucus | -5.1% | method_shift_ge_40pct; low_independent_support

## Largest terminal phylogenetic residuals

species | genome_pg | terminal_z | parent_expected_pg | terminal_shift | triggers
--- | --- | --- | --- | --- | ---
D. mavrokoilius | 17.05 | -2.65 | 20.11 | -15.2% | terminal_phylo_z_ge_2
D. valtos | 28.95 | 2.38 | 23.85 | 21.4% | fuscus_delta_ge_50pct; terminal_phylo_z_ge_2; method_shift_ge_40pct
D. fuscus | 16.36 | -1.98 | 20.06 | -18.4% | nearest_relative_jump_ge_25pct
D. intermedius | 23.18 | 1.94 | 20.51 | 13.0% | method_shift_ge_40pct
D. kanawha | 22.46 | 1.67 | 20.11 | 11.7% | nearest_relative_jump_ge_25pct; method_shift_ge_40pct
D. welteri | 17.07 | -1.48 | 19.92 | -14.3% | method_shift_ge_40pct
D. brimleyorum | 23.66 | 1.25 | 20.26 | 16.8% | method_shift_ge_40pct
D. marmoratus | 19.07 | -1.25 | 20.51 | -7.0% | none
D. ochrophaeus | 17.96 | -0.97 | 19.92 | -9.8% | nearest_relative_jump_ge_25pct; method_shift_ge_40pct
D. ocoee | 21.65 | 0.87 | 20.28 | 6.8% | method_shift_ge_40pct

## Largest nearest-relative jumps

species | genome_pg | nearest | nearest_genome_pg | nearest_jump | distance_mya
--- | --- | --- | --- | --- | ---
D. fuscus | 16.36 | D. bairdi | 24.42 | -33.0% | 15.20
D. kanawha | 22.46 | D. mavrokoilius | 17.05 | 31.7% | 5.45
D. ochrophaeus | 17.96 | D. bairdi | 24.42 | -26.5% | 17.43
D. mavrokoilius | 17.05 | D. kanawha | 22.46 | -24.1% | 5.45
D. intermedius | 23.18 | D. marmoratus | 19.07 | 21.5% | 5.01
D. auriculatus | 19.63 | D. bairdi | 24.42 | -19.6% | 14.51
D. valtos | 28.95 | D. bairdi | 24.42 | 18.6% | 8.43
D. marmoratus | 19.07 | D. intermedius | 23.18 | -17.7% | 5.01
D. bairdi | 24.42 | D. valtos | 28.95 | -15.7% | 8.43
D. perlapsus | 18.66 | D. ocoee | 21.65 | -13.8% | 11.51

## Largest current-to-primary-panel method shifts

species | current_pg | primary_panel_pg | shift_vs_current | shift_vs_original_balanced | support | triggers
--- | --- | --- | --- | --- | --- | ---
D. perlapsus | 10.73 | 18.66 | 74.0% | 0.0% | medium | method_shift_ge_40pct
D. kanawha | 13.25 | 22.46 | 69.5% | 4.5% | medium | nearest_relative_jump_ge_25pct; method_shift_ge_40pct
D. ochrophaeus | 11.06 | 17.96 | 62.4% | 0.0% | limited | nearest_relative_jump_ge_25pct; method_shift_ge_40pct
D. apalachicolae | 12.69 | 20.12 | 58.6% | 0.0% | medium | method_shift_ge_40pct
D. brimleyorum | 15.08 | 23.66 | 56.9% | 0.0% | medium | method_shift_ge_40pct
D. folkertsi | 13.79 | 21.60 | 56.6% | 0.0% | medium | method_shift_ge_40pct
D. tilleyi | 14.03 | 21.63 | 54.1% | 0.0% | medium | method_shift_ge_40pct
D. anicetus | 13.07 | 19.87 | 52.0% | 0.0% | low | method_shift_ge_40pct; low_independent_support
D. amphileucus | 14.28 | 21.51 | 50.7% | 0.0% | medium | method_shift_ge_40pct
D. bairdi | 16.29 | 24.42 | 49.9% | 0.0% | low | method_shift_ge_40pct; low_independent_support

## Figures

- `figures/genome_vs_fuscus_distance.png`
- `figures/terminal_phylo_outlier_rank.png`
- `figures/nearest_relative_jump_rank.png`
- `figures/tree_genome_audit.png`

## Audit rule

Prioritize species with multiple triggers, especially if the same species is
far from `D. fuscus`, has a high terminal phylogenetic residual, and also
shows a large current-to-balanced method shift. Those species should be
checked in the mask viewer before their genome-size values are treated as
biological signal.

# TE Age Spectrum Audit

## Purpose

This note audits whether the local RepeatMasker landscapes support a stable,
interpretable species-level TE tempo axis before the final cell/nucleus/genome
size layers are ready.

## Upstream inputs

- `results/landscapes/repeat_landscape_*.csv`
- `path_analysis/data/derived/te_path_features.csv`
- `path_analysis/data/derived/ltr_history_features.csv`

## Outputs

- `results/data/te_age_spectra/te_age_spectrum_metrics.csv`
- `results/data/te_age_spectra/te_age_spectrum_correlations.csv`

## Metric definition

Each landscape was collapsed to species-level shape summaries for `ALL`, `LTR`,
`LINE`, and `TIR` bins:

- recent mass: Kimura bins `< 5%`
- mid mass: Kimura bins `5-15%`
- old tail: Kimura bins `>= 20%`
- ancient tail: Kimura bins `>= 30%`
- recent/old logratio
- peak bin
- weighted mean and median bin
- peak mass fraction
- entropy-derived effective bin count

These are tempo descriptors, not direct insertion-age estimates.

## Coverage

- Landscapes processed: `34` species
- Species overlapping the paired-LTR history layer: `30`

## Main result

The landscape branch passes the basic interpretability screen. It recovers a
clear recent-vs-old species axis, and that axis strongly tracks the existing
RepeatMasker divergence summaries rather than behaving like random noise.

The strongest supported relationships are:

- `all_recent_mass_frac_0_5` vs `weighted_te_divergence_p90`: rho = -0.897, BH p = 1.71e-11, n = 34
- `all_weighted_mean_bin_pct` vs `weighted_te_divergence_p90`: rho = 0.871, BH p = 1.3e-10, n = 34
- `ltr_recent_mass_frac_0_5` vs `ltr_divergence_p90`: rho = -0.889, BH p = 3.05e-11, n = 34
- `line_weighted_mean_bin_pct` vs `line_divergence_p90`: rho = 0.865, BH p = 2.12e-10, n = 34
- `all_effective_bin_count` vs `order_pielou`: rho = 0.679, BH p = 3.02e-05, n = 34
- `all_recent_mass_frac_0_5` vs `ltr_history_median_k2p_distance`: rho = -0.008, BH p = 0.967, n = 30

## Interpretation

- The clearest stable tempo axis is `recent` versus `old retained tail`.
- That axis is strongest for the full TE landscape and for `LTR`/`LINE`
  sub-landscapes.
- The landscape-derived tempo summaries line up strongly with the existing
  divergence summaries, which is expected and validates the new metrics rather
  than replacing the existing divergence layer.
- The paired-LTR history layer is comparatively independent: it does not track
  these coarse landscape summaries strongly enough to treat them as the same
  signal.

## Species extremes

Most recent all-TE landscapes by `all_recent_mass_frac_0_5`:
- `fuscus`: all_recent_mass_frac_0_5=0.125, old_tail=0.453, mean_bin=18.65
- `orestes`: all_recent_mass_frac_0_5=0.116, old_tail=0.452, mean_bin=18.70
- `bairdi`: all_recent_mass_frac_0_5=0.110, old_tail=0.460, mean_bin=19.01
- `valtos`: all_recent_mass_frac_0_5=0.109, old_tail=0.457, mean_bin=18.88
- `santeetlah`: all_recent_mass_frac_0_5=0.107, old_tail=0.446, mean_bin=18.61

Oldest all-TE landscapes by `all_recent_mass_frac_0_5`:
- `organi`: all_recent_mass_frac_0_5=0.040, old_tail=0.652, mean_bin=24.12
- `wrighti`: all_recent_mass_frac_0_5=0.040, old_tail=0.659, mean_bin=24.37
- `aeneus`: all_recent_mass_frac_0_5=0.052, old_tail=0.584, mean_bin=22.29
- `aureatus`: all_recent_mass_frac_0_5=0.064, old_tail=0.530, mean_bin=20.93
- `gvnigeusgwotli`: all_recent_mass_frac_0_5=0.069, old_tail=0.521, mean_bin=20.74

## Bottom line

The TE landscape data support a usable species-level tempo axis. The main
descriptive contrast is not "single recent burst everywhere" but variation in
how much recent mass versus old retained tail each species carries. This is
worth keeping as a comparative layer, but it should be described as a
landscape-derived tempo summary rather than a substitute for paired-LTR age
estimation.

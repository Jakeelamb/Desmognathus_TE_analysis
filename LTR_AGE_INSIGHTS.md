# LTR Age Calibration Insights

This note interprets the calibrated paired-LTR age layer using the tracked
sequence-divergence outputs and the primary-literature substitution-rate set.

## Source Basis

- Primary substitution-rate sources:
  - Herrick and Sclavi 2014: https://doi.org/10.1016/j.crpv.2014.06.002
  - Crawford 2003: https://doi.org/10.1007/s00239-003-2513-7
- Canonical divergence workflow:
  - `results/reports/LTR_AGE_AUDIT.md`
  - `results/data/ltr_age/ltr_age_pairwise_divergence.csv`
  - `results/data/ltr_age/ltr_age_species_summary.csv`
- Calibration layer:
  - `LTR_SUBSTITUTION_RATE_CALIBRATION.md`
  - `results/data/ltr_age/ltr_age_calibration_summary.csv`
  - `results/data/ltr_age/ltr_age_species_summary_calibrated.csv`
- Signal/context checks:
  - `results/reports/LTR_DIVERGENCE_SIGNAL_AUDIT.md`
  - `results/data/ltr_age/ltr_divergence_phylogenetic_signal.csv`
  - `results/data/ltr_age/ltr_divergence_feature_correlations.csv`

## Main Stable Findings

- The defensible Desmognathus sensitivity window is `1.00e-9` to `1.53e-9`
  substitutions/site/year, with `1.227e-9` as the central representative rate.
- Across all successful paired-LTR estimates, the median inferred age stays old
  across that entire window: `35.282` Myr at the salamander lower anchor,
  `28.754` Myr at the central rate, and `23.060` Myr at the upper bound.
- The recommended high-confidence subset is slightly younger overall, but not
  qualitatively different: `33.101`, `26.977`, and `21.634` Myr across the same
  low, central, and high rates. At the central rate this is about `6.18%`
  younger than the all-pairs median.
- The central-rate age distribution is broad rather than concentrated in a
  narrow recent burst window: the pairwise `P10` is `11.416` Myr and the `P90`
  is `46.544` Myr.
- Because age is a simple linear rescaling of K2P distance, the substitution
  rate changes the absolute time axis but not the relative ordering of species
  or elements. The rate choice therefore affects calibration scale, not which
  species look older or younger.

## Biological Interpretation

- The paired-LTR layer does not support a story of predominantly recent LTR
  insertions. Even the fastest recommended calibration still puts the median
  paired-LTR age at about `23` Myr.
- The age signal is dominated by `Gypsy` elements, because `1283 / 1291`
  successful sequence-based estimates come from that superfamily. The current
  absolute-age interpretation is therefore strongest for Gypsy-rich LTR history,
  not for LTRs broadly in proportion to all superfamilies.
- Species differences are real but modest compared with the full calibration
  window. At the central rate, the oldest species-level median is `D.aureatus`
  at `36.591` Myr and the youngest is `D.aeneus` at `21.866` Myr, a spread of
  `14.726` Myr.
- The oldest central-rate median species are `D.aureatus`, `D.ocoee`, and
  `D.monticola`. The youngest are `D.aeneus`, `D.conanti`, and `D.balsameus`.
- The high-confidence subset can shift some species medians sharply in either
  direction, so species-level high-confidence summaries should be treated
  cautiously when the number of high-confidence pairs is small.

## What The Calibration Does Not Support

- It does not support claiming a single measured neutral substitution rate for
  Desmognathus. The rates are borrowed from salamander and frog nuclear
  literature, so age estimates must remain sensitivity-based.
- It does not support a strong phylogenetically conserved age axis. Species-level
  median K2P divergence has weak phylogenetic signal (`Blomberg's K = 0.359`,
  `p = 0.464`; `Pagel's lambda = 0.000073`, `p = 1`) in
  `results/data/ltr_age/ltr_divergence_phylogenetic_signal.csv`.
- It does not support a simple mapping between LTR age and the current TE/path
  predictor set. The strongest tested association was only a weak positive trend
  with ectopic ratio (`rho = 0.291`, `p = 0.119`), while the remaining tested
  correlations were weaker.

## Working Interpretation

- The safest interpretation is that much of the detectable paired-LTR signal in
  Desmognathus reflects relatively old, lineage-specific LTR history rather than
  a uniform recent burst shared across clades.
- The age scale is best reported as a sensitivity range rather than a point
  estimate, with the central `1.227e-9` value used only as a representative
  summary in tables or figures.
- Any results text should explicitly say that the calibrated ages are derived
  from sequence-based 5'/3' LTR divergence combined with externally sourced
  nuclear substitution rates.

# LTR Substitution-Rate Calibration

This note records the primary-literature substitution-rate calibration used to convert
paired-LTR K2P divergence into absolute age sensitivities.

## Bottom Line

- The recommended salamander-compatible calibration window is `1.00e-09` to `1.53e-09` substitutions/site/year.
- The recommended central sensitivity value is `1.227e-09` substitutions/site/year (`frog_cmyc_midpoint`).
- Across all successful paired-LTR estimates, the median inferred age is `35.282` Myr at the lower bound, `28.754` Myr at the central value, and `23.060` Myr at the upper bound.
- Across the recommended high-confidence subset, the corresponding medians are `33.101`, `26.977`, and `21.634` Myr.
- The faster frog tyrosinase estimate at `1.69e-09` is retained only as exploratory sensitivity, while the very fast `3.35e-09` tyrosinase estimate is excluded from the main Desmognathus sensitivity set.

## Primary Sources

- Herrick and Sclavi (2014) used salamander `rag1` synonymous substitution rates and reported that Plethodontidae rates bottom out at `0.001 dS/Mya`, which maps to `1.00e-9` substitutions/site/year. This is the best lineage-matched lower anchor.
- Crawford (2003) reported frog nuclear synonymous-site rates of `0.924e-9` to `1.53e-9` for `c-myc`, `1.03e-9` for `slug`, and `1.69e-9` to `3.35e-9` for `tyrosinase`.
- Crawford also explicitly noted that the higher ranid tyrosinase estimates could be biased upward if the calibration underestimated divergence ages, so they are not used as the main Desmognathus defaults.

## Recommended Use

- Main analysis text: report LTR ages as a sensitivity range using `1.00e-09` to `1.53e-09` substitutions/site/year.
- Supplementary tables/figures: use `frog_cmyc_midpoint` (`1.227e-09`) as the representative central estimate.
- Do not present a single absolute age as if the substitution rate were directly measured in Desmognathus.

## Central-Rate Species Extremes

- Oldest median species at the central rate: `D.aureatus` = `36.591` Myr
- Oldest median species at the central rate: `D.ocoee` = `33.318` Myr
- Oldest median species at the central rate: `D.monticola` = `32.293` Myr
- Youngest median species at the central rate: `D.aeneus` = `21.866` Myr
- Youngest median species at the central rate: `D.conanti` = `23.431` Myr
- Youngest median species at the central rate: `D.balsameus` = `26.230` Myr

## Files

- `results/data/ltr_age/ltr_substitution_rate_candidates.csv`
- `results/data/ltr_age/ltr_age_calibration_summary.csv`
- `results/data/ltr_age/ltr_age_species_summary_calibrated.csv`
- `results/data/ltr_age/ltr_age_pairwise_divergence_calibrated.csv`

## Source Links

- Herrick and Sclavi 2014: https://doi.org/10.1016/j.crpv.2014.06.002
- Crawford 2003: https://doi.org/10.1007/s00239-003-2513-7

# Figure 10: nuclear–cell scaling and paired allocation

Exact frozen inputs are copied with SHA-256 hashes to inputs/nuclear_allometry. Identity joins use the species crosswalk. No prior figure or frozen release is altered.

Primary areas are released pooled species medians. Equal-animal alternative: arithmetic mean of within-animal medians, separately for each trait. Nuclear fraction is computed for each actual pair, then summarized by animal median and equal-animal species mean; it is not a ratio of aggregate medians.

phylolm fits log10(nucleus area) ~ log10(cell area) by ML with Pagel lambda in [1e-7,1]. Slopes are unstandardized. Isometry test: t=(slope-1)/SE with n-2 degrees of freedom. Model intervals use t critical values conditional on fitted lambda. Four fixed sensitivities use released/equal-animal summaries and all24/excluding aeneus+wrighti (22).

Seed20260924; 1,000 draws resample animals within species, then paired cells within each sampled animal. Equal-animal means of sampled animal medians are fitted across all24 species, re-estimating lambda. All failures and warnings are retained; every fit must succeed before export. This is selected-measurement sensitivity, not a full errors-in-variables model or combined model confidence interval. Seven single-animal species cannot identify between-animal uncertainty.

Conditional response-on-predictor regression is directional; predictor error and axis choice can affect the exponent. Do not interpret the slope as a demonstrated mechanistic scaling law. Nuclear fraction is descriptive and is not tested against its own component variables as an independent outcome.

Outputs: allometry_model_summary.csv has four models and slope=0/slope=1 tests. paired_animal_summaries.csv and equal_animal_species_data.csv contain plotted data. paired_hierarchical_refits.csv contains all fits. paired_hierarchical_traits.csv.gz stores all24,000 draws. displayed_slope_intervals.csv distinguishes t confidence intervals from selected-measurement ranges.

Primary conditional slope 0.4500 [0.1701,0.7300], P(slope=1) 0.000504209; lambda 0.0000001.
Equal-animal full24 slope 0.4832; bootstrap median 0.4650, central95% [0.3056,0.6331].

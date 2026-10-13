# Figure 5 refinement

## Scope and source mapping

This is a display-only refinement in the isolated polished run. No observations were excluded and no regressions or resampling procedures were rerun.

- New 5a: prior Figure 5a, genome size versus cell area. Frozen `plotted_species_data.csv`, `model_summary.csv` and `fitted_relationships.csv` supply all observations, coefficients, intervals and predictions.
- New 5b: prior Figure 5b, genome size versus nucleus area, using the same frozen exports.
- New 5c: prior Figure 10a, conditional nucleus-on-cell area allometry. Frozen `released_species_data.csv`, `allometry_model_summary.csv` and `fitted_relationship.csv` supply the observations, regression, confidence band and slope-one reference. This replaces the redundant reverse-regression association display in prior 5c without refitting or changing either original result.
- New 5d: prior Figure 10b, paired nuclear fractions. Frozen animal and equal-animal species summaries and resampling intervals retain all 24 species, 42 animals and 1,152 accepted paired cells. Animal-count labels are present in the larger isolated export. Complete support counts remain in its data export.

## Visual and scientific decisions

- A four-panel 180 × 240 mm collage integrates genome relationships with a distinct nuclear-allocation result. Isolated scatterplots label every species; overview labels identify fixed reference species only. Every observation and uncertainty interval remains plotted.
- Shared Arial styling has an explicit 7.5 pt minimum. Reproductive strategies use the common color and shape mapping. The three species without a frozen egg-placement assignment remain gray diamonds; their known aquatic-larval development is not converted into an invented egg-placement category.
- Paired-fraction species follow the common phylogenetic display order. Open circles denote animal medians and diamonds denote species arithmetic means of animal medians, independently of the strategy shape mapping in the scatterplots. The caption makes the distinction explicit.
- Genome labels say “Genome size” as requested. The legend retains the IOD-derived, assembly-anchored measurement definition and unresolved calibration limitation.
- The allometry caption distinguishes a conditional regression slope from an errors-in-variables structural exponent and distinguishes two-dimensional areas from volumes. The slope-one line is a proportional reference, not a competing newly fitted model.
- The paired fraction is computed on actual cell/nucleus pairs. It is not a ratio of species medians and is not offered as an independent causal variable.
- Earlier review notes requested complete species coverage, tighter informative ranges and clearer genome-size labels; all are retained here. No heatmap was introduced.

## Validation

The script checks 24 species in every scatterplot and allocation panel, 42 animal summaries, 1,152 paired cells, three explicit missing strategy assignments, and identity of the nuclear-allometry and genome-association morphology point data. The two model pipelines retain their original standardization and interval definitions. Raster, vector and CMYK TIFF outputs are exported by the shared publication helper. Visual inspection is recorded in `qa/figures_5_7_visual_checks.md`.

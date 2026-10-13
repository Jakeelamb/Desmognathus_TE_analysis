# Figure 7 refinement

## Scope and source mapping

This is a display-only refinement in the isolated polished run. All fitted lines, coefficients and intervals are read directly from frozen exports.

- New 7a–c: prior 7b–d, species maximum SVL versus genome size, cell area and nucleus area. Frozen prior Figure 7 observations, regression summaries and prediction files supply the full24 and subset22 fits.
- New 7d: prior Figure 9c, the 12 conditional coefficients comparing released/equal-animal summaries in full24/subset22 panels. The three equations and focal predictors are stated without causal arrows.
- The prior body-size-by-species dot plot is not repeated in this main collage. Body size remains plotted directly in every scatterplot and retained numerically for all 24 species. The original figure remains intact in the earlier run.

## Visual and scientific decisions

- A four-panel 180 × 235 mm layout makes the comparison between body relationships and cellular robustness visible together. Three scatterplots share the same body-size breaks and log scale.
- All 24 species remain in every scatterplot. Full24 lines are solid; subset22 lines are dashed and are restricted to their own observed body-size range. Subset22 excludes exactly D. aeneus and D. wrighti. The full24 uncertainty ribbon is displayed and both models' complete interval statistics remain exported.
- Shared Arial styling has an explicit 7.5 pt minimum. Reproductive strategy colors and shapes are consistent with other polished figures; missing egg-placement assignments remain visibly unassigned. Isolated scatterplots label all 24 species.
- Coefficients use the frozen full24 log10 means and SDs. The coefficient forest contains conditional regressions, not claims of identified directional effects. Formula headings and a focal-predictor subtitle identify the quantity plotted.
- The caption distinguishes species maximum SVL from the size of a microscopy animal. It also explains that the two-species exclusion jointly changes body size, developmental mode and lineage representation; changes in significance do not test a difference between slopes or isolate causation.
- The two equal-animal definitions are deliberately distinguished: Figure 7d morphology uses the median of animal medians, while Figure 5d nuclear fractions use the mean of animal medians. Neither is silently substituted for the other.
- Sampling robustness does not resolve staining calibration or the fact that erythrocytes were selected for large area. Those limitations remain in the full legend.

## Validation

The script verifies 24 species, positive frozen predictions, and exactly 12 displayed coefficients: two aggregations × two cohorts × three conditional relations. Existing uncertainty bands, P values and sample counts are unchanged. Raster, vector and CMYK TIFF outputs use the common publication helper. Visual inspection is recorded in `qa/figures_5_7_visual_checks.md`.

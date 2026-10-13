# Animal-aware robustness of the body-adjusted path analysis

Exploratory sensitivity analysis; frozen objects and prior figures are preserved. These summaries concern the reviewed selected large-cell morphology panel and quality-matched IOD panel, not random samples of all blood cells.

## Estimands
Released morphology: pooled within-species medians. Equal-animal morphology: median across per-animal pooled-object medians for each trait. Cell/nucleus traits always use paired rows. This changes the species estimand and is a sensitivity comparison, not a correction to a known truth.
Released IOD: equal-image mean of within-image nuclear-IOD medians. Equal-animal IOD: average image medians within each animal, then average animal summaries. Genome size uses the same 16.1/0.978 pg/1C D. fuscus anchor, recomputing the common denominator per aggregation/draw; anchor and staining-calibration uncertainty remain excluded.

## Fits
All traits log10 transformed and standardized using the original released full24 means/SDs. Pagel-lambda PGLS on the frozen Path24 dated tree. Individual fitted coefficients use 95% t intervals with n minus regression coefficient count degrees of freedom. Retain exactly the original four body-adjusted candidates, all assuming SVL upstream. Weights are relative to this set, not identified causal directions. Fit full24 and the existing 22-species sensitivity excluding aeneus+wrighti under both aggregations.

## Hierarchical measurement sensitivity
250 draws fixed in advance; seed 2026092409. Within each species and assay, sample observed animals with replacement, sample each selected animal's images with replacement, then sample rows within each selected image with replacement. Retain original sample size at every level. Morphology cell/nucleus rows stay paired. Combine per-animal morphology medians by median; IOD image medians by mean within animals then mean across animals.
Assays are resampled independently even when animal IDs overlap. Cross-assay covariance, including shared objects, is unmodeled. Ranges describe conditional measurement sensitivity rather than full uncertainty. Single-animal species cannot estimate between-animal variability; resampling their images/objects does not manufacture replication. Maximum SVL is fixed without measurement uncertainty.
Every measurement draw retains 24 species. Separately, 24 leave-one-species-out fits use released summaries, fixed scaling, and the pruned dated tree. This assesses species leverage rather than specimen uncertainty. No species is deleted from the primary analysis based on results.

## Verification
Checks: exact identity joins, complete positive traits, animal/image nesting, support totals, source-data replay, Figure8 coefficient/ranking reproduction, candidate-weight normalization, finite CICc, fixed fuscus anchor, 250 complete draws, 24 complete deletion fits, source hashes. Warnings/failures retained in CSVs.

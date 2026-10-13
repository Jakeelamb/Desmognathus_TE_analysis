# Figure 6 refinement

This display-only consolidation reads frozen files under inputs/source_run. It does not fit, resample, select new species or edit the source run. Original Figures 6, 8 and 13 remain available. The original three-trait comparison supplies context, but its weights and robustness rates are not transplanted into the new eight-class comparison.

## Panel mapping

| Panel | Frozen source and rows | Decision |
|---|---|---|
| a: weighted H5 hypothesis | Figure 13 analysis/structural_coefficients.csv; full24, cell_body_branch, non-intercept rows | Show a plausible genome-first orientation using its actual three structural coefficients; direction is assumed. |
| b: conditional regressions | Figure 8 analysis/regression_coefficients.csv; full24, the relations below | Keep adjusted evidence separate from panel a because these equations have different conditioning sets. |
| c: eight-class comparison | Figure 13 analysis/model_comparison.csv; primary test order, both cohorts; candidate_edges.csv | Every candidate is shown with skeleton and weight; full24 and subset22 weights are separately normalized. |
| d: causal equivalence | Figure 13 analysis/equivalent_orientations.csv and equivalence_classes.csv; cell_body_branch | All four verified H5 orientations are drawn, using the complete four-variable graph. |

## Exact coefficient origins

All coefficients use log10 traits standardized with fixed full24 means and sample SDs.

| Display | Source relation and term | Estimate | 95% t interval | Lambda |
|---|---|---:|---|---:|
| a: G → C | Figure 13 H5, c ~ g, term g | 0.7166644026006085 | 0.408300113950521 to 1.025028691250696 | 0.0000001 |
| a: C → N | Figure 13 H5, n ~ c, term c | 0.579301282090048 | 0.218897830928528 to 0.939704733251568 | 0.0000001 |
| a: C → B | Figure 13 H5, b ~ c, term c | 0.3614728141580796 | 0.094502320913937 to 0.628443307402223 | 1 |
| b: G → C given B | Figure 8 genome_to_cell_body_adjusted, term g | 0.5453794365491734 | 0.191994729375368 to 0.898764143722979 | 0.0000001 |
| b: C → N given G+B | Figure 8 genome_and_cell_to_nucleus_body_adjusted, term c | 0.6523756966402472 | 0.074373167056630 to 1.230378226223865 | 0.0000001 |
| b: G → N given C+B | Figure 8 same relation, term g | 0.05058352519534526 | -0.498921826900186 to 0.600088877290877 | 0.0000001 |

C → B is 0.36, not the coefficient of the reversed B → C equation. Panel-b C → N is 0.65, not the 0.68 coefficient from the equation omitting genome. Full-precision coefficients, intervals, conditioning sets and coordinates are exported under display/.

Full24 H5/H7 weights are 0.5240300752926639 and 0.3836081092600889, separated by 0.6238551911985937 CICc. Subset22 H8/H7/H5 weights are 0.3873962849634908, 0.29540152926361846 and 0.2834157322571316. All 16 primary-order weights are exported under display/.

## Scientific boundaries

- H5 is representative, not uniquely supported; H7 is close. Arrows are one plausible orientation.
- H1's old four-model weight of 0.944 is not shown beside the expanded comparison; its expanded full24 weight is approximately 0.065.
- Panel-a arrow width is proportional to the absolute coefficient. Panel-d widths are schematic.
- Coefficient intervals condition on fitted lambda and omit measurement, tree and covariance-parameter uncertainty.
- Restricted-model robustness rates are not evidence validating the expanded eight-class comparison.
- Removing aeneus and wrighti jointly changes body-size range, development and taxon composition.
- Graphical equivalence need not give numerically identical scores after reversing regressions and re-estimating covariance.
- The grey genome–nucleus interval is not evidence of negligible causal effect or complete mediation.
- Staining calibration and upper-tail erythrocyte selection remain interpretation limits.

## Publication design and verification

- Collage is 180 × 235 mm, with all four isolated panels.
- Arial throughout; smallest declared text 7.5 pt. PDF embeds ArialMT and Arial-BoldMT.
- Vector PDF, 300 dpi PNG and 300 dpi 8-bit CMYK TIFF.
- Candidate skeletons use a diamond arrangement so a genome–nucleus edge cannot cross the cell node and imply false edges.
- Visually inspected collage and all four isolated PNGs; text and coefficient labels are legible and uncropped.
- R assertions verify three H5 edges, three adjusted coefficients, all 16 weights, normalization in both cohorts, and all four H5 roots/orientations.
- Build: run scripts/figure_6.R from the isolated run with /home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ and that environment's bin directory on PATH.
- pdffonts confirmed embedded Arial. pdfinfo confirmed one page of approximately 180 × 235 mm. ImageMagick confirmed the 2125 × 2775 CMYK TIFF.
- Prior notes asking for Genome size instead of Conditional genome size are respected. Conditional associations refers to adjusted regression, not the genome-size variable.


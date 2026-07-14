# Microscopy and image-IOD release audit — final 18 species, v1

## Release verdict

| Layer | Verdict | Approved interpretation |
| --- | --- | --- |
| Object trace and cell–nucleus linkage | **Approved with frozen-source provenance** | Exact source image, tile, cell mask, nucleus mask, object IDs, and one-to-one geometry are retained. |
| Segmentation accuracy | **Not approved as publication-grade validation** | The current two-tile check is diagnostic only. |
| Cell and nucleus morphology | **Sensitivity-only** | Median area of the 50 highest composite-ranked eligible linked erythrocytes per species; not typical cell size and not literal largest 50. |
| Nuclear IOD | **Sensitivity-only** | Relative nuclear IOD index, with fuscus set to 1 within each QC subset. |
| Absolute genome size from images | **Not approved** | No picogram/C-value claim enters the corrected release. |

This audit does not delete or overwrite the historical 21-species top-50 analysis.  It builds a fail-closed final-18 surface for the eventual path-analysis sensitivity suite.

## 1. Segmentation validation

The current manually labeled ground truth contains 12 large tiles (10 local-train, 2 test) from only *D. ochrophaeus* and *D. folkertsi*.  Neither species belongs to the focal 18-species analysis.  For the locally trained nucleus model, the two test tiles come from the same source images and species used for training; that split tests new tiles within two images, not new images or specimens.  The production cell model was externally pretrained `cpsam` and did not use the ten local training tiles, but its local evaluation is still limited to the same two test images and two excluded species.

The production run plan and March 7 run manifest show that the frozen cell masks were generated with the default Cellpose model (`cpsam`): no custom `--cellpose-model` argument was supplied.  The Desmognathus fine-tuned checkpoint was created on March 8 and is not the source of the frozen cell masks.  This audit scores the exact archived production masks, including their hashes and March 7 modification times.  The cached `cpsam` binary predates the run, but the original run did not store its model hash, so model-binary identity remains reconstructive rather than contemporaneously frozen.

- **Cell:** foreground IoU 0.634, instance F1@0.50 0.624, F1@0.75 0.613; object-count bias +62.3%, foreground-area bias +13.3%; 2 test tiles, 0/18 focal species represented.
- **Nucleus:** foreground IoU 0.720, instance F1@0.50 0.687, F1@0.75 0.656; object-count bias -39.5%, foreground-area bias -4.2%; 2 test tiles, 0/18 focal species represented.

The nucleus result is materially lower than the YOLO training dashboard's patch-validation mAP and is the more relevant full-tile diagnostic.  It still cannot estimate error in the focal panel.  The publication gate remains: annotate an image/specimen-held-out, taxonomically and technically stratified set; report mask overlap, boundary/area bias, split/merge rates, pairing error, and downstream species-summary bias.

The nucleus test manifest records 138 masks across the two rows, but the label arrays actually scored contain 243 instances.  Thus its per-row `n_masks` field is not a reliable nucleus-object count and must not be used as the validation denominator.

Figure: [`results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png)

## 2. What the frozen top 50 actually estimate

The frozen final-18 table has 900 linked pairs.  Of these, 324 are explicit manual keeps, 2 are manual maybes, and 574 are model-ranked rows without an individual keep decision.  All frozen rows satisfy the recorded one-to-one/physical linkage gates, but a grid-selected model rank is not equivalent to an independent manual mask validation.

The ranking is composite: pair-quality score, model keep probability, the earlier top-50 score, cell area, and nucleus darkness all contribute.  Consequently:

- only 17–38 of each species' 50 frozen rows are also in its literal largest 50;
- the frozen rows occupy different cell-area quantiles across species;
- candidate-pool size correlates strongly with selection depth (Spearman rho = 0.750); and
- three species have one contributing image/specimen: anicetus, bairdi, gvnigeusgwotli.

The honest label is **median area of the 50 highest composite-ranked eligible linked erythrocytes**, or more briefly **quality-screened upper-tail morphology**.  It is not an unbiased estimate of species-average or species-median erythrocyte morphology.  A recent salamander hematology study sampled 50 erythrocytes randomly per individual, while the closest genome/cell/nucleus comparative study used 50 nuclei per individual and species medians across 1–4 individuals ([Liu et al. 2023](https://doi.org/10.7717/peerj.15446); [Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)).

Species ranks are fairly but not perfectly stable against a literal-largest-100 sensitivity: cell-area rho = 0.961; nucleus-area rho = 0.973.  The full estimator table, not one preferred summary, must enter the downstream sensitivity analysis.

Figures:

- [`results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png)
- [`results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png)
- [`results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png)
- [`results/figures/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.png)

## 3. Uncertainty and replication

The hierarchical bootstrap resamples specimens and then cells within specimen, but it remains conditional on the frozen selection and segmentation.  For a one-specimen species it is explicitly labeled `conditional_object_only_single_specimen`; it is not a population-level species interval.  None of these intervals include stain batch, acquisition, segmentation-model, reviewer, or selection-rule uncertainty.

## 4. Nuclear IOD is not currently an absolute genome-size assay

The pixel equation is implemented correctly: object IOD equals `sum(log10(I_bg / I_pixel))`, and in the exported object table `IOD = area_px × mean_OD` to numerical precision (maximum relative identity error 2.22e-16).  That is necessary but not sufficient for genome-size densitometry.

The repository does not provide a slide-linked record of fixation, acid hydrolysis, Schiff/Feulgen batch, staining time, co-processed DNA standard, microscope/camera settings, exposure/linearity/saturation tests, or the provenance and 1C/2C interpretation of the historical `fuscus = 16.36 pg` constant.  The exact methods benchmark requires controlled Feulgen staining, optical-density conversion, and a same-batch/same-slide standard ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)).  The closest *Desmognathus* precedent used *Xenopus laevis* erythrocyte nuclei as an internal standard ([Sessions and Kezer 1986](https://doi.org/10.1007/BF00494802)).

The current [Animal Genome Size Database *D. fuscus* record](https://genomesize.com/result_species.php?id=553) also points to heterogeneous published Feulgen estimates rather than documenting the exact 16.36-pg constant.  Recovering a plausible species-level value would still not establish comparability among independently stained and imaged slides.

There is also built-in measurement dependence: species median IOD correlates with nucleus area (rho = 0.680, p = 0.0019) and mean OD (rho = 0.810, p = 4.595e-05), because IOD contains nuclear area algebraically.  A path `image-IOD genome size -> nucleus area` therefore cannot be treated as an independent measurement test.

The corrected release withdraws the picogram conversion and retains only relative-IOD sensitivity columns.  Absolute genome size and causal genome-to-nucleus claims require an independently calibrated assay or recovered same-batch Feulgen-standard evidence.

Figures:

- [`results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png)
- [`results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png)
- [`results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png`](../../results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png)

## 5. Files for the path-analysis audit

- [`results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv): review and specimen/image support.
- [`results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv): morphology estimator sensitivity.
- [`results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv): species-rank stability.
- [`results/data/corrected/microscopy/microscopy_selected50_hierarchical_bootstrap_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_selected50_hierarchical_bootstrap_analysis18_v1.csv): conditional hierarchical intervals.
- [`results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv): relative-IOD sensitivity only.
- [`results/data/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.csv): image-level IOD QC.
- [`results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv): segmentation validation and leakage audit.
- [`results/data/corrected/microscopy/microscopy_source_registry_analysis18_v1.csv`](../../results/data/corrected/microscopy/microscopy_source_registry_analysis18_v1.csv): immutable source hashes and runtime versions.

## Bottom line

The microscopy layer is now auditable but not wholly approved.  Linkage provenance is strong.  Morphology is usable only as a declared upper-tail sensitivity estimand with estimator and specimen-support sensitivity.  Segmentation generalization to the focal species remains unvalidated.  Image IOD is not an absolute genome-size measurement under the evidence currently stored in the repository.

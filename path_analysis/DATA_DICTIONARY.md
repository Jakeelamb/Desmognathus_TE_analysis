# Data Dictionary

This document maps the conceptual nodes in the path-analysis plan to the current observed variables.

## Provenance Namespaces

The genomic and microscopy evidence streams are joined at the species level;
they are not matched specimens.

| Namespace / field | Meaning |
|---|---|
| `te_sra_accession` | SRA resource used for the species-level TE analysis. |
| `te_assembly_accession` | Assembly resource used for assembly-dependent TE, LTR, or ectopic analysis. |
| `te_resource_source_id` | Stable resource row in `data/templates/source_manifest.csv`. |
| `specimen_id` | Independently collected microscopy specimen identifier. It must never be populated from a genomic voucher or BioSample. |
| `genome_size_pg` | Historical microscopy image-IOD rescaling. It is not a measurement from `te_assembly_accession` and is not approved as an absolute C-value. |

Until the `te_*` fields are propagated into rebuilt analysis tables, the active
selection remains in `input_data/lookup_table.txt` and the full *fuscus*
resource decisions remain in `data/templates/source_manifest.csv`. No genomic
accession should be added to a `cellprofiler_*` artifact.

## Core Variables

| Concept | Current column | Current source | Status | Notes |
|---|---|---|---|---|
| Historical image-IOD rescaling | `genome_size_pg` | `path_analysis/data/external/derived/cellprofiler_final_species_results.csv` | Preserved; not approved as absolute genome size | The fixed `fuscus = 16.36 pg` source/convention and between-slide Feulgen comparability are not established. Do not use this column in the corrected release. |
| Relative nuclear-IOD index | `relative_iod_to_fuscus_*` | `results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv` | Sensitivity-only | Fuscus equals 1 within all-selected, image-QC-pass, or high-QC subsets. This is not a C-value. |
| Historical image-IOD interval | `genome_size_se_pg`, `genome_size_ci_low_pg`, `genome_size_ci_high_pg` | historical import above | Conditional object bootstrap only | Does not include slide, stain, calibration, camera, selection, or segmentation uncertainty. |
| Genome QC status | `genome_result_status`, `genome_flag_summary` | same as above | Available | Use for stable-only or caution-exclusion sensitivity runs. |
| Nucleus size | `morph_nucleus_area_um2` plus corrected estimator table | historical import plus `results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv` | Upper-tail sensitivity trait | Frozen rows are composite quality-ranked, not literal largest 50 or typical cells; 574/900 focal rows lack individual keep labels. |
| Cell size | `morph_cell_area_um2` plus corrected estimator table | same as above | Upper-tail sensitivity trait | Use estimator variants and specimen-support labels in every comparative sensitivity. |
| Cytoplasm area | `morph_cytoplasm_area_um2` | same as above | Available | Optional downstream morphology proxy. |
| N:C ratio | `morph_nc_ratio` | same as above | Available | Derived quantity. Do not include it in the same DAG as both cell area and nucleus area. |
| TE composition, order level | `order_*` columns | `results/data/dnaPipeTE_order_breakdown.csv` plus corrected mass ledger | Relative composition; denominator audited | Closed percentages conditional on retained order-classified mass. Use `dnapipete_mass_accounting_analysis18_v2.csv` to preserve unresolved mass and accessions. |
| TE composition, superfamily level | `superfamily_*` columns | `results/data/dnaPipeTE_superfamily_breakdown.csv` plus corrected mass ledger | Relative composition; denominator audited | Higher-dimensional closed block with 17.507%–21.586% unresolved mass in the final panel. Prefer declared log-ratios/ILR for secondary analyses. |
| Repeat-aligned load | not promoted to path input | `results/data/corrected/dnapipete/dnapipete_absolute_load_sensitivity_analysis18_v1.csv` | Sensitivity-only | Fraction of the historically configured 1.5-Gb quantification sample mapped to annotated repeat contigs. Runtime log, container digest, library checksum, and sampling uncertainty are missing. Genomic and microscopy specimens remain independent. |
| TE diversity | `order_shannon`, `order_simpson`, `order_pielou` | historical `results/data/diversity_order_stats.csv`; corrected `results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv` | Approved descriptive, conditional on classified mass | `order_simpson` is Gini-Simpson `1 - sum(p_i^2)`, not a finite-count correction. Shannon uses natural logs. Preserve the unresolved-bin sensitivity and use `order_pielou` only with its denominator stated. |
| TE order CLR PCA | `PC1`, `PC2` in corrected scores table | `results/data/corrected/diversity_pca/te_pca_scores_analysis18_v1.csv` | Approved descriptive only | Primary order-classified analysis has 18 species and 10 nonzero features; PC1/PC2 explain 56.7%/28.9%. Do not promote axes to causal predictors without phylogenetic sensitivity. |
| TE superfamily CLR PCA | supplementary scores/loadings | same corrected diversity/PCA directory | Sensitivity-only | 27 retained features, 32 zero replacements, and unstable PC2 under alternative replacement/leave-one-out analysis. |
| TE clade PERMANOVA | historical clade test | historical PCA/PERMANOVA outputs | Not approved for inference | Unrestricted species-label permutations violate phylogenetic exchangeability. Use a phylogenetic multivariate model for an inferential claim. |
| Focal phylogeny | corrected 18-tip tree | `results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk` | Structurally approved; focal source provenance open | Exact prune, no substitution, rounding-only terminal padding. `fuscus` is retained and `planiceps` excluded. |
| Published phylogeny uncertainty | main plus 200 time trees | corrected phylogeny directory; Stewart-Wiens 2025 Supplementary Files S3/S4 | Approved sensitivity set | Exact 18-taxon matches. Use across final comparative/path models; bootstrap trees do not represent reticulation uncertainty. |
| TE divergence / turnover | `weighted_te_divergence_p90` | aggregated from `results/data/divergence/divergence_summary_statistics_by_species.csv` | Available | Current species-level TE turnover summary used in the compact TE feature panel. |
| RepeatMasker deletion-gap statistic | `weighted_te_deletions_p90` | historical summary plus corrected sensitivity branch | Sensitivity-only | Alignment gaps relative to repeat consensus, not a genomic DNA-loss or ectopic-recombination rate. Corrected values are in `results/data/corrected/divergence/`. |
| Terminal:internal LTR depth proxy | corrected robust species metrics | `results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv` | Exploratory only | Mapping/deletion-footprint proxy; not validated as an ectopic-recombination rate. Prefer median/geometric sensitivity and retain zero depth. |
| Ectopic support | `ectopic_n_rows_total`, `ectopic_n_elements`, `ectopic_n_complete_known`, `ectopic_complete_fraction` | same as above | Available | `ectopic_n_rows_total` counts all stored rows, `ectopic_n_elements` counts usable non-missing ratio rows, and `ectopic_complete_fraction` is calculated only over known `yes`/`no` completion states. |
| Body size | `body_size_proxy_mm` | `path_analysis/data/derived/organismal_traits_curated.csv` | Available | Current analysis-facing size proxy with confidence and provenance tracking. |
| Development mode | `development_mode` | same as above | Available | Currently the cleanest metamorphosis/direct-development axis. |
| Lifestyle / aquaticity | `aquaticity_index`, `microhabitat_class` | same as above | Available | Curated ecological covariates now staged for organismal-augmented models. |

## Cell And Nucleus Measurement Formulas

The path-analysis layer treats the imported verified CellProfiler/YOLO tables as
the authoritative current measurement snapshots. When those snapshots are
rebuilt by `path_analysis/scripts/pull_cellprofiler_estimates.py`, the following
formulas and units should hold:

| Quantity | Formula / unit | Notes |
|---|---|---|
| Cell area | `cell_area_um2`, square micrometers | Weighted median per species in `morph_cell_area_um2`; requires a verified linked cell mask. |
| Nucleus area | `nuc_area_um2`, square micrometers | Weighted median per species in `morph_nucleus_area_um2`; requires a verified linked nucleus mask. |
| Cytoplasm area | `cell_area_um2 - nuc_area_um2`, square micrometers | Imported as `morph_cytoplasm_area_um2` when upstream exports it. |
| N:C ratio | `nuc_area_um2 / cell_area_um2` | Imported as `morph_nc_ratio`; do not model it alongside both numerator and denominator. |
| Integrated optical density | `iod` / `nuc_iod` | Exactly `area_px * mean_od` in the audited export; retain only as relative sensitivity. This creates algebraic dependence with nucleus area. |
| Mean optical density | `iod / area_px` or upstream `mean_od` when exported | Requires upstream image-calibration metadata and thresholding state. |

Per-image traceability belongs in the imported sidecar tables
`cellprofiler_linked_genome_image_trace.csv`,
`cellprofiler_species_morphology_image_trace.csv`,
`cellprofiler_genome_sensitivity.csv`, and
`cellprofiler_image_iod_quality_summary.csv`, including source image, mask, tile
manifest, run manifest, support tier, and hash fields where available.

The current top-50 rerun also writes compact analysis summaries under
`path_analysis/results/`:

- `top50_size_estimate_spread_summary.csv`
- `top50_size_spearman_correlations.csv`
- `top50_model_ranking_summary.csv`
- `top50_size_analysis_summary.md`

## Deferred Variables

| Concept | Current state | Why deferred |
|---|---|---|
| Independent genome-size assay | Not available in this repo | Required before an absolute C-value or `genome -> nucleus -> cell` causal claim. |
| Developmental rate | Not yet staged in this repo | Same. |
| Richer reproductive life-history rates | Not yet staged in a comparative-ready table | Useful later, but current organismal layer focuses on size, development mode, and broad habitat/lifestyle. |
| Finer ecomorph / shape axes | Only partly staged | Broad habitat coding exists, but shape-residual and richer external morphology axes are still incomplete. |

## Preferred First-Pass Derived Variables

These are the planned observed variables used by the current R scaffold:

| Derived variable | Meaning | Construction |
|---|---|---|
| `ltr_balance` | Interpretable TE-composition summary | `log((order_LTR + p) / (order_LINE + p))`, where `p` is a small pseudocount |
| `te_evenness` | TE distribution evenness | Current order-level Pielou evenness |
| `ectopic_index` | Species-level ectopic recombination proxy | `log10(ectopic_mean_ratio)` after positivity checks |
| historical `gs` | Scaled image-IOD rescaling | Z-score of `log10(genome_size_pg)`; preserved historical scaffold only, not a corrected absolute-genome variable |
| `ns` | Scaled nucleus-size variable | Z-score of `log10(morph_nucleus_area_um2)` |
| `cs` | Scaled cell-size variable | Z-score of `log10(morph_cell_area_um2)` |

## Datasets Produced By The Builder

The builder writes the following path-ready tables to `path_analysis/data/derived/`:

| File | Intended use |
|---|---|
| `master_species_table.csv` | Union of all currently available species with availability flags and merged variables |
| `dataset_te_genome.csv` | Core TE-mechanism models |
| `dataset_te_genome_ectopic.csv` | TE + ectopic mechanism models |
| `dataset_genome_morphology.csv` | Nucleotypic planning models |
| `dataset_te_genome_morphology.csv` | Integrated TE -> genome -> morphology planning models |
| `dataset_te_genome_ectopic_morphology.csv` | Most restrictive integrated subset |
| `dataset_overlap_summary.csv` | Species counts and species lists for each staged dataset |

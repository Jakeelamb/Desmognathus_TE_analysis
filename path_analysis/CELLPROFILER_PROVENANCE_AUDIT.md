# CellProfiler Provenance Audit

This document records the current provenance status of the imported
CellProfiler layer used by `path_analysis/`.

It answers a narrower question than the broader methods notes:

Can the merged `cellprofiler_test -> path_analysis` bridge now be audited from
the preserved upstream run outputs without falling back to circular downstream
tables?

## Bottom Line

Yes for computational traceability. No for publication validation of the
segmentation models or absolute genome-size interpretation.

As of July 8, 2026:

- the imported genome bundle is no longer reconstructed from downstream
  `path_analysis` tables
- the current imported bridge snapshot uses the exact curated top-50 verified
  species dataset under
  `cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/`
- `cellprofiler_final_species_results.csv` now carries the verified
  OD-QC-primary genome estimate plus bootstrap intervals and support labels
- `cellprofiler_species_morphology_summary.csv` now carries the verified
  per-species cell and nucleus area estimates plus bootstrap intervals
- the traceability audits are currently header-only for gaps

## Historical March Rebuild

The earlier imported CellProfiler state was refreshed on March 10, 2026 with
the following active run tags:

- raw genome and raw cell source run: `full_dataset_v1`
- fresh mixed linkage run: `mixed_cellpose_yolo_pubrebuild_20260310T063533Z`

What was rerun for this refresh:

- brightfield-only mixed-linkage tile preparation
- YOLO nucleus measurement on the rebuilt tile manifest
- cell-to-nucleus linkage and morphology summary generation
- `path_analysis` import, master-table rebuild, panel rebuild, path-model
  reruns, and source-traceability audit

What was intentionally reused:

- `output/runs/full_dataset_v1/cell_size_segmentation`
- `output/runs/full_dataset_v1/nucleus_iod/brightfield`

Why the raw cell stage was reused:

- a fresh brightfield Cellpose rerun was started with GPU enabled, but the
  first image (`Process_337_raw_green.ome`) reached about `700.5 s` for the
  first `2` of `183` scored tiles, which is not overnight-feasible for a full
  60-image rerun on the current hardware
- that earlier shareable snapshot represented a clean rebuild from the
  preserved raw cell/genome runs plus a fresh mixed-linkage and downstream
  analysis rerun, not a de novo raw cell rerun from microscope images

## What Is Now Reconstructable

The imported CellProfiler layer now has an auditable local chain from:

1. verified reviewed linked-pair outputs under
   `cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/`
2. selected-pair and image-level trace tables
   `path_analysis/data/external/derived/cellprofiler_linked_genome_image_trace.csv`
   and
   `path_analysis/data/external/derived/cellprofiler_species_morphology_image_trace.csv`
3. species-by-state summary
   `path_analysis/data/external/derived/cellprofiler_genome_state_summary.csv`
4. imported species bundle
   `path_analysis/data/external/derived/cellprofiler_final_species_results.csv`
5. merged path-analysis tables such as
   `path_analysis/data/derived/master_species_table.csv` and
   `path_analysis/data/derived/path_input_master.csv`

That is the correct direction of dependency for an auditable analysis-facing
workflow.

## Historical Import Rule

The bridge:

- imports verified selected linked cell+nucleus pairs
- uses weighted species medians for cell and nucleus area
- uses OD-QC-pass linked nucleus-IOD rows for the primary genome estimate when
  available
- retains all-selected and high-OD-QC genome sensitivity rows in
  `cellprofiler_genome_sensitivity.csv`
- scales verified linked nucleus IOD to `D. fuscus = 16.36 pg`
- records support counts and support tier in the imported bundle

This import is explicitly linked-IOD-derived. It is a frozen historical bridge
snapshot, not a validated absolute genome-size assay. The corrected release
does not promote its picogram column.

## July 9 Publication-Release Audit

The final-18 audit is
`plans/publication-readiness-deep-audit/microscopy_release_audit_analysis18_v1.md`.
Its machine-readable products are under
`results/data/corrected/microscopy/` and its validation figures are under
`results/figures/corrected/microscopy/`.

The audit scores the exact archived production cell masks and the current YOLO
nucleus model on the only two manually labeled test tiles. Cell instance F1 is
0.624 at IoU 0.50 and nucleus instance F1 is 0.687. Neither test set contains a
focal analysis species, and the nucleus train/test tiles come from the same two
source images. These are diagnostic results, not cross-species validation.

For the final 18 species, 324/900 frozen rows are explicit manual keeps, two are
maybes, and 574 are model-ranked/unlabeled. The corrected morphology tables
therefore label the trait as a quality-screened upper-tail sensitivity estimand.
The corrected IOD table reports only a relative nuclear-IOD index with fuscus =
1 within each QC subset; it does not report picograms.

## Machine-Readable Audit Outputs

The canonical audit files are:

- `path_analysis/data/external/derived/cellprofiler_source_discovery.json`
- `path_analysis/data/external/derived/cellprofiler_traceability_audit_summary.csv`
- `path_analysis/data/external/derived/cellprofiler_traceability_audit_gaps.csv`
- `path_analysis/data/external/derived/cellprofiler_final_species_results_reconstruction.json`
- `path_analysis/data/external/derived/cellprofiler_genome_sensitivity.csv`
- `path_analysis/data/external/derived/cellprofiler_image_iod_quality_summary.csv`
Refresh them with:

```bash
scripts/run_in_dusky.sh python path_analysis/scripts/pull_cellprofiler_estimates.py \
  --verified-species-dir /home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest
```

## Remaining Warnings

The current traceability gap file is header-only. ROI ZIP paths are absent, but
this is not a gap for the verified mask-based workflow because saved cell and
nucleus masks are the matching artifacts.

## Analysis Use Guidance

What is defensible now:

- exact cell/nucleus object traceability and one-to-one linkage provenance
- per-species cell and nucleus upper-tail sensitivity estimates, with the
  estimator and specimen-support boundaries stated
- relative nuclear-IOD sensitivity indices, not absolute genome sizes
- exact top-50 size/spread summaries in
  `path_analysis/results/top50_size_analysis_summary.md`
- exploratory TE-to-relative-IOD comparative models, with the measurement
  proxy and QC subset named explicitly
- analysis supplements that cite the machine-readable audit outputs and the
  verified import rule

What should still remain sensitivity-only:

- any strong causal interpretation of `genome size -> nucleus size -> cell size`
- any claim that the current imported genome bundle is an independent
  non-IOD genome-size assay
- any claim that the current two-tile segmentation benchmark validates error
  rates across the focal species

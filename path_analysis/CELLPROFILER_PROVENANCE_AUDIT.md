# CellProfiler Provenance Audit

This document records the current provenance status of the imported
CellProfiler layer used by `path_analysis/`.

It answers a narrower question than the broader manuscript methods:

Can the merged `cellprofiler_test -> path_analysis` bridge now be audited from
the preserved upstream run outputs without falling back to circular downstream
tables?

## Bottom Line

Yes, with explicit warnings.

As of March 10, 2026:

- the imported genome bundle is no longer reconstructed from
  `path_analysis/data/derived/master_species_table.csv`
- when the legacy upstream species bundle is absent, the bridge now rebuilds
  `cellprofiler_final_species_results.csv` directly from the preserved raw
  nucleus-IOD run outputs
- the bridge also writes `cellprofiler_genome_state_summary.csv`, which makes
  the per-species brightfield/pmount support and state-selection rule explicit
- the publication-readiness audit currently returns overall status `warn`, not
  `fail`, because the remaining gaps are retained traceability warnings rather
  than circular provenance failures

## Latest Overnight Rebuild

The current imported CellProfiler state was refreshed on March 10, 2026 with
the following active run tags:

- raw genome and raw cell source run: `full_dataset_v1`
- fresh mixed linkage run: `mixed_cellpose_yolo_pubrebuild_20260310T063533Z`

What was rerun for this freeze:

- brightfield-only mixed-linkage tile preparation
- YOLO nucleus measurement on the rebuilt tile manifest
- cell-to-nucleus linkage and morphology summary generation
- `path_analysis` import, master-table rebuild, panel rebuild, path-model
  reruns, source-traceability audit, publication-readiness audit, and paper
  freeze manifest refresh

What was intentionally reused:

- `output/runs/full_dataset_v1/cell_size_segmentation`
- `output/runs/full_dataset_v1/nucleus_iod/brightfield`

Why the raw cell stage was reused:

- a fresh brightfield Cellpose rerun was started with GPU enabled, but the
  first image (`Process_337_raw_green.ome`) reached about `700.5 s` for the
  first `2` of `183` scored tiles, which is not overnight-feasible for a full
  60-image rerun on the current hardware
- the current publication freeze therefore represents a clean rebuild from the
  preserved raw cell/genome runs plus a fresh mixed-linkage and downstream
  analysis rerun, not a de novo raw cell rerun from microscope images

## What Is Now Reconstructable

The imported CellProfiler layer now has an auditable local chain from:

1. mixed Cellpose+YOLO linkage outputs under
   `cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/`
2. image-level trace table
   `path_analysis/data/external/derived/cellprofiler_linked_genome_image_trace.csv`
3. species-by-state summary
   `path_analysis/data/external/derived/cellprofiler_genome_state_summary.csv`
4. imported species bundle
   `path_analysis/data/external/derived/cellprofiler_final_species_results.csv`
5. merged path-analysis tables such as
   `path_analysis/data/derived/master_species_table.csv` and
   `path_analysis/data/derived/path_input_master.csv`

That is the correct direction of dependency for an auditable manuscript-facing
workflow.

## Reconstruction Rule

When the legacy upstream file
`cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv`
is missing, the bridge now:

- groups strict-core linked YOLO nuclei into image-level linked nucleus-IOD summaries
- prefers `analysis_ready_image = TRUE` images when a species has them, otherwise
  falls back to all strict-core linked images
- aggregates those image medians to specimen-level summaries
- scales specimen-level mean image-median linked nucleus IOD to
  `D. fuscus = 16.36 pg`
- records support counts and support tier in the imported bundle

This reconstruction is explicitly linked-IOD-derived and remains provisional for any
formal `genome -> nucleus -> cell` interpretation.

## Machine-Readable Audit Outputs

The canonical audit files are:

- `path_analysis/data/external/derived/cellprofiler_source_discovery.json`
- `path_analysis/data/external/derived/cellprofiler_traceability_audit_summary.csv`
- `path_analysis/data/external/derived/cellprofiler_traceability_audit_gaps.csv`
- `path_analysis/data/external/derived/cellprofiler_final_species_results_reconstruction.json`
- `path_analysis/data/derived/publication_readiness_checks.csv`
- `path_analysis/data/derived/publication_readiness_summary.json`

Refresh them with:

```bash
python3 path_analysis/scripts/pull_cellprofiler_estimates.py
python3 path_analysis/scripts/audit_publication_readiness.py
```

## Remaining Warnings

The current warnings are concentrated in upstream trace retention, not in the
merge logic itself:

- the morphology linkage layer still lacks ROI-ZIP trace paths

These warnings are tracked automatically in
`cellprofiler_traceability_audit_gaps.csv`.

## Publication Use Guidance

What is defensible now:

- TE-to-genome comparative models using the imported genome layer, with the
  linked-IOD-derived and warning-laden status made explicit
- manuscript supplements that cite the machine-readable audit outputs and the
  linked-run reconstruction rule

What should still remain sensitivity-only:

- any strong causal interpretation of `genome size -> nucleus size -> cell size`
- any claim that the current imported genome bundle is the final independent
  genome-size estimate rather than an audited interim reconstruction

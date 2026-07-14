# Path Analysis Workspace

This folder is the active workspace for the phylogenetic path-analysis layer. It keeps the current dataset assembly code, model scaffold, source manifests, and small derived audit tables in one place.

If returning after a context switch or changing focus within the project, start with `PATH_ANALYSIS_STATUS.md`.

## What Is Here

- `PATH_ANALYSIS_STATUS.md`
  The current high-level map of the workspace, pipeline, panel definitions, generated outputs, and claim boundaries. This should be the first file to read when restarting or shifting to a different task within the same analysis.
- `README.md`
  This overview and the current staged plan.
- `STUDY_SPECIES_PANELS.md`
  The authoritative TE34 / linked-cell21 / integrated-path18 contract. Read
  this before interpreting a missing species as an analytical exclusion.
- `DATA_DICTIONARY.md`
  Variable definitions, preferred observed proxies, and current caveats.
- `TE_DATA_AUDIT.md`
  Current audit status of the canonical TE source tables and derived path-analysis TE products, including validated checks and remaining warnings.
- `TE_METHODS_LANGUAGE.md`
  Methods wording for the TE layer, including Simpson and ectopic denominator semantics.
- `TE_DIVERSITY_CANONICALIZATION.md`
  Non-destructive reconstruction note for the TE diversity summary tables, including the audited path from threshold-grid outputs to the current canonical diversity snapshots.
- `TE_PROVENANCE_AUDIT.md`
  Upstream provenance map for the canonical TE summary tables, including the now-resolved diversity-summary reconstruction path.
- `CELLPROFILER_PROVENANCE_AUDIT.md`
  Human-readable audit note for the imported CellProfiler genome and morphology layer, including the verified-species import rule used by the refresh bridge.
- `INPUT_PREPARATION_PLAN.md`
  Detailed plan for literature trait mining, taxonomy crosswalks, and TE feature engineering.
- `SOURCE_TRACKING.md`
  Provenance rules and source-manifest workflow for reproducible traceability.
- `CANDIDATE_MODELS.md`
  The initial DAG families to compare with `phylopath`.
- `scripts/build_master_dataset.py`
  Builds overlap-ready species tables from the repo plus the imported CellProfiler snapshots in `path_analysis/data/external/derived/`.
- `scripts/pull_cellprofiler_estimates.py`
  Rebuilds the imported CellProfiler snapshots and traceability audits from the active `cellprofiler_test` run outputs, preferring the verified species dataset when present and falling back to legacy linked YOLO reconstruction otherwise.
- `scripts/path_model_scaffold.R`
  Prepares transformed analysis inputs and defines the current candidate model sets.
- `scripts/prepare_te_features.py`
  Builds TE composition, diversity, turnover, and ectopic proxies into a reusable feature table.
- `scripts/build_canonical_diversity_tables.py`
  Rebuilds scratch candidates for the canonical TE diversity summary tables and audits them against the current `results/data/` snapshots without overwriting those files.
- `scripts/prepare_external_morphometrics.py`
  Summarizes raw external morphometric appendices into source-linked species tables.
- `scripts/prepare_nc_biodiversity_traits.py`
  Parses local NC Biodiversity Project HTML snapshots into a source-linked species trait table.
- `scripts/prepare_curated_organismal_traits.py`
  Collapses external trait sources into an analysis-facing organismal trait table with retained provenance.
- `scripts/build_phylogenetic_trait_imputation.R`
  Builds a sensitivity-only phylogenetic nearest-neighbor imputation layer for missing or low-confidence organismal traits without overwriting the curated observed table.
- `scripts/build_te_model_panel.py`
  Cuts the full TE feature table down to the compact predictor panel intended for comparative models.
- `scripts/build_analysis_panels.py`
  Writes panel files for each model family plus primary and sensitivity subsets with confidence flags and species-level eligibility audits.
- `scripts/refresh_source_file_inventory.py`
  Rebuilds the raw-file hash inventory for stored external inputs.
- `scripts/audit_source_traceability.py`
  Builds a source-file registry plus gap reports so every source id used in the path-input tables can be audited.
- `TE_MODEL_INPUTS.md`
  Notes on which TE predictors are recommended for the main and sensitivity path-analysis families.
- `ANALYSIS_DATASETS.md`
  Defines the panel files that should be used for primary versus sensitivity path analyses.
- `data/derived/`
  Generated master tables, staged analysis datasets, and species-level readiness audits.
- `data/templates/`
  Templates for taxonomy crosswalks, source manifests, trait registries, and literature extraction.
- `data/external/`
  Raw and derived external source files with provenance tracking.
- `results/`
  Placeholder for model-ranking tables, coefficient summaries, and figures.

## Current Data Situation

These counts reflect the current files already present in the workspace:

- Vetted genomic-resource species (TE34): `34`
- Current linked-cell species (Cell21): `21`
- TE + linked-cell evidence intersection (path18): `18`
- TE + genome + ectopic overlap: `16`
- Genome + linked morphology overlap: `21`
- TE + genome + linked morphology overlap: `18`
- TE + genome + ectopic + linked morphology overlap: `16`

## Important Caveat

The historical CellProfiler bridge uses the exact frozen top-50 dataset under
`cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/`.
The publication-release interpretation is now narrower than the bridge column
names: cell and nucleus areas are quality-screened upper-tail sensitivity
traits, and image IOD is a relative nuclear-intensity proxy rather than an
approved absolute genome-size assay. See
`../plans/publication-readiness-deep-audit/microscopy_release_audit_analysis18_v1.md`
and `../results/data/corrected/microscopy/`.

For now:

- `TE -> relative nuclear IOD` is exploratory proxy analysis only.
- terminal:internal LTR depth is an exploratory deletion-footprint/mapping
  proxy, not an ectopic-recombination rate.
- `genome size -> nucleus size -> cell size` remains scaffolded until an
  independent genome-size assay is available.

## Corrected Final-18 Audit Branch

The non-destructive corrected branch is the integrated path18 layer only. It
does not replace the larger TE34 or Cell21 descriptive analyses, historical
tables, or promote relative IOD to genome size. See `STUDY_SPECIES_PANELS.md`.

```bash
scripts/run_in_dusky.sh python scripts/processing/build_corrected_path_inputs.py
scripts/run_in_dusky.sh Rscript scripts/processing/audit_corrected_path_models.R --phase all
scripts/run_in_dusky.sh Rscript scripts/processing/simulate_corrected_path_calibration.R
scripts/run_in_dusky.sh python scripts/processing/summarize_corrected_path_audit.py
scripts/run_in_dusky.sh python scripts/processing/build_publication_audit_notebook.py
scripts/run_in_dusky.sh jupyter nbconvert --to notebook --execute notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb --inplace --ExecutePreprocessor.timeout=600
scripts/run_in_dusky.sh python scripts/processing/audit_publication_notebook.py
```

Review these first:

- `../plans/publication-readiness-deep-audit/corrected_path_analysis_audit_analysis18_v1.md`
- `../results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv`
- `../notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb`
- `../notebooks/Desmognathus_publication_audit_analysis18_v1.manifest.json`

The model branch exports all rankings, basis-set components, best-model edges,
tree distributions, species-omission results, and simulation replicates. Every
row marks absolute-genome and publication-causal claims false.

## Recommended Analysis Order

1. Build the consolidated species tables.
2. Use the corrected final-18 TE and phylogeny inputs.
3. Run data-estimator, IOD-QC, leave-one-species-out, and 200-tree sensitivity.
4. Keep every image-IOD path result explicitly exploratory; do not promote a
   genome-size causal claim.

## How To Use

Run project commands through the Dusky wrapper:

```bash
scripts/run_in_dusky.sh python verify_setup.py --skip-data
```

The path-analysis dependencies are listed in `Dusky.yml`. Refresh the
environment from that file rather than installing packages inside analysis
scripts:

```bash
conda env update -f Dusky.yml
```

The wrapper prepends `$CONDA_PREFIX/bin` before execution, so Python and R
resolve from the active `Dusky` Conda environment rather than from the host
PATH. The scaffold script also prefers the active Conda R library when
`CONDA_PREFIX` is set, so it does not accidentally load incompatible user-level
packages.

Refresh the imported CellProfiler layer first:

```bash
scripts/run_in_dusky.sh python path_analysis/scripts/pull_cellprofiler_estimates.py
```

Then build the merged datasets:

```bash
scripts/run_in_dusky.sh python path_analysis/scripts/build_master_dataset.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_te_features.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_amphibio_traits.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_external_morphometrics.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_nc_biodiversity_traits.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_curated_organismal_traits.py
scripts/run_in_dusky.sh Rscript path_analysis/scripts/build_phylogenetic_trait_imputation.R
scripts/run_in_dusky.sh python path_analysis/scripts/build_te_model_panel.py
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_ltr_history_features.py
scripts/run_in_dusky.sh python path_analysis/scripts/build_path_input_master.py
scripts/run_in_dusky.sh python path_analysis/scripts/build_analysis_panels.py
scripts/run_in_dusky.sh python path_analysis/scripts/audit_source_traceability.py
```

The main analysis-panel products are:

- `path_analysis/data/derived/analysis_panel_summary.csv`
  Panel-level counts and species lists.
- `path_analysis/data/derived/analysis_species_readiness.csv`
  Species-level eligibility and exclusion reasons for every panel.
- `path_analysis/data/derived/source_traceability_gaps.csv`
  Traceability audit report; this should stay empty for reproducible analysis use.

The phylogenetic-imputation products are:

- `path_analysis/data/derived/phylogenetic_trait_imputation_long.csv`
  Long-form audit table showing, for each species-trait combination, whether the value is source-backed, missing and phylogenetically inferred, or low-confidence and either supported or contradicted by the phylogeny.
- `path_analysis/data/derived/organismal_traits_phylo_inference.csv`
  Wide-form species table with explicit `phylo_` columns merged into `path_input_master.csv`.
- `path_analysis/data/derived/phylogenetic_trait_imputation_summary.csv`
  Leave-one-out cross-validation summary for each imputed trait.
- `path_analysis/data/derived/phylofill_panel_comparison.csv`
  Observed-only versus `primary_phylofill` panel comparison showing whether phylogenetic fills actually change the current model subsets.

Important rule:

- the curated observed table remains the primary data layer
- `phylo_` columns are sensitivity-only and should never be silently substituted for source-backed values
- the current `phylofill` panel comparison is a no-op for the TE/genome model families: no species are added to the present overlap sets

Inspect the overlap summary:

```bash
sed -n '1,120p' path_analysis/data/derived/dataset_overlap_summary.csv
```

Preview the model families without requiring `phylopath`:

```bash
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --summary-only
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --summary-only
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_strict_body --summary-only
```

Once the Dusky environment includes `phylopath`, run a family:

```bash
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome
```

To run the same family on a panel-defined sensitivity subset without overwriting existing outputs:

```bash
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_organismal --panel te_genome_organismal_primary_strict_body
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ltr_history --panel te_genome_ltr_history_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic_organismal --panel te_genome_ectopic_organismal_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --panel te_genome_ectopic_primary_strict_body
```

When `--panel` is supplied, the scaffold reads `data/derived/panels/<panel>.csv` and writes results with the panel name as the filename prefix.

Current model-result boundary:

- Panel builders and source-traceability audits are current.
- Re-run panel-specific model rankings before making winner claims from the
  refreshed panels.
- The source-linked note for the paired-LTR historical-axis integration is
  `LTR_HISTORY_PATH_INTEGRATION.md`.

## External Inputs

The dataset builder now consumes imported snapshots under `path_analysis/data/external/derived/`:

- `cellprofiler_final_species_results.csv`
- `cellprofiler_genome_state_summary.csv`
- `cellprofiler_species_morphology_summary.csv`
- `cellprofiler_linked_genome_image_trace.csv`
- `cellprofiler_traceability_audit_summary.csv`
- `cellprofiler_traceability_audit_gaps.csv`
- `cellprofiler_source_discovery.json`

Those snapshots are rebuilt from the active `cellprofiler_test` outputs:

- `output/runs/mixed_cellpose_yolo_full_dataset_v1/linkage/`

`build_master_dataset.py` requires these imported snapshots. Refresh the
snapshots with `pull_cellprofiler_estimates.py` before rebuilding merged
path-analysis tables.

The pull step now imports the verified species dataset directly when available,
writes species-level morphology, genome, genome-sensitivity, image-QC, and
traceability sidecars, and records that import explicitly in the audit outputs
rather than silently reusing downstream path-analysis tables.

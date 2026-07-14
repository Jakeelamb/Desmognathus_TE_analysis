# Path Analysis Status

This file is the stable, human-readable map of the `path_analysis/` workspace.
It is meant to be the first thing to read when returning to this project, even
if the session continues in a different direction.

## What This Workspace Is For

The goal is a reproducible comparative path-analysis framework for genome
evolution in `Desmognathus`.

The main conceptual chain is:

`TE features -> genome size -> nucleus size -> cell size`

There is now also an organismal extension layer for:

- body size
- development mode
- lifestyle / aquaticity
- broad habitat class

## Current Bottom Line

The historical data-build is reproducible, but its image-IOD picogram column
and causal model rankings are not publication-approved. Corrected final-18
audits now exist for TE composition/diversity, terminal:internal depth,
phylogeny, microscopy, and path-model sensitivity. The corrected path branch
uses relative nuclear IOD only, never reads the historical `genome_size_pg`
column, and is explicitly exploratory.

The corrected implementation completed 958 family-level fits: 84 data/tree
specifications, 804 fits over the published main tree plus 200 time trees, and
70 leave-one-species-out fits. All completed without a fit failure. Stable
selection does not promote the analysis to a genome-size or causal claim:
relative IOD is uncalibrated and algebraically contains nuclear area,
segmentation lacks held-out final-panel species validation, upper-tail
morphometry has uneven biological support, and trait measurement error is not
jointly propagated.

Separately, the genome24 three-trait notebook reports the current conditional
*D. fuscus*-anchored genome-size estimates in **picograms**, never on a
relative-IOD axis. Its 250 measurement bootstraps propagate the audited image
and morphology sampling uncertainty, but the analysis remains exploratory
because the pg values are still derived from nuclear IOD and the three proposed
causal orientations are Markov-equivalent.

## Workspace Structure

The important files are:

- `../README.md`
  Repo-level workflow overview and cleanup boundary for ignored generated outputs.
- `README.md`
  General workspace overview and command examples.
- `PATH_ANALYSIS_STATUS.md`
  This file. Read this first when returning later.
- `CANDIDATE_MODELS.md`
  The current family definitions and the logic behind them.
- `ANALYSIS_DATASETS.md`
  Definitions of the panel files and sensitivity subsets.
- `DATA_DICTIONARY.md`
  Variable meanings and caveats.
- `TE_DATA_AUDIT.md`
  Current audit status of the repo-local TE source tables and derived path-analysis TE products.
- `TE_METHODS_LANGUAGE.md`
  Methods wording for the TE layer and denominator caveats.
- `TE_DIVERSITY_CANONICALIZATION.md`
  Documents the now-resolved TE diversity summary generation path and the non-destructive candidate reconstruction.
- `TE_PROVENANCE_AUDIT.md`
  Upstream provenance map for the TE summary tables used by the current path-analysis layer.

Important scripts:

- `scripts/build_master_dataset.py`
  Merges TE, genome, morphology, and organismal layers into species-level tables.
- `scripts/prepare_curated_organismal_traits.py`
  Builds the source-backed organismal trait table.
- `scripts/build_phylogenetic_trait_imputation.R`
  Builds the sensitivity-only phylogenetic imputation layer.
- `scripts/build_path_input_master.py`
  Builds the merged analysis input table.
- `scripts/build_analysis_panels.py`
  Writes panel CSVs and the species readiness audit.
- `scripts/path_model_scaffold.R`
  Prepares transformed model inputs and runs the `phylopath` families.
- `../scripts/processing/build_corrected_path_inputs.py`
  Builds the six-morphology-by-three-IOD final-18 sensitivity cube without the
  historical picogram column.
- `../scripts/processing/audit_corrected_path_models.R`
  Runs the fail-closed corrected candidate sets, 200-tree sensitivity, and
  leave-one-species-out audit.
- `../scripts/processing/simulate_corrected_path_calibration.R`
  Measures actual-tree null false selection, signal recovery, coefficient
  bias, and approximate interval coverage at `n=18`.
- `../scripts/processing/summarize_corrected_path_audit.py`
  Freezes review tables, publication gates, figures, and the corrected report.
- `../scripts/processing/build_publication_audit_notebook.py`
  Builds the executable collaborator workbench under `../notebooks/`.
- `../scripts/processing/build_research_review_notebooks.py`
  Builds eight independent analysis-domain notebooks from frozen corrected
  outputs under `../notebooks/research_review/`.
- `scripts/build_cell_nucleus_genome_path_notebook.py`
  Freezes the exact genome24 traits, all 25 labeled DAGs, 11 equivalence
  classes, and 250 paired measurement-bootstrap trait panels.
- `scripts/run_cell_nucleus_genome_phylogenetic_path_analysis.R`
  Fits the ten testable genome24 classes and the declared measurement, tree,
  species, evolutionary-model, and actual-tree simulation sensitivities.
- `scripts/build_cell_nucleus_genome_path_presentation.py`
  Fails closed on quick/incomplete path runs and writes eight PNG/PDF figures
  plus the canonical Notebook 07 source.
- `../scripts/processing/audit_research_review_notebooks.py`
  Fails closed on missing/unexecuted/error notebooks or any external-process
  execution surface, records optional microscopy viewer availability, and
  writes the bundle integrity manifest.

Important derived data:

- `data/derived/organismal_traits_curated.csv`
  The source-backed organismal trait table. This is the factual organismal layer.
- `data/derived/organismal_traits_phylo_inference.csv`
  Explicit `phylo_` trait columns for sensitivity-only inference.
- `data/derived/path_input_master.csv`
  The merged species table used downstream.
- `data/derived/analysis_panel_summary.csv`
  Panel counts and species lists.
- `data/derived/analysis_species_readiness.csv`
  Species-by-panel inclusion and exclusion reasons.
- `data/derived/phylofill_panel_comparison.csv`
  Shows whether phylogenetic backfilling changes current panel membership.
- `results/*.csv`
  Generated model rankings and edge tables for each family/panel combination.
  These are not tracked source artifacts and should be regenerated from the
  documented inputs when needed.
- `results/top50_size_analysis_summary.md`
  Compact readout of the exact top-50 cell, nucleus, genome, spread, correlation,
  and top model ranking summaries.

## How The Pipeline Works

The logic is:

1. Build species-level TE, genome, and morphology tables.
2. Build the source-backed organismal table.
3. Build the phylogenetic sensitivity layer in separate `phylo_` columns.
4. Merge everything into `path_input_master.csv`.
5. Build panel subsets for medium-plus, strict-body, ectopic, morphology, and
   organismal families.
6. Run `path_model_scaffold.R` on the chosen panel.

The normal rebuild sequence is:

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

To run models:

```bash
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --panel te_genome_ectopic_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_organismal --panel te_genome_organismal_primary_mediumplus
scripts/run_in_dusky.sh Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic_organismal --panel te_genome_ectopic_organismal_primary_mediumplus
```

## Current Data State

Current panel sizes:

- `te_genome_primary_mediumplus = 18`
- `te_genome_primary_strict_body = 9`
- `te_genome_organismal_primary_mediumplus = 18`
- `te_genome_organismal_primary_strict_body = 9`
- `te_genome_ltr_history_primary_mediumplus = 15`
- `te_genome_ltr_history_primary_strict_body = 8`
- `te_genome_ectopic_primary_mediumplus = 16`
- `te_genome_ectopic_primary_strict_body = 8`
- `te_genome_ectopic_organismal_primary_mediumplus = 16`
- `te_genome_ectopic_organismal_primary_strict_body = 8`
- `te_genome_morphology_primary_mediumplus = 18`
- `te_genome_morphology_primary_strict_body = 9`

Current organismal coverage:

- `body_size_proxy_mm`: 36 species
- `development_mode`: 36 species
- `aquaticity_index`: 36 species
- `microhabitat_class`: 36 species
- `elevation_mid_m`: 9 species

## Truth Hierarchy

This part is important.

There are now two different organismal layers:

1. `organismal_traits_curated.csv`
   This is the source-backed trait table and should be treated as the factual
   organismal layer.
2. `organismal_traits_phylo_inference.csv`
   This is a sensitivity-only inference layer and should never silently replace
   the observed table.

Rules:

- observed values remain primary
- `phylo_` values are explicitly inferred
- `phylo_inferred_for_missing` means "filled only for sensitivity analysis"
- `low_confidence_observed_supported_by_phylogeny` means "the weak observed
  value is phylogenetically corroborated", not "the source is suddenly strong"

The current phylogenetic sensitivity layer is analytically a no-op for the main
model subsets:

- `te_genome`
- `te_genome_organismal`
- `te_genome_ltr_history`
- `te_genome_ectopic`
- `te_genome_ectopic_organismal`
- `te_genome_morphology`

All observed-only and `primary_phylofill` panel species sets are currently
identical.

## Tree Inputs

There is one collaborator-supplied focal tree source for historical path analysis:

- `input_data/phylogeny/desmo900dated_test.tre`

The path-model scaffold and phylogenetic trait-imputation builder both read this
local, currently Git-ignored input tree. The imputation builder normalizes tip labels internally so it
does not depend on the generated `results/phylogeny/processed_phylogeny.nwk`
artifact.

The publication-readiness branch now adds:

- `results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk`
  — exact final-panel prune of the focal tree with rounding-only terminal
  padding; structurally approved, but focal citation/calibration provenance is
  still required.
- `results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk`
  — published optimal time-tree sensitivity.
- `results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex`
  — 200 published time-calibrated bootstrap trees for branch-time and path-model
  sensitivity.

The tracked `results/phylogeny/processed_phylogeny.nwk` is not a second tree
hypothesis: it has the same common-tip topology as the focal source and an
undocumented exact 0.8 branch-length scale. It is not approved for inference.

## Model Result Boundary

The current panel builders and source-traceability audits are the source of
truth. Model ranking CSVs under `path_analysis/results/` or `results/` are
generated outputs; regenerate them from the current panels before making winner
claims.

Current support rules:

- LTR-history panels require at least one high-confidence paired-LTR element.
- Morphology panels carry genome support-status counts so the CellProfiler
  support tier is visible.
- Phylogenetic trait backfilling remains sensitivity-only and currently changes
  no panel membership.

Corrected release artifacts:

- `../results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv`
- `../results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv`
- `../results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv`
- `../results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv`
- `../plans/publication-readiness-deep-audit/corrected_path_analysis_audit_analysis18_v1.md`
- `../notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb`
- `../notebooks/research_review/01_phylogeny_tree_trimming.ipynb`
- `../notebooks/research_review/02_data_tables_and_provenance.ipynb`
- `../notebooks/research_review/03_repeat_analysis_te34.ipynb`
- `../notebooks/research_review/04_ltr_deletion_footprint.ipynb`
- `../notebooks/research_review/05_cell_modeling_and_measurement.ipynb`
- `../notebooks/research_review/06_genome_size_estimation.ipynb`
- `../notebooks/research_review/07_cell_nucleus_genome_path_analysis.ipynb`
- `../notebooks/research_review/08_integrated_phylogenetic_path_analysis.ipynb`

The integrated TE-evenness and IOD-morphology chain is a robust exploratory
association pattern across the audited specifications. The terminal:internal
family has no globally supported model in any full-panel specification and
must not be described as an ectopic-recombination mechanism.

## What To Trust Most Right Now

Most trustworthy:

- panel definitions and species counts
- source-backed organismal table
- source traceability audits
- CellProfiler image/mask/tile traceability snapshots and exact archived-mask
  benchmark hashes
- corrected final-18 morphology-estimator and relative-IOD sensitivity tables

Not publication-approved:

- all historical generated model rankings
- the image-IOD-to-picogram conversion
- segmentation generalization beyond the two excluded-species test tiles
- morphology-linked causal interpretation without an independent genome assay
- phylogenetic imputation as anything beyond sensitivity analysis
- any interpretation that relies on the strict-body morphology subset

## If You Need To Re-Enter Later

Read these in order:

1. `PATH_ANALYSIS_STATUS.md`
2. `README.md`
3. `data/derived/analysis_panel_summary.csv`
4. `data/derived/phylofill_panel_comparison.csv`
5. regenerated ranking CSVs for:
   - `te_genome`
   - `te_genome_ectopic`
   - `te_genome_organismal`
   - `te_genome_ectopic_organismal`

## The Next Logical Step

If returning to this analysis later, the next step is not more data rescue.

It is:

- build a compact cross-family summary table
- write the results interpretation
- decide how to frame the shift from ectopic-dominant models to body-size-plus-TE
  models once organismal covariates are included

If more data work is ever needed later, the most useful remaining targets are:

- an independent non-IOD genome-size assay
- expanded linked morphology coverage
- possibly a cleaner tree for phylogenetic imputation coverage

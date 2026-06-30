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

The data-build phase is the trustworthy layer right now: source-backed
organismal traits, CellProfiler traceability snapshots, TE feature tables, panel
definitions, and species readiness audits. Panel-specific model rankings should
be regenerated from the current panels before making winner claims.

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

There is one canonical tree source for path analysis:

- `input_data/phylogeny/desmo900dated_test.tre`

The path-model scaffold and phylogenetic trait-imputation builder both read this
tracked input tree. The imputation builder normalizes tip labels internally so it
does not depend on the generated `results/phylogeny/processed_phylogeny.nwk`
artifact.

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

## What To Trust Most Right Now

Most trustworthy:

- panel definitions and species counts
- source-backed organismal table
- source traceability audits
- CellProfiler image/mask/tile traceability snapshots

More provisional:

- generated model rankings until rerun from current panels
- morphology-linked causal interpretation
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

- independent final genome-size estimates
- expanded linked morphology coverage
- possibly a cleaner tree for phylogenetic imputation coverage

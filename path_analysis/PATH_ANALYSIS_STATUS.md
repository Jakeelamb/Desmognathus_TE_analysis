# Path Analysis Status

This file is the stable, human-readable map of the `path_analysis/` workspace.
It is meant to be the first thing to read when returning to this project, even
if the session continues in a different direction.

## What This Workspace Is For

The goal is a publication-grade comparative path-analysis framework for genome
evolution in `Desmognathus`.

The main conceptual chain is:

`TE features -> genome size -> nucleus size -> cell size`

There is now also an organismal extension layer for:

- body size
- development mode
- lifestyle / aquaticity
- broad habitat class

## Current Bottom Line

The data-build phase is in good shape.

The important current result is:

- TE-only models are stable.
- Ectopic-only models are stable when organismal covariates are absent.
- Once organismal covariates are allowed, `body_size` matters.
- Once `body_size` and `ectopic_index` compete in the same candidate family,
  the winning model drops `ectopic_index`.
- Phylogenetic trait backfilling currently changes no model panel membership.

So the project has shifted from "collect and rescue traits" to
"interpret what the fitted families are now saying."

## Workspace Structure

The important files are:

- `README.md`
  General workspace overview and command examples.
- `PATH_ANALYSIS_STATUS.md`
  This file. Read this first when returning later.
- `SESSION_HANDOFF.md`
  Secondary restart notes and chronological context, not the main analysis map.
- `CANDIDATE_MODELS.md`
  The current family definitions and the logic behind them.
- `ANALYSIS_DATASETS.md`
  Definitions of the panel files and sensitivity subsets.
- `DATA_DICTIONARY.md`
  Variable meanings and caveats.
- `TE_DATA_AUDIT.md`
  Current audit status of the frozen TE source tables and derived path-analysis TE products.
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
  Model rankings and edge tables for each family/panel combination.

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
python3 path_analysis/scripts/build_master_dataset.py
python3 path_analysis/scripts/prepare_te_features.py
python3 path_analysis/scripts/prepare_amphibio_traits.py
python3 path_analysis/scripts/prepare_external_morphometrics.py
python3 path_analysis/scripts/prepare_nc_biodiversity_traits.py
python3 path_analysis/scripts/prepare_curated_organismal_traits.py
Rscript path_analysis/scripts/build_phylogenetic_trait_imputation.R
python3 path_analysis/scripts/build_te_model_panel.py
python3 path_analysis/scripts/build_path_input_master.py
python3 path_analysis/scripts/build_analysis_panels.py
python3 path_analysis/scripts/audit_source_traceability.py
```

To run models:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate Dusky

Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_mediumplus
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --panel te_genome_ectopic_primary_mediumplus
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_organismal --panel te_genome_organismal_primary_mediumplus
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic_organismal --panel te_genome_ectopic_organismal_primary_mediumplus
```

## Current Data State

Current panel sizes:

- `te_genome_primary_mediumplus = 27`
- `te_genome_primary_strict_body = 16`
- `te_genome_organismal_primary_mediumplus = 27`
- `te_genome_organismal_primary_strict_body = 16`
- `te_genome_ectopic_primary_mediumplus = 24`
- `te_genome_ectopic_primary_strict_body = 15`
- `te_genome_ectopic_organismal_primary_mediumplus = 24`
- `te_genome_ectopic_organismal_primary_strict_body = 15`
- `te_genome_morphology_primary_mediumplus = 18`
- `te_genome_morphology_primary_strict_body = 9`

Current organismal coverage:

- `body_size_proxy_mm`: 37 species
- `development_mode`: 37 species
- `aquaticity_index`: 37 species
- `microhabitat_class`: 37 species
- `elevation_mid_m`: 10 species

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
- `te_genome_ectopic`
- `te_genome_ectopic_organismal`
- `te_genome_morphology`

All observed-only and `primary_phylofill` panel species sets are currently
identical.

## Tree Inputs

There are two tree-related entry points in use:

- `input_data/phylogeny/desmo900dated_test.tre`
  Used by the path-model scaffold for the fitted comparative models.
- `results/phylogeny/processed_phylogeny.nwk`
  Used by the phylogenetic trait-imputation builder.

That distinction matters because the imputation tree currently excludes some
species that still exist in the broader modeling tree.

## Current Model Families And Winners

### 1. TE only

Files:

- `results/te_genome_primary_mediumplus_model_ranking.csv`
- `results/te_genome_primary_strict_body_model_ranking.csv`

Winner:

- `mediated_evenness`

Interpretation:

- the `ltr_balance -> te_evenness` path is stable
- the `te_evenness -> gs` edge is weak and changes sign across subsets

### 2. TE + ectopic

Files:

- `results/te_genome_ectopic_primary_mediumplus_model_ranking.csv`
- `results/te_genome_ectopic_primary_strict_body_model_ranking.csv`

Winner:

- `ectopic_only`

Interpretation:

- `ectopic_index -> gs` is stable and negative
- this is the strongest ectopic-only story before organismal competition

### 3. TE + organismal

Files:

- `results/te_genome_organismal_primary_mediumplus_model_ranking.csv`
- `results/te_genome_organismal_primary_strict_body_model_ranking.csv`

Winner:

- `body_size_additive`

Interpretation:

- `body_size -> gs` is positive
- `aquaticity` does not rank well
- body size is the only organismal covariate that currently matters enough to
  compete with the TE-only baseline

### 4. TE + ectopic + organismal

Files:

- `results/te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv`
- `results/te_genome_ectopic_organismal_primary_strict_body_model_ranking.csv`

Winner:

- `te_body_size_baseline`

Interpretation:

- once `body_size` is allowed to compete with `ectopic_index`, the winning model
  drops `ectopic_index`
- the current best-supported combined story is TE-evenness plus body size, not
  ectopic plus body size

This is currently the most important substantive modeling result added after the
trait rescue phase.

### 5. TE + genome + morphology

Files:

- `results/te_genome_morphology_primary_mediumplus_model_ranking.csv`
- `results/te_genome_morphology_primary_strict_body_model_ranking.csv`

Winner:

- medium-plus: `te_evenness_path`
- strict-body: do not interpret as stable; `n = 9` is too small

Interpretation:

- keep morphology families as planning/sensitivity analyses
- do not treat the strict-body morphology ranking as robust

## What To Trust Most Right Now

Most trustworthy:

- panel definitions and species counts
- source-backed organismal table
- TE-only and TE+ectopic panel stability
- body-size signal in the organismal-augmented families

More provisional:

- morphology-linked causal interpretation
- phylogenetic imputation as anything beyond sensitivity analysis
- any interpretation that relies on the strict-body morphology subset

## If You Need To Re-Enter Later

Read these in order:

1. `PATH_ANALYSIS_STATUS.md`
2. `SESSION_HANDOFF.md`
3. `data/derived/analysis_panel_summary.csv`
4. `data/derived/phylofill_panel_comparison.csv`
5. the ranking CSVs for:
   - `te_genome`
   - `te_genome_ectopic`
   - `te_genome_organismal`
   - `te_genome_ectopic_organismal`

## The Next Logical Step

If returning to this analysis later, the next step is not more data rescue.

It is:

- build a compact cross-family summary table
- write the results-language interpretation
- decide how to frame the shift from ectopic-dominant models to body-size-plus-TE
  models once organismal covariates are included

If more data work is ever needed later, the most useful remaining targets are:

- independent final genome-size estimates
- expanded linked morphology coverage
- possibly a cleaner tree for phylogenetic imputation coverage

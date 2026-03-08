# Path Analysis Workspace

This folder is a temporary, consolidated workspace for the phylogenetic path analysis section of the Desmognathus chapter and manuscript. It keeps the current planning documents, dataset assembly code, and model scaffold in one place so the work does not drift across the legacy analysis directories.

If returning after a context switch or changing focus within the project, start with `PATH_ANALYSIS_STATUS.md`.

## What Is Here

- `PATH_ANALYSIS_STATUS.md`
  The current high-level map of the workspace, pipeline, panel definitions, and best-supported model results. This should be the first file to read when restarting or shifting to a different task within the same analysis.
- `README.md`
  This overview and the current staged plan.
- `SESSION_HANDOFF.md`
  Secondary restart document with chronological context.
- `DATA_DICTIONARY.md`
  Variable definitions, preferred observed proxies, and current caveats.
- `TE_DATA_AUDIT.md`
  Current audit status of the frozen TE source tables and derived path-analysis TE products, including validated checks and remaining warnings.
- `TE_DIVERSITY_CANONICALIZATION.md`
  Non-destructive reconstruction note for the TE diversity summary tables, including the audited path from threshold-grid outputs to the current canonical diversity snapshots.
- `TE_PROVENANCE_AUDIT.md`
  Upstream provenance map for the frozen TE summary tables, including the now-resolved diversity-summary reconstruction path.
- `INPUT_PREPARATION_PLAN.md`
  Detailed plan for literature trait mining, taxonomy crosswalks, and TE feature engineering.
- `SOURCE_TRACKING.md`
  Provenance rules and source-manifest workflow for publication-grade traceability.
- `CANDIDATE_MODELS.md`
  The initial DAG families to compare with `phylopath`.
- `scripts/build_master_dataset.py`
  Builds overlap-ready species tables from `Desmognathus_TE` plus `~/Projects/cellprofiler_test`.
- `scripts/path_model_scaffold.R`
  Prepares transformed analysis inputs and defines the current candidate model sets.
- `scripts/prepare_te_features.py`
  Freezes TE composition, diversity, turnover, and ectopic proxies into a reusable feature table.
- `scripts/build_canonical_diversity_tables.py`
  Rebuilds scratch candidates for the canonical TE diversity summary tables and audits them against the frozen `results/data/` snapshots without overwriting those files.
- `scripts/prepare_external_morphometrics.py`
  Summarizes raw external morphometric appendices into source-linked species tables.
- `scripts/prepare_nc_biodiversity_traits.py`
  Parses local NC Biodiversity Project HTML snapshots into a source-linked species trait table.
- `scripts/prepare_curated_organismal_traits.py`
  Collapses external trait sources into a manuscript-facing organismal trait table with retained provenance.
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

- TE species in this repo: `34`
- Genome-estimate species in `cellprofiler_test`: `31`
- TE + genome overlap: `27`
- TE + genome + ectopic overlap: `24`
- Genome + linked morphology overlap: `21`
- TE + genome + linked morphology overlap: `18`
- TE + genome + ectopic + linked morphology overlap: `16`

## Important Caveat

The current `genome_size_pg` values imported from `cellprofiler_test` come from the current primary CellProfiler bundle and are still provisional for path-analysis purposes. In particular, the present `primary_genome_pg` values are area-derived estimates, so any model that simultaneously treats current genome size and nucleus size as separate causal variables should be interpreted as planning and sensitivity work, not the final manuscript result.

For now:

- `TE -> genome size` models are the cleanest first target.
- `TE -> genome size <- ectopic / DNA loss` mechanism models are the next layer.
- `genome size -> nucleus size -> cell size` models are scaffolded here, but should wait for the independent final genome-size estimates before formal interpretation.

## Recommended Analysis Order

1. Build the consolidated species tables.
2. Start with the `te_genome` family.
3. Add ectopic recombination in `te_genome_ectopic`.
4. Hold the morphology chain as a staged plan until the final genome estimates are independent enough to defend `genome -> nucleus` causality.

## How To Use

Activate the project R environment first:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate Dusky
```

Install the path-analysis packages into that environment:

```bash
conda install -n Dusky -c conda-forge -c bioconda -y r-igraph r-tidygraph r-graphlayouts r-ggraph r-ggm bioconductor-graph r-phylolm
Rscript -e 'envlib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library"); .libPaths(envlib); options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("phylopath", lib = envlib)'
```

The scaffold script prefers the active conda R library when `CONDA_PREFIX` is set, so it does not accidentally load incompatible user-level packages.

Build the merged datasets:

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

The main analysis-panel products are:

- `path_analysis/data/derived/analysis_panel_summary.csv`
  Panel-level counts and species lists.
- `path_analysis/data/derived/analysis_species_readiness.csv`
  Species-level eligibility and exclusion reasons for every panel.
- `path_analysis/data/derived/source_traceability_gaps.csv`
  Traceability audit report; this should stay empty for publication use.

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
Rscript path_analysis/scripts/path_model_scaffold.R --summary-only
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --summary-only
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_strict_body --summary-only
```

Once `phylopath` is installed, run a family:

```bash
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome
```

To run the same family on a panel-defined sensitivity subset without overwriting the legacy outputs:

```bash
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome --panel te_genome_primary_mediumplus
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_organismal --panel te_genome_organismal_primary_strict_body
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic_organismal --panel te_genome_ectopic_organismal_primary_mediumplus
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --panel te_genome_ectopic_primary_strict_body
```

When `--panel` is supplied, the scaffold reads `data/derived/panels/<panel>.csv` and writes results with the panel name as the filename prefix.

Current organismal-family result:

- `te_genome_organismal_primary_mediumplus` currently prefers `body_size_additive` over the TE-only baseline, while `aquaticity`-heavy models rank poorly.
- `te_genome_ectopic_organismal_primary_mediumplus` and `..._strict_body` currently prefer `te_body_size_baseline`, which means the combined family does not retain `ectopic_index` as a winning predictor once body size is allowed to compete.

## External Inputs

The dataset builder currently expects these `cellprofiler_test` outputs by default:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv`
- `/home/jake/Projects/cellprofiler_test/output/publication_analysis/species_morphology_summary.csv`

These defaults can be overridden on the command line if the bundle location changes later.

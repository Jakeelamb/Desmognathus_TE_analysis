# Path Analysis Workspace

This folder is a temporary, consolidated workspace for the phylogenetic path analysis section of the Desmognathus chapter and manuscript. It keeps the current planning documents, dataset assembly code, and model scaffold in one place so the work does not drift across the legacy analysis directories.

## What Is Here

- `README.md`
  This overview and the current staged plan.
- `DATA_DICTIONARY.md`
  Variable definitions, preferred observed proxies, and current caveats.
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
- `scripts/prepare_external_morphometrics.py`
  Summarizes raw external morphometric appendices into source-linked species tables.
- `scripts/prepare_nc_biodiversity_traits.py`
  Parses local NC Biodiversity Project HTML snapshots into a source-linked species trait table.
- `scripts/prepare_curated_organismal_traits.py`
  Collapses external trait sources into a manuscript-facing organismal trait table with retained provenance.
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

Inspect the overlap summary:

```bash
sed -n '1,120p' path_analysis/data/derived/dataset_overlap_summary.csv
```

Preview the model families without requiring `phylopath`:

```bash
Rscript path_analysis/scripts/path_model_scaffold.R --summary-only
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome_ectopic --summary-only
```

Once `phylopath` is installed, run a family:

```bash
Rscript path_analysis/scripts/path_model_scaffold.R --family te_genome
```

## External Inputs

The dataset builder currently expects these `cellprofiler_test` outputs by default:

- `/home/jake/Projects/cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv`
- `/home/jake/Projects/cellprofiler_test/output/publication_analysis/species_morphology_summary.csv`

These defaults can be overridden on the command line if the bundle location changes later.

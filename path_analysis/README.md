# Path Analysis Workspace

This folder is a temporary, consolidated workspace for the phylogenetic path analysis section of the Desmognathus chapter and manuscript. It keeps the current planning documents, dataset assembly code, and model scaffold in one place so the work does not drift across the legacy analysis directories.

## What Is Here

- `README.md`
  This overview and the current staged plan.
- `DATA_DICTIONARY.md`
  Variable definitions, preferred observed proxies, and current caveats.
- `CANDIDATE_MODELS.md`
  The initial DAG families to compare with `phylopath`.
- `scripts/build_master_dataset.py`
  Builds overlap-ready species tables from `Desmognathus_TE` plus `~/Projects/cellprofiler_test`.
- `scripts/path_model_scaffold.R`
  Prepares transformed analysis inputs and defines the current candidate model sets.
- `data/derived/`
  Generated master tables and staged analysis datasets.
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
```

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

# Desmognathus TE Analysis

Comprehensive analysis of transposable element evolution across *Desmognathus*
salamanders. The current panel contract is **TE34** for vetted genomic-resource
descriptions, **Cell21** for linked-cell descriptions, and **path18** only for
their exact tree-aligned intersection. It includes TE classification, divergence,
phylogenetic comparative methods, paired-LTR divergence, a terminal:internal
LTR deletion-footprint proxy, and diversity/ordination metrics.

Start here for the active TE workflow. Use `path_analysis/README.md` for the phylogenetic path-analysis workspace and `TODO.md` for the current analysis queue.

## Publication-Audit Release Boundary

There is no single denominator for every analysis. Corrected, non-destructive
audits are under `results/data/corrected/`, `results/figures/corrected/`, and
`plans/publication-readiness-deep-audit/`. The primary collaborator-review
surface is the eight executed notebooks documented in
[`notebooks/research_review/README.md`](notebooks/research_review/README.md).
Its TE and LTR figures are regenerated from audited current data with the
historical R/ggplot grammar. Its hashed frozen-input registry makes the static
review bundle re-executable without the ignored upstream result tree; only the
large interactive microscopy galleries remain optional and machine-local.
The review directory contains one canonical notebook per analysis domain,
including a dedicated genome-size/IOD notebook. The combined analysis18
workbench remains a secondary cross-domain provenance artifact, not the meeting
entry point.

Current permitted interpretation:

- TE diversity and order-level CLR/phylogenetic PCA are approved for
  descriptive use with the recorded mass and tree sensitivities.
- LTR terminal:internal depth is a deletion-footprint/mapping proxy, not a
  measured ectopic-recombination rate.
- linked cell/nucleus morphology is an upper-tail sensitivity analysis; the
  segmentation models are not validated on held-out final-panel species.
- quality-matched nuclear IOD is a relative image-intensity proxy, not an
  absolute C-value. The separate *D. fuscus*-anchored rescaling is retained only
  as a conditional genome-size estimate for the descriptive Shannon scatter; it
  is not promoted as a direct C-value or causal path-analysis input.
- corrected phylogenetic path models are exploratory association-model
  sensitivity only; they are not approved causal genome-size results.

See `results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv`
for the machine-readable claim boundary.

## Quick Start

```bash
# Verify setup
scripts/run_in_dusky.sh python verify_setup.py

# Open the prepared, executed collaborator-review bundle
scripts/run_in_dusky.sh jupyter lab notebooks/research_review

# Run a core processing script
scripts/run_in_dusky.sh python scripts/processing/dnaPipe.py
```

## Project Structure

```
./
├── input_data/                    # Raw input data (not tracked in git)
│   ├── dnaPipeTE/                 # dnaPipeTE classification files
│   ├── repeatmasker/              # RepeatMasker .align files
│   ├── phylogeny/                 # Phylogenetic tree files
│   ├── ectopic_recombination/     # LTR domain data
│   └── lookup_table.txt           # Species ID mapping
│
├── results/                       # Analysis outputs (upstream products ignored)
│   ├── data/                      # Processed CSV files
│   ├── figures/                   # Generated visualizations
│   ├── data/research_review/      # Tracked compact review tables/manifest
│   ├── figures/research_review/   # Tracked current-data PNG/PDF figure bundle
│   └── reports/                   # Generated prose reports
│
├── interim/                       # Intermediate processing files
├── path_analysis/                 # Path-analysis workspace and derived tables
│
├── scripts/
│   ├── config.py                  # Centralized path configuration
│   ├── processing/                # Data processing scripts
│   │   ├── dnaPipe.py             # dnaPipeTE data processing
│   │   ├── repeatmask.py          # RepeatMasker data processing
│   │   ├── ec.py                  # Ectopic recombination analysis
│   │   ├── divergence.py          # Divergence calculations
│   │   ├── diversity_stats.py     # Canonical diversity writer
│   │   ├── pca.R                  # Canonical compositional PCA
│   │   ├── pca_utils.R            # Shared compositional PCA utilities
│   │   ├── phylogenetic_pca_analysis.R
│   │   ├── clean_tree_phylo.R     # Phylogeny cleaning
│   │   └── analyze_phylogenetic_signal.R
│   └── visualization/             # Plotting scripts
│       ├── divergence.R
│       ├── hierarchical_donut_TE_diversity.R
│       └── plot_*.R
│
├── paths.yaml                     # Path configuration
├── DATA_MANIFEST.yml              # Share-facing local/tracked/generated data contract
├── Dusky.yml                      # Conda environment specification
├── verify_setup.py                # Setup verification script
└── README.md
```

## Environment Setup

### Using Conda (Recommended)

```bash
# Install Miniconda if needed. For audited environments, verify the installer
# checksum from Anaconda's release page before running it.
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda3

# Create the Dusky environment
source $HOME/miniconda3/etc/profile.d/conda.sh
conda env create -f Dusky.yml

# Verify installation
scripts/run_in_dusky.sh python verify_setup.py
```

Use `scripts/run_in_dusky.sh <command>` for repo commands even if you normally
activate Conda interactively. It prepends `$CONDA_PREFIX/bin` before execution,
which prevents shell PATH leakage from resolving system Python or R instead of
the `Dusky` environment.

### Environment Contents

The Dusky environment includes:
- **Python 3.8** with pandas, numpy, matplotlib, seaborn, biopython, pysam
- **R 4.x** with tidyverse, ggplot2, ape, phytools, factoextra
- **Bioinformatics tools**: samtools, bowtie2, bedtools, RepeatMasker, trf, tesorter

## Input Data

Input data is not tracked in git due to size. Required files:

| Directory | Contents | Source |
|-----------|----------|--------|
| `input_data/dnaPipeTE/` | `SRX*_reads_per_component_and_annotation` | dnaPipeTE output |
| `input_data/repeatmasker/` | `SRX*_Trinity.align` | RepeatMasker output |
| `input_data/phylogeny/` | `desmo900dated_test.tre` | Phylogenetic tree |
| `input_data/ectopic_recombination/` | `GCA_*_tabout.csv`, coverage files | TEsorter output |
| `input_data/lookup_table.txt` | Species-SRA-Genome mapping | Manual |

See `DATA_MANIFEST.yml` for the share-facing distinction between local
required inputs, imported CellProfiler snapshots, tracked small external
snapshots, and generated outputs.

## Processing Workflows

### 1. dnaPipeTE Processing

Classifies TEs from dnaPipeTE output into Class/Order/Superfamily.

```bash
scripts/run_in_dusky.sh python scripts/processing/dnaPipe.py
```

**Outputs:**
- `results/data/dnaPipeTE_merged_classifications.csv`
- `results/data/dnaPipeTE_class_breakdown.csv`
- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`

### 2. RepeatMasker Processing

Parses RepeatMasker alignment files and merges with dnaPipeTE classifications.

```bash
scripts/run_in_dusky.sh python scripts/processing/repeatmask.py
```

**Outputs:**
- `results/data/merged_repeatmasker_data.csv`
- `results/data/repeatmasker_detailed_classification_combined.csv`

### 3. Historical Terminal:Internal Proxy Workflow

The historical script calculates LTR terminal:internal depth ratios. These
ratios are not a validated ectopic-recombination rate and must not be described
as direct solo-LTR formation. They are retained for provenance.

```bash
scripts/run_in_dusky.sh python scripts/processing/ec.py
```

**Outputs:**
- `results/data/ectopic_recombination_master.csv`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

The analysis visualization step for this branch is:

```bash
scripts/run_in_dusky.sh Rscript scripts/visualization/plot_ectopic_recombination.R
```

That script reads the canonical filtered ectopic table directly and writes:
- `results/figures/ectopic_recombination/ectopic_ratio_violin_log10.png`
- `results/data/ectopic_recombination_species_summary.csv`
- `results/data/ectopic_recombination_species_tests.txt`

The release-facing corrected branch is:

```bash
scripts/run_in_dusky.sh python scripts/processing/build_corrected_ectopic_recombination.py
```

It restricts analysis to the final panel, retains zero-depth positions, requires
exact ≥5/6-domain elements, reports coverage and robust/influence summaries,
and writes only proxy-labeled products under `results/data/corrected/ectopic/`.

### 4. Divergence Analysis

Calculates sequence divergence metrics grouped by TE classification.

```bash
scripts/run_in_dusky.sh python scripts/processing/divergence.py
```

**Outputs:**
- `interim/divergence/class/*.csv`
- `interim/divergence/order/*.csv`
- `interim/divergence/superfamily/*.csv`

### 5. Diversity Statistics

Calculates Shannon, Simpson, and Pielou's evenness indices.

```bash
scripts/run_in_dusky.sh python scripts/processing/diversity_stats.py
```

**Outputs:**
- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`
- `results/data/comparison_diversity_order_stats_granular_0_5pct.csv`
- `results/data/comparison_diversity_superfamily_stats_granular_0_5pct.csv`
- `results/data/long_format_diversity_order_stats.csv`
- `results/data/long_format_diversity_superfamily_stats.csv`

For non-destructive verification of the existing canonical diversity snapshots,
use:

```bash
scripts/run_in_dusky.sh python path_analysis/scripts/build_canonical_diversity_tables.py
```

That script writes scratch candidates and an audit report under
`path_analysis/data/derived/` without overwriting the current `results/data/`
files. See `path_analysis/TE_DIVERSITY_CANONICALIZATION.md`.

### 6. Phylogeny Cleaning

Cleans and prepares phylogenetic tree for analysis.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/clean_tree_phylo.R
```

**Outputs:**
- `results/data/desmo900dated_test_cleaned_phylo.tre`
- `results/figures/phylogeny/rectangular_phylogeny.png`

### 7. PCA Analysis

Runs the canonical compositional PCA workflow on the current TE breakdown
tables. This is intended as a supplementary ordination layer, not the primary
comparative predictor block. See `TE_PCA_METHODS.md` for the exact filtering,
zero-replacement, and CLR rules.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/pca.R
```

**Outputs:**
- `results/tables/pca/te_pca_analysis_manifest.csv`
- `results/tables/pca/*_clr_matrix.csv`
- `results/tables/pca/*_scores.csv`
- `results/tables/pca/*_loadings.csv`
- `results/figures/pca/*_scores_pc1_pc2.png`
- `results/figures/pca/*_scree_plot.png`

### 8. Phylogenetic PCA

Runs supplementary phylogenetic PCA from the exact CLR matrices written by the
standard PCA workflow.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/phylogenetic_pca_analysis.R
```

**Outputs:**
- `results/tables/pca/te_ppca_analysis_manifest.csv`

## Generated Outputs

Large rebuild products stay out of git. Scripts write tables and figures under
`results/`, intermediate files under `interim/`, and generated prose reports
under `results/reports/`. The tracked repo should contain source code,
configuration, curated input templates, and small canonical derived tables, not
local generated report artifacts.

- `results/tables/pca/*_ppca_scores.csv`
- `results/tables/pca/*_ppca_loadings.csv`
- `results/figures/pca/*_ppca_scores_pc1_pc2.png`
- `results/figures/pca/*_ppca_phylomorphospace.png`

### 9. Trait Evolution Modeling

Estimates phylogenetic signal for TE traits using Pagel's lambda, Blomberg's K with a fixed RNG seed, PIC dispersion, and root-state summaries.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/trait_evolution.R
```

**Outputs:**
- `results/data/trait_evolution/evolutionary_model_comparison.csv`
- `results/figures/trait_evolution/ancestral_*.png`

### 10. LTR Sequence-Divergence and Age Audit

Audits the paired-LTR substrate and, when local genome assemblies are present,
computes true 5'/3' LTR sequence divergence directly from assembly coordinates.
The default validated output is divergence rather than absolute age in years,
because no substitution rate is imposed automatically. When regenerated locally,
`results/reports/LTR_AGE_AUDIT.md` records workflow status.
`LTR_SUBSTITUTION_RATE_CALIBRATION.md` records the primary-literature calibration
window, and `LTR_AGE_INSIGHTS.md` records the current interpretation note.

```bash
# Resolve or download the accession-linked genome assemblies when needed
scripts/run_in_dusky.sh python scripts/processing/fetch_genome_assemblies.py --download

# Then run the local divergence audit/estimation
scripts/run_in_dusky.sh python scripts/processing/ltr_age_estimation.py
```

**Outputs:**
- `results/data/ltr_age/ltr_age_readiness_by_species.csv`
- `results/data/ltr_age/ltr_age_readiness_overview.csv`
- `results/data/ltr_age/ltr_age_candidate_inventory.csv`
- `results/data/ltr_age/ltr_age_pairwise_divergence.csv`
- `results/data/ltr_age/ltr_age_species_summary.csv`
- `results/data/ltr_age/ltr_age_scaffold_extraction_summary.csv`
- `results/reports/LTR_AGE_AUDIT.md`

To convert those divergence outputs into literature-backed age sensitivities
without changing the canonical divergence table:

```bash
scripts/run_in_dusky.sh python scripts/processing/build_ltr_substitution_rate_calibration.py
```

**Calibration outputs:**
- `results/data/ltr_age/ltr_substitution_rate_candidates.csv`
- `results/data/ltr_age/ltr_age_calibration_summary.csv`
- `results/data/ltr_age/ltr_age_species_summary_calibrated.csv`
- `results/data/ltr_age/ltr_age_pairwise_divergence_calibrated.csv`
- `LTR_SUBSTITUTION_RATE_CALIBRATION.md`

Import the rebuilt LTR-history layer into the path-analysis workspace with:

```bash
scripts/run_in_dusky.sh python path_analysis/scripts/prepare_ltr_history_features.py
```

### 11. PGLS Regression

Phylogenetic Generalized Least Squares regression for exploratory compositional
screening between TE orders and superfamilies using `caper::pgls()` on
CLR-transformed compositions with ML lambda estimation and BH-corrected
p-values. See `COMPARATIVE_INFERENCE_AUDIT.md` for scope and interpretation
limits.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/pgls_analysis.R
```

**Outputs:**
- `results/data/pgls/pgls_order_pairwise.csv`
- `results/figures/pgls/pgls_volcano_plot.png`

### 12. PERMANOVA Group Comparisons

Formal statistical tests for TE compositional differences between curated
phylogenetic clades using `vegan::adonis2()` with Bray-Curtis and
CLR-Euclidean distances, beta dispersion tests, and PCoA ordination. See
`COMPARATIVE_INFERENCE_AUDIT.md` for the current audited interpretation.

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/permanova_analysis.R
```

**Outputs:**
- `results/data/permanova/permanova_summary.csv`
- `results/figures/permanova/pcoa_*_bray.png`

## Configuration

Path configuration is centralized in the repository-root `paths.yaml`. Python
scripts use `scripts/config.py`, and R scripts use `scripts/R/path_config_utils.R`
to resolve paths relative to the repository root.

```python
# Python usage
from config import paths, PROJECT_ROOT

input_dir = paths.input_data.dnaPipeTE
output_dir = paths.results.data
```

```r
# R usage
source("scripts/processing/pca_utils.R")
config <- load_config()
data_dir <- config$results$data
```

## Git Management

Large data files are excluded from git tracking:
- `input_data/` - Raw input data
- `results/` - Generated outputs
- `interim/` - Intermediate files

Tracked files include scripts, configuration, documentation, source manifests,
and small canonical derived/audit tables under `path_analysis/`. Large raw
inputs and local generated rebuild products stay out of git.

## Troubleshooting

### Conda not found
```bash
source $HOME/miniconda3/etc/profile.d/conda.sh
```

### Import errors
Run commands through the Dusky wrapper so `$CONDA_PREFIX/bin` is first on PATH:
```bash
scripts/run_in_dusky.sh python -c "import pandas, matplotlib, scipy"
```

### Verify setup
```bash
scripts/run_in_dusky.sh python verify_setup.py
```

## License

MIT License

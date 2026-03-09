# Desmognathus TE Analysis

Comprehensive analysis of transposable element evolution across 34 Desmognathus salamander species. Includes genome-wide TE classification, divergence quantification, phylogenetic comparative methods (PGLS, PERMANOVA, BM/OU modeling), a sequence-based paired-LTR divergence branch for insertion-age inference, ectopic recombination analysis, and diversity metrics — spanning 12 analysis stages with 30+ processing and visualization scripts in Python and R.

For the current paper-facing frozen snapshot, start with `PAPER_FREEZE_MANIFEST.md`.

## Quick Start

```bash
# Activate the conda environment
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate Dusky

# Verify setup
python verify_setup.py

# Rebuild the tracked paper freeze manifest
python scripts/processing/build_paper_freeze_manifest.py

# Run a processing script
python scripts/processing/dnaPipe.py
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
├── results/                       # Analysis outputs (not tracked in git)
│   ├── data/                      # Processed CSV files
│   └── figures/                   # Generated visualizations
│
├── interim/                       # Intermediate processing files
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
├── config/paths.yaml              # Path configuration
├── Dusky.yml                      # Conda environment specification
├── verify_setup.py                # Setup verification script
└── README.md
```

## Environment Setup

### Using Conda (Recommended)

```bash
# Install Miniconda if not already installed
curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda3

# Create the Dusky environment
source $HOME/miniconda3/etc/profile.d/conda.sh
conda env create -f Dusky.yml

# Activate the environment
conda activate Dusky

# Verify installation
python verify_setup.py
```

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

## Processing Workflows

### 1. dnaPipeTE Processing

Classifies TEs from dnaPipeTE output into Class/Order/Superfamily.

```bash
python scripts/processing/dnaPipe.py
```

**Outputs:**
- `results/data/dnaPipeTE_merged_classifications.csv`
- `results/data/dnaPipeTE_class_breakdown.csv`
- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`

### 2. RepeatMasker Processing

Parses RepeatMasker alignment files and merges with dnaPipeTE classifications.

```bash
python scripts/processing/repeatmask.py
```

**Outputs:**
- `results/data/merged_repeatmasker_data.csv`
- `results/data/repeatmasker_detailed_classification_combined.csv`

### 3. Ectopic Recombination Analysis

Analyzes LTR depth ratios to identify potential ectopic recombination.

```bash
python scripts/processing/ec.py
```

**Outputs:**
- `results/data/ectopic_recombination_master.csv`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

The paper-facing visualization step for this branch is:

```bash
Rscript scripts/visualization/plot_ectopic_recombination.R
```

That script reads the canonical filtered ectopic table directly and writes:
- `results/figures/ectopic_recombination/ectopic_ratio_violin_log10.png`
- `results/data/ectopic_recombination_species_summary.csv`
- `results/data/ectopic_recombination_species_tests.txt`

### 4. Divergence Analysis

Calculates sequence divergence metrics grouped by TE classification.

```bash
python scripts/processing/divergence.py
```

**Outputs:**
- `interim/divergence/class/*.csv`
- `interim/divergence/order/*.csv`
- `interim/divergence/superfamily/*.csv`

### 5. Diversity Statistics

Calculates Shannon, Simpson, and Pielou's evenness indices.

```bash
python scripts/processing/diversity_stats.py
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
python path_analysis/scripts/build_canonical_diversity_tables.py
```

That script writes scratch candidates and an audit report under
`path_analysis/data/derived/` without overwriting the current `results/data/`
files. See `path_analysis/TE_DIVERSITY_CANONICALIZATION.md`.

### 6. Phylogeny Cleaning

Cleans and prepares phylogenetic tree for analysis.

```bash
Rscript scripts/processing/clean_tree_phylo.R
```

**Outputs:**
- `results/data/desmo900dated_test_cleaned_phylo.tre`
- `results/figures/phylogeny/rectangular_phylogeny.png`

### 7. PCA Analysis

Runs the canonical compositional PCA workflow on the frozen TE breakdown
tables. This is intended as a supplementary ordination layer, not the primary
comparative predictor block. See `TE_PCA_METHODS.md` for the exact filtering,
zero-replacement, and CLR rules.

```bash
Rscript scripts/processing/pca.R
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
Rscript scripts/processing/phylogenetic_pca_analysis.R
```

**Outputs:**
- `results/tables/pca/te_ppca_analysis_manifest.csv`

## Paper Freeze

The repo keeps `results/` out of git, so the current paper snapshot is captured
by the tracked checksum inventories and freeze note rather than by versioning
all generated outputs directly.

Rebuild the manifest with:

```bash
python scripts/processing/build_paper_freeze_manifest.py
```

Tracked freeze artifacts:
- `PAPER_FREEZE_MANIFEST.md`
- `paper_freeze/key_file_manifest.csv`
- `paper_freeze/results_inventory.csv`
- `paper_freeze/results_summary.csv`

To build the manuscript-facing summary table and figure plan from the frozen
primary results:

```bash
python scripts/processing/build_manuscript_assets.py
```
- `results/tables/pca/*_ppca_scores.csv`
- `results/tables/pca/*_ppca_loadings.csv`
- `results/figures/pca/*_ppca_scores_pc1_pc2.png`
- `results/figures/pca/*_ppca_phylomorphospace.png`

### 9. Trait Evolution Modeling

Compares Brownian Motion vs Ornstein-Uhlenbeck models for TE trait evolution using `geiger::fitContinuous()` with AICc model selection and ancestral state reconstruction via `phytools::fastAnc()`.

```bash
Rscript scripts/processing/trait_evolution.R
```

**Outputs:**
- `results/data/trait_evolution/evolutionary_model_comparison.csv`
- `results/figures/trait_evolution/ancestral_*.png`

### 10. LTR Sequence-Divergence and Age Audit

Audits the paired-LTR substrate and, when local genome assemblies are present,
computes true 5'/3' LTR sequence divergence directly from assembly coordinates.
The default validated output is divergence rather than absolute age in years,
because no substitution rate is imposed automatically. See `LTR_AGE_AUDIT.md`
for workflow status, `LTR_DIVERGENCE_SIGNAL_AUDIT.md` for the current
biological signal summary, and `LTR_SUBSTITUTION_RATE_CALIBRATION.md` for the
primary-literature calibration window. The current interpretation note is
tracked in `LTR_AGE_INSIGHTS.md`.

```bash
# Resolve or download the accession-linked genome assemblies when needed
python scripts/processing/fetch_genome_assemblies.py --download

# Then run the local divergence audit/estimation
conda activate Dusky
python scripts/processing/ltr_age_estimation.py
```

**Outputs:**
- `results/data/ltr_age/ltr_age_readiness_by_species.csv`
- `results/data/ltr_age/ltr_age_readiness_overview.csv`
- `results/data/ltr_age/ltr_age_candidate_inventory.csv`
- `results/data/ltr_age/ltr_age_pairwise_divergence.csv`
- `results/data/ltr_age/ltr_age_species_summary.csv`
- `results/data/ltr_age/ltr_age_scaffold_extraction_summary.csv`
- `LTR_AGE_AUDIT.md`

To convert those divergence outputs into literature-backed age sensitivities
without changing the canonical divergence table:

```bash
python scripts/processing/build_ltr_substitution_rate_calibration.py
```

**Calibration outputs:**
- `results/data/ltr_age/ltr_substitution_rate_candidates.csv`
- `results/data/ltr_age/ltr_age_calibration_summary.csv`
- `results/data/ltr_age/ltr_age_species_summary_calibrated.csv`
- `results/data/ltr_age/ltr_age_pairwise_divergence_calibrated.csv`
- `LTR_SUBSTITUTION_RATE_CALIBRATION.md`

### 11. PGLS Regression

Phylogenetic Generalized Least Squares regression for exploratory compositional
screening between TE orders and superfamilies using `caper::pgls()` on
CLR-transformed compositions with ML lambda estimation and BH-corrected
p-values. See `COMPARATIVE_INFERENCE_AUDIT.md` for scope and interpretation
limits.

```bash
Rscript scripts/processing/pgls_analysis.R
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
Rscript scripts/processing/permanova_analysis.R
```

**Outputs:**
- `results/data/permanova/permanova_summary.csv`
- `results/figures/permanova/pcoa_*_bray.png`

## Configuration

Path configuration is centralized in the repository-root `paths.yaml`, with
`config/paths.yaml` retained as a compatibility mirror for older callers.
Python scripts use `scripts/config.py`, and the rebuilt R helper layer resolves
either location relative to the repository root.

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

Only scripts, configuration, and documentation are tracked.

## Troubleshooting

### Conda not found
```bash
source $HOME/miniconda3/etc/profile.d/conda.sh
```

### Import errors
Ensure you're in the Dusky environment:
```bash
conda activate Dusky
```

### Verify setup
```bash
python verify_setup.py
```

## License

MIT License

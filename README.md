# Desmognathus TE Analysis

Comprehensive analysis of transposable element evolution across 34 Desmognathus salamander species. Includes genome-wide TE classification, divergence quantification, phylogenetic comparative methods (PGLS, PERMANOVA, BM/OU modeling), LTR insertion age estimation, ectopic recombination analysis, and diversity metrics — spanning 12 analysis stages with 30+ processing and visualization scripts in Python and R.

## Quick Start

```bash
# Activate the conda environment
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate Dusky

# Verify setup
python verify_setup.py

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
├── archive/
│   └── legacy_diversity/          # Archived exploratory diversity scripts
│
├── scripts/
│   ├── config.py                  # Centralized path configuration
│   ├── processing/                # Data processing scripts
│   │   ├── dnaPipe.py             # dnaPipeTE data processing
│   │   ├── repeatmask.py          # RepeatMasker data processing
│   │   ├── ec.py                  # Ectopic recombination analysis
│   │   ├── divergence.py          # Divergence calculations
│   │   ├── diversity_stats.py     # Canonical diversity writer
│   │   ├── pca.R                  # PCA analysis
│   │   ├── pca_utils.R            # Shared PCA utilities
│   │   ├── phylogenetic_pca_analysis.R
│   │   ├── clean_tree_phylo.R     # Phylogeny cleaning
│   │   └── analyze_phylogenetic_signal.R
│   └── visualization/             # Plotting scripts
│       ├── divergence.R
│       ├── hierarchical_donut_TE_diversity.R
│       └── plot_*.R
│
├── paths.yaml                     # Path configuration
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
- `results/figures/rectangular_phylogeny.png`

### 7. PCA Analysis

Performs PCA on TE composition data.

```bash
Rscript scripts/processing/pca.R
```

**Outputs:**
- `results/figures/*_pca_scatter_plot.png`
- `results/figures/*_scree_plot.png`

### 8. Phylogenetic PCA

PCA with phylogenetic correction and phylomorphospace visualization.

```bash
Rscript scripts/processing/phylogenetic_pca_analysis.R
```

**Outputs:**
- `results/figures/*_pPCA_phylomorphospace_plot.png`

### 9. Trait Evolution Modeling

Compares Brownian Motion vs Ornstein-Uhlenbeck models for TE trait evolution using `geiger::fitContinuous()` with AICc model selection and ancestral state reconstruction via `phytools::fastAnc()`.

```bash
Rscript scripts/processing/trait_evolution.R
```

**Outputs:**
- `results/data/trait_evolution/evolutionary_model_comparison.csv`
- `results/figures/trait_evolution/ancestral_*.png`

### 10. LTR Insertion Age Estimation

Estimates LTR retrotransposon insertion times from intra-element (5' vs 3' LTR) divergence, converted to age via substitution rate.

```bash
python scripts/processing/ltr_age_estimation.py
```

**Outputs:**
- `results/data/ltr_age/ltr_insertion_ages.csv`
- `results/figures/ltr_age/ltr_age_by_species.png`

### 11. PGLS Regression

Phylogenetic Generalized Least Squares regression for phylogenetically-corrected pairwise correlations between TE orders and superfamilies using `caper::pgls()` with ML lambda estimation and BH-corrected p-values.

```bash
Rscript scripts/processing/pgls_analysis.R
```

**Outputs:**
- `results/data/pgls/pgls_order_pairwise.csv`
- `results/figures/pgls/pgls_volcano_plot.png`

### 12. PERMANOVA Group Comparisons

Formal statistical tests for TE compositional differences between phylogenetic clades using `vegan::adonis2()` with Bray-Curtis and CLR-Euclidean distances, beta dispersion tests, and PCoA ordination.

```bash
Rscript scripts/processing/permanova_analysis.R
```

**Outputs:**
- `results/data/permanova/permanova_summary.csv`
- `results/figures/permanova/pcoa_*_bray.png`

## Configuration

Path configuration is centralized in `paths.yaml`. Python scripts use `scripts/config.py` and R scripts use `scripts/processing/pca_utils.R` to load paths consistently.

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

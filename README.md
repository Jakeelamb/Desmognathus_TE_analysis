# Desmognathus TE Analysis

Analysis of transposable elements (TEs) in Desmognathus salamander genomes.

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
├── scripts/
│   ├── config.py                  # Centralized path configuration
│   ├── processing/                # Data processing scripts
│   │   ├── dnaPipe.py             # dnaPipeTE data processing
│   │   ├── repeatmask.py          # RepeatMasker data processing
│   │   ├── ec.py                  # Ectopic recombination analysis
│   │   ├── divergence.py          # Divergence calculations
│   │   ├── diversity.py           # Diversity metrics
│   │   ├── diversity_stats.py     # Diversity statistics
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

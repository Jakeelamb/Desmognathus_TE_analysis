# Desmognathus TE Analysis - Complete Usage Guide

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Quick Start](#quick-start)
3. [Processing Pipeline](#processing-pipeline)
4. [Visualization Pipeline](#visualization-pipeline)
5. [Troubleshooting](#troubleshooting)

---

## Environment Setup

### Current System Status
- Python: System Python 3.13.7 available at `/usr/bin/python3`
- R: R 4.5.1 installed with required packages (tidyverse, ggplot2, viridis, dplyr, gridExtra)
- Conda: Not currently installed or not in PATH

### Option 1: Install Python Packages via pip (Recommended if Conda unavailable)

```bash
# Create a virtual environment
cd /home/jake/Projects/Desmognathus_TE
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate

# Install required packages
pip install pandas numpy matplotlib seaborn biopython tqdm pyyaml
```

After installation, always activate the environment before running scripts:
```bash
source /home/jake/Projects/Desmognathus_TE/venv/bin/activate
```

### Option 2: Install/Setup Conda Environment (As per original design)

If you want to use the original Conda setup:

```bash
# Install miniconda if not installed
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh

# Create the Dusky environment
conda env create -f Dusky.yml

# Activate the environment
conda activate Dusky
```

### Check R Packages

The following R packages should be installed (most are already available):
- tidyverse, ggplot2, viridis, dplyr, gridExtra
- factoextra, compositions, ggrepel, cluster, readr, tidyr
- ape, phytools, geiger (for phylogenetic analyses)

To install missing R packages:
```R
install.packages(c("factoextra", "compositions", "ggrepel", "ape", "phytools", "geiger"))
```

---

## Quick Start

### Verify Your Installation

```bash
# Test Python installation
python3 -c "import pandas, numpy, matplotlib, seaborn, tqdm; print('Python packages OK')"

# Test R installation
Rscript -e "library(tidyverse); library(ggplot2); print('R packages OK')"
```

---

## Processing Pipeline

All processing scripts are in `scripts/processing/`. Run these in order:

### 1. Ectopic Recombination Analysis

**Script:** `ec.py`

**Purpose:** Process ectopic recombination data

**Output:**
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
- `results/data/ectopic_recombination_master.csv`

**Command:**
```bash
cd /home/jake/Projects/Desmognathus_TE
python3 scripts/processing/ec.py
```

---

### 2. dnaPipeTE Processing

**Script:** `dnaPipe.py`

**Purpose:** Process dnaPipeTE data to classify transposable elements

**Output:**
- `results/data/dnaPipeTE_class_breakdown.csv`
- `results/data/dnaPipeTE_merged_classifications.csv`
- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`

**Command:**
```bash
python3 scripts/processing/dnaPipe.py
```

**Note:** This script uses multiprocessing for efficiency. It processes data from `input_data/dnaPipeTE/` and uses the lookup table at `input_data/lookup_table.txt`.

---

### 3. RepeatMasker Processing

**Script:** `repeatmask.py`

**Purpose:** Process RepeatMasker output files (.align files)

**Output:**
- `results/data/repeatmasker_detailed_classification_combined.csv`

**Command:**
```bash
python3 scripts/processing/repeatmask.py
```

**Input:** Reads .align files from `input_data/repeatmasker/`

---

### 4. Phylogenetic Tree Cleaning

**Script:** `clean_tree_phylo.R`

**Purpose:** Clean and prepare phylogenetic tree for analysis

**Output:**
- `results/data/desmo900dated_test_cleaned.tre`
- `results/figures/rectangular_phylogeny.png`

**Command:**
```bash
Rscript scripts/processing/clean_tree_phylo.R
```

**Input:** `input_data/phylogeny/desmo900dated_test.tre`

---

### 5. Divergence Analysis

**Script:** `divergence.py`

**Purpose:** Calculate divergence statistics from RepeatMasker output

**Output:**
- `results/data/divergence/divergence_summary_statistics_by_species.csv`
- `interim/divergence/` (intermediate files per TE type)

**Command:**
```bash
python3 scripts/processing/divergence.py
```

---

### 6. Diversity Statistics

**Script:** `diversity_stats.py`

**Purpose:** Calculate diversity metrics (Shannon, Simpson, Pielou) for TE
order and superfamily compositions and write the canonical repo-level diversity
summary tables.

**Output:**
- `results/data/comparison_diversity_order_stats_granular_0_5pct.csv`
- `results/data/comparison_diversity_superfamily_stats_granular_0_5pct.csv`
- `results/data/long_format_diversity_order_stats.csv`
- `results/data/long_format_diversity_superfamily_stats.csv`
- `results/data/diversity_superfamily_stats.csv`
- `results/data/diversity_order_stats.csv`

**Command:**
```bash
python3 scripts/processing/diversity_stats.py
```

---

### 7. PCA Analysis

**Script:** `pca.R`

**Purpose:** Perform principal component analysis on TE composition

**Output:**
- `results/figures/Superfamily_Diversity_CLR_pca_scatter_plot.png`
- `results/figures/Order_Diversity_CLR_pca_scatter_plot.png`
- Cluster analysis outputs (elbow plots, silhouette plots)
- `results/data/cluster_assignments_superfamily.csv`
- `results/data/cluster_assignments_order.csv`

**Command:**
```bash
Rscript scripts/processing/pca.R
```

**Requirements:** Needs output from dnaPipe.py (the breakdown CSV files)

---

### 8. Phylogenetic PCA Analysis

**Script:** `phylogenetic_pca_analysis.R`

**Purpose:** Perform phylogenetic PCA with phylogeny overlay (phylomorphospace)

**Output:**
- `results/figures/Superfamily_Diversity_CLR_pPCA_phylomorphospace_plot.png`
- `results/figures/Order_Diversity_CLR_pPCA_phylomorphospace_plot.png`

**Command:**
```bash
Rscript scripts/processing/phylogenetic_pca_analysis.R
```

---

### 9. Phylogenetic Signal Analysis

**Script:** `analyze_phylogenetic_signal.R`

**Purpose:** Analyze phylogenetic signal in TE composition

**Command:**
```bash
Rscript scripts/processing/analyze_phylogenetic_signal.R
```

---

## Visualization Pipeline

All visualization scripts are in `scripts/visualization/`:

### 1. Divergence Visualization

**Script:** `divergence.R`

**Output:**
- Boxplots: `results/figures/divergence/boxplots/{level}/{te_type}_{metric}_0.9.png`
- Phylogeny bubble plots: `results/figures/divergence_phylogeny/{level}/{te_type}_median_score_phylo_0.9.png`

**Command:**
```bash
Rscript scripts/visualization/divergence.R
```

---

### 2. Diversity Stats on Phylogeny

**Script:** `diversity_stats_phylogeny.R`

**Command:**
```bash
Rscript scripts/visualization/diversity_stats_phylogeny.R
```

---

### 3. TE Diversity Donut Plots

**Script:** `donut_TE_diversity.R` or `hierarchical_donut_TE_diversity.R`

**Output:**
- `results/figures/TE_diversity_hierarchical_donuts/{species_id}_hierarchical_donut.png`

**Command:**
```bash
Rscript scripts/visualization/hierarchical_donut_TE_diversity.R
```

---

### 4. TE Landscape Plots

**Script:** `plot_te_landscape_analysis.R` or `plot_all_te_landscapes.R`

**Output:**
- Various heatmaps and barplots

**Command:**
```bash
Rscript scripts/visualization/plot_te_landscape_analysis.R
```

---

### 5. Phylogenetic Signal Plots

**Script:** `plot_phylogenetic_signal.R`

**Command:**
```bash
Rscript scripts/visualization/plot_phylogenetic_signal.R
```

---

### 6. Simple Phylogeny Plot

**Script:** `plot_simple_phylogeny.R`

**Command:**
```bash
Rscript scripts/visualization/plot_simple_phylogeny.R
```

---

### 7. Ectopic Recombination Plots

**Script:** `plot_ectopic_recombination.R`

**Command:**
```bash
Rscript scripts/visualization/plot_ectopic_recombination.R
```

---

## Complete Workflow (Recommended Order)

```bash
# 1. Activate environment
source venv/bin/activate  # or: conda activate Dusky

# 2. Process raw data
python3 scripts/processing/dnaPipe.py
python3 scripts/processing/repeatmask.py
python3 scripts/processing/ec.py

# 3. Process phylogeny
Rscript scripts/processing/clean_tree_phylo.R

# 4. Calculate metrics
python3 scripts/processing/divergence.py
python3 scripts/processing/diversity_stats.py

# 5. Run analyses
Rscript scripts/processing/pca.R
Rscript scripts/processing/phylogenetic_pca_analysis.R
Rscript scripts/processing/analyze_phylogenetic_signal.R

# 6. Generate visualizations
Rscript scripts/visualization/divergence.R
Rscript scripts/visualization/diversity_stats_phylogeny.R
Rscript scripts/visualization/hierarchical_donut_TE_diversity.R
Rscript scripts/visualization/plot_te_landscape_analysis.R
Rscript scripts/visualization/plot_phylogenetic_signal.R
Rscript scripts/visualization/plot_ectopic_recombination.R
```

---

## Troubleshooting

### Python Package Issues

If you get `ModuleNotFoundError`:
```bash
# Ensure virtual environment is activated
source venv/bin/activate

# Or install missing package
pip install <package_name>
```

### R Package Issues

If you get package loading errors:
```R
# In R console
install.packages("<package_name>")
```

### Path Issues

All scripts use absolute paths based on `/home/jake/Projects/Desmognathus_TE`. If you've moved the project, you'll need to update paths in:
- `scripts/processing/dnaPipe.py` (line 25)
- `scripts/processing/repeatmask.py` (line 21)
- `scripts/processing/pca.R` (lines 13-14)
- Other scripts as needed

### Missing Input Data

Ensure you have data in:
- `input_data/dnaPipeTE/`
- `input_data/repeatmasker/` (.align files)
- `input_data/phylogeny/desmo900dated_test.tre`
- `input_data/lookup_table.txt`

### Memory Issues

Some scripts (dnaPipe.py, repeatmask.py) use multiprocessing. If you encounter memory issues:
- Close other applications
- Check script for `Pool()` initialization and reduce worker count
- Process files in smaller batches

---

## Quick Reference Card

| Task | Command |
|------|---------|
| Activate environment | `source venv/bin/activate` |
| Process dnaPipeTE | `python3 scripts/processing/dnaPipe.py` |
| Process RepeatMasker | `python3 scripts/processing/repeatmask.py` |
| Calculate diversity | `python3 scripts/processing/diversity_stats.py` |
| Run PCA | `Rscript scripts/processing/pca.R` |
| Clean phylogeny | `Rscript scripts/processing/clean_tree_phylo.R` |
| Divergence analysis | `python3 scripts/processing/divergence.py` |
| Make all plots | See "Complete Workflow" section |

---

**Last Updated:** October 19, 2025

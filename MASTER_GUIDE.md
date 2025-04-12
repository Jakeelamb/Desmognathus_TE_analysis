# Desmognathus TE Analysis Master Guide

## Quick Start

1. Activate the Dusky environment:
```bash
conda activate Dusky
```

2. Run the complete analysis pipeline:
```bash
./scripts/batch_analyze_te_landscape.sh --all
```

## Directory Structure

```
desmognathus_te/
├── data/
│   ├── raw/
│   │   ├── repeatmasker/    # RepeatMasker .align files
│   │   ├── lookup/          # Species metadata
│   │   │   └── lookup_table.txt  # Species-SRX-Genome mapping
│   │   └── Phylogeny/       # Phylogenetic tree files
│   │       └── desmo900dated_test.tre  # Raw phylogenetic tree (2.5-30.6 MYA)
│   ├── interim/             # Intermediate files
│   └── processed/           # Processed data
│       └── desmognathus_processed.tre  # Processed phylogenetic tree
├── results/
│   ├── figures/             # All visualizations
│   │   ├── processed_phylogeny_rectangular.png  # Time-calibrated tree
│   │   ├── processed_phylogeny_circular.png     # Circular tree layout
│   │   └── processed_phylogeny_clades.png       # Tree with major clades
│   ├── landscapes/          # TE landscape data
│   └── tables/              # Analysis tables
└── scripts/                 # Analysis scripts
```

## Analysis Workflows

### 1. TE Landscape Analysis

**Purpose**: Analyze TE landscapes for individual samples

**Input Files**:
- `data/raw/repeatmasker/{SRX_ID}_Trinity.align`
- `data/interim/{SRX_ID}_reads_per_component_and_annotation_processed`
- `data/raw/lookup/lookup_table.txt`  # Species mapping

**Output Files**:
- `results/landscapes/repeat_landscape_{SRX_ID}.csv`
- `results/figures/landscape/overall_landscape_{SRX_ID}.png`
- `results/figures/landscape/class_landscape_{SRX_ID}_{Class}.png`
- `results/figures/landscape/order_landscape_{SRX_ID}_{Order}.png`

**Commands**:
```bash
# Single sample
./scripts/analyze_te_landscape.sh SRX19953421

# Multiple samples
./scripts/batch_analyze_te_landscape.sh --samples SRX19952657,SRX19952891

# All samples
./scripts/batch_analyze_te_landscape.sh --all
```

### 2. Combined Visualization

**Purpose**: Generate combined visualizations across all samples

**Input Files**:
- All files in `results/landscapes/`
- `data/raw/lookup/lookup_table.txt`  # Species mapping

**Output Files**:
- `results/figures/combined_te_landscape_all_samples.png`
- `results/figures/heatmap_te_landscape_all_samples.png`
- `results/figures/te_order_distribution.png`

**Command**:
```bash
Rscript scripts/R/visualization/visualize_all_landscapes.R
```

### 3. Phylogenetic Analysis

**Purpose**: Process and visualize phylogenetic relationships

**Input Files**:
- `data/raw/Phylogeny/desmo900dated_test.tre`  # Raw phylogenetic tree (2.5-30.6 MYA)
- `data/raw/lookup/lookup_table.txt`  # Species mapping

**Output Files**:
- `data/processed/desmognathus_processed.tre`  # Processed tree file
- `results/figures/processed_phylogeny_rectangular.png`  # Time-calibrated tree
- `results/figures/processed_phylogeny_circular.png`     # Circular tree layout
- `results/figures/processed_phylogeny_clades.png`       # Tree with major clades

**Major Clades**:
- Quadramaculatus Group
- Fuscus Group
- Ochrophaeus Group

**Commands**:
```bash
# Process phylogenetic tree and generate visualizations
Rscript scripts/R/analysis/process_phylogeny.R
```

### 4. TE Superfamily Analysis

**Purpose**: Analyze TE superfamily distributions

**Input Files**:
- Superfamily breakdown CSV
- Processed landscape data
- `data/raw/lookup/lookup_table.txt`  # Species mapping

**Output Files**:
- `results/tables/superfamily_proportions.csv`
- `results/figures/superfamily_distribution.png`

**Commands**:
```bash
# Generate proportions
python scripts/python/preprocessing/generate_superfamily_proportions.py

# Run diversity analysis
python scripts/python/analysis/run_diversity_analysis.py

# Run PCA analysis
Rscript scripts/R/analysis/te_pca_analysis.R
```

## Dependencies

### Python Packages
- pandas
- numpy
- matplotlib
- seaborn
- PyYAML
- BioPython

### R Packages
- tidyverse
- ggplot2
- ape
- ggtree
- viridis
- phytools

## Troubleshooting

1. **Missing Files**: Ensure all input files are in the correct locations
2. **Environment Issues**: Always activate the Dusky environment
3. **Permission Errors**: Make scripts executable with `chmod +x script.sh`
4. **Package Errors**: Install missing packages using conda or R's install.packages()

## Notes

- All scripts must be run from the project root directory
- The Dusky environment must be active for all analyses
- Check `SCRIPTS_DOCUMENTATION.md` for detailed script information
- Check `README.md` for project overview and setup instructions 
# Desmognathus TE Analysis Scripts Documentation

This document provides a comprehensive overview of all scripts in the Desmognathus TE Analysis project, including their purposes, inputs, outputs, and how they fit into the overall workflow.

## Table of Contents

- [Directory Structure](#directory-structure)
- [Workflow Overview](#workflow-overview)
- [Core Scripts](#core-scripts)
  - [TE Landscape Analysis Scripts](#te-landscape-analysis-scripts)
  - [Batch Processing Scripts](#batch-processing-scripts)
  - [Visualization Scripts](#visualization-scripts)
  - [Phylogenetic Analysis Scripts](#phylogenetic-analysis-scripts)
- [Utility Scripts](#utility-scripts)
- [Dependencies](#dependencies)
- [Script Reference Table](#script-reference-table)

## Directory Structure

```
desmognathus_te/
├── config/                 # Configuration files
│   └── paths.yaml          # Path configuration for data, results, and scripts
├── data/                   # Data directory
│   ├── raw/                # Raw data files
│   │   ├── repeatmasker/   # RepeatMasker .align files
│   │   ├── lookup/         # Species metadata lookup tables
│   │   └── Phylogeny/      # Phylogenetic tree files
│   ├── interim/            # Intermediate data files
│   │   └── pivot_tables/   # Pivot tables for analyses
│   └── processed/          # Processed data ready for analysis
│       ├── diversity/      # Diversity analysis data
│       └── te_superfamily/ # TE superfamily data
├── results/                # Results and outputs
│   ├── figures/            # Generated figures and visualizations
│   │   ├── diversity/      # Diversity analysis figures
│   │   ├── landscape/      # TE landscape figures
│   │   ├── phylo_landscape/# Phylogenetic landscape figures
│   │   └── pca/            # PCA analysis figures
│   ├── landscapes/         # Processed landscape data files
│   └── tables/             # Generated data tables
├── scripts/                # Analysis scripts
│   ├── python/             # Python scripts 
│   │   ├── preprocessing/  # Data preprocessing scripts
│   │   ├── analysis/       # Analysis scripts
│   │   ├── utils/          # Utility functions and modules
│   │   └── tests/          # Test scripts
│   └── R/                  # R scripts
│       ├── analysis/       # Analysis scripts
│       └── visualization/  # Visualization scripts
```

## Workflow Overview

The Desmognathus TE Analysis workflow consists of several key steps:

1. **Data Preprocessing**
   - Parsing RepeatMasker `.align` files
   - Joining with classification metadata
   - Calculating Kimura divergence values

2. **TE Landscape Analysis**
   - Generating repeat landscapes for individual samples
   - Creating visualization plots
   - Aggregating results across samples

3. **Comparative Analysis**
   - Comparing TE landscapes across species
   - Mapping TE data to phylogenetic trees
   - Generating summary statistics

4. **Visualization**
   - Creating individual sample plots
   - Generating combined visualizations
   - Producing phylogenetic visualizations with TE data

## Core Scripts

### TE Landscape Analysis Scripts

#### `scripts/analyze_te_landscape.sh`

**Purpose**: Analyzes TE landscape data for a single SRX sample.

**Inputs**:
- SRX ID (command-line argument)
- RepeatMasker `.align` file: `input_data/repeatmasker/{SRX_ID}_Trinity.align`
- Canonical classification table: `results/data/dnaPipeTE_merged_classifications.csv`

**Outputs**:
- CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`
- Plots: Various visualization plots in `results/figures/landscape/`

**Usage**: `./scripts/analyze_te_landscape.sh SRX19953421`

**Dependencies**:
- Dusky conda environment
- Python script: `scripts/processing/parse_repeatmasker_landscape.py`
- R script: `scripts/R/visualization/plot_te_landscape.R`

#### `scripts/python/preprocessing/parse_repeatmasker_landscape.py`

**Purpose**: Backward-compatible wrapper for the canonical RepeatMasker landscape parser.

**Inputs**:
- SRX ID (command-line argument)
- RepeatMasker `.align` file: `input_data/repeatmasker/{SRX_ID}_Trinity.align`
- Canonical classification table: `results/data/dnaPipeTE_merged_classifications.csv`

**Outputs**:
- CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`

**Usage**: `python scripts/python/preprocessing/parse_repeatmasker_landscape.py SRX19953421`

**Note**: The canonical implementation now lives at `scripts/processing/parse_repeatmasker_landscape.py`. The preprocessing path remains as a compatibility wrapper so older commands still work.

### Batch Processing Scripts

#### `scripts/batch_analyze_te_landscape.sh`

**Purpose**: Batch analyzes TE landscape data for multiple SRX samples.

**Options**:
- `--all`: Process all available samples
- `--sample-list FILE`: Process samples listed in FILE (one SRX ID per line)
- `--samples SRX1,SRX2`: Process specific samples (comma-separated)
- `--help`: Show help message

**Inputs**:
- List of SRX IDs from command-line options
- For each SRX ID:
  - RepeatMasker `.align` file: `input_data/repeatmasker/{SRX_ID}_Trinity.align`
  - Canonical classification table: `results/data/dnaPipeTE_merged_classifications.csv`
  - Species lookup table: `input_data/lookup_table.txt`

**Outputs**:
- For each SRX ID:
  - CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`
  - Individual plots in `results/figures/landscape/`
- Combined visualization plots:
  - `results/figures/landscape/te_order_distribution.png`
- Summary file: `results/landscapes/landscape_analysis_summary.csv`

**Usage**: `./scripts/batch_analyze_te_landscape.sh --samples SRX19952657,SRX19952891`

**Dependencies**:
- Dusky conda environment
- `scripts/analyze_te_landscape.sh`
- GNU Parallel when available; otherwise the script falls back to serial execution
- `scripts/visualization/plot_all_te_landscapes.R`

### Visualization Scripts

#### `scripts/R/visualization/plot_te_landscape.R`

**Purpose**: Creates visualizations for a single SRX sample's TE landscape data.

**Inputs**:
- SRX ID (command-line argument)
- Processed landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`

**Outputs**:
- Overall landscape plot: `results/figures/landscape/overall_landscape_{SRX_ID}.png`
- Class-specific landscape plots: `results/figures/landscape/class_landscape_{SRX_ID}_{Class}.png`
- Order-specific landscape plots: `results/figures/landscape/order_landscape_{SRX_ID}_{Order}.png`

**Usage**: `Rscript scripts/R/visualization/plot_te_landscape.R SRX19953421`

#### `scripts/R/visualization/visualize_all_landscapes.R`

**Purpose**: Creates combined visualizations for all processed landscape data.

**Inputs**:
- All processed landscape CSV files in `results/landscapes/`
- Species lookup table: `data/raw/lookup/lookup_table.txt`

**Outputs**:
- Combined line plot: `results/figures/combined_te_landscape_all_samples.png`
- Heatmap visualization: `results/figures/heatmap_te_landscape_all_samples.png`
- TE order distribution: `results/figures/te_order_distribution.png`

**Usage**: `Rscript scripts/R/visualization/visualize_all_landscapes.R`

**Key Features**:
- Filters out "Unknown" species
- Uses improved data visualization techniques
- Groups data by species for comparative analysis

### Phylogenetic Analysis Scripts

#### `scripts/R/visualization/te_phylo_landscape.R`

**Purpose**: Creates visualizations that integrate TE landscape data with a phylogenetic tree of Desmognathus species.

**Inputs**:
- Phylogenetic tree file: `input_data/phylogeny/desmo900dated_test.tre`
- Species lookup table: `input_data/lookup_table.txt`
- Processed landscape files: `results/landscapes/repeat_landscape_*.csv`

**Outputs**:
- Basic tree visualization: `results/figures/phylo_landscape/basic_phylogeny.png`
- Phylogeny with TE landscape heatmap: `results/figures/phylo_landscape/phylogeny_with_landscape.png`
- Phylogeny with TE class distribution: `results/figures/phylo_landscape/phylogeny_with_classes.png`
- Runtime note: no figures are produced unless `results/landscapes/repeat_landscape_*.csv` is populated locally

**Usage**: `Rscript scripts/R/visualization/te_phylo_landscape.R`

**Key Features**:
- Processes tree tips to remove suffixes
- Prunes species not found in the lookup table
- Maps SRX values to tree tip labels
- Integrates TE landscape data with phylogenetic relationships

## Utility Scripts

#### `scripts/python/utils/create_gitkeep.py`

**Purpose**: Creates `.gitkeep` files in empty directories to maintain directory structure in git.

**Usage**: `python scripts/python/utils/create_gitkeep.py`

#### `scripts/run_tests.sh`

**Purpose**: Runs all tests for the project.

**Usage**: `./scripts/run_tests.sh`

#### `scripts/utilities/split_fasta.py`

**Purpose**: Splits a FASTA file into smaller chunks.

**Usage**: `python scripts/utilities/split_fasta.py input.fasta chunk_size output_prefix`

## Dependencies

All scripts require the "Dusky" conda environment to be activated before running:

```bash
conda activate Dusky
```

### Python Dependencies
- pandas
- numpy
- matplotlib
- seaborn
- PyYAML
- BioPython

### R Dependencies
- tidyverse
- ggplot2
- ape (for phylogenetic trees)
- ggtree (for phylogenetic visualization)
- viridis (for color palettes)
- phytools (for tree manipulation) 

## Script Reference Table

| Script Name | Language | Purpose | Main Inputs | Main Outputs |
|-------------|----------|---------|-------------|--------------|
| **TE Landscape Analysis** |
| `analyze_te_landscape.sh` | Bash | Process single SRX sample | SRX ID, .align file, classification file | CSV landscape data, visualization plots |
| `parse_repeatmasker_landscape.py` | Python | Parse RepeatMasker files | .align file, classification file | CSV landscape data |
| `batch_analyze_te_landscape.sh` | Bash | Process multiple SRX samples | List of SRX IDs | Multiple CSV files, combined plots, summary |
| **Visualization** |
| `plot_te_landscape.R` | R | Single sample visualizations | SRX ID, landscape CSV | Various PNG plots |
| `visualize_all_landscapes.R` | R | Combined visualizations | All landscape CSVs, lookup table | Combined plots, heatmaps |
| `te_phylo_landscape.R` | R | Phylogenetic visualizations | Tree file, landscape CSVs, lookup table | Tree plots with TE data |
| `plot_te_landscape_analysis.R` | R | TE landscape summary plots | `results/data/dnaPipeTE_*_breakdown.csv` | Summary visualizations |
| **Data Processing** |
| `process_te_landscape.py` | Python | Complete TE analysis workflow | Raw data files | Processed data, plots |
| `generate_superfamily_proportions.py` | Python | Calculate TE proportions | TE breakdown CSV | Proportion tables |
| `run_diversity_analysis.py` | Python | Diversity metrics calculation | Superfamily proportions | Diversity metrics, plots |
| **Phylogenetic Analysis** |
| `pca.R` | R | Canonical compositional PCA on frozen TE breakdown tables | `results/data/dnaPipeTE_*_breakdown.csv` | PCA manifests, CLR matrices, loadings, scores, plots |
| `phylogenetic_pca_analysis.R` | R | Supplementary phylogenetic PCA from saved CLR matrices | `results/tables/pca/*_clr_matrix.csv`, phylogeny | pPCA manifests, scores, loadings, phylomorphospace plots |
| `analyze_phylogenetic_correlogram.R` | R | Moran's I correlograms across phylogenetic distance bins | `results/data/dnaPipeTE_superfamily_breakdown.csv`, phylogeny | `results/tables/phylogenetic_signal/phylogenetic_correlogram_moran_per_trait.csv`, correlogram plot |
| **Utility Scripts** |
| `create_gitkeep.py` | Python | Maintain git directory structure | None | .gitkeep files |
| `split_fasta.py` | Python | Split FASTA files | FASTA file, chunk size | Multiple FASTA chunks |
| `run_tests.sh` | Bash | Run project tests | Test files | Test reports | 

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
- RepeatMasker `.align` file: `data/raw/repeatmasker/{SRX_ID}_Trinity.align`
- Classification file: `data/interim/{SRX_ID}_reads_per_component_and_annotation_processed`

**Outputs**:
- CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`
- Plots: Various visualization plots in `results/figures/landscape/`

**Usage**: `./scripts/analyze_te_landscape.sh SRX19953421`

**Dependencies**:
- Dusky conda environment
- Python script: `scripts/python/preprocessing/parse_repeatmasker_landscape.py`
- R script: `scripts/R/visualization/plot_te_landscape.R`

#### `scripts/python/preprocessing/parse_repeatmasker_landscape.py`

**Purpose**: Parses RepeatMasker `.align` files and joins with classification metadata to generate a dataset for TE repeat landscape analysis.

**Inputs**:
- SRX ID (command-line argument)
- RepeatMasker `.align` file: `data/raw/repeatmasker/{SRX_ID}_Trinity.align`
- Classification file: `data/interim/{SRX_ID}_reads_per_component_and_annotation_processed`

**Outputs**:
- CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`

**Usage**: `python scripts/python/preprocessing/parse_repeatmasker_landscape.py SRX19953421`

**Key Functions**:
- `parse_align_file()`: Extracts contig names, aligned base pairs, and Kimura distances
- `load_classification_file()`: Loads TE classification metadata
- `merge_align_classification()`: Merges alignment data with classification info
- `aggregate_landscape()`: Aggregates data into Kimura distance bins

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
  - RepeatMasker `.align` file: `data/raw/repeatmasker/{SRX_ID}_Trinity.align`
  - Classification file: `data/interim/{SRX_ID}_reads_per_component_and_annotation_processed`

**Outputs**:
- For each SRX ID:
  - CSV landscape data: `results/landscapes/repeat_landscape_{SRX_ID}.csv`
  - Individual plots in `results/figures/landscape/`
- Combined visualization plots:
  - `results/figures/combined_te_landscape_all_samples.png`
  - `results/figures/heatmap_te_landscape_all_samples.png`
  - `results/figures/te_order_distribution.png`
- Summary file: `results/landscapes/landscape_analysis_summary.csv`

**Usage**: `./scripts/batch_analyze_te_landscape.sh --samples SRX19952657,SRX19952891`

**Dependencies**:
- Dusky conda environment
- `scripts/analyze_te_landscape.sh`
- `scripts/R/visualization/visualize_all_landscapes.R`

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
- Phylogenetic tree file: `data/raw/Phylogeny/desmo900dated_test.tre`
- Species lookup table: `data/raw/lookup/lookup_table.txt`
- Processed landscape files: `results/landscapes/repeat_landscape_*.csv`

**Outputs**:
- Basic tree visualization: `results/figures/phylo_landscape/desmognathus_phylogeny.png`
- Phylogeny with TE landscape heatmap: `results/figures/phylo_landscape/phylogeny_with_te_landscape_heatmap.png`
- Phylogeny with TE class distribution: `results/figures/phylo_landscape/phylogeny_with_te_class_distribution.png`
- Circular phylogeny: `results/figures/phylo_landscape/desmognathus_circular_phylogeny.png`

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
| `te_landscape_plots.R` | R | TE landscape summary plots | Processed data tables | Summary visualizations |
| **Data Processing** |
| `process_te_landscape.py` | Python | Complete TE analysis workflow | Raw data files | Processed data, plots |
| `generate_superfamily_proportions.py` | Python | Calculate TE proportions | TE breakdown CSV | Proportion tables |
| `run_diversity_analysis.py` | Python | Diversity metrics calculation | Superfamily proportions | Diversity metrics, plots |
| **Phylogenetic Analysis** |
| `te_pca_analysis.R` | R | PCA on TE distributions | Superfamily proportions | PCA plots, tables |
| **Utility Scripts** |
| `create_gitkeep.py` | Python | Maintain git directory structure | None | .gitkeep files |
| `split_fasta.py` | Python | Split FASTA files | FASTA file, chunk size | Multiple FASTA chunks |
| `run_tests.sh` | Bash | Run project tests | Test files | Test reports | 
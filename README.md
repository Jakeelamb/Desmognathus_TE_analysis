# Desmognathus TE Analysis

This repository contains scripts and tools for analyzing transposable elements (TEs) in Desmognathus salamander genomes.

## Important Note

**All scripts must be run within the "Dusky" conda environment.**

Always ensure the Dusky environment is active before running any scripts:

```bash
conda activate Dusky
```

## Project Structure

```
desmognathus_te/
├── config/               # Configuration files
│   └── paths.yaml        # Path configuration for data, results, and scripts
├── data/                 # Data directory
│   ├── raw/              # Raw data files (fastq, assemblies, etc.)
│   ├── interim/          # Intermediate data files
│   │   └── pivot_tables/ # Pivot tables for analyses
│   └── processed/        # Processed data ready for analysis
│       ├── diversity/    # Diversity analysis data
│       └── te_superfamily/ # TE superfamily data
├── results/              # Results and outputs
│   ├── figures/          # Generated figures and visualizations
│   │   ├── diversity/    # Diversity analysis figures
│   │   ├── landscape/    # TE landscape figures
│   │   └── pca/          # PCA analysis figures
│   └── tables/           # Generated data tables
│       ├── diversity/    # Diversity analysis tables
│       └── pca/          # PCA analysis tables
└── scripts/              # Analysis scripts
    ├── python/           # Python scripts
    │   ├── preprocessing/ # Data preprocessing scripts
    │   ├── analysis/     # Analysis scripts
    │   ├── utils/        # Utility functions and modules
    │   └── visualization/ # Visualization scripts
    └── R/                # R scripts
        ├── analysis/     # Analysis scripts
        └── visualization/ # Visualization scripts
```

## Workflows

The project supports several analysis workflows:

### TE Landscape Analysis

The TE landscape analysis workflow includes the following steps:

1. **Generate superfamily proportions**: Calculate the proportions of each TE superfamily.
2. **Diversity analysis**: Calculate diversity metrics for TE superfamilies.
3. **PCA analysis**: Perform principal component analysis on TE superfamily proportions.

Run the complete workflow with:

```bash
conda activate Dusky
python scripts/process_te_landscape.py
```

#### Workflow Options

- `--skip-proportions`: Skip generating superfamily proportions
- `--skip-diversity`: Skip diversity analysis
- `--skip-pca`: Skip PCA analysis
- `--min-species-presence INT`: Minimum number of species a TE must be present in (default: 3)
- `--verbose`: Enable verbose logging

## Individual Scripts

### Generate Superfamily Proportions

Generates the proportions of each TE superfamily from a breakdown CSV file:

```bash
conda activate Dusky
python scripts/python/preprocessing/generate_superfamily_proportions.py [OPTIONS]
```

Options:
- `--input PATH`: Path to superfamily breakdown CSV
- `--output PATH`: Path to save the output file
- `--diversity-copy`: Copy output to diversity directory for PCA analysis
- `--verbose`: Enable verbose logging

### Run Diversity Analysis

Calculates diversity metrics (Shannon, Simpson, etc.) for TE superfamily proportions:

```bash
conda activate Dusky
python scripts/python/analysis/run_diversity_analysis.py [OPTIONS]
```

Options:
- `--input PATH`: Path to superfamily proportions CSV
- `--output-dir PATH`: Directory to save output files
- `--no-plots`: Skip generating plots
- `--verbose`: Enable verbose logging

### Run PCA Analysis

Performs principal component analysis on TE superfamily proportions:

```bash
conda activate Dusky
Rscript scripts/R/analysis/te_pca_analysis.R [OPTIONS]
```

Options:
- `--input PATH`: Path to superfamily proportions CSV
- `--output PATH`: Directory to save output files
- `--min-species-presence INT`: Minimum number of species a TE must be present in (default: 3)
- `--verbose`: Enable verbose logging

### Generate TE Landscape Visualizations

Creates visualizations for TE landscape analysis:

```bash
conda activate Dusky
Rscript scripts/R/visualization/te_landscape_plots.R [OPTIONS]
```

Options:
- `--verbose`: Enable verbose logging

## Project Cleanup

A cleanup script is provided to remove legacy directories after migrating to the new structure:

```bash
conda activate Dusky
./cleanup_legacy_dirs.sh
```

This script will:
1. Remove redundant directories within the new structure
2. Optionally remove legacy directories (Data, Projects, Results, Output, old_scripts)

**Note:** Make sure you have a backup before running this script, as it permanently deletes files.

## Dependencies

### Python Dependencies

- pandas
- numpy
- matplotlib
- seaborn
- PyYAML

### R Dependencies

- tidyverse
- FactoMineR
- factoextra
- yaml

## Installation

Clone the repository and install dependencies:

```bash
# Clone the repository
git clone https://github.com/yourusername/desmognathus_te.git
cd desmognathus_te

# Create and activate the Dusky conda environment
conda create -n Dusky python=3.12
conda activate Dusky

# Install Python dependencies
conda install -y pandas numpy matplotlib seaborn pyyaml

# Install R and required packages
conda install -y -c conda-forge r-base r-tidyverse r-factoextra r-factominer r-yaml
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Git Management

This repository is configured to track only essential code files while ignoring large data files:

- All scripts and configuration files are tracked
- Raw data, processed outputs, and results files are ignored
- Directory structure is maintained using `.gitkeep` files

To generate `.gitkeep` files for proper directory tracking:

```bash
conda activate Dusky
python scripts/python/utils/create_gitkeep.py
```

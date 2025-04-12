# Desmognathus TE Project - New Directory Structure

This document outlines the new directory structure for the Desmognathus TE project and provides a guide for transitioning from the old structure to the new one.

## New Directory Structure

```
Desmognathus_TE/
├── config/                          # Configuration files
│   └── paths.yaml                   # Central path configuration
├── data/                            # All data files (lowercase)
│   ├── raw/                         # Original, unchangeable data
│   │   ├── dnaPipeTE_outputs/       # Raw outputs from dnaPipeTE
│   │   ├── phylogenetic_tree/       # Original phylogenetic tree files
│   │   └── lookup_table.txt         # Species mapping table
│   ├── processed/                   # Processed, analysis-ready data
│   │   ├── te_classifications/      # Classified TE files
│   │   ├── te_superfamily/          # Superfamily data
│   │   └── diversity/               # Diversity indices
│   ├── interim/                     # Intermediate data
│   │   └── pivot_tables/            # Pivot tables for visualization
│   └── external/                    # External data sources
├── results/                         # Analysis results (lowercase)
│   ├── figures/                     # Visualizations and plots
│   │   ├── te_landscape/            # TE landscape plots
│   │   ├── pca/                     # PCA analysis plots
│   │   ├── phylogeny/               # Phylogenetic tree visualizations
│   │   └── ectopic_recombination/   # Ectopic recombination plots
│   ├── tables/                      # Tabular results
│   ├── models/                      # Model outputs
│   └── reports/                     # Analysis reports
├── scripts/                         # All scripts
│   ├── python/                      # Python scripts
│   │   ├── preprocessing/           # Data preparation
│   │   ├── analysis/                # Analysis scripts
│   │   ├── visualization/           # Visualization scripts
│   │   └── utils/                   # Utility functions
│   │       └── paths.py             # Path configuration utility
│   ├── R/                           # R scripts
│   │   ├── analysis/                # Analysis scripts
│   │   └── visualization/           # Visualization scripts
│   └── bash/                        # Shell scripts
│       ├── pipeline/                # Full pipelines
│       └── utils/                   # Utility shell scripts
└── notebooks/                       # Jupyter notebooks
```

## Migration Plan

The migration to the new directory structure will be done incrementally to avoid breaking existing workflows:

### Phase 1: Setup (Complete)
- [x] Create the new directory structure
- [x] Create the paths.yaml configuration file
- [x] Create the Python path utility module

### Phase 2: Script Migration (In Progress)
- [ ] Move Python scripts to scripts/python/
- [ ] Move R scripts to scripts/R/
- [ ] Move shell scripts to scripts/bash/
- [ ] Update path references in scripts to use the new configuration

### Phase 3: Data Migration
- [ ] Copy raw data to data/raw/
- [ ] Copy processed data to data/processed/
- [ ] Update workflows to use the new data paths

### Phase 4: Results Migration
- [ ] Move figures to results/figures/
- [ ] Move tables to results/tables/
- [ ] Move reports to results/reports/

### Phase 5: Cleanup
- [ ] Remove obsolete files and directories
- [ ] Update documentation
- [ ] Final testing of all workflows

## Using the New Structure

### Path Configuration

All paths are now defined in `config/paths.yaml`. This file serves as a central location for all path references. The Python utility module `scripts/python/utils/paths.py` provides functions to access these paths programmatically:

```python
from scripts.python.utils.paths import get_path, get_and_ensure_dir

# Get a path from the configuration
raw_data_path = get_path("data.raw")

# Get a path and ensure the directory exists
output_dir = get_and_ensure_dir("results.figures.pca")
```

## Notes on Backward Compatibility

During the transition, both the old and new directory structures will be maintained. Legacy paths are defined in the `legacy` section of the paths.yaml configuration file.

## Tasks for Completion

- [ ] Update all scripts to use the new path configuration
- [ ] Create symbolic links for backward compatibility
- [ ] Document each script's purpose and dependencies
- [ ] Create a requirements.txt file for Python dependencies
- [ ] Create a central pipeline script that calls individual components 
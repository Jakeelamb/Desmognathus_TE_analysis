#!/bin/bash

# Run full pipeline for Desmognathus TE analysis
set -e  # Exit immediately if a command exits with a non-zero status

# Activate conda environment
echo "Activating Dusky conda environment..."
conda activate Dusky || { echo "Failed to activate Dusky environment"; exit 1; }

# Step 1: Generate individual landscape files
echo "Step 1: Generating individual landscape files..."
./scripts/batch_analyze_te_landscape.sh --all --cpus 13

# Step 2: Generate combined visualizations
echo "Step 2: Generating combined visualizations..."
Rscript scripts/R/visualization/visualize_all_landscapes.R

# Step 3: Generate phylogenetic visualizations
echo "Step 3: Generating phylogenetic visualizations..."
Rscript scripts/R/visualization/te_phylo_landscape.R

echo "Pipeline completed successfully!"
echo "Results available in:"
echo "- Individual landscapes: results/landscapes/"
echo "- Combined visualizations: results/figures/landscape/"
echo "- Phylogenetic visualizations: results/figures/phylo_landscape/" 
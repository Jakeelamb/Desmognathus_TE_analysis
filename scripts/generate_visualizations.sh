#!/bin/bash
#
# generate_visualizations.sh
#
# Description: Generates all visualizations for the Desmognathus TE project
#

set -e  # Exit on error

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate Dusky || { echo "Failed to activate Dusky environment"; exit 1; }

# Create necessary directories
mkdir -p results/figures/phylo_landscape
mkdir -p results/figures/landscape

echo "Step 1: Generating basic phylogeny visualization..."
Rscript scripts/R/visualization/simple_phylogeny.R

echo "Step 2: Generating TE order distribution..."
Rscript scripts/R/visualization/visualize_all_landscapes.R

echo "Step 3: Generating phylogeny with TE landscape and class visualizations..."
# Create a modified copy of the original script that skips the circular plot
cp scripts/R/visualization/te_phylo_landscape.R scripts/R/visualization/te_phylo_landscape_modified.R
sed -i 's/# Save circular phylogeny/# Circular phylogeny is skipped\n  # ggsave(/' scripts/R/visualization/te_phylo_landscape_modified.R

# Run the modified script
Rscript scripts/R/visualization/te_phylo_landscape_modified.R

# Clean up
rm scripts/R/visualization/te_phylo_landscape_modified.R

echo "All visualizations completed successfully."
echo "Results available in:"
echo "- Basic phylogeny: results/figures/phylo_landscape/basic_phylogeny.png"
echo "- TE landscape distribution: results/figures/landscape/te_order_distribution.png"
echo "- Phylogeny with TE landscape: results/figures/phylo_landscape/phylogeny_with_landscape.png"
echo "- Phylogeny with TE classes: results/figures/phylo_landscape/phylogeny_with_classes.png" 
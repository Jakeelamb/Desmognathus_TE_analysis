#!/bin/bash
#
# generate_visualizations.sh
#
# Description: Generates all visualizations for the Desmognathus TE project
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Activate conda environment
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda command not found."
    exit 1
fi

if [ "${CONDA_DEFAULT_ENV:-}" != "Dusky" ]; then
    eval "$(conda shell.bash hook)"
    conda activate Dusky || { echo "Failed to activate Dusky environment"; exit 1; }
fi

# Create necessary directories
mkdir -p "${PROJECT_ROOT}/results/figures/phylo_landscape"
mkdir -p "${PROJECT_ROOT}/results/figures/landscape"

echo "Step 1: Generating basic phylogeny visualization..."
Rscript "${PROJECT_ROOT}/scripts/R/visualization/simple_phylogeny.R"

echo "Step 2: Generating TE order distribution..."
Rscript "${PROJECT_ROOT}/scripts/visualization/plot_all_te_landscapes.R"

echo "Step 3: Generating phylogeny with TE landscape and class visualizations..."
Rscript "${PROJECT_ROOT}/scripts/R/visualization/te_phylo_landscape.R"

echo "All visualizations completed successfully."
echo "Results available in:"
echo "- Basic phylogeny: ${PROJECT_ROOT}/results/figures/phylogeny/rectangular_phylogeny.png"
echo "- TE landscape distribution: ${PROJECT_ROOT}/results/figures/landscape/te_order_distribution.png"
echo "- Phylogeny with TE landscape: ${PROJECT_ROOT}/results/figures/phylo_landscape/phylogeny_with_landscape.png"
echo "- Phylogeny with TE classes: ${PROJECT_ROOT}/results/figures/phylo_landscape/phylogeny_with_classes.png"

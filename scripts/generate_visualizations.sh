#!/bin/bash
#
# generate_visualizations.sh
#
# Description: Generates all visualizations for the Desmognathus TE project
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RUN_IN_DUSKY="${PROJECT_ROOT}/scripts/run_in_dusky.sh"

# Create necessary directories
mkdir -p "${PROJECT_ROOT}/results/figures/phylo_landscape"
mkdir -p "${PROJECT_ROOT}/results/figures/landscape"

echo "Step 1: Generating basic phylogeny visualization..."
"${RUN_IN_DUSKY}" Rscript "${PROJECT_ROOT}/scripts/visualization/plot_simple_phylogeny.R"

echo "Step 2: Generating TE order distribution..."
"${RUN_IN_DUSKY}" Rscript "${PROJECT_ROOT}/scripts/visualization/plot_all_te_landscapes.R"

echo "Step 3: Generating phylogeny with TE landscape and class visualizations..."
"${RUN_IN_DUSKY}" Rscript "${PROJECT_ROOT}/scripts/visualization/plot_phylogeny_with_te_landscape.R"

echo "All visualizations completed successfully."
echo "Results available in:"
echo "- Basic phylogeny: ${PROJECT_ROOT}/results/figures/phylogeny/rectangular_phylogeny.png"
echo "- TE landscape distribution: ${PROJECT_ROOT}/results/figures/landscape/te_order_distribution.png"
echo "- Phylogeny with TE landscape: ${PROJECT_ROOT}/results/figures/phylo_landscape/phylogeny_with_landscape.png"
echo "- Phylogeny with TE classes: ${PROJECT_ROOT}/results/figures/phylo_landscape/phylogeny_with_classes.png"

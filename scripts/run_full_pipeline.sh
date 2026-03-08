#!/bin/bash

# Run the current TE landscape pipeline from canonical repo paths.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda command not found."
    exit 1
fi

if [ "${CONDA_DEFAULT_ENV:-}" != "Dusky" ]; then
    echo "Activating Dusky conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate Dusky || { echo "Failed to activate Dusky environment"; exit 1; }
fi

if [ "$#" -eq 0 ]; then
    BATCH_ARGS=(--all --cpus 4)
else
    BATCH_ARGS=("$@")
fi

echo "Step 1: Generating landscape CSVs and per-sample plots..."
"${PROJECT_ROOT}/scripts/batch_analyze_te_landscape.sh" "${BATCH_ARGS[@]}"

echo "Step 2: Regenerating aggregate landscape visualizations..."
Rscript "${PROJECT_ROOT}/scripts/visualization/plot_all_te_landscapes.R"

echo "Step 3: Regenerating phylogeny-linked landscape visualizations..."
Rscript "${PROJECT_ROOT}/scripts/R/visualization/te_phylo_landscape.R"

echo "Pipeline completed successfully."
echo "Results available in:"
echo "- Individual landscapes: ${PROJECT_ROOT}/results/landscapes/"
echo "- Combined visualizations: ${PROJECT_ROOT}/results/figures/landscape/"
echo "- Phylogenetic visualizations: ${PROJECT_ROOT}/results/figures/phylo_landscape/"

#!/bin/bash

# Run the current TE landscape pipeline from canonical repo paths.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
RUN_IN_DUSKY="${PROJECT_ROOT}/scripts/run_in_dusky.sh"

if [ "$#" -eq 0 ]; then
    BATCH_ARGS=(--all --cpus 4)
else
    BATCH_ARGS=("$@")
fi

echo "Step 1: Generating landscape CSVs and per-sample plots..."
"${RUN_IN_DUSKY}" bash "${PROJECT_ROOT}/scripts/batch_analyze_te_landscape.sh" "${BATCH_ARGS[@]}"

echo "Step 2: Regenerating aggregate landscape visualizations..."
"${RUN_IN_DUSKY}" Rscript "${PROJECT_ROOT}/scripts/visualization/plot_all_te_landscapes.R"

echo "Step 3: Regenerating phylogeny-linked landscape visualizations..."
"${RUN_IN_DUSKY}" Rscript "${PROJECT_ROOT}/scripts/visualization/plot_phylogeny_with_te_landscape.R"

echo "Pipeline completed successfully."
echo "Results available in:"
echo "- Individual landscapes: ${PROJECT_ROOT}/results/landscapes/"
echo "- Combined visualizations: ${PROJECT_ROOT}/results/figures/landscape/"
echo "- Phylogenetic visualizations: ${PROJECT_ROOT}/results/figures/phylo_landscape/"

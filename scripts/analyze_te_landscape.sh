#!/bin/bash
#
# analyze_te_landscape.sh
#
# Description: Analyzes TE landscape data for a single SRX sample
#
# Usage: ./analyze_te_landscape.sh SRX_ID
#
# Example: ./analyze_te_landscape.sh SRX19953421
#

# Exit on error
set -e

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda command not found."
    echo "Please install conda and try again."
    exit 1
fi

# Check if the Dusky environment exists
if ! conda env list | grep -q "Dusky"; then
    echo "Error: Dusky conda environment not found."
    echo "Please create the Dusky environment as specified in the README."
    exit 1
fi

# Ensure conda environment is activated
if [ "$CONDA_DEFAULT_ENV" != "Dusky" ]; then
    echo "Activating Dusky conda environment..."
    eval "$(conda shell.bash hook)"
    conda activate Dusky
fi

# Check for command line argument
if [ $# -lt 1 ]; then
    echo "Error: No SRX ID provided"
    echo "Usage: $0 SRX_ID"
    exit 1
fi

SRX_ID=$1
echo "Analyzing TE landscape for sample: $SRX_ID"

# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PARSER_SCRIPT="${SCRIPT_DIR}/processing/parse_repeatmasker_landscape.py"
PLOT_SCRIPT="${SCRIPT_DIR}/visualization/plot_te_landscape.R"

# Create necessary directories
mkdir -p "${PROJECT_ROOT}/results/landscapes"
mkdir -p "${PROJECT_ROOT}/results/figures/landscape"

# Check if input files exist
ALIGN_FILE="${PROJECT_ROOT}/input_data/repeatmasker/${SRX_ID}_Trinity.align"

if [ ! -f "$ALIGN_FILE" ]; then
    echo "Error: Align file not found: $ALIGN_FILE"
    exit 1
fi

echo "Input files:"
echo "- Alignment file: $ALIGN_FILE"
echo "- Classification source: ${PROJECT_ROOT}/results/data/dnaPipeTE_merged_classifications.csv (preferred)"

# Make scripts executable
chmod +x "$PARSER_SCRIPT" "$PLOT_SCRIPT"

# Step 1: Parse the .align file
echo "Step 1: Parsing RepeatMasker .align file..."
python3 "$PARSER_SCRIPT" "$SRX_ID"
if [ $? -ne 0 ]; then
    echo "Error: Parsing failed"
    exit 1
fi

# Step 2: Create visualizations
echo "Step 2: Generating visualizations..."
Rscript "$PLOT_SCRIPT" "$SRX_ID"
if [ $? -ne 0 ]; then
    echo "Error: Visualization failed"
    exit 1
fi

echo "Analysis completed successfully."
echo "Output files:"
echo "- CSV data: ${PROJECT_ROOT}/results/landscapes/repeat_landscape_${SRX_ID}.csv"
echo "- Plots: ${PROJECT_ROOT}/results/figures/landscape/"

# List generated files
echo "Generated plots:"
ls -l "${PROJECT_ROOT}/results/figures/landscape/"*"${SRX_ID}"*

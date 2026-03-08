#!/bin/bash

echo "Warning: scripts/process_ltr.sh is an upstream raw-data utility that still expects the historical Data/ and Output/ layout."
echo "The canonical paper-facing ectopic workflow is scripts/processing/ec.py plus scripts/visualization/plot_ectopic_recombination.R."

set -e  # Exit immediately if a command exits with a non-zero status

# Set temporary directory
export TMPDIR="./tmp"
mkdir -p $TMPDIR

# --- Configuration ---
DATA_DIR="Data"
OUTPUT_DIR="Output"
RESULTS_DIR="Results"
LTR_DIR="${DATA_DIR}/LTRs"
LOG_FILE="${DATA_DIR}/Log_file.txt"

# --- Input Validation ---
if [ $# -lt 3 ]; then
    echo "Usage: $0 <read_file1> <read_file2> <species_name>"
    exit 1
fi

READ_FILE1="$1"
READ_FILE2="$2"
SPECIES="$3"

# Check if input files exist
[[ -f $READ_FILE1 ]] || { echo "Read file 1 not found: $READ_FILE1"; exit 1; }
[[ -f $READ_FILE2 ]] || { echo "Read file 2 not found: $READ_FILE2"; exit 1; }
[[ -d $LTR_DIR ]] || { echo "LTR directory not found: $LTR_DIR"; exit 1; }

# Create needed directories
mkdir -p "${OUTPUT_DIR}/Depth/${SPECIES}" || { echo "Failed to create depth directory"; exit 1; }

# Initialize log file
echo "Processing started at $(date)" > "$LOG_FILE"
echo "Species: $SPECIES" >> "$LOG_FILE"
echo "Read files: $READ_FILE1, $READ_FILE2" >> "$LOG_FILE"

# Function to process a single LTR pair
process_ltr() {
    local LTR_file=$1

    # Create a unique temporary directory for this process
    local PROCESS_TMPDIR="${TMPDIR}/$(basename ${LTR_file})_${RANDOM}"
    mkdir -p $PROCESS_TMPDIR
    
    echo "Processing file: $LTR_file" >> "$LOG_FILE"
    
    # Define output files
    local SAM_FILE="${PROCESS_TMPDIR}/$(basename ${LTR_file}).sam"
    local SORTED_BAM="${LTR_file}.sorted.bam"
    local DEPTH_FILE="${OUTPUT_DIR}/Depth/${SPECIES}/$(basename ${LTR_file}).depth.txt"
    
    # Build indices and map reads
    echo "Building index for $LTR_file" >> "$LOG_FILE"
    bowtie2-build $LTR_file ${PROCESS_TMPDIR}/$(basename ${LTR_file})_index >> "$LOG_FILE" 2>&1
    
    echo "Mapping reads to $LTR_file" >> "$LOG_FILE"
    bowtie2 -1 $READ_FILE1 -2 $READ_FILE2 -x ${PROCESS_TMPDIR}/$(basename ${LTR_file})_index -L 20 --very-sensitive-local > $SAM_FILE 2>> "$LOG_FILE"

    # Process SAM files
    echo "Processing SAM to BAM for $LTR_file" >> "$LOG_FILE"
    samtools view -bS $SAM_FILE | samtools sort -T $PROCESS_TMPDIR -o $SORTED_BAM
    samtools index $SORTED_BAM
    
    echo "Calculating depth for $LTR_file" >> "$LOG_FILE"
    samtools depth -aa $SORTED_BAM > $DEPTH_FILE

    # Clean up temporary files
    rm -rf $PROCESS_TMPDIR

    echo "Processing complete for $LTR_file at $(date)" >> "$LOG_FILE"
}

export -f process_ltr
export READ_FILE1 READ_FILE2 LOG_FILE OUTPUT_DIR SPECIES

# Process all LTRs in parallel
total_LTRs=$(ls $LTR_DIR/*.fa | wc -l)
echo "Starting pipeline with $total_LTRs LTR files at $(date)" >> "$LOG_FILE"

# Process each LTR file (using find + parallel for parallel processing)
find "$LTR_DIR" -name "*.fa" | xargs -n 1 -P 4 -I {} bash -c 'process_ltr "$@"' _ {}

echo "Pipeline completed at $(date)" >> "$LOG_FILE"
echo "Output depth files are in ${OUTPUT_DIR}/Depth/${SPECIES}/"

# Run the T:I ratio calculation script if it exists
if [ -f "scripts/calculate_ti_ratio.sh" ]; then
    echo "Calculating T:I ratios..."
    bash scripts/calculate_ti_ratio.sh
    echo "T:I ratio calculation complete."
else
    echo "T:I ratio calculation script not found. Skipping."
fi

# Generate plots if the R script exists
if [ -f "scripts/visualization/ectopic_recomb_plots.R" ]; then
    echo "Generating plots..."
    Rscript scripts/visualization/ectopic_recomb_plots.R
    echo "Plots generated."
else
    echo "R visualization script not found. Skipping plot generation."
fi 

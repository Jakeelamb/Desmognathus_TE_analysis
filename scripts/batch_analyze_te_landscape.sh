#!/bin/bash
#
# batch_analyze_te_landscape.sh
#
# Description: Batch analyzes TE landscape data for multiple SRX samples using GNU Parallel
#
# Usage: ./batch_analyze_te_landscape.sh [options]
#
# Options:
#   --all               Process all available samples
#   --sample-list FILE  Process samples listed in FILE (one SRX ID per line)
#   --samples SRX1,SRX2 Process specific samples (comma-separated)
#   --cpus INT          Number of parallel jobs to run (default: 4)
#   --help              Show this help message
#

# Exit on error
set -e

# --- Configuration ---
NUM_CPUS=4 # Default number of parallel jobs

# Function to show usage
show_usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --all               Process all available samples"
    echo "  --sample-list FILE  Process samples listed in FILE (one SRX ID per line)"
    echo "  --samples SRX1,SRX2 Process specific samples (comma-separated)"
    echo "  --cpus INT          Number of parallel jobs to run (default: $NUM_CPUS)"
    echo "  --help              Show this help message"
    echo ""
    echo "Example: $0 --all --cpus 8"
}

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

# Check for GNU Parallel
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU Parallel is not installed."
    echo "Please install it (e.g., 'conda install -c conda-forge parallel' or 'sudo apt-get install parallel')."
    exit 1
fi

# Define paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")" # Assumes script is in scripts/
ANALYSIS_SCRIPT="${PROJECT_ROOT}/scripts/analyze_te_landscape.sh"
ALIGN_DIR="${PROJECT_ROOT}/data/raw/repeatmasker"
CLASS_DIR="${PROJECT_ROOT}/data/interim"
OUTPUT_DIR="${PROJECT_ROOT}/results/landscapes"
FIGURE_DIR="${PROJECT_ROOT}/results/figures/landscape"
LOG_DIR="${PROJECT_ROOT}/results/logs" # Directory for logs
PARALLEL_LOG="${LOG_DIR}/batch_parallel.log"
SUMMARY_FILE="${OUTPUT_DIR}/landscape_analysis_summary.csv"

# Make sure the analysis script is executable
chmod +x "$ANALYSIS_SCRIPT"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FIGURE_DIR"
mkdir -p "$LOG_DIR"

# Function to get a list of all available samples
get_all_samples() {
    # Find all .align files and extract SRX IDs
    for file in "$ALIGN_DIR"/*.align; do
        # Check if file exists and is readable before processing
        if [ -f "$file" ]; then
            basename "$file" | sed 's/_Trinity.align//'
        fi
    done
}

# Function to check if both align and classification files exist
check_sample_files() {
    local srx_id="$1"
    local align_file="${ALIGN_DIR}/${srx_id}_Trinity.align"
    local class_file="${CLASS_DIR}/${srx_id}_reads_per_component_and_annotation_processed"

    if [ -f "$align_file" ] && [ -f "$class_file" ]; then
        return 0  # Both files exist
    else
        return 1  # One or both files missing
    fi
}

# Initialize variables
process_all=false
sample_list_file=""
samples_arg=""
declare -a sample_ids_to_process # Array to hold IDs to actually process

# Parse command line arguments
if [ $# -eq 0 ]; then
    show_usage
    exit 1
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --all)
            process_all=true
            shift
            ;;
        --sample-list)
            if [ -f "$2" ]; then
                sample_list_file="$2"
                shift 2
            else
                echo "Error: Sample list file '$2' not found"
                exit 1
            fi
            ;;
        --samples)
            samples_arg="$2"
            shift 2
            ;;
        --cpus)
            NUM_CPUS="$2"
            if ! [[ "$NUM_CPUS" =~ ^[0-9]+$ ]] || [ "$NUM_CPUS" -lt 1 ]; then
                echo "Error: --cpus requires a positive integer."
                exit 1
            fi
            shift 2
            ;;
        --help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# --- Determine and Filter Sample List ---
declare -a initial_sample_list # Array to hold initial list
total_initial_samples=0

if [ "$process_all" = true ]; then
    echo "Gathering all available samples..."
    mapfile -t initial_sample_list < <(get_all_samples)
    total_initial_samples="${#initial_sample_list[@]}"
    echo "Found $total_initial_samples potential samples."

elif [ -n "$sample_list_file" ]; then
    echo "Reading samples from list: $sample_list_file"
    mapfile -t initial_sample_list < "$sample_list_file"
    # Remove empty lines or comments
    initial_sample_list=("${initial_sample_list[@]/#\#*/}") # Remove comments
    initial_sample_list=("${initial_sample_list[@]//$'\n'/}") # Remove empty lines (may need adjustment)
    # A more robust way to filter empty lines:
    temp_list=()
    for item in "${initial_sample_list[@]}"; do
        if [[ -n "$item" ]]; then
            temp_list+=("$item")
        fi
    done
    initial_sample_list=("${temp_list[@]}")
    total_initial_samples="${#initial_sample_list[@]}"
    echo "Read $total_initial_samples samples from file."

elif [ -n "$samples_arg" ]; then
    echo "Processing specific samples: $samples_arg"
    IFS=',' read -ra initial_sample_list <<< "$samples_arg"
    total_initial_samples="${#initial_sample_list[@]}"
    echo "Processing $total_initial_samples specified samples."
else
    echo "Error: No processing option specified"
    show_usage
    exit 1
fi

# Filter the initial list to find samples with required files
echo "Checking required files for each sample..."
skipped_list=()
for srx_id in "${initial_sample_list[@]}"; do
    if check_sample_files "$srx_id"; then
        sample_ids_to_process+=("$srx_id")
    else
        echo "Skipping $srx_id - missing required files."
        skipped_list+=("$srx_id")
    fi
done

total_to_process="${#sample_ids_to_process[@]}"
total_skipped="${#skipped_list[@]}"

if [ "$total_to_process" -eq 0 ]; then
    echo "No valid samples found to process. Exiting."
    exit 0
fi

echo "$total_to_process samples will be processed."
if [ "$total_skipped" -gt 0 ]; then
    echo "$total_skipped samples were skipped."
fi

# --- Run Parallel Processing ---
echo "----------------------------------------"
echo "Starting parallel processing with $NUM_CPUS CPUs..."
echo "Log file: $PARALLEL_LOG"

# Export the path to the analysis script so parallel jobs can find it
export ANALYSIS_SCRIPT

# Run GNU Parallel
# We pass the list of samples via stdin
# --jobs specifies the number of parallel jobs
# --joblog records the status of each job
# --eta provides an estimated completion time
# {} is the placeholder for each input line (sample ID)
printf "%s\\n" "${sample_ids_to_process[@]}" | \
    parallel --jobs "$NUM_CPUS" --joblog "$PARALLEL_LOG" --eta --bar \
    "$ANALYSIS_SCRIPT" {}

echo "Parallel processing finished."
echo "----------------------------------------"

# --- Process Parallel Log and Generate Summary ---
echo "Processing job log and generating summary..."
successful_samples=0
failed_samples=0
failed_list=""

# Create CSV header for summary
echo "SRX_ID,Status,Date_Processed" > "$SUMMARY_FILE"

# Append skipped samples to summary
current_date=$(date '+%Y-%m-%d %H:%M:%S')
for srx_id in "${skipped_list[@]}"; do
    echo "$srx_id,Skipped,$current_date" >> "$SUMMARY_FILE"
done

# Process the parallel job log
# Skip the header line and read tab-separated fields
{
    # Skip header line
    read
    while IFS=$'\t' read -r seq host starttime jobruntime send recv exitval signal command; do
        # Extract SRX ID from the command (it's the last argument)
        srx_id=$(echo "$command" | awk '{print $NF}')
        
        # Check exit value (0 = success)
        if [ "$exitval" -eq 0 ]; then
            ((successful_samples++))
            status="Success"
        else
            ((failed_samples++))
            failed_list="$failed_list $srx_id"
            status="Failed (Exit: $exitval)"
        fi
        
        # Add to summary file
        echo "$srx_id,$status,$current_date" >> "$SUMMARY_FILE"
    done
} < "$PARALLEL_LOG"

# Sort summary file by SRX_ID (optional, keeps it tidy)
{ head -n 1 "$SUMMARY_FILE" && tail -n +2 "$SUMMARY_FILE" | sort -t, -k1; } > "${SUMMARY_FILE}.tmp" && mv "${SUMMARY_FILE}.tmp" "$SUMMARY_FILE"

# Print summary
echo "----------------------------------------"
echo "Batch Processing Summary"
echo "----------------------------------------"
echo "Total samples attempted: $total_initial_samples"
echo "Samples skipped (missing files): $total_skipped"
echo "Samples processed: $total_to_process"
echo "Successfully processed: $successful_samples"
echo "Failed: $failed_samples"
if [ "$failed_samples" -gt 0 ]; then
    echo "Failed samples:$failed_list"
fi
echo "----------------------------------------"
echo "Parallel job log: $PARALLEL_LOG"
echo "Summary saved to: $SUMMARY_FILE"
echo "Output files in: $OUTPUT_DIR"
echo "Plots in: $FIGURE_DIR"

# --- Run Combined Visualization ---
if [ "$successful_samples" -gt 0 ]; then
    echo "----------------------------------------"
    echo "Generating combined visualizations..."

    VISUALIZATION_SCRIPT="${PROJECT_ROOT}/scripts/R/visualization/visualize_all_landscapes.R"

    if [ -f "$VISUALIZATION_SCRIPT" ]; then
        chmod +x "$VISUALIZATION_SCRIPT"
        Rscript "$VISUALIZATION_SCRIPT"
        if [ $? -eq 0 ]; then
            echo "Combined visualizations completed successfully."
        else
            echo "Error: Combined visualization script failed."
        fi
    else
        echo "Warning: Visualization script not found: $VISUALIZATION_SCRIPT"
    fi
else
    echo "Skipping combined visualization as no samples were processed successfully."
fi

echo "----------------------------------------"
echo "Script finished." 
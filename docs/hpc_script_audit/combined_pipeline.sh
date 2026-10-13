#!/bin/bash
#SBATCH --job-name=EctopicRecomb
#SBATCH --partition=week-long-cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o EctopicRecomb_%j.out
#SBATCH -e EctopicRecomb_%j.err

set -e  # Exit immediately if a command exits with a non-zero status

# Set custom temporary directory
export TMPDIR="/nfs/home/jlamb/Projects/Necturus/Data/tmp"
mkdir -p $TMPDIR

# Load environment
source ~/.bashrc
conda activate map_trim || { echo "Failed to activate conda environment"; exit 1; }

cd "/nfs/home/jlamb/Projects/Necturus/" || { echo "Failed to change directory"; exit 1; }

# Input parameters
READ_FILE1=$1
READ_FILE2=$2
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$
READ_FILE=$

species=$3
LTR_dir="Data/LTRs"
Log_file="Data//Log_file.txt"

# Check if input files exist
[[ -f $READ_FILE1 ]] || { echo "Read file 1 not found: $READ_FILE1"; exit 1; }
[[ -f $READ_FILE2 ]] || { echo "Read file 2 not found: $READ_FILE2"; exit 1; }
[[ -d $LTR_dir ]] || { echo "LTR directory not found: $LTR_dir"; exit 1; }

# Function to process a single LTR pair
process_ltr() {
    local LTR_file=$1

    # Create a unique temporary directory for this process
    local PROCESS_TMPDIR="${TMPDIR}/$(basename ${LTR_file})_${RANDOM}"
    mkdir -p $PROCESS_TMPDIR
    
    echo "Processing file: $LTR_file" >> "$Log_file"
    
    # Define output files
    local SAM_FILE="${PROCESS_TMPDIR}/$(basename ${LTR_file}).sam"
    local SORTED_BAM="${LTR_file}.sorted.bam"
    local DEPTH_FILE="${LTR_file}.depth.txt"
    
    # Build indices and map reads
    bowtie2-build $LTR_file ${PROCESS_TMPDIR}/$(basename ${LTR_file})_index >> "$Log_file" 2>&1
    bowtie2 -1 $READ_FILE1 -2 $READ_FILE2 -x ${PROCESS_TMPDIR}/$(basename ${LTR_file})_index -L 20 --very-sensitive-local > $SAM_FILE 2>> "$Log_file"

    # Process SAM files
    samtools view -bS $SAM_FILE | samtools sort -T $PROCESS_TMPDIR -o $SORTED_BAM
    samtools index $SORTED_BAM
    samtools depth -aa $SORTED_BAM > $DEPTH_FILE

    # Clean up temporary files
    rm -rf $PROCESS_TMPDIR

    echo "Processing complete for $LTR_file at $(date)" >> "$Log_file"
}

export -f process_ltr
export READ_FILE1 READ_FILE2 Log_file

# Process all LTRs in parallel
total_LTRs=$(ls $LTR_dir/*.fa | wc -l)
echo "Starting pipeline with $total_LTRs LTR files at $(date)" >> "$Log_file"

find "$LTR_dir" -name "*.fa" | parallel -j 4 process_ltr

echo "Pipeline completed at $(date)" >> "$Log_file"

#!/bin/bash

# --- Configuration ---
DATA_DIR="Data"
READS_DIR="${DATA_DIR}/Reads"
INDEX_DIR="${DATA_DIR}/Indexs"
LOOKUP_FILE="${DATA_DIR}/Lookup_table.txt"
RESULTS_DIR="Results/Alignment"
NUM_CPUS=4 # Adjust based on your system

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR" || { echo "Failed to create results directory: $RESULTS_DIR"; exit 1; }

# --- Setup ---
echo "Starting read mapping at $(date)"

# --- Input Checks ---
if [ ! -f "$LOOKUP_FILE" ]; then
    echo "Lookup table not found: $LOOKUP_FILE"
    exit 1
fi

if [ ! -d "$READS_DIR" ]; then
    echo "Reads directory not found: $READS_DIR"
    exit 1
fi

if [ ! -d "$INDEX_DIR" ]; then
    echo "Index directory not found: $INDEX_DIR"
    exit 1
fi

# --- Main Processing Loop ---
# Read the lookup table and process each line
while IFS=$'\t' read -r species sra genome || [ -n "$genome" ]; do # Handle lines without trailing newline
    # Skip the header line
    if [ "$species" = "Species" ]; then
        continue
    fi

    # Define input and output files
    read_file1="${READS_DIR}/${sra}_nuclear_1.fastq"
    read_file2="${READS_DIR}/${sra}_nuclear_2.fastq"
    index_prefix="${INDEX_DIR}/${genome}" # Assumes index base name is the genome identifier
    output_bam="${RESULTS_DIR}/${species}_${genome}_aligned.bam" # More informative BAM name

    echo "--- Processing: $species ($sra mapped to $genome) ---"
    echo "Read file 1: $read_file1"
    echo "Read file 2: $read_file2"
    echo "Index prefix: $index_prefix"
    echo "Output BAM: $output_bam"

    # Check required files for this entry
    if [ ! -f "$read_file1" ] || [ ! -f "$read_file2" ]; then
        echo "Warning: Read files not found for SRA $sra. Skipping $species."
        continue
    fi
    # Check if *any* index file exists with this prefix (bowtie2 needs multiple .bt2 files)
    if ! ls "${index_prefix}".*.bt2 1> /dev/null 2>&1; then
         echo "Warning: Bowtie2 index files not found with prefix $index_prefix. Skipping $species."
         echo "Ensure index was built using 'bowtie2-build <ref.fasta> $index_prefix'"
         continue
    fi

    echo "Starting mapping for $species at $(date)"

    # Run bowtie2 and pipe to samtools to create BAM file
    bowtie2 -x "$index_prefix" \
            -1 "${read_file1}" \
            -2 "${read_file2}" \
            -p $NUM_CPUS | \
    samtools view -bS - > "$output_bam"

    # Check bowtie2 exit status
    if [ $? -ne 0 ]; then
        echo "Error: bowtie2 mapping failed for $species. Check logs."
        # Optional: Decide whether to exit or continue with next species
        # exit 1
    else
        echo "Finished mapping $species at $(date)"
    fi
    echo "---------------------------------------------------"

done < "$LOOKUP_FILE"

echo "Read mapping completed at $(date)" 
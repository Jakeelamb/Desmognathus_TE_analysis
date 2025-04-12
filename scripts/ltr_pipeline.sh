#!/bin/bash

# --- Configuration ---
GENOME_DIR="Genomes"
LOOKUP_FILE="Data/Lookup_table.txt"
ANNOTATION_DIR="Annotations"
OUTPUT_DIR="Output/LTR_Analysis"
NUM_CPUS=4  # Adjust based on your system

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create output directory: $OUTPUT_DIR"; exit 1; }

# --- Input Validation ---
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_genome_fasta>"
    exit 1
fi

GENOME_FILE_PATH="$1"
if [ ! -f "$GENOME_FILE_PATH" ]; then
    echo "Genome file not found: $GENOME_FILE_PATH"
    exit 1
fi

# Extract filename without extension (e.g., GCA_030264955.1)
filename=$(basename "$GENOME_FILE_PATH" .fna)
gca_input="$filename"

# Construct expected LTR annotation file path
LTR_ANNOT_FASTA="${ANNOTATION_DIR}/${filename}_LTR_annoted.fasta"
if [ ! -f "$LTR_ANNOT_FASTA" ]; then
    echo "LTR annotation fasta file not found: $LTR_ANNOT_FASTA"
    echo "Please ensure it exists in the '$ANNOTATION_DIR' directory and follows the naming convention '<filename>_LTR_annoted.fasta'."
    exit 1
fi

if [ ! -f "$LOOKUP_FILE" ]; then
    echo "Lookup table not found: $LOOKUP_FILE"
    exit 1
fi

# --- Lookup ---
scientific_name=""
srx=""

# Read the lookup file and find the matching entry
while IFS=$'\t' read -r name srx_value gca_value || [ -n "$gca_value" ]; do
    trimmed_gca_value=$(echo "$gca_value" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    if [ "$trimmed_gca_value" = "$gca_input" ]; then
        scientific_name="$name"
        srx="$srx_value"
        break
    fi
done < "$LOOKUP_FILE"

# Check if a match was found
if [ -z "$scientific_name" ]; then
    echo "Warning: No match found for GCA '$gca_input' in lookup file '$LOOKUP_FILE'. Continuing without scientific name."
fi

echo "--- Processing Genome: $filename ---"
echo "Scientific Name: ${scientific_name:-N/A}"
echo "SRX: ${srx:-N/A}"
echo "Genome File: $GENOME_FILE_PATH"
echo "LTR Annotation Fasta: $LTR_ANNOT_FASTA"
echo "Output Directory: $OUTPUT_DIR"
echo "-------------------------------------"

# Define output file paths
GFF3_OUT="${OUTPUT_DIR}/${filename}.gff3"
SUFFIXERATOR_IDX="${OUTPUT_DIR}/${filename}"
LTRHARVEST_OUT="${OUTPUT_DIR}/${filename}.ltrharvest.out"
LTRHARVEST_OUTINNER="${OUTPUT_DIR}/${filename}.ltrharvest.outinner"
LTRHARVEST_GFF3_UNSORTED="${OUTPUT_DIR}/${filename}.ltrharvest.gff3"
SORTED_GFF3="${OUTPUT_DIR}/${filename}_sorted.gff3"
LTRDIGEST_PREFIX="${OUTPUT_DIR}/${filename}.ltrdigest"

# --- LTR Annotation Pipeline ---

# 1. Run DANTE to generate GFF3 annotation from LTR fasta
echo "Step 1: Running DANTE..."
dante -q "$LTR_ANNOT_FASTA" -D Metazoa_v3.1 -o "$GFF3_OUT" -c $NUM_CPUS || { echo "DANTE failed"; exit 1; }
echo "DANTE finished."

# 2. Create suffixerator index for the genome
echo "Step 2: Creating genome index with gt suffixerator..."
gt suffixerator -db "$GENOME_FILE_PATH" -indexname "$SUFFIXERATOR_IDX" -tis -suf -lcp -des -ssp -sds -dna || { echo "gt suffixerator failed"; exit 1; }
echo "Suffixerator finished."

# 3. Run LTRharvest using the DANTE GFF as hints
echo "Step 3: Running LTRharvest..."
gt ltrharvest -index "$SUFFIXERATOR_IDX" -gff3 "$GFF3_OUT" -out "$LTRHARVEST_OUT" -outinner "$LTRHARVEST_OUTINNER" -outfile "$LTRHARVEST_GFF3_UNSORTED" || { echo "gt ltrharvest failed"; exit 1; }
echo "LTRharvest finished."

# 4. Sort the LTRharvest GFF3 output
echo "Step 4: Sorting LTRharvest GFF3 output..."
gt gff3 -sort "$LTRHARVEST_GFF3_UNSORTED" > "$SORTED_GFF3" || { echo "gt gff3 sort failed on LTRHarvest output"; exit 1; }
echo "GFF sorting finished."

# 5. Run LTRdigest
echo "Step 5: Running LTRdigest..."
gt ltrdigest -outfileprefix "$LTRDIGEST_PREFIX" "$SORTED_GFF3" "$SUFFIXERATOR_IDX" || { echo "gt ltrdigest failed"; exit 1; }
echo "LTRdigest finished."

echo "--- LTR Pipeline completed for $filename ---" 
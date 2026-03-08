#!/bin/bash

set -euo pipefail

echo "Warning: scripts/calculate_ti_ratio.sh is an upstream raw-data utility that still expects the historical Data/ and Output/ layout."
echo "The canonical paper-facing ectopic workflow is scripts/processing/ec.py plus scripts/visualization/plot_ectopic_recombination.R."

# --- Configuration ---
DATA_DIR="Data"
LOOKUP_FILE="${DATA_DIR}/Lookup_table.txt"
LTR_ANALYSIS_DIR="Output/LTR_Analysis"
DEPTH_DIR_BASE="Output/Depth" # Base directory where depth files reside (e.g., Output/Depth/SpeciesName/GenomeName/)
RESULTS_DIR="Results/Coverage"
OUTPUT_FILE="${RESULTS_DIR}/T_I_ratios.tsv"

# --- Setup ---
echo "Starting T:I ratio calculation at $(date)"

# Create results directory
mkdir -p "$RESULTS_DIR" || { echo "Failed to create results directory: $RESULTS_DIR"; exit 1; }

# --- Input File Checks ---
if [ ! -f "$LOOKUP_FILE" ]; then
    echo "Error: Lookup table not found: $LOOKUP_FILE"
    exit 1
fi

# --- Initialize Output File ---
# Use Tab Separated Values (TSV)
echo -e "Species\tGenome\tLTR_ID\tContig\tTotal_length\tLeft_LTR_length\tInternal_length\tRight_LTR_length\tLeft_LTR_MeanCov\tLeft_LTR_MinCov\tLeft_LTR_MaxCov\tInternal_MeanCov\tInternal_MinCov\tInternal_MaxCov\tRight_LTR_MeanCov\tRight_LTR_MinCov\tRight_LTR_MaxCov\tCombined_LTR_MeanCov" > "$OUTPUT_FILE"

# --- Main Processing Loop --- 
# Read the lookup table
while IFS=$'\t' read -r species sra genome || [ -n "$genome" ]; do
    # Skip the header line
    [[ "$species" == "Species" ]] && continue

    echo "--- Processing species: $species, genome: $genome ---"

    # Define expected input paths for this genome
    ltr_gff3_file="${LTR_ANALYSIS_DIR}/${genome}.ltrdigest.gff3"
    depth_dir="${DEPTH_DIR_BASE}/${species}/${genome}"

    # Check required inputs for this species/genome
    if [ ! -f "$ltr_gff3_file" ]; then
        echo "Warning: LTR GFF3 file not found: $ltr_gff3_file. Skipping $species."
        continue
    fi
    if [ ! -d "$depth_dir" ]; then
        echo "Warning: Depth directory not found: $depth_dir. Skipping $species."
        continue
    fi

    # Find all unique LTR IDs from the depth files in the directory
    # Assumes depth files are named like <ltr_id>.fasta.depth.txt
    find "$depth_dir" -maxdepth 1 -name '*.fasta.depth.txt' -printf '%f\n' | sed 's/\.fasta\.depth\.txt$//' | while read -r ltr_id; do
        
        echo "  Processing LTR ID: $ltr_id"
        depth_file="${depth_dir}/${ltr_id}.fasta.depth.txt"

        if [ ! -f "$depth_file" ]; then
             echo "    Warning: Depth file expected but not found: $depth_file. Skipping LTR ID $ltr_id."
             continue
        fi

        # --- Extract coordinates from LTR GFF3 --- 
        # Use awk to parse the GFF3 for this specific LTR ID
        coords=$(awk -F'\t' -v id="$ltr_id" '
            BEGIN { OFS="\t"; contig=""; ltr5s=-1; ltr5e=-1; ltr3s=-1; ltr3e=-1; }
            $9 ~ ("Parent=" id "($|;)") && $3 == "repeat_region" {
                if ($9 ~ /type=LTR5/) { ltr5s=$4; ltr5e=$5; contig=$1; }
                if ($9 ~ /type=LTR3/) { ltr3s=$4; ltr3e=$5; contig=$1; }
            }
            END {
                if (contig != "" && ltr5s != -1 && ltr3s != -1) {
                    print contig, ltr5s, ltr5e, ltr3s, ltr3e;
                } else {
                    # Print dummy values if not found
                    print "NA", -1, -1, -1, -1;
                }
            }
        ' "$ltr_gff3_file")

        read -r contig left_LTR_start left_LTR_end right_LTR_start right_LTR_end <<< "$coords"

        if [ "$contig" == "NA" ] || [ "$left_LTR_start" -eq -1 ] || [ "$right_LTR_start" -eq -1 ]; then
            echo "    Warning: Could not find 5' and 3' LTR coordinates for Parent ID '$ltr_id' in $ltr_gff3_file. Skipping."
            continue
        fi
        
        # Calculate internal coordinates (handle adjacent LTRs)
        internal_start=$((left_LTR_end + 1))
        internal_end=$((right_LTR_start - 1))
        if [ "$internal_start" -gt "$internal_end" ]; then
            internal_start=-1 # Mark as no internal sequence
            internal_end=-1
        fi

        echo "    Coordinates found: Contig=$contig, LTR5=[$left_LTR_start-$left_LTR_end], Internal=[$internal_start-$internal_end], LTR3=[$right_LTR_start-$right_LTR_end]"

        # --- Function to calculate coverage stats --- 
        calculate_coverage_stats() {
            local file=$1
            local start=$2
            local end=$3
            local description=$4 # For error messages

            # Handle invalid ranges (e.g., no internal sequence)
            if [ "$start" -lt 0 ] || [ "$end" -lt 0 ] || [ "$start" -gt "$end" ]; then
                echo "0 0 0" # Mean, Min, Max = 0 if range is invalid
                return
            fi

            # awk processes the depth file (contig pos depth)
            awk -F'\t' -v s="$start" -v e="$end" '
                BEGIN { sum=0; count=0; min=-1; max=-1; }
                $2 >= s && $2 <= e {
                    sum += $3;
                    count++;
                    if (min == -1 || $3 < min) { min = $3; }
                    if (max == -1 || $3 > max) { max = $3; }
                }
                END {
                    if (count > 0) {
                        printf "%.2f\t%d\t%d\n", sum / count, min, max;
                    } else {
                        # No coverage data in the range
                        printf "0.00\t0\t0\n"; 
                    }
                }
            ' "$file"
        }

        # --- Calculate coverage for each region --- 
        read -r left_mean left_min left_max <<< $(calculate_coverage_stats "$depth_file" "$left_LTR_start" "$left_LTR_end" "Left LTR")
        read -r int_mean int_min int_max <<< $(calculate_coverage_stats "$depth_file" "$internal_start" "$internal_end" "Internal")
        read -r right_mean right_min right_max <<< $(calculate_coverage_stats "$depth_file" "$right_LTR_start" "$right_LTR_end" "Right LTR")

        # --- Calculate lengths and combined stats --- 
        left_len=$(echo "$left_LTR_end - $left_LTR_start + 1" | bc)
        right_len=$(echo "$right_LTR_end - $right_LTR_start + 1" | bc)
        if [ "$internal_start" -eq -1 ]; then
            int_len=0
        else
            int_len=$(echo "$internal_end - $internal_start + 1" | bc)
        fi
        total_len=$(echo "$left_len + $int_len + $right_len" | bc)

        # Combined LTR mean coverage (weighted average? or simple average? Using simple average of means for now)
        # Avoid division by zero if one LTR has zero coverage/length?
        # Let's just average the means if both are non-zero
        combined_ltr_mean=0
        count_ltrs=0
        if [ $(echo "$left_mean > 0" | bc) -eq 1 ]; then 
            combined_ltr_mean=$(echo "$combined_ltr_mean + $left_mean" | bc)
            count_ltrs=$((count_ltrs + 1))
        fi
        if [ $(echo "$right_mean > 0" | bc) -eq 1 ]; then
            combined_ltr_mean=$(echo "$combined_ltr_mean + $right_mean" | bc)
            count_ltrs=$((count_ltrs + 1))
        fi
        if [ $count_ltrs -gt 0 ]; then
            combined_ltr_mean=$(echo "scale=2; $combined_ltr_mean / $count_ltrs" | bc)
        fi

        # --- Write results to output file ---
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n" \
            "$species" "$genome" "$ltr_id" "$contig" "$total_len" \
            "$left_len" "$int_len" "$right_len" \
            "$left_mean" "$left_min" "$left_max" \
            "$int_mean" "$int_min" "$int_max" \
            "$right_mean" "$right_min" "$right_max" \
            "$combined_ltr_mean" >> "$OUTPUT_FILE"

    done # End loop through LTR IDs

done < "$LOOKUP_FILE" # End loop through species

echo "T:I ratio calculation completed at $(date)"
echo "Output written to $OUTPUT_FILE"

exit 0 

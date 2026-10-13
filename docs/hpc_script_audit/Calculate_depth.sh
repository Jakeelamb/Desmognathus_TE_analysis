#!/bin/bash
#SBATCH --job-name=Calculate_depth
#SBATCH --output=Calculate_depth.out
#SBATCH --error=Calculate_depth.err
#SBATCH --partition=short-cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

set -euo pipefail

echo "Current dir: $(pwd)"
cd "/nfs/home/jlamb/Projects/Ectopic_recombination/"
echo "Changing directory to: /nfs/home/jlamb/Projects/Ectopic_recombination/"

# Load environment
source ~/.bashrc
conda activate nextflow_ltr

# Start logging
echo "Job started at $(date)"

# Initialize the output file with a header
echo -e "Species\tLTR_name\tTotal_length\tLeft_LTR_length\tInternal_sequence_length\tRight_LTR_length\tLTR_mean_coverage\tLeft_LTR_min_coverage\tLeft_LTR_max_coverage\tInternal_sequence_mean_coverage\tInternal_sequence_min_coverage\tInternal_sequence_max_coverage\tRight_LTR_mean_coverage\tRight_LTR_min_coverage\tRight_LTR_max_coverage" > Data/coverage.txt

# Read the lookup table and process each line
while IFS=$'\t' read -r species _ genome; do
    # Skip the header line
    [[ "$species" == "Species" ]] && continue

    echo "Processing species: $species"

    # Define input and output files
    output_coverage="Data/${species}_coverage.txt"
    Genome="Data/${species}/${genome}.fna"
    LTR_dir="Data/${species}/LTRs"
    CSV="Data/${species}/${genome}_tabout.csv"

    #for each LTR, calculate the depth of coverage for the leading LTR, internal sequence, and lagging LTR
    for LTR in "$LTR_dir"/*.fa.depth.txt; do

        #extract the first 20 characters of the LTR name
        LTR_name=$(basename "$LTR" .fa.depth.txt | cut -c1-20)
        echo "Processing LTR: $LTR_name"

        # Extract coordinates from the CSV file
        read -r Total_length left_LTR_start_coord left_LTR_end_coord right_LTR_start_coord right_LTR_end_coord < <(awk -F',' -v ltr="$LTR_name" '$1 == ltr {print $3,$5,$6,$8,$9}' "$CSV")

        # Calculate internal sequence coordinates
        internal_sequence_start_coord=$((left_LTR_end_coord + 1))
        internal_sequence_end_coord=$((right_LTR_start_coord - 1))

        echo "Coordinates for $LTR_name:"
        echo "Total length: $Total_length"
        echo "Left LTR: $left_LTR_start_coord - $left_LTR_end_coord"
        echo "Internal: $internal_sequence_start_coord - $internal_sequence_end_coord"
        echo "Right LTR: $right_LTR_start_coord - $right_LTR_end_coord"

        # Use a function to calculate coverage statistics
        calculate_coverage_stats() {
            awk -F',' -v start="$1" -v end="$2" '
                $1 >= start && $1 <= end {
                    sum += $3
                    count++
                    if (NR == 1 || $3 < min) min = $3
                    if (NR == 1 || $3 > max) max = $3
                }
                END {
                    print sum / count, min, max
                }
            ' "$LTR"
        }

        # Calculate coverage statistics for each region
        read -r left_LTR_mean_coverage left_LTR_min_coverage left_LTR_max_coverage < <(calculate_coverage_stats "$left_LTR_start_coord" "$left_LTR_end_coord")
        read -r internal_sequence_mean_coverage internal_sequence_min_coverage internal_sequence_max_coverage < <(calculate_coverage_stats "$internal_sequence_start_coord" "$internal_sequence_end_coord")
        read -r right_LTR_mean_coverage right_LTR_min_coverage right_LTR_max_coverage < <(calculate_coverage_stats "$right_LTR_start_coord" "$right_LTR_end_coord")

        # Sum the Mean left and right LTR coverage
        LTR_mean_coverage=$(echo "$left_LTR_mean_coverage + $right_LTR_mean_coverage" | bc)

        #additional variables to output
        left_LTR_length=$(echo "$left_LTR_end_coord - $left_LTR_start_coord" | bc)
        right_LTR_length=$(echo "$right_LTR_end_coord - $right_LTR_start_coord" | bc)
        internal_sequence_length=$(echo "$internal_sequence_end_coord - $internal_sequence_start_coord" | bc)

        # Write the results to a file
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%.2f\t%s\t%s\t%.2f\t%s\t%s\n" \
            "$species" "$LTR_name" "$Total_length" "$left_LTR_length" "$internal_sequence_length" "$right_LTR_length" \
            "$LTR_mean_coverage" "$left_LTR_min_coverage" "$left_LTR_max_coverage" \
            "$internal_sequence_mean_coverage" "$internal_sequence_min_coverage" "$internal_sequence_max_coverage" \
            "$right_LTR_mean_coverage" "$right_LTR_min_coverage" "$right_LTR_max_coverage" >> Data/coverage.txt

    done
done < Data/Lookup_table.txt

echo "Job finished at $(date)"

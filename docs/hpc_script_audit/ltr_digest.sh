#!/bin/bash
#SBATCH --job-name=gt_ltr
#SBATCH --partition=day-long-highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH -o LTR_%j.o
#SBATCH -e LTR_%j.e

# Corrected: removed the dot after the tilde
source ~/.bashrc
conda activate dante_ltr

# Initialize variables
file="$1"  # Added quotes for safety
# example = GCA_030264955.1.fna
# extract the filename without the extension, should look like this: GCA_030264955.1
filename=$(basename "$file" .fna)  # Added quotes
gca_input="$filename"
scientific_name=""
srx=""
gff3="${filename}.gff3"  # Added quotes

# Read the lookup file and find the matching entry
while IFS=$'\t' read -r name srx_value gca_value; do
    if [ "$gca_value" = "$gca_input" ]; then  # Changed == to = for POSIX compatibility
        scientific_name="$name"
        srx="$srx_value"
        break
    fi
done < /nfs/home/jlamb/Projects/RepeatMasker/Data/lookup.txt

# Check if a match was found
if [ -z "$scientific_name" ] || [ -z "$srx" ]; then
    echo "No match found for GCA: $gca_input"
    exit 1
else
    echo "Scientific Name: $scientific_name"
    echo "SRX: $srx"
    echo "GCA: $gca_input"
fi

echo "Starting to run LTR annotation pipeline"
echo "filename is: $filename"
path="/nfs/home/jlamb/Projects/dnaPipeTE_rnd_2/${srx}/dnaPipeTE_15gb_10gc_RMt15_smpl_2/Annotation/LTR_annoted.fasta"
echo "Path to the LTR annotation file is: $path"

# Change to the directory where the genome file is located
cd /nfs/home/jlamb/Projects/RepeatMasker/Data || exit 1

# Step 3: Run LTRdigest
gt ltrdigest -outfileprefix "$filename" "${filename}.gff3_sorted.gff3" "$filename"

#!/bin/bash
#SBATCH --job-name=bowtie2_mapping
#SBATCH --partition=day-long-highmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH -o bowtie2_mapping_%j.out
#SBATCH -e bowtie2_mapping_%j.err

echo "Current dir:" pwd
cd "/nfs/home/jlamb/Projects/Ectopic_recombination/"
echo "Changing directory to: /nfs/home/jlamb/Projects/Ectopic_recombination/"

# Load environment
source ~/.bashrc
conda activate nextflow_ltr

# Start logging
echo "Job started at $(date)"

# Read the lookup table and process each line
while IFS=$'\t' read -r species sra genome; do
    # Skip the header line
    if [ "$species" = "Species" ]; then
        continue
    fi

    # Define input and output files
    read_file1="Data/Reads/${sra}_nuclear_1.fastq"
    read_file2="Data/Reads/${sra}_nuclear_2.fastq"
    index_prefix="Data/Indexs/${genome}.fna"
    output_bam="Results/${species}_aligned.bam"
    echo "Species: $species"
    echo "read file 1 : $read_file1"
    echo "read file 2 : $read_file2"
    echo "index file : $index_prefix"
    echo "output location : $output_bam"
    echo "Starting mapping for $species ($sra) to $genome at $(date)"

    # Run bowtie2 and pipe to samtools to create BAM file
    bowtie2 -x "$index_prefix" \
            -1 "${read_file1}" \
            -2 "${read_file2}" \
            -p $SLURM_CPUS_PER_TASK | \
    samtools view -bS - > "$output_bam"
#
    echo "Finished mapping $species ($sra) to $genome at $(date)"

done < Data/Lookup_table.txt

echo "Job finished at $(date)"

#!/bin/bash
#SBATCH --account=csu82_alpine1
#SBATCH --qos mem
#SBATCH --partition amem
#SBATCH --job-name=RepeatModeler
#SBATCH -c 32
#SBATCH --time=168:00:00
#SBATCH -N 1 #nodes
#SBATCH -o RM_%j.o
#SBATCH -e RM_%j.e

cd /scratch/alpine/jlamb1@colostate.edu

module purge
module load singularity/3.7.4

#$1 should be the path to the fasta file
#$2 should be the species name from the dictionary

# Run the commands
singularity exec -H /scratch/alpine/jlamb1@colostate.edu/ tetools_1.87.sif BuildDatabase -name "$2" "$1"
singularity exec -H /scratch/alpine/jlamb1@colostate.edu/ tetools_1.87.sif RepeatModeler -database "$2" -threads 32 -LTRStruct -genomeSampleSizeMax 10000000000

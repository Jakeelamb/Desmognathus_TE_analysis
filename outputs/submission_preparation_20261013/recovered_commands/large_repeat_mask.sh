#!/bin/bash
#SBATCH --partition=week-long-cpu
#SBATCH --job-name=RepeatMasker
#SBATCH --cpus-per-task=128
#SBATCH -N 1 #nodes
#SBATCH -o rmask_%j.o

Dir="/nfs/home/jlamb/Projects/RepeatMasker/Data"
cd $Dir

JOB_ID=$SLURM_JOB_ID
TMP_DIR="${Dir}/temp_$JOB_ID"
mkdir -p "$TMP_DIR" && export TMPDIR="$TMP_DIR"
mkdir -p "${1}_rmask_out"

Fasta=$1
source ~/.bashrc
conda activate /nfs/home/jlamb/bin/miniconda3/envs/RepeatMasker

/nfs/home/jlamb/bin/RepeatMasker/RepeatMasker -s -pa 32 -lib /nfs/home/jlamb/TE_libs/dedupe_telib.fasta\
 -dir "${Fasta}_rmask_out" "$Fasta"

exit



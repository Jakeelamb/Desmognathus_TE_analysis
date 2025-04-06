#!/bin/bash
#SBATCH --partition week-long-highmem
#SBATCH --time=168:00:00
#SBATCH --job-name=cdhit
#SBATCH --cpus-per-task=64
#SBATCH -N 1 #nodes
#SBATCH -o CDHIT_%j.o


cd /nfs/home/jlamb/TE_libs
source ~/.bashrc
conda activate /nfs/home/jlamb/bin/miniconda3/envs/RepeatMasker
cd-hit-est -i master_dfam_rmodeler_rmueller_TELib.fasta -o master_dfam_rmodeler_rmueller_TELib_80.fasta -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000
exit

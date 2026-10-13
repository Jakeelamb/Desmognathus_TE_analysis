---
tags: [project, domain/genomics, organism/desmognathus]
---
This is a README.md file for the project to keep track versions, scripts, arguments, and data origins 

# Data acquisition & processing

35 low coverage shot gun reads of desmognathus were downloaded from  ncbi <https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=43045>

4 species were not included form the total 39 due to contamination ( D.lycos, D. fuscus, D. aeneus, D.ocoee)

Reads were downloaded using the sratoolkit
```
fasterq-dump --split-files
```

The raw reads were trimmed using fastp with default parameters
```
fastp --in1 $1 --in2 $2 --out1=${BASENAME}_trimmed_1.fq.gz --out2=${BASENAME}_trimmed_2.fq.gz -t 32 --detect_adapter_for_pe
```

The desmognathus fuscus mitochondrial genome was then downloaded from ncbi 
https://www.ncbi.nlm.nih.gov/nuccore/NC_006339.1?report=fasta

The Mitochondrial genome was indexed and the raw reads of each genome were mapped against the index. 

```
bwa index $mito_genome
bwa mem -t 8 $mito_genome $reads1 $reads2 %3E ${base}_mapped_reads.sam
samtools view -@ 8 -bS ${base}_mapped_reads.sam > ${base}_mapped_reads.bam
samtools sort -@ 8 ${base}_mapped_reads.bam -o sorted_${base}_mapped_reads.bam
samtools view -@ 8 -b -f 12 -F 256 sorted_${base}_mapped_reads.bam > ${base}_nuclear_reads.bam
samtools fastq ${base}_nuclear_reads.bam > ${base}_nuclear.fastq
```

Sort bam and split back to fastq
```

samtools sort -n -o sorted.bam input.bam

samtools fastq input.bam -1 output_R1.fastq -2 output_R2.fastq
```

The resulting reads are now have the adapters, low quality reads, and mitochondrial reads removed. 

# TE library creation & curation

34 Desmognathus genomes, were downloaded from ncbi <https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=43045>

Ambystoma reference chromosomes were downloaded from ncbi https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_002915635.3/
The chromosomes arms were concatenated using
```
cat p#.fasta q#.fasta > Chrom#.fasta
```

Desmognathus Fuscus scaffolded genome was downloaded from ncbi : 

The scaffolds were split into 3g
```
(base) jlamb1@Ovis:~/Projects/Fuscus$ cat split_fasta.py
from Bio import SeqIO

def split_fasta_safe(fasta_file, output_prefix, max_size=3*10**9):  # Adjust max_size for 3GB
  reader = SeqIO.parse(fasta_file, "fasta")
  current_file = None
  current_size = 0

  for record in reader:
    record_size = len(f">{record.id}\n{record.seq}\n")  # Calculate record size
    if not current_file or current_size + record_size > max_size:
      if current_file:
        current_file.close()
      counter = 1 if not current_file else counter + 1
      current_file = open(f"{output_prefix}_{counter}.fasta", "w")
      current_size = 0
    current_file.write(f">{record.id}\n{record.seq}\n")
    current_size += record_size

  if current_file:
    current_file.close()

# Replace with your actual file names
split_fasta_safe("aDesFus_merged.purged2.fa", "Fuscus")
```

All 34 Desmognathus gneomes, the Amystoma chromosomes, and the spit D.fuscus genome were individually ran through repeatmodeler using these parameters 

```
module load singularity/3.7.4

#$1 should be the path to the fasta file
#$2 should be the species name from the dictionary

# Run the commands
singularity exec -H /scratch/alpine/jlamb1@colostate.edu/ tetools_1.87.sif BuildDatabase -name "$2" "$1"
singularity exec -H /scratch/alpine/jlamb1@colostate.edu/ tetools_1.87.sif RepeatModeler -database "$2" -threads 32 -LTRStruct -genomeSampleSizeMax 10000000000
```

![](https://lh7-us.googleusercontent.com/CnIcg2IB7lRsRD7iioSgKi2Pe4uhSzNiTlHPIcF5yBcrg1g6Mmmc-O7gIvP9m0ImyHSl2uIQ1wsEoqpNhh0jKyWUquISE9J_IMGxawgbVVyIPoWnSXiOkvrjRvFElZ2INSGx3nncaAba9BxZXokBIyDBkQ=s2048)

The unkown TEs were ran through DeepTE, TESorter, and TEClass2 for further annotation. 
```
Deep TE
source activate /nfs/home/jlamb/bin/miniconda3/envs/py36
python3 /nfs/home/jlamb/bin/DeepTE/DeepTE.py -d "$Working_dir" -o "$Output_dir" -i "$ipt_seq" -sp M -m_dir /nfs/home/jlamb/bin/Metazoans_model
exit
```

```
TEsorter
TEsorter $1 -db rexdb-metazoa -p 64
```

```
TEClass2
TEClass/TEclass-2.1.3d/TEclassTest.pl $1
```
The TE library was filtered using Wharton 808080 rule 

```
CDHIT-est 
cd-hit-est -i TELib.fasta -o TELib_80.fasta -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000
```

The DFAM TE libraries were downloaded https://www.dfam.org/releases/Dfam_3.8/families/FamDB/

The Amphibian TE libraries were parsed from embl to fasta using the DFAM toolkit
```
python3 famdb.py -i DFAM families Amphibia -a -d -f fasta_name --include-class-in-name > amphibia_a_d.fasta
```

The Fasta Libraries were merged

```
cat amphibia_a_d.fasta TELib_80 > Master_TElib.fasta
```

The library was again filtered using Wharton 808080 rule 

```
CDHIT-est 
cd-hit-est -i Master_TELib.fasta -o Master_TElib_80.fasta -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000
```

# dnaPipeTE 

We did some param scanning tests to figure out which parameters yielded the best results. 

```
import os
import subprocess

SRXs = ["SRX21425398", "SRX21394251"]
genome_coverages = [0.01, 0.05, 0.1, 0.5, 0.75, 0.1, 0.15, 0.2, 0.25, 0.3]
rm_t_values = [0.1, 0.15, 0.2, 0.25, 0.3]
sample_numbers = [1, 2, 3, 4, 5]

def submit_job(command):
    subprocess.run(command, shell=True, check=True)

# Function to submit jobs with different genome coverage values
def gc_compare(srx, rm_t):
    for gc in genome_coverages:
        for _ in range(2):  # Run twice for each setting
            output_dir = f"{srx}_GC_{gc}_RmT_{rm_t}_smpl_num_2"
            command = f"sbatch /nfs/home/jlamb/scripts/paramscan_dnaPipeTE.sh {srx} {gc} {output_dir} {rm_t} 2"
            submit_job(command)

# Function to submit jobs with different rm_t values
def rm_t_compare(srx, gc):
    for rm_t in rm_t_values:
        for _ in range(2):  # Run twice for each setting
            output_dir = f"{srx}_GC_{gc}_RmT_{rm_t}_smpl_num_2"
            command = f"sbatch /nfs/home/jlamb/scripts/paramscan_dnaPipeTE.sh {srx} {gc} {output_dir} {rm_t} 2"
            submit_job(command)

# Function to submit jobs with different sample numbers
def sample_num_compare(srx, gc, rm_t):
    for sample_num in sample_numbers:
        for _ in range(2):  # Run twice for each setting
            output_dir = f"{srx}_GC_{gc}_RmT_{rm_t}_smpl_num_{sample_num}"
            command = f"sbatch /nfs/home/jlamb/scripts/paramscan_dnaPipeTE.sh {srx} {gc} {output_dir} {rm_t} {sample_num}"
            submit_job(command)

# Submit jobs for different genome coverage values
for srx in SRXs:
    gc_compare(srx, 0.15)

# Submit jobs for different rm_t values
for srx in SRXs:
    rm_t_compare(srx, 0.1)

# Submit jobs for different sample numbers
for srx in SRXs:
    sample_num_compare(srx, 0.1, 0.15)
```

```
#DNAPipeTE script

SRX=$1
genome_coverage=$2
output=$3
#output dir name should match this pattern /nfs/home/jlamb/dnaPipeTE_test/${SRX}_GC_${2}_RmT_${4}_smpl_num_${5}
rm_t=$4
sample_num=$5

cd "/nfs/home/jlamb/Projects/dnaPipeTE_test/"

module load singularity/4.1.2
module load python3

# Use singularity exec instead of shell to run the subsequent commands in the container
singularity exec --bind /nfs/home/jlamb/Projects/dnaPipeTE_test/$1:/mnt ~/bin/dnaPipeTE/dnapipte.img /bin/bash <<EOF

# Set locale environment variables to avoid locale warnings
export LANG=C.UTF-8
export LC_ALL=C.UTF-8

cd /opt/dnaPipeTE

python3 dnaPipeTE.py \
-input /mnt/${1}_nuclear_1.fastq \
-output /mnt/${3} \
-cpu 8 \
-genome_size 15000000000 \
-genome_coverage ${2} \
-RM_lib /nfs/home/jlamb/TE_libs/master_dfam_rmodeler_rmueller_TELib.fasta \
-RM_t ${4} \
-sample_number ${5} \

EOF
```

Each Desmognathus sample was processed using genome_size 15Gb, genome_coverage__, RM_t__, sample_number__ using dnaPipeTE 

# Figure generation 

## TEs 
ec
Figure showing variance of total composition 
![[Pasted image 20240613162415.png]]
https://www.mdpi.com/genes/genes-11-00170/article_deploy/html/images/genes-11-00170-g001-550.jpg
### Plan to make it
- [x] Convert Counts.txt to a single tsv file
```
(map_trim) [jlamb@login001 Scripts]$ cat Counts.txt_to_dataframe.py
import pandas as pd
import os
import sys
import re

def txt_to_dataframe(file_path):
    """
    Converts a TXT file into a pandas DataFrame with a single row.

    Args:
        file_path (str): The path to the TXT file.

    Returns:
        pandas.DataFrame: A DataFrame containing the data from the TXT file with the SRX value as the first column.
    """
    # Extract the SRX value from the file path using regex
    srx_value = re.search(r'SRX\d+', file_path).group()
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                category, count = line.split('\t')
                data[category] = [int(count)]

    df = pd.DataFrame(data)
    df.insert(0, 'Species', srx_value)
    return df

def process_file(file_path):
    # Check if Counts.txt file exists
    if not os.path.isfile(file_path):
        srx_value = re.search(r'SRX\d+', file_path).group()
        print(f"Counts.txt does not exist for {srx_value}")
        return

    df = txt_to_dataframe(file_path)

    # Check if "Desmo_counts.tsv" exists
    if os.path.isfile("Desmo_counts.tsv"):
        # If it exists, append the DataFrame to the file
        df.to_csv("Desmo_counts.tsv", mode='a', header=False, sep='\t', index=False)
    else:
        # If it doesn't exist, create the file and write the header and the DataFrame
        df.to_csv("Desmo_counts.tsv", sep='\t', index=False)

# Example usage
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the file path as an argument.")
        sys.exit(1)

    file_path = sys.argv[1]
    process_file(file_path)
```
Then to run all of the output together run this command
```
ls -d1 /nfs/home/jlamb/Projects/dnaPipeTE_rnd_2/SRX*/*/Counts.txt | xargs -n 1 -I {} python3 Counts.txt_to_dataframe.py {}
```
- [ ] Make Pie charts for each species
- [ ] 
Figure showing variance among the TEs
![[Pasted image 20240613162402.png]]
https://www.mdpi.com/genes/genes-11-00170/article_deploy/html/images/genes-11-00170-g002-550.jpg

Distribution of TE age divergence between the genomes 
![[Pasted image 20240613162511.png]]
https://www.mdpi.com/genes/genes-11-00170/article_deploy/html/images/genes-11-00170-g003-550.jpg

Scatter plot and analysis of linear regression of the relative genome proportions of shared TE families between the Rhodnius genomes in terms of genome percentage (log10 scale). Each dot represents a shared TE family, defined by the highest BLAST hit between the TE family reference contig of each species
![[Pasted image 20240613164007.png]]
https://www.mdpi.com/genes/genes-11-00170/article_deploy/html/images/genes-11-00170-g003-550.jpg

*this would be interesting to run against all species but there could be cool findings comparing most diverse species in genus (genome size wise, body size wise, different life history characteristics, cell size, etc)*


![[Pasted image 20240624100726.png]]
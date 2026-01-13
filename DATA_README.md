# Data Acquisition Guide

This document describes how to obtain the input data required to run the Desmognathus TE analysis pipeline.

## Overview

The analysis requires data from 35 Desmognathus salamander species. Raw data is not included in this repository due to size constraints (~40GB total).

## Required Data Files

### 1. dnaPipeTE Output Files

**Location:** `input_data/dnaPipeTE/`

**Source:** Run dnaPipeTE on RNA-seq data from NCBI SRA

**Files needed:** One file per species named `{SRX_ID}_reads_per_component_and_annotation`

**How to obtain:**
1. Download RNA-seq reads from NCBI SRA using the accessions in `input_data/lookup_table.txt`
2. Run dnaPipeTE on each sample
3. Copy the `reads_per_component_and_annotation` output files to `input_data/dnaPipeTE/`

Example SRA accessions (see lookup_table.txt for complete list):
- D. fuscus: SRX20497025
- D. monticola: SRX19955091
- D. ocoee: SRX19952863

### 2. RepeatMasker Alignment Files

**Location:** `input_data/repeatmasker/`

**Source:** Run RepeatMasker on genome assemblies from NCBI

**Files needed:** One `.align` file per species named `{SRX_ID}_Trinity.align`

**How to obtain:**
1. Download genome assemblies from NCBI using GCA accessions in `input_data/lookup_table.txt`
2. Run RepeatMasker with the `-align` flag
3. Copy the `.align` output files to `input_data/repeatmasker/`

Example genome accessions:
- D. fuscus: GCA_032353935.1
- D. monticola: GCA_031753675.1
- D. ocoee: GCA_030264995.1

### 3. Phylogenetic Tree

**Location:** `input_data/phylogeny/`

**Files needed:** `desmo900dated_test.tre` (Newick format)

**Source:** Time-calibrated phylogeny from published Desmognathus phylogenomic analysis

### 4. Ectopic Recombination Data

**Location:** `input_data/ectopic_recombination/`

**Source:** TEsorter and coverage analysis output

**Files needed:**
- `GCA_*_tabout.csv` - TEsorter domain annotations per genome
- `combined_sequences.fasta.rexdb-metazoa.cls.tsv` - Combined classification
- `coverage_results.tsv` - Read depth coverage data
- `*.fa.depth.txt` - Per-element depth files

**How to obtain:**
1. Run TEsorter on LTR sequences extracted from each genome
2. Calculate read depth coverage using samtools depth
3. Organize outputs in `input_data/ectopic_recombination/`

### 5. Species Lookup Table

**Location:** `input_data/lookup_table.txt`

**This file IS included in the repository.** It maps species names to SRA and genome accessions.

Format (tab-separated):
```
Species    SRA_Accension    Genome_Accension
D.fuscus   SRX20497025      GCA_032353935.1
...
```

## Directory Structure

After obtaining all data, your `input_data/` directory should look like:

```
input_data/
├── dnaPipeTE/
│   ├── SRX19952657_reads_per_component_and_annotation
│   ├── SRX19952691_reads_per_component_and_annotation
│   └── ... (35 files total)
├── repeatmasker/
│   ├── SRX19952657_Trinity.align
│   ├── SRX19952691_Trinity.align
│   └── ... (35 files total)
├── phylogeny/
│   └── desmo900dated_test.tre
├── ectopic_recombination/
│   ├── GCA_030264455.1_tabout.csv
│   ├── combined_sequences.fasta.rexdb-metazoa.cls.tsv
│   ├── coverage_results.tsv
│   └── ... (depth files)
└── lookup_table.txt
```

## Species List

| Species | SRA Accession | Genome Accession |
|---------|---------------|------------------|
| D. abditus | SRX19952928 | GCA_030264455.1 |
| D. adatsihi | SRX20502416 | GCA_032275005.1 |
| D. aeneus | SRX19953415 | GCA_030264635.1 |
| D. amphileucus | SRX20502537 | GCA_032360725.1 |
| D. anicetus | SRX21425396 | GCA_033119125.1 |
| D. apalachicolae | SRX19952657 | GCA_030264955.1 |
| D. aureatus | SRX20502577 | GCA_032274985.1 |
| D. auriculatus | SRX19953419 | GCA_030264895.1 |
| D. bairdi | SRX21425397 | GCA_034783875.1 |
| D. balsameus | SRX20502826 | GCA_032353675.1 |
| D. campi | SRX21425398 | GCA_033118305.1 |
| D. carolinensis | SRX19952691 | GCA_031753635.1 |
| D. catahoula | SRX21425399 | GCA_034783935.1 |
| D. cheaha | SRX20502857 | GCA_032353815.1 |
| D. conanti | SRX19953420 | GCA_030264475.1 |
| D. fuscus | SRX20497025 | GCA_032353935.1 |
| D. gvnigeusgwotli | SRX20502827 | GCA_032353775.1 |
| D. intermedius | SRX20502776 | GCA_032275315.1 |
| D. kanawha | SRX20502939 | GCA_032353855.1 |
| D. lycos | SRX21425400 | GCA_035049765.1 |
| D. marmoratus | SRX19954879 | GCA_031753695.1 |
| D. mavrokoilius | SRX20503049 | GCA_032354855.1 |
| D. monticola | SRX19955091 | GCA_031753675.1 |
| D. ocoee | SRX19952863 | GCA_030264995.1 |
| D. orestes | SRX19952890 | GCA_030264915.1 |
| D. organi | SRX19952929 | GCA_030180145.1 |
| D. pascagoula | SRX19953381 | GCA_030264975.1 |
| D. perlapsus | SRX20503048 | GCA_032275285.1 |
| D. santeetlah | SRX19955089 | GCA_031753655.1 |
| D. tilleyi | SRX21394251 | GCA_034782015.1 |
| D. valentinei | SRX19953382 | GCA_030264875.1 |
| D. valtos | SRX20503015 | GCA_032357565.1 |
| D. welteri | SRX19958882 | GCA_030264495.1 |
| D. wrighti | SRX19958881 | GCA_030265035.1 |

## Verification

After setting up your data, run the verification script:

```bash
conda activate Dusky
python verify_setup.py
```

This will check that all required files are in place.

## Data Availability

Raw sequencing data is available from NCBI SRA under BioProject [ADD_BIOPROJECT_ID].

Genome assemblies are available from NCBI GenBank under the accessions listed above.

## Contact

For questions about data acquisition, please open an issue on the repository.

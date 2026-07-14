# Data Acquisition Guide

This document describes how to obtain the input data required to run the Desmognathus TE analysis pipeline.

## Overview

The active TE/genome analysis scope is the 18-species
`te_genome_primary_mediumplus` panel. The lookup and staged directories retain
additional resources for provenance and historical reconstruction, but those
resources are not automatically analysis inputs. Large sequencing and
genome-analysis inputs are not included in this repository due to size
constraints (~40GB total). `path_analysis/data/external/raw/` may include small
source-backed literature or external-data snapshots when their license and
sensitivity boundary are documented.

## Required Data Files

### 1. dnaPipeTE Output Files

**Location:** `input_data/dnaPipeTE/`

**Source:** Run dnaPipeTE on nuclear-filtered genomic WGS reads from NCBI SRA

**Files needed:** One file per species named `{SRX_ID}_reads_per_component_and_annotation`

**How to obtain:**
1. Download genomic WGS reads from NCBI SRA using the accessions in `input_data/lookup_table.txt`
2. Reproduce the nuclear-read filtering and run dnaPipeTE on one read mate
3. Copy the `reads_per_component_and_annotation` output files to `input_data/dnaPipeTE/`

The historical upstream command is recoverable from Git commit
`cc64ffb134968250671ff0e3456255063a1aba97` (`dnaPipeTE.sh`). It used a
15,000,000,000-bp configured genome, 0.1× coverage, two Trinity samples, a 0.15
RepeatMasker annotation threshold, and
`/nfs/home/jlamb/TE_libs/dedupe_telib.fasta`. The historical container digest,
library checksum, runtime logs, and sampling replicates are not stored, so the
reconstructed absolute-load output remains sensitivity-only.

Suffixes added after an SRX accession identify computational retries, not new
specimens or biological replicates. Select exactly one completed output per SRX
using a run manifest; never sum retries.

Example SRA accessions (see lookup_table.txt for complete list):
- D. fuscus: SRX20497025
- D. monticola: SRX19955091
- D. ocoee: SRX19952863

### 2. RepeatMasker Alignment Files

**Location:** `input_data/repeatmasker/`

**Source:** RepeatMasker `-align` output for the dnaPipeTE `Trinity.fasta`
repeat-contig assembly, using the custom repeat library

**Files needed:** One `.align` file per species named `{SRX_ID}_Trinity.align`

**How to obtain:**
1. Retain the `Trinity.fasta` repeat-contig assembly from each dnaPipeTE run
2. Run RepeatMasker with the same custom library and the `-align` flag
3. Copy the `.align` output files to `input_data/repeatmasker/`

Despite the `Trinity.align` name, these are not alignments over the GCA genome
assemblies. Assembly-dependent LTR-pair and ectopic-recombination inputs are
separate resources linked through `input_data/lookup_table.txt`.

### 3. Phylogenetic Tree

**Location:** `input_data/phylogeny/`

**Files needed:** `desmo900dated_test.tre` (Newick format)

**Focal source:** Collaborator-supplied time-calibrated *Desmognathus* tree.
Its exact publication, archive, tree identifier, and calibration record still
need to be supplied; the release audit locks its SHA-256 and derives only an
exact 18-tip, rounding-corrected tree.

**Published uncertainty sensitivity:** Stewart and Wiens (2025),
DOI `10.1016/j.ympev.2024.108272`, Supplementary File S3 (optimal dated tree)
and Supplementary File S4 (200 time-calibrated bootstrap trees). Publisher
archives, extracted files, checksums, and final-panel derivatives are recorded
in `results/data/corrected/phylogeny/phylogeny_source_registry_v1.csv`.

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

This file is a local required input and is not tracked in the cleaned repository. It maps species names to SRA and genome accessions.

Format (tab-separated):
```
Species    SRA_Accension    Genome_Accension
D.fuscus   SRX20497025      GCA_032353935.1
...
```

### 6. Microscopy release-audit inputs

The raw microscopy data and model artifacts live in the sibling
`/home/jake/Projects/cellprofiler_test` workspace. The corrected final-18 audit
records absolute source paths and SHA-256 hashes in
`results/data/corrected/microscopy/microscopy_source_registry_analysis18_v1.csv`.
It scores the exact archived production cell masks, audits the locally trained
nucleus model, compares six morphology estimators, builds conditional
hierarchical intervals, and withdraws the unsupported picogram conversion.

Run:

```bash
scripts/run_in_dusky.sh python scripts/processing/audit_microscopy_release.py
```

The approved release variable is a relative nuclear-IOD sensitivity index, not
absolute genome size. Historical bridge outputs remain preserved.

### 7. Corrected path-audit products

Corrected path inputs and results live under
`results/data/corrected/path_analysis/`. They are generated only from the
declared final-18 TE features, corrected microscopy estimators, relative-IOD
subsets, and corrected tree products. The builder never reads the historical
`genome_size_pg` field.

```bash
scripts/run_in_dusky.sh python scripts/processing/build_corrected_path_inputs.py
scripts/run_in_dusky.sh Rscript scripts/processing/audit_corrected_path_models.R --phase all
scripts/run_in_dusky.sh Rscript scripts/processing/simulate_corrected_path_calibration.R
scripts/run_in_dusky.sh python scripts/processing/summarize_corrected_path_audit.py
scripts/run_in_dusky.sh python scripts/processing/build_publication_audit_notebook.py
scripts/run_in_dusky.sh jupyter nbconvert --to notebook --execute notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb --inplace --ExecutePreprocessor.timeout=600
scripts/run_in_dusky.sh python scripts/processing/audit_publication_notebook.py
```

The output bundle includes every candidate ranking, basis-set component,
best-model edge, 200-tree fit, species-omission fit, simulation replicate,
review figure, and release gate. It is explicitly exploratory because the
image trait is not an independently calibrated genome-size measurement.

### 8. Canonical research-review notebooks from frozen outputs

The primary collaborator-review surface is split into eight analysis-domain
notebooks under `notebooks/research_review/`. Building and executing them reads
frozen corrected inputs and runs only the lightweight frozen image-IOD summary;
it does not rerun RepeatMasker, RepeatModeler, dnaPipeTE, read mapping, model
training, or cell/nucleus segmentation.

```bash
scripts/run_in_dusky.sh python scripts/processing/build_corrected_repeat_landscape.py
scripts/run_in_dusky.sh python scripts/processing/build_research_review_notebooks.py
scripts/run_in_dusky.sh Rscript scripts/processing/build_audited_historical_style_figures.R
scripts/run_in_dusky.sh Rscript path_analysis/scripts/run_cell_nucleus_genome_phylogenetic_path_analysis.R
scripts/run_in_dusky.sh python path_analysis/scripts/build_cell_nucleus_genome_path_presentation.py
uv run --with pandas --with numpy --with matplotlib --with nbformat --with nbclient --with nbconvert --with jupyter-client --with ipykernel python path_analysis/scripts/publish_cell_nucleus_genome_path_notebook.py
for notebook in notebooks/research_review/[0-9][0-9]_*.ipynb; do
  scripts/run_in_dusky.sh jupyter nbconvert --to notebook --execute "$notebook" --inplace --ExecutePreprocessor.timeout=600
done
scripts/run_in_dusky.sh python path_analysis/scripts/build_frozen_genome_iod_notebook.py --finalize
scripts/run_in_dusky.sh python scripts/processing/audit_research_review_notebooks.py
scripts/run_in_dusky.sh jupyter lab notebooks/research_review
```

The compact RepeatMasker landscape uses inclusive hit query-coordinate bp and
exactly conserves the corrected hit-level manifest. Its denominator is
within-species corrected aligned hit bp, not genome or assembly span.

## Directory Structure

After obtaining all data, your `input_data/` directory should look like:

```
input_data/
├── dnaPipeTE/
│   ├── SRX19952657_reads_per_component_and_annotation
│   ├── SRX19952691_reads_per_component_and_annotation
│   └── ... (34 files total)
├── repeatmasker/
│   ├── SRX19952657_Trinity.align
│   ├── SRX19952691_Trinity.align
│   └── ... (34 files total)
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
scripts/run_in_dusky.sh python verify_setup.py
```

This will check that all required files are in place.

## Data Availability

Raw sequencing data should be fetched from NCBI SRA using the accessions above. This repository does not currently assert one canonical BioProject ID.

Genome assemblies are available from NCBI GenBank under the accessions listed above.

## Contact

For questions about data acquisition, please open an issue on the repository.

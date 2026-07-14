# Fuscus genomic-resource benchmark

**Generated comparison:** `SRX20497025 / GCA_032353935.1` versus `GCA_050004315.1`
**Interpretation:** diagnostic method-and-resource sensitivity only; no value below authorizes silent substitution.

## What is established now

| Metric | Current repository resource | 2025 published resource | Comparability |
|---|---:|---:|---|
| Assembly span | 0.653 Gb | 16.118 Gb | Same unit, radically different completeness |
| Repeat-aligned load | 66.010% of configured 1.5-Gb quantification sample | about 75% of the assembly | Same percentage scale; different specimens and methods |
| TE-classified load | 62.536% of configured 1.5-Gb quantification sample | about 75% whole-assembly repeat coverage | Current value excludes Other/Unknown/unresolved mass |
| LTR | 64.696% of classified dnaPipeTE orders | about 36% of the assembly | Different method and denominator |
| LINE | 16.360% of classified dnaPipeTE orders | about 15% of the assembly | Different method and denominator |
| Gypsy / Ty3 | 57.965% of represented superfamilies | about 31% of repeat coverage | Different method and denominator |

The old assembly spans only 4.05% of the 2025 assembly. That makes the modern resource particularly important for assembly-dependent LTR and structural analyses. It does not invalidate the current shallow-read dnaPipeTE comparison, which should remain a separately named branch. The recovered dnaPipeTE denominator makes the 66.010% repeat-aligned value useful as a sensitivity check, not a silent replacement for the published assembly annotation.

The current paired-LTR baseline is 40 estimated pairs, 9 high-confidence pairs, and median K2P 0.070580.

## Public-data blocker

The [2025 paper](https://doi.org/10.1093/g3journal/jkaf157) states that its EarlGrey annotation is deposited, but the live [Zenodo record](https://zenodo.org/records/15255946) checked on 2026-07-09 exposes only a 4.76-GB compressed FASTA and a 1.30-MB AGP. The promised filtered-repeat GFF is not present. Downloading the FASTA alone would not reproduce the published repeat classification and would consume about 21 GB when compressed and decompressed copies are retained.

## Minimum asset request

Request the EarlGrey `*_summaryFiles/` directory plus:

- `aDesFus1-2.1.filteredRepeats.gff`;
- `*.highLevelCount.txt`;
- `*_divergence_summary_table.tsv`;
- `*_combined_library.fasta`;
- the raw LTR_FINDER GFF;
- the EarlGrey 4.4.2 command/log and checksums.

The filtered GFF is enough for TE bp, count, length, class/family, and repeat-to-consensus divergence comparisons without the full FASTA. The raw LTR_FINDER GFF is needed to determine whether the repository's paired 5-prime/3-prime LTR K2P metric can be reproduced.

## Decision

Keep `SRX20497025 / GCA_032353935.1` as the current comparative branch. Add `GCA_050004315.1` as an accession-labeled validation branch once the annotation assets are available. Do not merge either genomic specimen identity into the independently collected microscopy records.

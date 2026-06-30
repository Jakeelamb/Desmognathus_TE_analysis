# TE Provenance Audit

This document audits the upstream provenance chain for the repo-local TE summary
tables currently used in `path_analysis/`.

This is not a rerun log. It is a transparency map answering a narrower
question:

How much of the TE pipeline is still locally inspectable from raw or
near-raw staged inputs through the summary tables that feed the comparative
and phylogenetic path-analysis workflows?

## Bottom Line

The upstream TE provenance chain is mostly transparent locally.

Strong local provenance exists for:

- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`
- `results/data/repeatmasker_detailed_classification_combined.csv`
- `results/data/divergence/divergence_summary_statistics_by_species.csv`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

The main historical provenance weakness was the exact writer for:

- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`

That gap is now resolved non-destructively. See
`path_analysis/TE_DIVERSITY_CANONICALIZATION.md` for the audited reconstruction
path and scratch candidate outputs.

## Local Upstream Chain

### 1. dnaPipeTE classification and breakdowns

Primary script:

- `scripts/processing/dnaPipe.py`

Expected inputs:

- `input_data/dnaPipeTE/`
- `input_data/lookup_table.txt`

Observed local evidence:

- `input_data/dnaPipeTE/` exists with `35` staged entries
- `input_data/lookup_table.txt` exists with `34` lookup rows
- all lookup SRA accessions are represented in the staged dnaPipeTE inputs
- the extra count is explained by duplicate accession-prefix material rather
  than missing lookup coverage

Observed outputs:

- `results/data/dnaPipeTE_merged_classifications.csv`
- `results/data/dnaPipeTE_class_breakdown.csv`
- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`

Assessment:

- local provenance is strong
- the raw staged inputs, the processing script, and the canonical breakdown
  outputs are all present

### 2. RepeatMasker merge layer

Primary script:

- `scripts/processing/repeatmask.py`

Expected inputs:

- `input_data/repeatmasker/*.align`
- `results/data/dnaPipeTE_merged_classifications.csv`
- `input_data/lookup_table.txt`

Observed local evidence:

- `input_data/repeatmasker/` exists with `38` `.align` files
- all lookup SRA accessions are represented
- four extra `.align` files are present but not in the current lookup:
  - `SRX19952891_Trinity.align`
  - `SRX19953421_Trinity.align`
  - `SRX19953983_Trinity.align`
  - `SRX19958874_Trinity.align`
- the script explicitly excludes only the four documented out-of-scope SRX ids
  not found in the lookup table; any new unmapped `.align` file is a hard
  failure

Observed outputs:

- `results/data/merged_repeatmasker_data.csv`
- `results/data/repeatmasker_detailed_classification_combined.csv`

Assessment:

- local provenance is moderately strong
- the raw staged `.align` inputs and merged output are present, and the extra
  unmapped `.align` files are treated as documented out-of-scope inputs until
  explicitly reconciled with the lookup table

### 3. Divergence summary layer

Primary script:

- `scripts/processing/divergence.py`

Expected inputs:

- `results/data/repeatmasker_detailed_classification_combined.csv`

Observed local evidence:

- `results/data/repeatmasker_detailed_classification_combined.csv` exists
- divergence interim partitions exist:
  - `interim/divergence/class/` with `8` files
  - `interim/divergence/order/` with `9` files
  - `interim/divergence/superfamily/` with `22` files

Observed outputs:

- `results/data/divergence/divergence_summary_statistics_by_species.csv`

Assessment:

- local provenance is strong
- both the merged RepeatMasker source table and the divergence interims are
  still on disk, so the summary table is locally inspectable without HPC reruns

### 4. Ectopic recombination layer

Primary script:

- `scripts/processing/ec.py`

Expected inputs:

- `input_data/ectopic_recombination/GCA_*_tabout.csv`
- `input_data/ectopic_recombination/*.fa.depth.txt`
- `input_data/ectopic_recombination/combined_sequences.fasta.rexdb-metazoa.cls.tsv`
- `input_data/lookup_table.txt`

Observed local evidence:

- `input_data/ectopic_recombination/` exists
- `35` `GCA_*_tabout.csv` files are present
- `5190` `*.fa.depth.txt` files are present
- `combined_sequences.fasta.rexdb-metazoa.cls.tsv` is present

Lookup reconciliation:

- four raw ectopic `GCA_*_tabout.csv` files are not represented in the current
  lookup table:
  - `GCA_030264435.1_tabout.csv`
  - `GCA_030264935.1_tabout.csv`
  - `GCA_030265015.1_tabout.csv`
  - `GCA_030265095.1_tabout.csv`
- three lookup genomes currently have no matching raw ectopic `GCA_*_tabout.csv`
  file:
  - `D.catahoula` / `GCA_034783935.1`
  - `D.kanawha` / `GCA_032353855.1`
  - `D.valtos` / `GCA_032357565.1`
- `D.lycos` does have a raw ectopic `GCA_*_tabout.csv` file, so its absence from
  the filtered ectopic summary reflects filtering/retention, not missing raw
  input

Observed outputs:

- `results/data/ectopic_recombination_master.csv`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

Assessment:

- local provenance is moderately strong
- the raw staged inputs and the canonical filtered output are present
- the main provenance caveat is accession reconciliation, not missing code

### 5. Diversity summary layer

Relevant scripts found locally:

- `scripts/processing/diversity_stats.py`

Observed outputs:

- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`
- `results/data/long_format_diversity_order_stats.csv`
- `results/data/long_format_diversity_superfamily_stats.csv`
- `results/data/comparison_diversity_order_stats_granular_0_5pct.csv`
- `results/data/comparison_diversity_superfamily_stats_granular_0_5pct.csv`

Resolved generation path:

- `scripts/processing/diversity_stats.py` writes the comparison-threshold files:
  - `results/data/comparison_diversity_order_stats_granular_0_5pct.csv`
  - `results/data/comparison_diversity_superfamily_stats_granular_0_5pct.csv`
- the canonical summary tables are the corresponding breakdown tables with the
  `0.0` threshold metrics merged back in as:
  - `Simpson_Diversity <- Simpson_0.0`
  - `Shannon_Diversity <- Shannon_0.0`
  - `Pielou_Evenness <- Pielou_0.0`
- `path_analysis/scripts/build_canonical_diversity_tables.py` now reconstructs
  scratch candidates for both canonical files and audits them against the
  current `results/data/` snapshots
- older standalone diversity helpers were removed after canonicalization
  cleanup because they implemented raw `1 - sum(p_i^2)` Simpson rather than the
  corrected definition used by the stored diversity summaries

Assessment:

- local provenance is now strong enough for research use
- the exact construction path is documented and reproducible from local files
  without rerunning upstream HPC work
- the remaining cleanup question is organizational rather than methodological:
  whether to later promote the audited candidate writer into `results/data/`
  generation and retire older ambiguous helper scripts

## Transparency Status By TE Input Used In Path Analysis

### High confidence local provenance

- `dnaPipeTE_order_breakdown.csv`
- `dnaPipeTE_superfamily_breakdown.csv`
- `divergence_summary_statistics_by_species.csv`
- `ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

### Moderate confidence local provenance

- `repeatmasker_detailed_classification_combined.csv`
  because extra unmapped `.align` files are present in the raw input directory

## Implications For Analysis Use

The current TE source tables are usable for downstream comparative and
phylogenetic path analysis, but the methods notes should reflect the true
confidence level of each upstream step.

Practical interpretation:

- the order-breakdown, diversity, divergence, and ectopic filtered files now
  all have a defensible local provenance chain
- the remaining TE-methods cleanup is optional housekeeping rather than a
  blocker for current analysis use

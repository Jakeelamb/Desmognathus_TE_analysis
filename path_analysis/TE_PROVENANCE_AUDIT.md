# TE Provenance Audit

This document audits the upstream provenance chain for the repo-local TE summary
tables currently used in `path_analysis/`.

This is not a rerun log. It is a transparency map answering a narrower
question:

How much of the TE pipeline is still locally inspectable from raw or
near-raw staged inputs through the summary tables that feed the comparative
and phylogenetic path-analysis workflows?

## Bottom Line

The upstream TE provenance chain is mostly transparent locally, but the
current RepeatMasker-derived divergence layer is a preserved pre-audit branch,
not a publication-ready result.

Strong local file provenance exists for:

- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`
- `results/data/repeatmasker_detailed_classification_combined.csv` (pre-audit classification state)
- `results/data/divergence/divergence_summary_statistics_by_species.csv` (derived from that pre-audit state)
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
- the extra count is `SRX19952890R2`, a computational retry label for the same
  *D. orestes* accession rather than a second specimen, SRA run, or missing
  lookup row
- both *D. orestes* outputs are full, similarly sized dnaPipeTE realizations;
  they must not be summed as independent evidence

Observed outputs:

- `results/data/dnaPipeTE_merged_classifications.csv`
- `results/data/dnaPipeTE_class_breakdown.csv`
- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`

Assessment:

- local provenance is strong
- the raw staged inputs, the processing script, and the canonical breakdown
  outputs are all present
- the historical breakdown percentages are closed relative compositions, not
  absolute genome-wide TE fractions
- a corrected, analysis-scoped mass ledger now exists at
  `results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv`
  with an adjacent manifest and mass-augmented order table
- v2 independently matches every observed SRX to the active lookup and records
  both `te_sra_accession` and `te_assembly_accession`; the earlier v1 audit
  snapshot is preserved but superseded for analysis
- that ledger contains exactly the declared 18-species TE/genome panel, retains
  total aligned-repeat bases plus unresolved mass at class, order, and
  superfamily levels, and explicitly excludes *D. orestes*
- order-level unresolved mass ranges from 3.548% to 6.179% (median 5.171%);
  superfamily-level unresolved mass ranges from 17.507% to 21.586% (median
  19.995%)
- recomputed retained-order percentages reproduce the historical table within
  7.105e-15 percentage points, showing that the new ledger changes denominator
  transparency rather than silently changing the retained composition
- Git history recovers the upstream run parameters from `dnaPipeTE.sh` at
  commit `cc64ffb134968250671ff0e3456255063a1aba97`: 15-Gb configured genome,
  0.1× coverage, two samples, `RM_t=0.15`, nuclear-filtered R1 input, and a
  custom `dedupe_telib.fasta`
- the configured quantification denominator is therefore 1.5 Gb, allowing a
  sensitivity-only repeat-aligned fraction at
  `results/data/corrected/dnapipete/dnapipete_absolute_load_sensitivity_analysis18_v1.csv`
- final-panel repeat-aligned fractions range from 60.619% to 68.196% (median
  65.744%); *D. fuscus* is 66.010% repeat-aligned and 62.536% assigned to a DNA-
  or retrotransposon class
- this branch is not confirmatory because runtime logs, container version or
  digest, custom-library checksum, replicate quantification samples, and
  sampling uncertainty are absent
- the corrected writer refuses more than one computational output per SRX in
  the final panel, preventing retry suffixes from being silently double-counted

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
- `results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv`
  for the declared 18-species final TE/genome panel

Assessment:

- local file provenance is strong, but the stored canonical combined table is
  scientifically pre-audit
- its generic `Class`, `Order`, and `Superfamily` fields inherited one
  dnaPipeTE annotation per contig rather than each hit's native
  `repeat_class`
- the corrected code retains dnaPipeTE fields under `dnapipete_*`, classifies
  each hit under `repeatmasker_*`, and keeps generic names only as
  backward-compatible aliases of the hit-level result
- the four extra `.align` files are documented out-of-scope inputs and are not
  analyzed merely because they are staged locally
- the corrected analysis-eligible table contains exactly 18 SRX IDs and 18
  species; its manifest records 5,132,397 hits, 1,155,346,511 inclusive aligned
  bp, readback alias equality, and SHA-256
- an all-34-resource rebuild exists only as an audit/conservation proof and is
  explicitly marked `eligible_for_path_analysis: false`
- full mismatch counts and unmapped labels are preserved under
  `plans/publication-readiness-deep-audit/`

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

- local provenance is strong enough to reconstruct the historical output
- the current divergence interims and species summary inherit the contig-level
  classification error and must remain labeled pre-audit
- regeneration must use the corrected hit-level table and a separate output
  location before any corrected result replaces the historical branch
- the corrected 18-species summary now exists at
  `results/data/corrected/divergence/divergence_summary_statistics_by_species_analysis18_v1.csv`
- its all-hit branch conserves 5,132,397 hits at each classification level; its
  legacy-comparable threshold branch conserves 4,715,646 hits with dnaPipeTE
  contig-threshold context
- 416,751 hits (8.12%; 7.86% of aligned bp) lack that context and are retained
  in the all-hit branch plus a per-species coverage sidecar rather than being
  silently dropped
- the corrected weighted species predictors are rank-stable relative to the
  historical values, but RepeatMasker percent deletions/insertions remain
  alignment-gap statistics relative to repeat consensus, not direct DNA-loss
  or ectopic-recombination rates

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
- `results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv`
- `results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv`

Assessment:

- file provenance is moderately strong, but mapping-method provenance is not
  sufficient for a mechanistic ectopic-recombination claim
- the corrected branch processes only the final panel, joins TEsorter by exact
  `sequence_start_end`, and retains zero-depth positions in regional means
- 567 exact 5/6-domain elements are available across 16 panel species;
  *D. kanawha* and *D. valtos* lack the required assembly `tabout` resources
- 565/567 elements pass 80% positive-position coverage in both LTRs and the
  internal region; low regional coverage is therefore not the dominant issue
- 97 element ratios change when 2,872 explicit zero-depth positions are
  retained; the maximum absolute ratio change is 1.706581
- the historical nonzero-only statistic is reproduced within 4.441e-16
- species arithmetic means are not robust: for *D. intermedius*, one element
  accounts for 75.0% of the ratio sum; the zero-aware arithmetic mean is 3.114,
  median 0.646, and geometric mean 0.800
- mapping command/reference, multimapper handling, MAPQ, secondary and
  supplementary alignments, duplicates, and direct solo:intact-LTR validation
  remain unavailable; this layer is not approved as an ectopic-recombination
  rate or confirmatory path predictor

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
- both production and audit writers now implement the scale-invariant
  Gini-Simpson definition `1 - sum(p_i^2)`; the stored values are not a
  finite-count-corrected Simpson estimator

Corrected final-panel audit products:

- `results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_composition_matrices_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_pca_clr_matrices_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_pca_scores_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_pca_loadings_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_pca_variance_analysis18_v1.csv`
- `results/data/corrected/diversity_pca/te_pca_stability_analysis18_v1.csv`
- `plans/publication-readiness-deep-audit/te_diversity_pca_corrected_analysis18_v1.md`

Assessment:

- the historical order and superfamily diversity values reproduce within
  `3.331e-16` and `4.441e-16`, respectively
- the corrected table names each estimand explicitly: Shannon entropy (natural
  log), Gini-Simpson, Simpson dominance, Hill q1/q2, Pielou evenness, and
  observed richness
- classified-conditional versus unresolved-bin ranks remain high (Spearman
  rho 0.965-0.975 at order level and 0.930-0.986 at superfamily level), but the
  two denominators remain separate because unresolved mass is a technical bin
- order-level CLR PCA is approved as a descriptive ordination: it has 18
  species, 10 nonzero features, PC1/PC2 explain 56.7%/28.9%, and the minimum
  leave-one-species-out score correlations are 0.990/0.953
- superfamily CLR PCA remains supplementary because 27 features and 32 zeros
  make the geometry materially dependent on zero replacement; PC2 minimum
  leave-one-out score correlation falls to 0.498 under feature-specific
  replacement
- the historical unrestricted clade PERMANOVA is not approved for inference;
  related species are not freely exchangeable permutation units
- fixed-tree phylogenetic PCA uses a full-rank 9-coordinate ILR projection of
  the order CLR matrix; pPC1/pPC2 explain 53.5%/23.6%, and ordinary-versus-pPCA
  species-score correlations are 0.983/0.913
- across all 200 Stewart-Wiens time-calibrated bootstrap trees, no axis swap is
  required; minimum ordinary-versus-pPCA score correlations are 0.988/0.910
  and loading correlations are 0.955/0.723
- order PCA is therefore approved as a descriptive ordination across published
  tree-time uncertainty, but not as a causal path variable

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
- `ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

### Preserved pre-audit provenance

- `repeatmasker_detailed_classification_combined.csv`
- `divergence_summary_statistics_by_species.csv`

The issue is not the presence of out-of-scope raw files; it is the confirmed
contig-versus-hit classification error in the historical merge.

## Implications For Analysis Use

The dnaPipeTE composition and ectopic source tables remain available for
downstream audit. RepeatMasker-derived divergence and alignment-gap summaries
must not enter a corrected comparative or phylogenetic analysis until rebuilt
from the hit-level table.

Practical interpretation:

- the dnaPipeTE order-breakdown, diversity, and ectopic filtered files retain
  locally inspectable provenance
- RepeatMasker divergence is a correctness blocker, not optional housekeeping
- current divergence-dependent path outputs should be retained as historical
  sensitivity artifacts and regenerated non-destructively

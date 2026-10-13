# Wang et al. (2025) Figure 2 TE-diversity source audit

**Research conclusion (checked 2026-08-12):** use the 84-row Figure 2 workbook deposited in Dryad as the authoritative published comparison source. Do **not** describe either `/home/jake/Projects/Desmognathus_Results/wang_ji_shannondiversity.csv` or the archive's cleaned 87-row derivative as the published Figure 2 dataset. Those historical files contain rounded values, omit *Desmognathus ochrophaeus*, add four caecilians that are absent from Figure 2, and replace the deposited duck name *Anas platyrhynchos* with `Anas anas`.

The authoritative primary sources are the [published article](https://doi.org/10.1038/s42003-025-08127-3) and its [Dryad dataset](https://doi.org/10.5061/dryad.zpc866tkv). Dryad states that `Fig_2.zip` contains the spreadsheet and R code used for Figure 2, and licenses the dataset under [CC0 1.0](https://spdx.org/licenses/CC0-1.0.html).

## Authoritative deposited artifacts

The public [Dryad v2 dataset API](https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.zpc866tkv) identifies the published dataset as version 4, internal version ID `360042`, published 2025-04-22. The [version file listing](https://datadryad.org/api/v2/versions/360042/files?per_page=100) identifies Figure 2 as file ID `4030762`.

| Artifact | Size | SHA-256 |
|---|---:|---|
| Dryad `Fig_2.zip` | 305,298 bytes | `aff6c3fdcd566a5ec1622ce3b256a38e9c337520a812ef75005a8fd896671dc9` |
| `Fig_2/Data_TE_diversity.xlsx` | 21,129 bytes | `367a8102fd18f9c17251c442b0a0112af6afc008b51bf4b110f72a98fcc05c70` |
| `Fig_2/Fig2_code.R` | 6,593 bytes | `a8b649514a942dd9054a7b93d91a810b94468305118515dbea30eb86ef638985` |
| `Fig_2/Fig2_April_17ps_final.pdf` | 287,109 bytes | `d40595bf1b6c9c2ce162c34bbe9f21f07a874ac466fe3735d1e1a22295a3e481` |
| `Fig_2/Fig2_R_output_April_17.pdf` | 31,695 bytes | `dc996e9f1cb0f2303941a2d6933c31f5527a2d6c8b5b155302028c1e2c39d355` |

The workbook has one sheet, 84 data rows, and 11 populated columns:

`English_name`, `species`, `Family`, `GS`, `SI`, `GSI`, `type`, `type1`, `type2`, `GSlog2`, and `uncommon_fish_type`.

`GS` is genome size in gigabases, `SI` is Shannon diversity, and `GSI` is Gini-Simpson diversity. `GSlog2` is a spreadsheet formula-derived base-2 transformation of `GS`; the deposited R code instead maps `GS` and applies a base-2 transformed x-axis. `uncommon_fish_type` stores plotting shapes for the six non-actinopterygian fishes.

## Exact Figure 2 sample

The workbook and the article's Methods agree on 84 species:

| Deposited Figure 2 group (`type1`) | Rows | Article description |
|---|---:|---|
| Fishes | 19 | 13 ray-finned, 3 cartilaginous, and 3 lobe-finned fishes |
| Anura | 10 | frogs |
| Caudata | 11 | salamanders |
| Non-avian Reptillia | 16 | 8 snakes, 3 lizards, 2 crocodilians, 2 turtles, and the tuatara |
| Aves | 14 | birds |
| Mammalia | 14 | mammals |
| **Total** | **84** | |

The spelling `Non-avian Reptillia` is present in the deposited workbook and code; it should be retained in a raw source snapshot but can be normalized to `Non-avian reptiles` in a presentation column.

The deposited R code reads workbook rows 1–13 and 20–84 into the primary plotting table, then reads rows 14–19 separately as the six non-actinopterygian fishes. All 84 rows are plotted. The fitted line in the fish facet uses only the 13 actinopterygians, while the other six fishes are displayed with separate symbols.

The eleven deposited salamanders are:

- *Desmognathus ochrophaeus*
- *Ranodon sibiricus*
- *Tylototriton verrucosus*
- *Batrachoseps nigriventris*
- *Ambystoma mexicanum*
- *Pachytriton brevipes*
- *Cynops orientalis*
- *Aneides flavipunctatus*
- *Andrias davidianus*
- *Paramesotriton hongkongensis*
- *Cryptobranchus alleganiensis*

## Metric contract and reproducibility boundary

The article defines the indices as:

- Shannon: `H' = -sum(P_i ln(P_i))`
- Gini-Simpson: `D = 1 - sum(P_i^2)`

Here, `P_i` is the proportion of sequences assigned to TE superfamily `i`. TE superfamilies are treated as community “species,” and only the ten most abundant superfamilies are included because rare-superfamily abundance was considered uncertain. `ln` establishes the natural-log Shannon convention. These definitions appear in the article's “Comparison of genomic TE community diversity across genome sizes” Methods section; the deposited workbook contains the final indices rather than the underlying ten abundances for all 84 species.

The article does not state explicitly whether the ten retained proportions were reclosed to sum to one after rarer superfamilies were excluded. Therefore, a new calculation for this study must document its denominator and any top-ten reclosure rather than claiming exact computational identity based only on the paper's prose.

The article describes Figure 2 fits and statistics as phylogenetic generalized least-squares (PGLS) results. The deposited `Fig2_code.R`, however, renders each smooth with `stat_smooth(method = lm)` and manually places the reported `F` and `P` labels. The Figure 2 deposit includes no phylogenetic trees or PGLS-fitting code. Consequently, the archive reproduces the plotted data and approximate rendered lines, but it does not independently reproduce the reported PGLS analyses. This project should use the deposited rows for descriptive cross-study context only unless the phylogenies and full PGLS workflow are separately reconstructed and validated.

## Reconciliation with the local 87-row file

| Local artifact | Rows | Columns | SHA-256 |
|---|---:|---:|---|
| `/home/jake/Projects/Desmognathus_Results/wang_ji_shannondiversity.csv` | 87 | 7 | `4ec3237a6ecd361c7e040b11dccfc9e2ac635cf10b739a4e95920015edccbb0e` |
| `/home/jake/Projects/Desmognathus_TE_archive/results/data/research_review/wang_ji_2025_shannon_context_v1.csv` | 87 | 9 | `6e40ecc7a18aacfdffc394cdb9ca6a11eba1adf5892b13382cd6d8dc93cd71eb` |
| Archive provenance JSON | n/a | n/a | `faff44098695d60e1141ca584c94401bbbe40c9741e95f614f5435fa83c9507f` |

The local raw columns are `English_name`, `species`, `Family`, `GS`, `SI`, `GSI`, and `Type` after trimming trailing spaces from headers. Its group counts are 13 Bony Fishes, 3 Lobe-finned fish, 3 Cartilaginous fish, 4 Caecilian, 10 Frogs, 10 Salamanders, 16 Non-avian reptiles, 14 Birds, and 14 Mammals.

After Unicode whitespace normalization, 82 species names overlap between the local and deposited tables. Reconciliation is exact at the sample-definition level:

- deposited but absent locally: *Desmognathus ochrophaeus* and *Anas platyrhynchos*;
- local but absent from the deposit: *Geotrypetes seraphini*, *Microcaecilia unicolor*, *Rhinatrema bivittatum*, *Ichthyophis bannanicus*, and `Anas anas`;
- the `Anas anas` row has the same genome size and rounded diversity values as deposited *A. platyrhynchos*, so it is a replacement/mislabel rather than an additional biological record;
- the net row difference is therefore the removal of one salamander and addition of four caecilians: `84 - 1 + 4 = 87`.

For all 82 name-matched rows, local `SI` and `GSI` values equal the deposited values rounded to two decimal places, within at most 0.005. The local `GS` value is identical for 79 rows and rounded to two decimals for the remaining three. The local table also discards the deposited `type`, `type1`, `type2`, `GSlog2`, and shape metadata and substitutes a nine-level `Type` grouping.

The archive's existing provenance JSON correctly discloses that its cleaned 87-row table is not an exact export of the paper's 84-species analysis. Its claim that “published values [are] unchanged” is only true relative to the already rounded local CSV, not relative to the full-precision Dryad workbook.

## Acceptance decision

**Accepted:** the exact 84-row `Data_TE_diversity.xlsx` values from Dryad, normalized losslessly into a source snapshot with original row numbers, original full-precision numeric values, dataset DOI, article DOI, license, archive/member hashes, and retrieval date.

**Rejected as publication source:** the 87-row historical CSV and its cleaned archive derivative. They may remain only as historical provenance and should not feed current figures, workbooks, tests, or reported sample counts.

**Permitted use:** descriptive cross-study panels that clearly distinguish Wang et al. values from this study's newly calculated *Desmognathus* values. The comparison should retain the paper's “top ten TE superfamilies” scope in axis/legend text and disclose that the studies used different genomic inputs and annotation pipelines.

**Not accepted without additional work:** pooled regression, cross-study hypothesis tests, reuse of the paper's PGLS statistics for a combined dataset, or a claim that a newly calculated *Desmognathus* index is method-identical when its top-ten denominator/reclosure rule has not been made explicit.

## Minimum provenance for a combined release table

Each row in a combined contextual table should retain at least:

- `record_id` and `source_row`;
- `data_origin` (`wang_2025_dryad` or `this_study`);
- `source_dataset_doi` and `source_article_doi`;
- `source_file`, `source_archive_sha256`, and `source_member_sha256`;
- `license`;
- `scientific_name`, `common_name`, and `family`;
- `comparison_group`;
- `genome_size_gb` and the genome-size provenance/scope;
- `shannon_top10_superfamilies` and `gini_simpson_top10_superfamilies`;
- `metric_taxonomic_unit`, `metric_scope`, `abundance_basis`, `log_base`, and `top10_reclosure_rule`;
- `figure_inclusion` and `notes`.

This keeps published observations verbatim, makes study-specific transformations visible, and prevents the descriptive combined table from being mistaken for a single method-homogeneous analysis.

# RepeatMasker corrected-feature sensitivity

This comparison isolates the effect of replacing inherited dnaPipeTE contig classes with each RepeatMasker hit's native class. Both branches retain the historical 0.9 dnaPipeTE contig-confidence threshold and the same dnaPipeTE order weights.

| Feature | Spearman rho | Median absolute change | Maximum absolute change |
|---|---:|---:|---:|
| `weighted_te_deletions_p90` | 1.0000 | 0.001074 | 0.005648 |
| `weighted_te_divergence_p90` | 0.9979 | 0.006316 | 0.022584 |
| `weighted_te_insertions_p90` | 1.0000 | 0.001166 | 0.004953 |

The weighted species predictors are numerically stable to the classification correction, even though many individual hits were assigned to the wrong order/superfamily in the historical detailed table. This happens because the downstream predictors use order-weighted medians, which dampen many row-level reassignments.

This stability does not convert `percent_deletions` into a DNA-loss rate. It remains a RepeatMasker alignment-gap statistic relative to a repeat consensus, and the corrected divergence branch remains a sensitivity layer rather than direct evidence of genomic deletion.

## Fuscus values

| Feature | Historical | Corrected | Change |
|---|---:|---:|---:|
| `weighted_te_deletions_p90` | 0.322902 | 0.328229 | +0.005328 |
| `weighted_te_divergence_p90` | 12.233229 | 12.246841 | +0.013612 |
| `weighted_te_insertions_p90` | 0.191154 | 0.192436 | +0.001282 |

## Inputs

- Historical: `results/data/divergence/divergence_summary_statistics_by_species.csv` (`d22265d56becc5f50d39271a099ec3a0f3d10cb4206faddb8569bd8acfd159a4`)
- Corrected: `results/data/corrected/divergence/divergence_summary_statistics_by_species_analysis18_v1.csv` (`ac4ffd2f2f5dbb5f4601fa2cb85d6bf0f86505dd88e41d6535130cf5ba731504`)
- Species panel: `path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv` (`9682bb74582570edae0f6b9b8672fdb9b76e5775a2d7ce45e622114151e399b8`)

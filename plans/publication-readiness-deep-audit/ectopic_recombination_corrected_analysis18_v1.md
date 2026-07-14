# Corrected ectopic-recombination depth audit: final panel

This branch rebuilds element identity and depth summaries for only the declared 18-species TE/genome panel. It preserves historical outputs unchanged.

## What was corrected

- TEsorter annotations are joined by exact `sequence_start_end`, not contig name alone.
- Zero-depth positions remain in terminal and internal regional means.
- Per-element left-LTR, right-LTR, terminal, and internal positive-position coverage is explicit.
- Species summaries include arithmetic mean, median, geometric mean, trimmed/winsorized means, bootstrap intervals, and leave-one-element-out influence.
- 5/6-domain elements are reported under all-element, >=80% coverage, TEsorter-complete, and combined sensitivity branches.

## Audit result

- Exact 5/6-domain elements: 567 across 16 of the 18 panel species.
- Species with no retained element-level estimate: 2.
- Missing assembly-level `tabout` resources: kanawha, valtos.
- Elements passing >=80% positive-position coverage in both LTRs and the internal region: 565/567 (99.6%).
- Maximum absolute change caused by retaining zeros instead of deleting them: 1.706581 in the element ratio.
- Maximum nonzero-only reproduction error against the historical stored ratio: 4.441e-16.
- The historical contig-only TEsorter join created 1 extra row and 1 coordinate-mismatched annotation; that row has fewer than five domains and did not enter the current 5/6-domain predictor.
- Strongest arithmetic-mean influence remains *D. intermedius*: removing one element can shift its species mean by 2.302.

## Interpretation boundary

Terminal:internal read depth is a screening proxy for excess LTR-like sequence relative to internal sequence. It is not yet a measured ectopic-recombination rate or a direct solo-LTR count. Mapping command, reference construction, multimapper handling, MAPQ, secondary/supplementary alignment policy, duplicate policy, and independent validation against structurally called solo/intact LTRs are absent. No corrected ectopic variable is approved for confirmatory path analysis until those are resolved.

## Validation figures

- `results/figures/corrected/ectopic/ectopic_element_log2_ratio_by_species_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_mean_vs_median_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_coverage_qc_analysis18_v1.png`
- `results/figures/corrected/ectopic/ectopic_mean_influence_analysis18_v1.png`

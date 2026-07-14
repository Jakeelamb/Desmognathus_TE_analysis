# RepeatMasker hit-classification audit

The current merged table contains 9,533,364 RepeatMasker hits. Its generic class/order/superfamily columns were inherited from one dnaPipeTE annotation per contig; this audit reclassified each row from its own native `repeat_class` without rewriting the stored table.

- Known-order disagreement: 18.148% of comparable hits and 11.426% of comparable aligned bp.
- Known-superfamily disagreement: 20.381% of comparable hits and 12.762% of comparable aligned bp.
- Native RepeatMasker labels not covered by the shared map: 137 hits (0.0014%).
- Aligned length uses inclusive RepeatMasker coordinates: `abs(query_end - query_start) + 1`.

The corrected merge code retains dnaPipeTE classifications under `dnapipete_*`, writes hit-derived fields under `repeatmasker_*`, and keeps generic `Class`, `Order`, and `Superfamily` only as backward-compatible aliases of the hit-derived result. Existing divergence outputs remain pre-audit until regenerated from that corrected table.

Unmapped-label detail is stored in `repeatmasker_unmapped_labels.csv` (19 labels).

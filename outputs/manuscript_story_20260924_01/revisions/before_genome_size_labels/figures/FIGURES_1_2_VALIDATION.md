# Figures 1 and 2: isolated candidate validation

Build: `scripts/figures_1_2.R`, with working directory equal to this isolated run.

Executable used: `/home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ`.
The initially suggested `desmo-gbe-figures` environment is not installed.

## Numerical scope

- Figure 1 retains all 34 TE species, with the canonical supplied-tree row order.
- Reproductive-strategy counts are 7 aquatic eggs, 24 terrestrial eggs with aquatic larvae, and 3 direct developers.
- Exactly 21 existing genome-size estimates are displayed, with 13 explicit `n/a` rows. Numeric labels are rounded to one decimal; exported analysis data retain the source values and intervals.
- Four displayed composition categories close to one for every species. `Other` sums the remaining classified orders; no unresolved mass is recoded into this category.
- Figure 2 exports 5,100 species/order/bin rows, 450 strategy/order/bin means, and 102 low-divergence fractions.
- Profile data retain the all-species aligned-hit-bp denominator, zero-fill missing bins 0–49, and do not reclose within orders. Group means weight each species equally.
- All observed values fit the fixed 0–3% profile and 0–6% low-divergence axes. No observations were clipped or removed.
- No inferential analysis was performed.

## Output checks

- Two collages and five isolated panels were exported as PDF, PNG and TIFF.
- Collage PDFs embed Arial regular, bold and italic fonts.
- All seven TIFF exports were confirmed CMYK at 300 pixels per inch.
- Both collages were visually inspected. Profile facet spacing was widened to separate boundary tick labels; the same spacing is used in both rows. The composition panel right margin retains the full 100% tick label.
- Isolated genome-size and composition panels include their own species labels. Isolated divergence and low-divergence panels include panel-specific keys and descriptions.
- Source and script checksums are recorded separately in each figure's `analysis/input_checksums.sha256`.
- The copied order-composition source matches the original SHA-256: `c7b00470d8e26b1ed6063a9d68b3c8403618df3298d485346c16818fb9c1b965`.

All writes for these figures were confined to this isolated run; original figures, data, review notes and shared `scripts/common.R` were not modified.

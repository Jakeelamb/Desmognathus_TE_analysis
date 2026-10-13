# Figures 5 and 7 visual checks

Both final collages were opened as rendered PNGs with `view_image`; isolated panels 5c, 5d and 7a were also inspected at larger size. The complete figure objects are exported as PDF, PNG and CMYK TIFF, not screenshots.

- Figure 5: four panels at 180 × 240 mm. All 24 species remain visible, paired fractions retain all 42 animals, the 24 fraction labels are readable at the explicit 7.5 pt size, and all uncertainty intervals are shown. A duplicated legend was removed. The fraction-panel order was aligned with the other figures, with D. anicetus at top and D. wrighti at bottom. Borderless white label backgrounds prevent fitted lines from crossing species text.
- Figure 7: four panels at 180 × 235 mm. Three SVL scatterplots use consistent horizontal breaks and log units. The 12-coefficient forest has separate equation strips, readable row labels and a visible zero reference. Species labels have borderless white backgrounds. The solid/dashed line legend distinguishes the full 24-species and 22-species sensitivity fits.
- Isolated scatterplots: all 24 species are labelled; the exported PDF text is checked programmatically for every species. Isolated 5d includes the full animal/cell support column and its own strategy color legend.
- Font inspection with `pdffonts` identified an initial fallback caused by decorative Unicode circle/diamond characters. Those characters were replaced with words. Final PDF font subsets are exclusively Arial regular, italic and bold, all embedded.
- `identify` confirms CMYK TIFFs at 300 dpi in the figure-render stage; the export check accepts later publication-resolution upgrades of 300 dpi or more.

Scientific content was checked independently of appearance: the verification script compares copied model and prediction tables exactly with their frozen sources, including all 12 conditional coefficients and intervals. There are no model-refitting or resampling calls in the plotting script. Complete statistics remain available in each figure's analysis directory.

Commands:

```bash
env PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH /home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ scripts/figures_5_7.R
env PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH /home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ qa/verify_figures_5_7.R
```

Run from `/home/jake/Projects/Desmognathus_TE/outputs/manuscript_polished_20260924_01`. Numerical, label, font and TIFF checks are saved in `figures_5_7_export_verification.txt`; the render log is `figures_5_7_render.log`.

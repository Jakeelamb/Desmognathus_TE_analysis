# Figures 3 and 5 isolated-candidate review

- Render command completed successfully with no warnings using `/home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ scripts/figures_3_5.R` from this run directory.
- Inspected both final overview PNGs and all five isolated-panel PNGs.
- Figure 3 contains 21 species in each panel and one full-sample fit per metric. The overview labels six reference species; the isolated exports label all 21. This is label selection only, not data selection.
- Figure 5 contains all 24 species per panel. The main genome-cell panel and each isolated export label all 24; the two smaller overview panels label five reference species. All measurement intervals are retained.
- Coefficient annotations occupy clear upper margins. Label tethers and multiple measurement intervals remain dense in the fully labelled reference exports; the overview limits labels where panels are smaller.
- All seven PDFs embed Arial regular/italic/bold fonts. Text in this script has a minimum size of 7 pt.
- Both overview TIFFs and all five isolated-panel TIFFs are CMYK at 300 pixels per inch.
- Figure 3 numerical checks match the existing full-sample F/P annotations after rounding. New coefficient and prediction CSVs are derived from that unchanged model specification.
- Figure 5 uses the existing S26 coefficients, covariances and t-based intervals without refitting, with predictor/outcome moments verified against the plotted 24-species data.
- Original figure candidates, original code, shared common.R and saved review notes were not edited.

These are publication-format review candidates; the biological interpretation boundaries remain in the legends.

Final consistency pass removed the genome-size axis asterisk, matching the other isolated manuscript figures: `Estimated genome size (pg/1C)`. All seven PDFs were regenerated and their extracted axis text checked. SHA256 checks confirm that all seven analysis tables/tree files are byte-identical to the pre-wording outputs (see `figures_3_5_before_axis_wording.sha256`).

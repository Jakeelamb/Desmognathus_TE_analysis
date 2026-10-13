# Figure 6 review

The 180 × 270 mm main collage, all four isolated panels, and the separate simulation calibration PNG were inspected after the final render. Text, arrows, model labels, uncertainty intervals, legends, and captions are legible and unclipped. The three DAG orientations preserve the source enumeration and have no numerical edge labels. The coefficient panels share the same horizontal scale, and the 18-species tree check explicitly states shared support for the nucleus-bridge class.

The final render completed with exit code 0 and no R warnings. The numerical assertions are recorded in `figure_6_numerical_checks.txt`; input checksums are stored with the figure's numerical exports.

Verification commands:

```bash
PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH /home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ scripts/figure_6.R
find figures/figure_6 -name '*.pdf' -exec pdffonts '{}' \;
magick identify -format '%f %w x %h %[colorspace] %x x %y %U\n' figures/figure_6/*.tif figures/figure_6/panels/*.tif figures/figure_6/validation/*.tif
```

Results: six PDFs contain only embedded Arial regular/bold fonts. All six TIFFs are CMYK, 300 × 300 pixels per inch. The collage is 2125 × 3188 pixels; all isolated exports and the support plot are 2125 pixels wide. No source data, original figure, existing review note, or original analysis was edited or refitted.

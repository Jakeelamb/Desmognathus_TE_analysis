# GBE figure release

The project house style follows the documented GBE/OUP requirements in
`../journal_requirements/GBE_REQUIREMENTS.md`:

- Arial or Helvetica-family sans serif text;
- labels at least 7 pt at final size (GBE absolute minimum: 6 pt);
- vector PDF master with embedded fonts;
- 300-dpi color TIF derivative when a raster revision file is required;
- accessible colors plus non-color encodings;
- no figure-internal title unless scientifically necessary;
- one legend and one alt-text record for every finalized figure.

GBE Articles allow eight combined main-text display items (figures plus tables),
not eight of each. The curated placement and scientific role of every current
figure are in `GBE_FIGURE_PLAN.md`; diagnostic and audit figures belong in the
supplementary-information PDF or DOI repository.

## Rebuild

```bash
make setup-r  # first run only
make figures
```

`FIGURE_MANIFEST.csv` maps each figure to its source supplemental datasets,
release status, dimensions, font, legend, alt text, and PDF/PNG/TIFF files.
The current host does not include Arial itself, so the reproducible build uses
Nimbus Sans, a Helvetica-compatible sans-serif font, and embeds it in every PDF.
`Path_DAG_equivalence_class_reference.{pdf,png,tif}` is a separate review
reference asset: it is intentionally outside `FIGURE_MANIFEST.csv`, is not a
numbered manuscript figure, and shows representative DAGs for the 10 testable
Path24 equivalence classes while leaving the saturated class unscored.

The visual grammar follows the supplied salamander-genome reference where it is
compatible with the journal contract: complete black panel boxes, no grids,
parenthesized panel tags, compact gutters, italic scientific names, and a
coherent coral/cyan/viridis-derived palette. Figure 1 includes the exact audited
24-tip Path24 time tree; branch time is shown as supplied while its publication
identifier, calibration method, and branch-length provenance remain pending.

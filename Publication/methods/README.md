# Methods release

The editable working draft is
`Materials_and_Methods_GBE_draft.docx`; its auditable source is
`Materials_and_Methods_GBE_draft.md`.

The DOCX has 1-inch margins, double-spaced body text, continuous line numbers,
page numbers, and seven embedded supporting figures. It is intentionally
explicit about the current claim boundary: image-derived nuclear IOD is a
relative phenotype, not an independently validated absolute genome size.

## Evidence and missing information

- `METHODS_EVIDENCE_MATRIX.csv` maps every planned Methods section to repository
  evidence, supplemental data, figure support, and manuscript disposition.
- `METHODS_REQUIRED_AUTHOR_INPUT.csv` is the blocking author checklist. It
  contains the wet-lab, specimen, imaging, upstream software-command, and
  authorship details that were not recoverable and must not be invented.

## Figure coverage

| Analysis | Supporting figure |
|---|---|
| species panels and genomic-resource QC | Figure 1 |
| TE diversity and CLR ordination | Figure 2 |
| RepeatMasker divergence landscape | Figure 3 |
| LTR terminal:internal deletion-footprint proxy | Figure 4 |
| reviewed cell and nucleus morphometry | Figure 5 |
| pairwise Pagel-lambda PGLS | Figure 6 |
| exploratory phylogenetic path analysis | Figure 7 |

Specimen collection, erythrocyte harvesting, slide preparation/staining, and
image acquisition are protocol sections rather than completed analyses. Their
text and any workflow schematic cannot be finalized until the corresponding
author-input rows are completed.

## Rebuild

```bash
uv run --with python-docx python scripts/publication/build_methods_docx.py
```

# Methods working materials

`METHODS_BULLET_SCAFFOLD.md` is the sole active factual Methods authority. It is
an internal evidence ledger: its labels, repository paths, and author-action
items are deliberately excluded from researcher-facing prose.

`MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` is the author-edited narrative
authority, and `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md` is its synchronized
plain-text mirror. The DOCX is the sole shareable Word draft. It uses conventional
scientific paragraphs and reports
limitations as ordinary prose. It contains no audit tags, internal dataset
codes, repository paths, or action-item placeholders. Verify the two files with:

```bash
make methods-word
```

The narrative remains a working draft rather than an independent factual
authority or submission-final manuscript. The verification command and test
suite compare visible prose and refuse to overwrite an existing author-edited
DOCX; explicit replacement from Markdown requires the builder's `--replace` flag.

## Evidence and missing information

- `METHODS_EVIDENCE_MATRIX.csv` is the current section-to-evidence crosswalk.
  Every repository path now points to the compact active analysis surface;
  unresolved upstream or author-supplied details remain explicitly marked
  partial or blocked rather than routed through retired analysis18 scripts.
- `METHODS_REQUIRED_AUTHOR_INPUT.csv` is the blocking author checklist. It
  contains the wet-lab, specimen, imaging, upstream software-command, and
  authorship details that were not recoverable and must not be invented.

## Figure coverage

| Analysis | Supporting figure |
|---|---|
| species panels and genomic-resource QC | Figure 1 |
| TE diversity context and Bonsai composition representation | Figure 2 |
| RepeatMasker divergence landscape | Figure 3 |
| LTR terminal:internal deletion-footprint proxy | Figure 4 |
| reviewed cell and nucleus morphometry | Figure 5 |
| conditional genome-size pairwise Pagel-lambda PGLS | Figure 6 |
| conditional genome-size exploratory phylogenetic path analysis | Figure 7 |

The researcher draft states only the specimen preparation, staining, imaging,
and computational details supported by the current evidence. Missing laboratory
or upstream-command details remain in `METHODS_REQUIRED_AUTHOR_INPUT.csv`; they
are not converted into bracketed warnings or invented prose.

## Submission-prose boundary

Historical prose and the redundant tagged Word export remain recoverable from
Git history and the pre-cleanup archive. The single researcher draft remains
subordinate to the evidence scaffold and author-input checklist. Final
submission prose still requires author review of the unresolved source-tree,
laboratory, segmentation, and upstream-command details.

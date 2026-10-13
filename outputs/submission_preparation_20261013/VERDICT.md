# Reconciled Methods; submission readiness remains open

Prepared 13 October 2026 against repository HEAD
`15dfb3c18755cc8659b626fa213087c62c6db02f` and the preserved dirty working tree.
No commit, scientific refit, upstream producer rerun, or external deposition
was performed.

The [Word Methods](../../Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx)
and [Markdown Methods](../../Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md)
now describe the completed analyses in the final eight-figure presentation.
This includes the strategy/SVL, diversity21, LTR19, paired allometry,
animal-weighting, expanded body-role, and exact21 TE-adjusted analyses.
The [evidence matrix](../../Publication/methods/METHODS_EVIDENCE_MATRIX.csv)
contains 26 sections and 78 distinct existing local source paths.

The imaging account incorporates the author's clarification: every retained
morphometry and IOD image is brightfield 40×; poor 40× oil-immersion reimages
were excluded. The unconfirmed all-animals ethanol-fixation statement was
removed. The actual genome-size denominator and the distinction between selected
t intervals and legacy normal-approximation intervals are explicit.

## Verification

| Check | Actual result | Receipt |
|---|---|---|
| Complete maintained test suite | 232 passed in 251.73 s; exit 0 | [command and log hash](checks/full-suite.json) |
| Final Methods/repository/DAG tests | 25 passed in 6.07 s; exit 0 | [final focused receipt](checks/methods-author-scope-final.json) |
| Publication validator | 43 tables; TE34=34, Path24=24, overlap21=21; exit 0 | [validator receipt](checks/release-validator.json) |
| Scoped Ruff | PASS; exit 0 | [lint receipt](checks/lint-author-scope-final.json) |
| Scoped diff whitespace check | PASS; exit 0 | [diff receipt](checks/diff-author-scope-final.json) |
| Word prose and package preservation | 61 paragraphs; 29 original paragraphs preserved; final scope correction preserves 60/61; only document.xml changed | [initial receipt](word_sync.initial.json), [final receipt](word_sync.json) |
| Bundled LibreOffice render | 8 pages; all inspected at original resolution, no clipping, missing glyphs, or orphan headings | [render command](checks/word-render-final.json), [rendered PDF](docx_render_final/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.pdf) |
| Protected final artwork and frozen inputs | All 614 baseline files unchanged; ZIP CRC and all 26 manifested entries match masters | [artifact integrity](artifact_integrity.json) |

The complete suite preceded a lint-equivalent set-comprehension edit, an
LF-only CSV normalization, and a four-word qualification limiting the
heart-sampling protocol to directly sampled specimens, with its guard assertion.
The final focused suite and lint checks cover the final files. Initial
wording/lint failures and the independent review's initial
CRLF failure remain retained, followed by successful correction receipts.
The builder now renders subsection headings as Word Heading 2 and still refuses
to overwrite an existing author document by default. Its visible behavior and
the refusal are tested. Scoped diffs are retained in this run; the builder's
pre-change copy was reconstructed by reversing only the coordinator's three
changes.

Independent review checked source claims, actual Word contents and pages,
and negative controls for canonical-source isolation and scientific caption
caveats. Its [report](INDEPENDENT_REVIEW.md) qualifies the reconciliation.
The old auxiliary DAG PDF differs from its currently tested renderer and is
excluded from the selected eight-figure submission inventory. It was not
regenerated; renderer acceptance is not a claim about that legacy asset.

## What remains before submission

1. Recover unknown-slide → staining-run → standard-slide records. Review the
   designated standards under the target QC pipeline, resolve the recorded
   approximately 41% standard-slide median discrepancy, and establish compatible
   acquisition settings before deciding whether to recalibrate. Current estimates
   use target D. fuscus images, not run-matched standards. See the
   [plain explanation and worksheet](CALIBRATION_RECORDS.md).
2. Supply or reconcile specimen/smear origins and preservation, applicable ethics
   records, acquisition settings, segmentation checkpoint lineage, tree citation
   and calibration, and actual upstream software/library execution records.
   Historical command recovery supplies configuration, not per-sample execution.
   See the [author-input checklist](../../Publication/methods/METHODS_REQUIRED_AUTHOR_INPUT.csv).
3. Assemble and review the Abstract, Introduction, Results, Discussion, primary
   references, title/authorship/affiliations, funding, data availability, and
   required disclosures. The supplied draft has placeholders and causal/rate
   claims that exceed the retained endpoints. The
   [manuscript reconciliation](MANUSCRIPT_RECONCILIATION.md) supplies a neutral
   title, revised aim, eight measured-endpoint Results units, and five proposed
   supplementary tables.
4. Assemble the manuscript and supplement PDFs, secure final author review,
   complete data/code deposition identifiers, and check the final files against
   the [official GBE instructions](https://academic.oup.com/gbe/pages/General_Instructions).

Artwork integrity, Methods coverage, and software checks pass. They do not
establish validated absolute genome sizes, an LTR deletion rate, causal
orientation, author approval, or a completed publication submission.

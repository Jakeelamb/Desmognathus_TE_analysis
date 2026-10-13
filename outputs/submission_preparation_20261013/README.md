# Submission preparation

This run reconciles the Methods with the eight figures in
`outputs/manuscript_polished_20260924_01`. It preserves the prior author draft
and records unresolved submission facts. It does not refit scientific models.

Completion requires that the Word and Markdown Methods agree, every final
figure's analytical procedure has source evidence, the maintained checks pass,
and the missing author records remain explicit until supplied. Manuscript
assembly and public deposition remain open until their actual artifacts and
identifiers are available.

## Workflow

1. Preserve the existing Methods and pin the final figure run.
2. Reconcile the scientific narrative, scaffold, and evidence matrix with the
   completed analyses. Keep unresolved author facts in the checklist.
3. Repair the two stale verification contracts while retaining their intended
   repository and scientific safeguards.
4. Synchronize the Word draft and visually inspect every rendered page.
5. Run source, document, and repository checks. Review the changed material
   independently and preserve the commands, exit codes, and logs.
6. Prepare a matched submission inventory and identify the next author inputs.

`baseline.json` records the prior hashes. `baseline/` contains exact copies of
the owned files before reconciliation. `decisions.tsv` records decisions and
their evidence. Checks and review receipts remain in this run.

## Scope and ownership

The Methods writer owns the narrative Markdown, factual scaffold, evidence
matrix, author checklist, and Methods README. A separate worker owns the two
stale repository and DAG-reference test contracts. The coordinator owns the
Word document, Methods tests, submission inventory, integration, and execution.
Provenance recovery is read-only. The final artwork and retained scientific
results are protected by the baseline hashes.

Throughput checkpoint: three bounded parallel tasks, one coordinator for Word
integration and verification. No scientific producer or HPC rerun is planned.

## Current review artifacts

- [Reconciled Word Methods](../../Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx)
  and its [Markdown companion](../../Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md).
- [Manuscript reconciliation](MANUSCRIPT_RECONCILIATION.md): revised aims,
  eight Results units, table placement, and claims needing revision.
- [Calibration explanation](CALIBRATION_RECORDS.md) and
  [staining-run worksheet](STAINING_RUN_MAP_TO_COMPLETE.csv).
- [Selected submission inventory](SUBMISSION_INVENTORY.csv): current eight
  collages, five supporting tables, workbook, Methods, captions, and package.
- [Writer source evidence](methods_writer_evidence.md),
  [provenance recovery](provenance_recovery.md), and [raw check receipts](checks/).

The full suite passed 232 tests. The final focused suite passed 25 tests after
an equivalent lint correction, normalization of the evidence CSV's line
endings, and qualification of the direct-sampling account. The publication
validator passed 43 tables and the 34/24/21 species contracts. All eight Word
pages were visually inspected. Only the Word body
part changed; unchanged paragraphs and all other package parts were preserved.

Submission readiness remains open. The selected artwork and completed analyses
do not supply the missing laboratory calibration, source-tree provenance,
specimen/ethics and upstream runtime records, or an assembled manuscript.
The older seven-figure release and its stale auxiliary DAG PDF are not included
in the selected inventory. Passing renderer tests does not certify that older
PDF as refreshed.

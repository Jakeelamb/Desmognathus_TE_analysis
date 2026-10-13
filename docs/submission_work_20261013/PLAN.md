# Methods work plan

This work follows the remaining-gap review at `d64f5a13981e1ff0edc2e30d66940110ac04be06`. The author narrowed the task to Methods. Abstract, Introduction, Results, Discussion and submission packaging are outside this task. It does not certify absolute genome sizes or authorize changing the scientific release.

## Acceptance

- Existing data, images, masks, tables, figures and frozen outputs retain their captured SHA-256 fingerprints. Supported Methods corrections are explicit exceptions and retain Word/Markdown agreement.
- A standards-only workflow reads the existing masks, preserves all 47 original identities, records linked-pair eligibility, and projects technical quality against the unchanged 40-object template and scalar calibration. Manual acceptance remains a separate field. Subset balance and full-template balance are reported separately.
- Author-supplied MVZ and RAP records are preserved exactly and joined to released animal inventories. Original and suggested names, collection dates, source namespaces and collecting permits remain distinguishable.
- A bounded records search logs locations and results for specimen, acquisition, tree, model and upstream-job provenance. Missing facts remain open instead of becoming inferred facts.
- Relevant code behavior, actual generated review artifacts, Methods source agreement and preservation checks are independently reviewed before landing.

## Phases

1. Capture the current scientific baseline and inspect the native microscopy interfaces.
2. Implement and exercise the standards-only existing-mask linkage/review/projection workflow in a fresh output directory. Do not rerun segmentation, fit new quality scaling, pool unreviewed standards or replace the released denominator.
3. Reconcile the supplied specimen records and search existing laboratory, model, tree and upstream-job records. Describe only supported facts in Methods.
4. Add an isolated, locked Pixi image-QC environment and local audit tasks. Preserve existing environment lock records and installed prefixes.
5. Independently inspect artifacts, test meaningful rejection controls, check the complete protected inventory, and land only the reviewed small files. Report author decisions and experimental evidence still needed.

## Ownership and execution

The coordinator owns integration, Pixi changes, Methods edits, canonical status records, scientific preservation and landing. Standards and specimen workers use separate detached worktrees. A records investigator writes a separate report. Manuscript work was stopped when the author clarified scope. No worker edits the current scientific release or another worker's files.

Existing Pixi environments and lock records are preserved when adding any isolated image-QC dependencies. Historical UV/Conda notes are not current execution instructions. All project commands use Pixi.

The standard review and specimen-record work can run independently after the baseline. Decisions and acceptance outcomes are recorded in `decisions.tsv` as they occur. Manual object review, confirmed specimen mappings, photometric validation and unavailable HPC assets remain external inputs; elapsed time does not supply them. The author asked to set aside the approval-record question; it is outside the remaining work.

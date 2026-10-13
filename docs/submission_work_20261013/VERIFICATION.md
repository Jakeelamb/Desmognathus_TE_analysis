# Methods audit verification — 13 October 2026

The Methods draft now describes the specimen-record comparison and the
designated-standard quality audit. Word and Markdown agree, and all nine
rendered Word pages were inspected. Work remained within Methods and local
reproducibility. The other manuscript sections were not drafted.

## Observed checks

| Check | Result |
| --- | --- |
| Maintained `make test` | 365 passed; 11 archive-dependent tests skipped |
| Dedicated native image-QC tests | 23 passed, including all 11 skipped native cases; one SciPy KS calculation fell back to its asymptotic method |
| `make lint` and `make methods-word` | Passed |
| Installed-package checks | Default, Bonsai and image-QC environments matched the lockfile; R, Cairo, Arial and the private R library passed |
| Real Bonsai handoff | Clean environment guard and nested default helper passed; retained comparison reproduced exactly |
| Registered specimen audit | Completed; 97 MVZ and 47 RAP source records preserved; 41/42 morphology and 49/50 IOD animals have unique numeric candidate joins |
| Registered standard audit | Completed; all 47 identities retained; 9 and 6 strict technical matches; no manually accepted matches |
| Isolated review-browser check | 47 cards; 141/141 images on each page; no console, page or request errors; exact fixture export |
| Preservation | 12,282 unique protected entries checked, including the four unchanged historical dependency snapshots; only the six authorized Methods files changed; no entries missing or added |
| Existing Pixi lock records | All 287 default and 290 Bonsai package records unchanged |

The [compact receipt](verification_receipts.json) records clocks, hashes and
commands. Full logs, failed attempts, the baseline inventory and execution
intermediates are retained privately under `.git/submission_gap_work_20261013`.
The browser decisions were test fixtures, were cleared, and were not imported
as author acceptance.

The checks exposed inherited Python startup settings, a hidden Pixi executable
and host-probe tools, and temporary Git fixtures inheriting the user's commit
hook. Maintained commands now use `--clean-env`; child Pixi calls use `PIXI_EXE`
and prefix-first host probe tools. Hook and signing isolation applies only to
temporary test repositories. Existing dependency sets and scientific outputs
were preserved.

## Remaining Methods evidence

- Both designated standards fail the existing image-quality balance thresholds.
  Manual review, photometric validation, a justified DNA-content reference and
  calibration uncertainty remain unresolved. Genome-size estimates retain
  their conditional assembly anchor and original denominator.
- The author reports animal 3463 as a likely transcription error. Its intended
  record remains unconfirmed; RAP3462 and the released animal counts remain
  separate. Source namespaces and the RAP3518 comment also need interpretation.
- Detailed preservation and acquisition settings, the legacy checkpoint-to-mask
  mapping, the focal tree's calibration record and complete HPC run/library
  linkage remain partial. Recovered model bytes and historical scripts do not
  supply those missing records.

The [author-input ledger](../../Publication/methods/METHODS_REQUIRED_AUTHOR_INPUT.csv)
tracks these details. This audit verifies the local implementation and preserved
artifacts; it does not establish a validated absolute genome-size scale or
certify submission readiness.

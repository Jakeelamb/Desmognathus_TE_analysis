# Scaffold and ledger edit notes

Edited only `Publication/methods/METHODS_BULLET_SCAFFOLD.md`, `METHODS_REQUIRED_AUTHOR_INPUT.csv`, `README.md`, and this report. The full narrative Markdown and Word document belong to the other editors.

The scaffold now uses sentences for panel definitions, historical provenance, PCA preparation, RepeatMasker and LTR rules, image acquisition, IOD aggregation, and selected figure procedures. Long clauses were split where needed. Factual labels, citations, commands, settings, biological limitations, and author actions remain.

The calibration account now records the author-confirmed single Feulgen run and the 53-row run map. Standards 413/414, specimens 32469/32470, were identified as 40× brightfield images without oil, using the green channel and the same fixation protocol as the other slides. These facts come from the retained author answers, not inferred acquisition metadata. The shared protocol does not establish exact fixation reagents or lifetime preservation.

The implemented denominator remains the target-species *D. fuscus* estimate from Processes 314/315/320/321. Accepted standard QC, the approximately 41% between-slide median difference, camera linearity, and an independently validated C-value remain unresolved. W06 remains `partial`; W07 retains its existing status. Only the evidence cells in W01, W03, W06, W07, W08, W09, and W16 changed.

The README describes the user-authorized prose edit and Word synchronization, links the editorial package, and distinguishes synchronization from rendered-page inspection. It retains the builder's default protection of an existing author-edited DOCX.

## Checks

The stdlib comparison ran through the existing hydrated environment with `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python`.

- PASS: CSV schema has the same six fields and all 16 rows in the same order. Item IDs, statuses, required details, and all other non-evidence fields are unchanged.
- PASS: The scaffold retains every baseline numeric token, inline-code token, and existing link target. New confirmation evidence and links were added.
- PASS: The retained run map has 53 rows, all assigned to `author_confirmed_single_run`.

The coordinator owns the integrated Methods tests, Word synchronization, rendering, and visual inspection. These local checks do not establish completed Word review or manuscript readiness.

An initial default `uv run python` attempted to build pandas and lxml before the edit script started. That owned invocation was interrupted with exit 130, and no related build process remained. Editing and comparisons then used the pre-existing hydrated environment without dependency provisioning.

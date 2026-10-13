# Designated-standard mask review

This workflow links the existing masks for Process 413/specimen 32469 and
Process 414/specimen 32470, computes their native shape/focus diagnostics, and
projects eligible candidates onto the preserved 40-object technical-quality
template. It keeps all 47 original nucleus identities. It does not segment
images, repair masks, fit a scaler, pool standards, replace the released
denominator, or validate absolute genome sizes.

The original masks, measurement records and native helpers remain external
archival inputs. `scripts/microscopy/standard_qc_source_pins.json` identifies
their required SHA-256 fingerprints. The registry itself is pinned by the
adapter. Source versions and template/scaler contents cannot be substituted
silently. The archived measurement records contain absolute source paths;
those paths must remain readable for the native producer. This workflow is
not a self-contained distribution of the microscopy archive.

## Run a fresh review

Use the repository's locked `image-qc` Pixi environment. `make setup` installs
the maintained environments; the image-QC task supplies the image packages.
Choose an output path whose parent exists and which has never been used.

```bash
pixi run -e image-qc --locked --no-install --clean-env standard-qc \
  --paper-root /home/jake/Projects/Desmognathus_TE \
  --cell-root /home/jake/Projects/cellprofiler_test \
  --archive-root /home/jake/Projects/Desmognathus_TE_archive \
  --output-dir /path/to/new/standard-review
```

The command refuses an existing destination, symlinked output ancestry,
input/output overlap, altered source fingerprints, unexpected standard
identities or an incomplete 47-object census. It invokes the genuine native
existing-mask linkage CLI with `--image-type pmount --window-radius 1`.
`pmount` preserves the historical metadata label; the author has identified
these retained images as 40× brightfield without oil, in the green channel.

Source-global coordinates remain unchanged for physical linkage and native
shape/focus computation. Nuclear pixel metrics come from the validated
`pixel_audit_v2` local-coordinate view, joined by the original nucleus identity.
Native linkage does not supply a cell-edge flag. The adapter measures the
recorded cell label on the source mask's borders, rather than copying the
historical scoring script's default `False` flag. Missing label evidence and
nonfinite required shape/focus/matching evidence cannot enter the strict pool.

## Review the queue

Open `review/index.html` from the new directory. Its paired raw/cell/nucleus
views use the native static reviewer, with target measurements blinded,
default decisions blank and automatic keep filling disabled. Every original
candidate remains visible, including technically excluded objects. The
unchanged full-label nuclear audit crops are retained alongside the paired
views for boundary checks. The native viewer saves manual choices in browser
storage; export its decisions CSV to retain them outside that browser.

The valid manual decisions are `keep`, `problem` and `unsure`. A blank decision
remains unreviewed. Exported decisions do not edit masks or source tables.
Run a new projection in another fresh directory with the export:

```bash
pixi run -e image-qc --locked --no-install --clean-env standard-qc \
  --paper-root /home/jake/Projects/Desmognathus_TE \
  --cell-root /home/jake/Projects/cellprofiler_test \
  --archive-root /home/jake/Projects/Desmognathus_TE_archive \
  --decisions-csv /path/to/exported_decisions.csv \
  --output-dir /path/to/new/reviewed-standard-projection
```

Repeat `--decisions-csv` for later review rounds. Later files replace earlier
decisions for the same identity. Duplicate identities within one export and
decisions referring to targets or unknown objects are rejected.

`problem` objects are excluded before technical matching. A proposed match
becomes an accepted matched object only when its actual manual decision is
`keep`. An unmatched or technically ineligible keep remains a manual keep;
it does not enter the accepted matched summary. `unsure` is not acceptance.

## Interpret the outputs

`standard_census.csv` contains all 47 original identities, linkage and
shape/focus fields, technical eligibility, gate reasons, source-boundary
diagnostics, proposed matches and manual decisions. Original coordinates are
retained. `proposed_technical_matches.csv` records the frozen assignment and
template IDs. `frozen_template.csv` and `frozen_scaling.json` preserve the exact
reference rows and parameters used. `standard_summary.json` keeps the two
standards separate.

The summary distinguishes input, linked, finite-evidence, strict, finite-match,
matched, manual-keep and accepted-matched counts. Native assigned-subset
balance is reported against only the assigned reference slots. Full-template
balance uses all 40 unchanged reference rows and separately lists missing
template IDs. A small balanced subset does not establish full-template
coverage. The historical 30-kept guard belongs to a species review; no new
30-per-image threshold is applied to these standards.

An empty eligible group produces an explicit zero-candidate status and
40-slot shortfall. A group with 16 candidates can fill at most 16 slots and
has at least a 24-slot shortfall. Any accepted IOD median is reported for that
standard alone, as a review diagnostic. No pooled denominator or genome-size
recalculation is produced. Photometric linearity, a justified DNA-content
reference and the calibration uncertainty policy remain separate work.

`linkage.log` and `linkage_command_receipt.json` retain the native producer's
command, clock and exit status. `receipt.json` records the successful adapter's
interpreter, packages, source fingerprints and final preservation check.
Successful technical execution does not provide human review decisions.

## Acceptance checks

Portable boundaries run in the default environment as part of the test suite.
The native controls require the image-QC environment and the pinned archive:

```bash
pixi run -e image-qc --locked --no-install --clean-env --executable env \
  STANDARD_QC_NATIVE=1 \
  STANDARD_QC_PAPER_ROOT=/home/jake/Projects/Desmognathus_TE \
  STANDARD_QC_CELL_ROOT=/home/jake/Projects/cellprofiler_test \
  STANDARD_QC_ARCHIVE_ROOT=/home/jake/Projects/Desmognathus_TE_archive \
  python -m pytest -q tests/test_standard_qc.py
```

These checks reject incorrect process/specimen identity, duplicate objects,
altered source hashes, corrupt template IDs/features/scalars, missing finite
evidence and unreviewed objects counted as accepted. Genuine native matching
is exercised with zero and 16 candidates, including its subset-balance
shortfall behavior. A skip reports that native archival acceptance was not
executed; it does not count as a native pass.

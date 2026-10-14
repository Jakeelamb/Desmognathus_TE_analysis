# Microscopy recovery tooling verification

Baseline: `8a5c00b33c0359b4a47f2df2ead145039554d623`. The source fingerprints,
command clocks, exit statuses and retained log hashes are recorded in
[receipts.json](receipts.json). Local diagnostics and full raw logs remain under
`.git/recovery_tooling_20261013/` and `.pixi/runs/`; they are not release datasets.

## Result and scope

The tooling is ready for the first joint boundary review. It creates immutable
original-pixel packets, imports separately authored contours, prepares blinded
queues for explicit proposal states, records named decisions and reports
measurement sensitivities. No real nucleus was accepted or repaired by these
checks. No absolute DNA-content calibration was calculated or promoted.

The real pilot contains all 47 designated standards and 35 released *D. fuscus*
target nuclei. All 82 mask areas match their source labels exactly. Recorded-
background IOD replay has a maximum absolute difference of
`4.895912297797622e-7`, below the `1e-6` tolerance. Backgrounds are 255 for 67
objects, 218 for five and 216 for ten. Source-edge contact occurs for one object;
two crops have clipped context padding. These are different conditions.

The original packet and baseline report contain zero authored proposals and
decisions. The real queue stayed blank during GUI tests. The packet hash remains
`fd9afdd0fc827d304553d31a003dc9cef6c46a0caaeed0680bbed106650baa94`.

## Executed checks

| Check | Observed result |
|---|---|
| `make check-env` | PASS; installed default prefix, locked packages, R, Cairo, Arial and private library |
| `make test` | 415 passed, 29 explicit skips |
| `make lint` | PASS |
| Recovery backend, core and queue tests | 132 passed, no skips |
| Final core source after comment trimming | 67 passed |
| Existing standard-QC tests with explicit native archive bindings | 23 passed, no skips |
| Independent native single-ROI/proposal controls | 16 expected exits matched, including eight refusals |
| Actual browser queue and export | Original and proposed mask states rendered; blank/reason/reviewer controls refused invalid exports |
| Actual ImageJ GUI export | Native `.roi` and multi-ROI ZIP consumed without changing source pixels |

The default skips comprise 17 native backend controls, 11 archive-dependent
standard-QC controls, and one module-level skip for the 67 core image controls.
The default environment lacks Pillow and tifffile. The explicit recovery and
native archive runs cover those controls; skipped tests are not counted as
native acceptance. The final optional-import change preserves all 32 test
function bodies, verified independently by AST comparison.

The independent review is retained in [independent_review.md](independent_review.md)
and [independent_receipt.json](independent_receipt.json). The GUI scope and pins
are in [gui_review.md](gui_review.md) and [gui_receipt.json](gui_receipt.json).

## Producer handoffs

ImageJ's own ROI decoder and mask rasterizer passed asymmetric, concave,
subpixel, composite-hole and boundary geometry controls. Generic polygon
fillers selected different pixels in the prototype and were not adopted.
Independent ImageJ measurements matched selected pixel counts and raw sums.
The Rust reducer preserves signed OD, the pixel floor of one and exact uint8
histogram medians.

An actual GUI-drawn 3-by-2 rectangle was saved as `nucleus.roi`. Two imports
wrapped its same 142 bytes into byte-identical, deterministic ZIP archives.
The genuine ImageJ raster contained six pixels; independently calculated IOD
matched the Rust result within `1e-12`. An import alone left all 82 fixture
cases unreviewed. A separately exported browser event bound to the proposed
state was then imported and selected in a diagnostic report. All such decisions
and corrections were labeled synthetic fixtures and remained separate from
the real packet.

Rejection controls include modified source/crop/mask pixels, wrong identities
or coordinates, rehashed identifying reviewer assets, unsupported ROI geometry,
malformed or oversized archives, invalid native ROI filenames, overlapping
background selections, stale bindings, duplicate events, symlinks, existing
destinations and output paths inside preserved source trees. Reviewer conflicts
remain unresolved rather than selecting the last event.

## Preserved failures and corrections

- The first setup rejected Conda's `21.0.10-internal` label. The locked runtime
  was valid; the semantic version check was corrected. The failed cache remains.
- A real ImageJ save dialog hung loading host GTK. The launcher now disables
  GTK dialogs and uses private home/preferences. Actual saving passed afterward.
- ImageJ saves one ROI as `.roi`. The importer now wraps the native bytes;
  a reviewer does not need an extra background ROI to obtain a ZIP.
- The old validator accepted identifying files and altered HTML after their
  hashes were rebound. Six reproduced controls now fail under the exact file/
  directory layout, pinned HTML, canonical protocol and crosswalk checks.
- The first default test collections failed on optional image imports. Explicit
  module skips now preserve the original dependency set; image tests run in the
  declared image environments.
- One independent evidence projection expected a nonexistent `receipt.status`
  field. It was repaired from retained `exit=1` evidence without repeating the
  successful native producers.

Cua native accessibility was degraded and browser binding refused a profile
grant. GUI checks used the isolated desktop's rendered browser/DOM and native
keyboard/menu routes. This is not accessibility-tree acceptance. Jake's physical
screen and input focus were not used.

## Preservation

Post-check hashing compares 12,282 protected entries containing 434,825,677
regular-file bytes against the initial inventory. Data, Methods, figures,
scientific manifests, existing outputs and the four legacy dependency snapshots
have zero additions, removals, changes or exceptions.

The old default, Bonsai and image-QC environment blocks and all referenced full
package records remain identical: 287, 290 and 301 records respectively. An
in-memory package-hash corruption is detected by all three comparisons. The
new `recovery` environment and setup cache are separate.

Biological contour correctness, historical camera response and chemical
stoichiometry, independent reference DNA content and original VSI/ETS metadata
remain unverified. The next work is human boundary review, followed by the
predeclared background and acquisition experiments in the
[measurement protocol](../research/feulgen_recovery_measurement_evidence.md).

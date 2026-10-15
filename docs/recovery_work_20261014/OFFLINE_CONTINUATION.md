# Offline source-reader verification

Jake resumed the [paused work](../handoffs/20261014_SOURCE_RECOVERY_PAUSED.md)
without the easystore image drive. The local checks below ran on 15 October 2026
UTC, from checkpoint `0254d00c17830c7c1a4a319e51d032b8613f7033`.

Independent review passed the supported TIFF-page and stat-only reader behavior.
It found no blocking code defect. The reader, tests, Pixi manifest and lock stayed
at their checkpoint hashes. This continuation changes documentation only.
It does not qualify absolute genome-size calibration or publication readiness.

## Code and artifact checks

| Check | Result |
| --- | --- |
| `make check-env` | Passed |
| `make lint` | Passed |
| `make test` | 429 passed, 58 skipped, one warning |
| Explicit four-file recovery suite | 175 passed, no skips, one warning |
| Independent supported-TIFF controls | Passed seam, image-edge, one-pixel, metadata and read-only memory-map cases |
| Missing or truncated selected tiles | Refused with failure receipts and no crop output |
| Source changed after decoding | Final receipt failed even though a crop and completed decoder record existed |
| Cached target comparisons | All 35 local TIFF crops equal the packet's actual green PNG pixels |
| Fresh packet equivalence | All 82 case states and metadata equal, with 246 corresponding PNGs equal in bytes and pixels |
| Original packet preservation | Its 498-file tree still matches its saved baseline |
| Protected scientific artifacts | All 12,282 entries, totaling 434,825,677 file bytes, match the retained baseline |

The default quiet test output does not state skip reasons. Its 58 skips are not
claimed as image-source coverage. The recovery suite ran every selected test.
Its duplicate-ROI ZIP warning comes from the deliberately malformed archive
fixture.

The independent verifier used genuine big-endian tiled RGB TIFFs with rectangular
16 × 32 tiles and an image size of 43 × 77 pixels. It observed selected byte
ranges, decoder positions and complete assignment of each output pixel. Missing
offset/count controls failed after earlier seam tiles had decoded, rather than
writing partially initialized output. These were synthetic controls, separate
from biological data and acceptance records.

The 35 cached comparisons cover 12 Process 314 crops, eight Process 315 crops,
five Process 320 crops and ten Process 321 crops. This check independently reads
saved crop and packet bytes and audits the original command, tile-range and log
records. It does not decode the unplugged source exports again. Current source
availability, timestamps and payload identity were not rechecked. The original
source reads used stat-only fingerprints, not full large-file hashes.

Only a complete final `receipt.json` qualifies a reader output. The operator
guide now distinguishes that final result from a completed decoder command,
and distinguishes successful header-description retention from early failures.

## Repeat the selected local tests

The maintained Make checks use locked Pixi. The selected recovery invocation is:

```bash
pixi run -e recovery --clean-env --locked --no-install --executable python -m pytest \
  tests/test_recovery_backend.py tests/test_microscopy_recovery.py \
  tests/test_recovery_review.py tests/test_source_reader.py -q
```

Use pytest's default temporary directory or a fresh `--basetemp` outside the
checkout. A temporary directory inside `.git` triggers the recovery output guard.
Do not add `.` as another collection target; archived upstream and baseline tests
are outside this selected suite.

Both invocation mistakes occurred during this continuation and their failed logs
remain preserved. The broad collection attempt exited 2. The in-checkout
temporary-directory attempt exited 1. The corrected four-file run in a fresh
`/tmp` directory exited 0 with 175 passes. No guard was relaxed or archived
dependency installed. A separate help-recording attempt also failed before
reader execution because it tried to launch nested Pixi inside a clean Pixi
environment. The corrected check used the already selected Python executable.

## Recovered-context review

Separate `gpt-6.1-sol` and `gpt-6-luna` reviewers inspected the two saved wider
contexts without the source crosswalk, prior measurements or each other's
findings. The coordinator also inspected the unchanged raw images and overlays.
This was an image review, not a biological acceptance decision. The reviewers
are separate models in the same family, not a cross-family panel.

| Case | Findings | Diagnostic next step |
| --- | --- | --- |
| `N815f697d43b9efed92c4` | Both reviewers find the original mask partial and consider an exploratory contour reasonable. Sol sees the apparent dark body fully in-frame; Luna considers its exact boundary unresolved. The upper/right or lower/right transitions remain soft. | A separately authored contour trial may follow the dark body while retaining boundary uncertainty. Do not use the larger pale outer perimeter as the nucleus. |
| `N3291d5e0a6bc3fe2af23` | Both reviewers hold a nuclear contour. Diffuse interior material and irregular lower-left structure do not supply a defensible closed boundary. | Keep the case unresolved. Display enlargement alone cannot recover missing structural detail. |

One reviewer described N815's cyan fragment as being at the top. The coordinator
does not adopt that spatial description. The unchanged mask occupies
`[96, 65, 141, 83)` within a 237 × 179 crop and covers the lower part of the
apparent dark body. The reviewers' agreement about a partial mask does not erase
their difference about boundary completeness.

No ROI, review acceptance, mask replacement or IOD calculation was created.
The recovered contexts remain separate from the original packet. Integrating a
wider contour requires explicit crop, mask-coordinate and proposal-state bindings;
an old acceptance or state hash cannot transfer to a wider crop.

## Evidence and limits

Complete receipts, exact argv, UTC clocks, exits, raw logs and file hashes remain
on this computer under `.git/drive_recovery_resume_20261014/`:

* `verifier/` contains independent code, cached-artifact and native-TIFF evidence.
* `checks/` contains the maintained checks, corrected suite and failed attempts.
* `visual-sol/` and `visual-luna/` retain independent image findings and input hashes.
* `root/` contains preservation checks, CLI help and coordinator image findings.

These ignored records and cached images are not included in a portable Git clone.
The small [verification summary](offline_verification.json) pins the local evidence
used here without adding raw images, native tools or datasets to Git.

Full native-source re-extraction and current drive-state checks remain unavailable
until the original files are connected. TIFF pixel-decoding limits do not impose
a general metadata-allocation bound on arbitrary malformed TIFFs. Physical camera
response, biological boundary truth and independent DNA-content calibration remain
unresolved. No scientific output, manuscript Methods, figure, mask or retained
dependency snapshot changed during this continuation.

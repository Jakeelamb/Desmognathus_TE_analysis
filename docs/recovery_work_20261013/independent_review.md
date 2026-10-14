# Independent final microscopy tooling review

**PASS for the reviewed native ROI import, proposal review queue, explicit authored-state handoff and preservation requirements.** This verdict does not establish biological mask validity or absolute genome-size calibration. It does not substitute for the coordinator's final maintained repository checks.

The identified base revision is `8a5c00b33c0359b4a47f2df2ead145039554d623` plus the working-tree fingerprints in [independent_receipt.json](independent_receipt.json). Final `recovery.py` SHA256 is `9aafedba21482498e23c02b95eef52b839d18740662991332d1a95dddbc15273`; the first phase executed `07fd3f1d846c861a6c045cf2426d6da3568e032176b7cfb35608d796cc4b75f3`. I independently compared their parsed ASTs, excluding source-location attributes: exact equality. The final phase ran the remaining controls on the final source. The final test module fingerprint is `0208cc42f34a25abc11676c9c1e80c58f5edb9e64f477f01e35df21f84006f8a`. Its 32 function bodies are AST-identical to the original module; only optional image dependency loading differs. The check is recorded in `final-test-source-check.json`.

## Executed behavior

I inspected the actual single-ROI implementation, ImageJ bridge, packet validation, proposal selection, authored-event validation, report aggregation and shipped page logic. Native toolchain identity was `5c35cb7aaba13706e1a11c0c819a65c841b964d735aa0ea9ab28a243e2c1aae5`.

Using the existing synthetic-only `gui_fixture_bridge.py`, **16 CLI controls matched all expected exits: eight successes and eight refusals**. The raw GUI-produced `nucleus.roi` was imported twice. Both archives contained exactly one stored `nucleus.roi`, a fixed ZIP timestamp and the original 142 bytes; archive hashes and proposal states were identical. ImageJ rasterization matched an independently decoded 3-by-2 integer rectangle exactly. The Rust IOD matched an independent signed `math.fsum(log10(background/max(pixel,1)))` calculation within `1e-12`.

Import alone produced zero decisions, an empty approved subset and 82 missing reviews. The selected proposal queue contained exactly 167 served files, exactly six allowed metadata keys per case, all 82 unchanged raw crop hashes and the selected proposal's mask/state. Its HTML and protocol matched the canonical sources. An event for that proposal failed unless the explicit proposal directory was supplied. A named synthetic event then produced one kept case and 81 missing cases.

Negative controls rejected another filename, an empty ROI, a file exceeding 1 MiB, malformed bytes, stale authorship, an existing destination, an unknown proposal event, and a copied packet containing a rehashed identifying TSV. Early refusals created no output. The malformed attempt retained its failed receipt and published no proposal. Both original synthetic and real packets retained blank proposal/decision arrays.

Exact argument lists, per-command clocks, exit statuses and stdout/stderr are in the `commands` entries in [independent_receipt.json](independent_receipt.json). The executions were:

```bash
pixi run -e recovery --clean-env --locked --no-install python .git/recovery_tooling_20261013/final_review/verify_final.py
pixi run -e recovery --clean-env --locked --no-install python .git/recovery_tooling_20261013/final_review/complete_retained.py
```

These are single-use evidence scripts with fixed fresh output paths. The first exited 1 because its evidence projection incorrectly expected `receipt.status`; all 12 completed CLI controls had already matched their expected exits. The implementation's receipt correctly records `exit=1`. The second repaired the projection from retained artifacts and executed only the four remaining controls, exiting 0. No successful native producer was repeated to manufacture evidence. Raw logs and wrapper clocks/exits: `integration.log`, `execution.json`, `completion.log`, `completion-execution.json`.

## Preservation and environment checks

I reran the existing independent preservation and environment probes through the locked **default** Pixi environment, changing only their assigned evidence output paths. Their original script SHA256s and that output-only substitution are recorded in `preservation-probe-pin.json` and `environment-probe-pin.json`.

- **12,282 protected entries / 434,825,677 regular-file bytes unchanged**: no additions, removals, modifications or exceptions. This includes retained Methods, data, figures, outputs and the four legacy dependency snapshots. Evidence: `preservation.json`.
- Old default, Bonsai and image-QC environment blocks and all **287 / 290 / 301** referenced full package records were byte-identical. An in-memory shared-package SHA mutation was detected by all three comparisons. Evidence: `environment-records.json`.

## GUI artifact inspection and limits

I read the actual latest browser producer receipt, its single TEST event, the consumed event artifact and import receipt, and visually inspected both 800% screenshots. They show the selected six-pixel rectangular overlay and its removal. The event hash `533c4ee4eef5a02f94f98001ec5359ffdcc3db4af36d92b511ea9e6545d65392` matches the consumed export, and its proposal state matches my native import. The real queue's separate browser receipt records zero authored events. These producer artifacts are under `../gui_final/`; I did not operate the GUI independently during this final delta review.

The Cua receipts explicitly report degraded AT-SPI access, and the consumer-profile grant was refused. Browser work used the isolated desktop's DOM route and actual rendered screenshots; this is not accessibility-tree acceptance. The documented host GTK file-dialog hang remains relevant; current `prepare-editor` disables GTK dialogs and uses private preferences. I inspected that implementation but did not independently relaunch it. The default environment lacks PIL/tifffile, so the new core image test module explicitly skips there; native/image fixtures require the recovery environment. PyYAML was absent; the environment probe used an exact byte-block comparison without installing a package.

No infrastructure evidence establishes a correct biological contour, historic camera linearity, chemical stoichiometry, independent reference C-value, missing original VSI/ETS metadata, or the benefits of reimaging. Hashes provide integrity checks, not cryptographic proof of reviewer identity. Broader scientific recovery still requires the planned human review and acquisition/reference experiments. I did not duplicate the coordinator's full test/lint runs or infer their results from a worker summary.

Only assigned private evidence and synthetic fixture outputs were written by this reviewer. No implementation, acceptance check, retained dataset, Method, figure, or real-case decision/proposal was edited.

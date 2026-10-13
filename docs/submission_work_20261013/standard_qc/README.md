# Standard-slide audit

This snapshot comes from the registered `standard-qc` Pixi task, using existing
masks and the pinned native helpers. It contains compact audit tables and a
3 MB review gallery. Raw mosaics, model weights and original analysis outputs
were not copied or changed.

| Standard | Original nuclei | Linked nuclei | Strict technical matches | Accepted matches |
| --- | ---: | ---: | ---: | ---: |
| Process_413 / MVZ32469 | 31 | 19 | 9 | 0 |
| Process_414 / MVZ32470 | 16 | 10 | 6 | 0 |

Neither group passes the assigned-subset or full 40-nucleus template balance
thresholds. All 47 identities remain in the census and review. Manual decisions
are blank; no accepted standard median or calibration denominator was produced.
These technical results do not validate absolute DNA content or explain the
approximately 41% difference between the original slide medians.

Open [the paired-mask review](review/index.html) and
[the full-label nuclear crops](review/nuclear_audit.html). Keep, problem and
unsure decisions must come from actual review. A technical match does not imply
acceptance. Exported decisions can be passed to a fresh audit run.

`SNAPSHOT.json` pins the included files and identifies the private execution
directory. `receipt.json` preserves the successful task's runtime and source
checks. Its original paths and complete linkage intermediates refer to that
execution directory; they are not a portable raw-data package. The sibling
microscopy project and archive are still required to rerun the image audit.

See [the runnable workflow](../../STANDARD_QC.md) for commands, source
requirements, rejection controls and the separation of subset balance,
full-template coverage and manual acceptance. Existing scientific results and
the conditional genome-size conversion remain unchanged.

# What we need to resolve the genome-size calibration

All retained cell, nucleus, and IOD measurements use brightfield 40× images.
The author confirmed this on 13 October 2026. Poor 40× oil-immersion reimages
were excluded. The older draft's 60× statement does not describe this release.

Nuclear IOD is the measured stain signal. The current calculation converts it
to conditional genome size by dividing a species' equal-image IOD estimate by
the released D. fuscus estimate and multiplying by 16.462167689 pg/1C. That
number comes from a 16.1-Gbp assembly estimate and the Gbp-to-pg conversion.
It is an assembly-derived anchor, not an independent DNA-content measurement.

The current denominator uses target images Processes 314, 315, 320, and 321.
The designated standard slides are different: Process 413/specimen 32469 and
Process 414/specimen 32470. The retained provenance records that they were
stained in the experimental runs with unknowns, but does not identify which
unknown belongs to which run or standard. It also records an approximately
41% difference between the standard-slide median IOD values. This preparation
has not remeasured that difference or established its cause.

The missing link is a record like:

> target slide or Process ID → staining run/date → standard slide or Process ID

A notebook page, slide inventory, staining worksheet, or existing spreadsheet
can supply that link. An image scan date alone cannot establish a staining run.
[The worksheet](STAINING_RUN_MAP_TO_COMPLETE.csv) contains the 51 retained target
images and the two designated standards. Its run and record-location columns
are intentionally blank. A mapping at slide level is sufficient if multiple
images share a slide; the source record should establish that relationship.

After recovering the map, we need to review and measure the standards using
the same QC rules as the targets, establish compatible acquisition settings,
and resolve the standard discrepancy before deciding whether calibration
changes are warranted. Staining or acquisition differences could affect
between-species IOD ratios; the record does not yet establish that they did.
A single global conversion constant cannot correct a run-specific effect.

The [current source manifest](../../analyses/02_morphology/provenance/iod_source_manifest.json)
and [recovery report](provenance_recovery.md) retain the exact denominator,
anchor, standard identities, and unresolved checks. No estimates or fitted
models were changed during this preparation.

The Pyron-supplied blood smears also need a confirmed collection and
preservation history. The broad statement that every animal was ethanol-fixed
was removed from the Methods pending that information.

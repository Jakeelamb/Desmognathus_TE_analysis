# Recovery pilot handoff

Read the [workflow](../MICROSCOPY_RECOVERY.md),
[research](../research/feulgen_recovery_measurement_evidence.md), and
[verification record](VERIFICATION.md).

The first joint-review packet is local and ignored:
`.pixi/runs/microscopy-recovery/pilot-v1/packet`. It contains all 47 standards
and 35 same-species target nuclei. Serve only its `reviewer` folder. Keep the
private source crosswalk outside the served root. All cases start unreviewed.

1. Review complete biological boundaries with the mask hidden and shown.
   Record visible defects and uncertainty before seeing measurement changes.
2. Prepare independent contour proposals on copied raw pixels. Keep original
   masks, authors and proposal states separate.
3. Review proposed states using an explicit queue policy. Original-state
   acceptance does not transfer to a correction.
4. Compare geometry-only measurements at the recorded background first.
   Report candidate-background sensitivity separately.
5. Use paired reimaging and controlled camera/illumination evidence if the
   pilot identifies an optical limitation. Retain source edge cases explicitly.
6. Validate independent reference DNA content and nuclear DNA class before
   making an absolute calibration claim.

Do not choose a contour, background, focal plane or exclusion rule because it
makes the two standards agree. Do not replace the retained release during this
pilot. The new toolchain supports these experiments; installation and numerical
parity do not establish their scientific outcome.

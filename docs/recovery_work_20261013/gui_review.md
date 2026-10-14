# Final proposal-review GUI acceptance

PASS within the assigned GUI scope. Checked at base revision
`8a5c00b33c0359b4a47f2df2ead145039554d623` with final uncommitted source
pins in `after.json`: core `9aafedba21482498e23c02b95eef52b839d18740662991332d1a95dddbc15273`,
review page `0b8e6c8b6d0e7001e897c1e8cab8097292724145215bd89e7bd4ba69fda540d1`.
Chrome was `155.0.8059.79`. This check did not modify implementation or existing
acceptance checks.

## Actual observations

- Real 82-case queue opened with blank reviewer, no selected decision and zero
  authored events. Native crop dimensions were 143 by 122; 400 percent display
  dimensions were 572 by 488. Mask visibility toggled correctly.
- The synthetic proposal queue contained 82 cases. Case `N3575a86c80f658661c3a`
  showed the six-pixel imported ROI mask at state
  `172039ef16c65cff8db6740fe4cebe6d7a0614511b330636cdf174b354ad19af`.
  Actual HTTP raw and mask bytes matched the selected queue files; the mask hash
  matched the proposal and the raw hash matched the unchanged crop.
- Actual browser export produced exactly one explicit TEST-only event with
  the proposal state, named reviewer, timezone-bearing time, decision and reason.
  The other 81 blank cases were omitted. Exported fields stayed disabled after
  navigating away and back.
- Negative controls refused empty export, missing reason and missing reviewer,
  with zero events exported before the completed TEST review.
- The focused maintained UI tests passed: `18 passed in 1.38s`.
- All 498 files in the real pilot packet matched the before inventory. The packet
  retained 82 cases, zero authored decisions and zero authored proposals.

## Exact execution

The browser recipe is `proposal_browser.py`; its actual invocation was the
isolated `mcp__agent_desktop__browser_script` tool with session
`microscopy-recovery-fina-c299b46b3a7c` and:

```python
from pathlib import Path
exec(compile(Path('/home/jake/Projects/Desmognathus_TE/.git/recovery_tooling_20261013/gui_final/proposal_browser.py').read_text(), '/home/jake/Projects/Desmognathus_TE/.git/recovery_tooling_20261013/gui_final/proposal_browser.py', 'exec'))
```

The selected queue was served read-only on loopback:

```bash
pixi run -e recovery --clean-env --locked --no-install --executable python -u -m http.server 8876 --bind 127.0.0.1 --directory .pixi/runs/independent-recovery-final/selected-queue/reviewer
```

Focused tests were run through the locked Pixi runtime and capture helper:

```bash
pixi run -e recovery --clean-env --locked --no-install --executable python .git/recovery_tooling_20261013/gui_final/capture_cli.py ui-tests /home/jake/Projects/Desmognathus_TE/.pixi/envs/recovery/bin/python -B -m pytest tests/test_recovery_review.py -p no:cacheprovider --basetemp=.pixi/runs/microscopy-recovery/gui-final/ui-tests -q
```

Exact Cua commands, timestamps, raw stdout/stderr and exit codes are retained in
the `cua-*.receipt.json` files. The successful browser recipe's clocks and step
observations are in `proposal-browser-2026-10-14T03-15-08.537268+00-00.receipt.json`.

## Evidence

- `receipt.json`: final scope, source pins, limits and evidence hashes.
- `proposal-inputs.json`: selected queue, proposal state and six-pixel membership.
- `proposal-test-event.json`: genuine browser-downloaded TEST event.
- `real-native.png`, `real-400-mask-hidden.png`: real queue native/display checks.
- `proposal-mask-800.png`, `proposal-hidden-800.png`: selected proposal rendering.
- `proposal-test-exported.png`, `proposal-native-exported.png`: fixed TEST event.
- `before.json`, `after.json`: real-packet preservation and final source identities.
- `ui-tests.receipt.json`: exact focused test command, clocks, output and exit.

## Tool failures and limits

The private Cua default X11 service returned no windows. Enabling its Wayland
backend discovered the exact private Chrome window and captured it, but AT-SPI
was unavailable (`org.a11y.atspi.Registry` had no owner). The typed browser bind
refused `browser_consent_required` with
`consumer_profile_endpoint_requires_grant`. The declared fallback was the
isolated agent-desktop browser API, verified with native Cua screenshots. No
physical desktop input or personal profile was used.

This continuation inspected the retained prior GUI ROI lineage and verified its
imported final artifact; it did not redraw or resave that ROI. The browser event
still requires the coordinator's independent CLI import/report check with the
proposal supplied. Biological mask correctness, camera response, preparation
equivalence and absolute genome-size calibration remain unproven by GUI tests.

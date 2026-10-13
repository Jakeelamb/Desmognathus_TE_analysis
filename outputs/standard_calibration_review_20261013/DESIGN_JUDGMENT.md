# Design judgment — 2026-10-13

**Recommend B as the implementation base; design review passes, runtime acceptance remains unqualified.** Both proposals respect the bounded retained-mask audit. B uses native imported helper functions with fewer execution mechanisms; its two additional import dependencies are an explicit, reviewable cost. No design parity defect found.

Scores: 2 = explicit and adequate for this bounded design, 1 = partial, 0 = absent. Scores do not certify execution.

| Requirement | A | B |
|---|---:|---:|
| Actual target/standard masks and reconstructed IOD parity | 2 | 2 |
| Explicit global-standard / historical-target coordinate provenance and parity | 2 | 2 |
| Protection against unrelated scientific/source changes | 2 | 2 |
| Exact existing algorithm reuse | 2 | 2 |
| Small one-file CLI and clean dependency boundary | 1 | 2 |
| Falsification, limits, and no claimed manual acceptance | 2 | 2 |
| **Total /12** | **11** | **12** |

A's unchanged AST-extracted bodies preserve intended formulas, but dependency/constant reconstruction adds a custom loader to maintain. B's full-tile `measure_repaired_iod` follows the documented identical logarithmic formula, including negative OD pixels; native pixel-helper calls preserve historical target fallback and float32 percentile behavior. Both correctly adapt only copied global standard coordinates, preserve original rows, and distinguish reconstruction from the unperformed eligibility/template/manual stages.

Two grafts from A:

1. Fix eight-metric target parity at absolute `1e-12` before execution; retain every delta and fail without tolerance relaxation.
2. Require standard crop coverage against complete label pixels, retain padding/clipping diagnostics, and prohibit display normalization from entering measurements.

Inspected final artifacts: `design_a/DESIGN.md`, `design_b/DESIGN.md`, `SOURCE_INVENTORY.md`, `MEASUREMENT_FLOW.md`, and `AUTHOR_IMAGE_CONFIRMATION.json`, all beneath this directory. The author confirmation resolves 40× brightfield/no-oil/green-channel/shared-fixation-protocol identification; it does not recover camera response or historical preservation details.

Exact read-only commands: `cat` on those five paths; `git status --short`; `git rev-parse HEAD` in the paper and sibling checkouts; `sha256sum` on the two designs and two grounding documents. All exited 0. Extensive existing paper changes were present and left untouched. Observed heads: paper `15dfb3c18755cc8659b626fa213087c62c6db02f`; sibling `d7c37fdaf2869245964b6710ca76696306c2daba`.

SHA256 pins:

```text
A 228c64bf3989a08ebf3ae585fa75bda4d708a0c90c9dd28e100d2a204d6fda09
B 4cd4b4e723a35e3d0eb05e21dd6178c38455e08b536526657be5d96beb536f0a
inventory 05b3855be42f2a39eefe01ec5fe3563e94ace634ff782fd673784e52d5f38d58
flow 3444131cec4591eb6a5003bbe1e185a04d8fa64dcb5b36d97b42fc9d631cd665
```

Missing evidence: neither artifact implements the proposed CLI. Candidate probe claims were not independently rerun in this design-only assignment. No runtime negative control was executed; B must demonstrate its missing-label, removed-pixel, altered-IOD, coordinate/hash, and existing-output refusals. Native imports/dependencies, full 47-row reconstruction, source-preservation guards and gallery correctness remain unproven. Four target controls alone cannot establish all-target parity. Frozen-template QC, accepted standard nuclei, acquisition linearity and absolute calibration remain outside this verdict. No broken tool observed.

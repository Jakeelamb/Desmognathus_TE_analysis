# Independent Methods editorial verification

**Verdict: PASS for the authorized editorial revision. No remaining material meaning change was found.** This is not a finding of publication readiness or validated absolute genome-size calibration.

Reviewed revision: dirty working tree at `15dfb3c18755cc8659b626fa213087c62c6db02f`, identified by the six hashes below, after the coordinator's final LF-normalization freeze. The six saved baseline files matched `BASELINE.json`. I read the baseline and actual revised prose, author-confirmation JSONs, relevant procedure implementations, and all eight rendered Word pages; writer and coordinator summaries were not the acceptance evidence.

## Preserved scientific meaning

- Panels remain TE34/Path24/overlap21 = 34/24/21, with no imputation or assumed animal matching. The two *orestes* technical runs are separately normalized and equally averaged; maximum SVL and the 52/50-mm selection rule remain intact.
- Top-ten positive superfamily selection, deterministic ties, reclosure, Shannon/Gini-Simpson equations, shared24 CLR-PCA, zero policy, centering, signs, command flags, settings, and reported software versions remain intact.
- LTR eligibility, inclusive within-species IQR fences calculated before coverage filtering, regional 80% gate, 408/381/407/380 counts, 30-species summaries, 2,000 element draws, and the 19-species/239-element regression remain intact. Elements remain neither independent model observations nor regression weights.
- Morphometry retains 1,152 pairs/42 animals/24 species, 22×50/44/8 allocation, reviewed upper-tail selection, 0.12/0.0144 scales, and conditional paired-object intervals. IOD retains 805 nuclei/51 images/50 animals/24 species, equal-image means of nucleus medians, 2,000 image-then-nucleus draws, excluded matching traits, and SMD 0.10/KS 0.25 flags.
- Full24 versus subset22, fixed full24 log10 scaling versus original within-fit scaling, maximum-SVL interpretation, paired fractions rather than ratios of aggregate summaries, median-of-animal-medians versus arithmetic-mean-of-animal-medians estimands, 1,000 animal/pair refits, and 250 animal/image/row draws remain distinct.
- The original 25 DAGs/11 classes/10 testable classes and 250 measurement/200 tree/800 simulation suite remain separate from the expanded eight-class comparison without measurement/tree bootstrap. Test orders, class weights, orientation ambiguity, saturated-model exclusion, interval methods, and omitted uncertainty remain explicit.
- Exact21 TE regressions retain full34 logit-TE/full24 cellular scaling, six fixed equations per endpoint, matched-cohort comparisons, 16 residual degrees of freedom, and exploratory/unadjusted status. All final Figure 1–8 procedure coverage remains present.

New facts agree with the retained author answers: one Feulgen run for both standards and all retained unknowns; standards 413/414, specimens 32469/32470, at 40× brightfield without oil, green channel, with the same preparation fixation protocol. The 53-row map identifies 51 unknown images and two standards. This does not establish lifetime preservation, exact fixation execution, camera response, or accepted standard nuclei. M16's stale run-map claim and the duplicate-tip removal wording were corrected; duplicate labels are rejected by the implementation.

Calibration still uses reviewed target Processes 314/315/320/321, not standards 413/414. The fixed 16.1/0.978 = 16.462167689 anchor remains conditional; the approximately 41% standard discrepancy, standard acceptance, and absolute C-value remain unresolved. No 47-mask audit was added as a final manuscript analysis.

## Executed checks and evidence

All commands used cwd `/home/jake/Projects/Desmognathus_TE`; Python used `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13`.

| Exact command / observation | Observed result |
|---|---|
| `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python -m pytest -q -p no:cacheprovider tests/test_methods_release.py` | Post-LF independent run, 18:25:13–18:25:14 UTC: exit 0, **12 passed in 0.37s**. Unchanged checks include overwrite refusal and out-of-sync rejection. |
| `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python scripts/publication/build_methods_document.py --check` | 18:24:23–18:24:24 UTC: exit 0, visible Word prose matches Markdown. Word was unchanged by the later CSV-only normalization; the post-LF suite also checks equality. |
| `git diff --check -- Publication/methods` | Post-LF independent run: exit 0, no output. |
| Independent inline Python comparison through the same UV invocation, `python -` (stdin retained in verifier tool transcript) | 18:26:44 UTC: exit 0. All 61 paragraphs match through an independent Markdown/XML comparison; only `word/document.xml` changed; six subscript-i, one superscript-2, and one superscript-minus-7 runs match baseline. No numeric-token or heading loss. Author-input IDs/statuses/requirements unchanged; only seven evidence cells changed. Only M16's status/source/disposition cells changed in the evidence matrix. |
| Negative control on a disposable copy of the actual final Word artifact | Changed “805 reviewed nuclei” to “804 reviewed nuclei”; the maintained checker raised `Author-edited Word document and Markdown companion are out of sync`. Original files were untouched. |
| Direct image inspection of `docx_render/page-1.png` through `page-8.png` | All eight pages readable; no observed clipping, missing glyphs, or broken equations. PDF SHA-256: `dfc286befdb2348bec6cef1af65b37bf5c386b3135dba7105eb71bd55598f8cd`. |
| Independent inline SHA-256 comparison against `outputs/submission_preparation_20261013/baseline.json` and `outputs/standard_calibration_review_20261013/FINAL_REVIEW_PINS.json` | 18:27:59 UTC: exit 0; 614 frozen figure-package files and 362 review files match, zero changes. |

Evidence paths are relative to `outputs/methods_deslop_20261013/` unless stated otherwise: `baseline/Publication/methods/`, `BASELINE.json`, `docx_render/`, `ROOT_EXECUTION_RECEIPTS.json`, and this report. The independent inline command payloads/raw results are in the verifier execution transcript, not separate evidence files.

The first freeze's diff check exited 2 because CSV serialization introduced CRLF; root normalized LF and the final check passed. One verifier inline assertion initially used assumed evidence-matrix column names and exited 1; correcting the observation to actual `evidence_status`/`manuscript_disposition` fields produced the successful audit above. No maintained acceptance checks or implementation were changed by this verifier. No broken tool remains.

## Final active file identities

| File under `Publication/methods/` | Bytes | SHA-256 |
|---|---:|---|
| `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md` | 28758 | `c2d4accc934f537c5c595d441ec83dca6da4fb550b02cf2d558532dc2ec375c9` |
| `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` | 32979 | `d84e2eeafcb1b6ecb6b1765d9594b6c4778e034de4ab2ac5d7800ee975749f40` |
| `METHODS_BULLET_SCAFFOLD.md` | 38508 | `441fb17ca89b4685f76fb049e1da44e7c80656f64547305d949e4c9152267b04` |
| `METHODS_EVIDENCE_MATRIX.csv` | 14475 | `588a858cbcb8f673f05c1b221c6d1e4788b3e3f4b014dc701432211f75a2e300` |
| `METHODS_REQUIRED_AUTHOR_INPUT.csv` | 8717 | `fe599b5457346afee17b1deae78af430d336d1a7e6ae2ef26b7646445dd76b1a` |
| `README.md` | 5857 | `dae84086e9be26d4124238052cc5e1a56ed1871cede13edb8e6229a5a4f639bf` |

Limits: this verifies editorial fidelity, synchronization, and observed rendering. It does not rerun scientific pipelines, establish the historical upstream runtime, certify calibration or causation, resolve ethics/deposition/author records, or guarantee any AI-detector outcome. Numeric-token checks alone cannot establish meaning; the semantic verdict relies on direct reading and scoped source inspection as well. No broad repository test suite was required or executed for this prose-only pass.

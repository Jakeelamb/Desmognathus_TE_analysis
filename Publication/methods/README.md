# Methods working materials

`METHODS_BULLET_SCAFFOLD.md` is the factual Methods authority. Its evidence labels, source links, and author-action items form the internal ledger. `METHODS_EVIDENCE_MATRIX.csv` maps procedures to retained sources and final figures. `METHODS_REQUIRED_AUTHOR_INPUT.csv` records unresolved specimen, laboratory, acquisition, tree, and upstream-run details.

`MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` is the author-edited Word draft. `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md` is its plain-text mirror and the source for the authorized prose edit and synchronization. The narrative preserves supported author prose and describes completed analyses in scientific paragraphs. It has no audit tags, internal dataset codes, repository paths, or author-action placeholders. The checklist retains missing facts.

The reconciliation covers the retained [eight-figure manuscript package](../../outputs/manuscript_polished_20260924_01/FINDINGS_AND_CAPTIONS.md). Its frozen [source run](../../outputs/manuscript_polished_20260924_01/inputs/source_run/README.md) contains the completed exploratory extensions. Source-run figure numbers differ from final manuscript numbers.

## Final eight-figure coverage

| Final figure | Methods coverage |
|---|---|
| 1 | TE34 repeat architecture, species phylogeny, overlap21 conditional genome-size availability |
| 2 | Repeat divergence; within-LTR low-divergence score; strategy and SVL regressions in 34 species; descriptive cross-pipeline QC |
| 3 | Top-ten-superfamily diversity and untransformed genome-size PGLS in 21 species |
| 4 | Primary 380-element LTR proxy in 30 species; log–log genome-size comparison in 19 species |
| 5 | Selected cell/nucleus morphometry; conditional allometry; paired nuclear fractions; 1,000 animal/pair resamples and allometry refits |
| 6 | Expanded eight body-role equivalence classes; full24/subset22 and fixed/reverse independence-test order; separate conditional regressions |
| 7 | Maximum-SVL associations; full24/subset22; animal-aware body-adjusted sensitivity |
| 8 | Exact21 TE-adjusted cellular regressions; fixed full24 cellular and full34 logit-TE scaling; alternative hit denominator |

The older controlled release under `Publication/figures` has seven main figures with different numbering, plus auxiliary composition displays. Its shared24 CLR-PCA and Bonsai analyses remain supporting context for the final eight-figure narrative. The original three-trait path suite retains 25 DAGs, 11 equivalence classes, ten testable classes, 250 measurement draws, 200 published bootstrap trees, and 800 simulations. Those sensitivities apply to the original suite. The later expanded eight-class comparison has no measurement or tree bootstrap. The saturated three-trait diagram is illustrative and was not ranked.

The animal-aware body analysis uses medians of animal morphology medians. Paired allocation and the allometry alternative use arithmetic means of animal medians. These are distinct estimands. Selected S26/final Figure5 pairwise coefficient intervals use t quantiles with 22 residual degrees of freedom. Supporting detailed S28 regression and S31 path-edge tables use estimate ± 1.96 SE. Retained extension intervals and fitted bands use residual-df t quantiles.

Conditional genome size uses the released *D. fuscus* equal-image species estimate from target Processes 314/315/320/321 as its denominator. The designated standards 413/414 were not used as that denominator. The author has confirmed one Feulgen run containing both standards and all retained unknowns, and identified the standard images as 40× brightfield without oil, using the green channel and the same fixation protocol as the other slides. The [calibration decision](../../outputs/standard_calibration_review_20261013/CALIBRATION_DECISION.md) retains the conditional assembly anchor because accepted standard QC and the approximately 41% between-slide median difference remain unresolved. The confirmations do not establish exact fixation reagents, lifetime preservation, photometric linearity, or an independent C-value. Tree and selected-cell sampling limitations also remain explicit.

## Word synchronization and submission review

Check Word synchronization with:

```bash
make methods-word
```

The builder refuses to replace an existing author-edited DOCX by default. The user authorized the current prose edit and synchronization from reviewed Markdown. The synchronization script checks the original Word file against the saved baseline, replaces changed paragraphs, and preserves the other Word package components. Word synchronization checks that the document matches its Markdown source. Rendering and visual inspection check the resulting pages. Detailed procedure-to-source evidence is retained in [the reconciliation report](../../outputs/submission_preparation_20261013/methods_writer_evidence.md). The [editorial package](../../outputs/methods_deslop_20261013/PLAN.md) preserves the baseline and review plan; [ledger edit notes](../../outputs/methods_deslop_20261013/NOTES_EDIT_NOTES.md) record this pass.

W15 records completed exploratory integration with author review pending. Recovered dnaPipeTE and two RepeatMasker command variants establish historical configurations; W10 and W11 retain the unresolved runtime versions, library bytes, and applied variants. W06 records the author-confirmed run map and remains partial because standard acceptance and absolute calibration are unresolved. W07 retains the corrected microscopy account and adds the green-channel standard confirmation. Source-tree calibration, camera settings, segmentation lineage, specimen/ethics, and upstream runtime details still require records before submission. Historical and redundant Word materials remain recoverable from Git history and the pre-cleanup archive.

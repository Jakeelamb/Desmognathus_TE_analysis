# Methods working materials

`METHODS_BULLET_SCAFFOLD.md` is the active factual Methods authority. Its evidence labels, source links, and author-action items belong to the internal ledger and are excluded from researcher prose. `METHODS_EVIDENCE_MATRIX.csv` maps procedures to retained sources and final figures. `METHODS_REQUIRED_AUTHOR_INPUT.csv` records unresolved specimen, laboratory, acquisition, tree, and upstream-run details.

`MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` is the author-edited Word draft. `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md` is its plain-text mirror and the reviewable source for the authorized reconciliation. The narrative preserves supported author prose and describes completed analyses as scientific paragraphs. It has no audit tags, internal dataset codes, repository paths, or author-action placeholders. Missing facts remain in the checklist and are not invented.

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

The older controlled release under `Publication/figures` has seven main figures with different numbering, plus auxiliary composition displays. Its shared24 CLR-PCA and Bonsai analyses remain supporting context. They do not define the final eight-figure narrative. The original three-trait path suite retains 25 DAGs, 11 equivalence classes, ten testable classes, 250 measurement draws, 200 published bootstrap trees, and 800 simulations. Those sensitivities do not qualify the later expanded eight-class comparison, which has no measurement or tree bootstrap. The saturated three-trait diagram is illustrative and was not ranked.

The animal-aware body analysis uses medians of animal morphology medians. Paired allocation and the allometry alternative use arithmetic means of animal medians. These are distinct estimands. Selected S26/final Figure5 pairwise coefficient intervals use t quantiles with 22 residual degrees of freedom. Supporting detailed S28 regression and S31 path-edge tables use estimate ± 1.96 SE; retained extension intervals and fitted bands use residual-df t quantiles. Conditional genome size uses the released D. fuscus equal-image species estimate as its denominator, not recovered run-matched standards. Anchor, staining-run, tree, and selected-cell sampling limitations remain explicit.

## Word synchronization and submission review

Verify visible prose with:

```bash
make methods-word
```

The builder refuses to replace an existing author-edited DOCX by default. An explicitly authorized synchronization from reviewed Markdown uses its `--replace` flag. Word synchronization, document rendering, and final artifact checks are separate from factual reconciliation. Detailed procedure-to-source evidence is retained in [the reconciliation report](../../outputs/submission_preparation_20261013/methods_writer_evidence.md).

W15 now records completed exploratory integration with author review pending. Historical dnaPipeTE and two RepeatMasker command variants have been recovered, so W10 and W11 distinguish known configurations from unverified actual versions, library bytes, and applied variants. Remaining source-tree calibration, target-to-standard staining-run matching, standard QC, acquisition, segmentation lineage, specimen/ethics, and upstream runtime details still require author records before submission. Historical and redundant Word materials remain recoverable from Git history and the pre-cleanup archive.

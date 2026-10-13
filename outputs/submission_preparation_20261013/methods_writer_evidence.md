# Methods reconciliation evidence

Reconciled 2026-10-13 against the final eight-figure package at `outputs/manuscript_polished_20260924_01`, its frozen `inputs/source_run`, and retained compact analysis implementations. Repository HEAD at dispatch was `15dfb3c18755cc8659b626fa213087c62c6db02f`; extensive existing author changes were preserved. No scientific analysis, refit, figure build, Word replacement, or scientific-data edit was performed by this worker.

## Changed files and ownership

- `Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md`: preserved the existing supported specimen/upstream/measurement prose; added retained diversity, low-divergence strategy, LTR19, paired-allometry, body/SVL, animal, expanded-class, and exact21 TE procedures. Removed obsolete unfitted-integration claims. Clarified implemented calibration denominator and uncertainty.
- `Publication/methods/METHODS_BULLET_SCAFFOLD.md`: updated evidence date; retained supporting PCA/Bonsai context separately from final figures; added source-linked procedures for all extensions, distinguished selected t and supporting normal-approximation intervals/scaling, and retained unresolved author/upstream facts.
- `Publication/methods/METHODS_EVIDENCE_MATRIX.csv`: retained M01–M20 with current dispositions and figure mappings; added M21–M26 for low-divergence strategy, LTR19, paired allometry, body-size models, animal robustness, and exact21 TE integration. All local primary-source paths exist.
- `Publication/methods/METHODS_REQUIRED_AUTHOR_INPUT.csv`: W15 reflects completed exploratory analyses with review pending. W10/W11 distinguish recovered historical configurations from unknown applied runtime variants, versions, and library bytes. W06 distinguishes designated standards from the implemented target-species denominator.
- `Publication/methods/README.md`: distinguishes eight-figure narrative from older seven-figure controlled release; summarizes coverage and remaining review boundaries.
- This report.

The parent owns Word synchronization, test changes, and final acceptance. Existing deleted/changed DOCX files in the working tree were not altered by this worker.

## Procedure-to-source crosswalk

Paths below are repository-relative. `S` abbreviates `outputs/manuscript_polished_20260924_01/inputs/source_run/` only within this report.

| Final figure / procedure | Actual retained source and code pins | Reconciled description and limits |
|---|---|---|
| 1: TE34 architecture and overlap21 conditional size | `S/scripts/figures_1_2.R`; `S/figures/figure_1/analysis/species_tracks.csv`; final `FINDINGS_AND_CAPTIONS.md` Fig.1 | Classified aligned-base proportions; species phylogeny distinct from composition/Bonsai tree; 34 rows, 21 genome estimates; no whole-genome occupancy claim. |
| Supporting CLR-PCA | `analyses/01_transposable_elements/recompute_pca.py`; existing scaffold and canonical S09–S11 tables | Shared24 positive superfamilies across34; retry averaging, closure, CLR, centered/no variance scaling, deterministic signs; supporting description separate from final eight figures. |
| 2: low-divergence definition and denominator identity | `S/scripts/te_integration.R:60–119`; `S/figures/figure_12/analysis/README.md`; `S/figures/figure_12/analysis/repeat_score_species34.csv` | LTR bp bins0–4/allLTR bp across all available bins; allhit alternative = withinLTR × totalLTR fraction; fixed threshold, no endpoint pseudocount, no genome-load/age/rate claim. |
| 2: strategy/SVL and QC | `S/scripts/te_integration.R:193–239`; same README | Four ML equations, lambda fitted separately; raw logit group contrasts vs aquatic; approximate chi-square LR2df;34 species, groups7/24/3; curves within observed ranges. Two predefined Spearman QC checks retained, unadjusted/post-hoc and cross-pipeline limitations explicit. Seed20260925 is set but no model resampling added. |
| 3: diversity vs genome21 | `S/scripts/figures_3_5.R:1–58`; `S/figures/figure_3/analysis/model_summary.csv` | Untransformed outcomes/genomepg; lambda[0,1]; t_(19) intervals, covariance bands, one full-sample fit per related metric; measurement intervals displayed but not propagated. |
| 4: primary380 cohort and interval reproducibility | `analyses/01_transposable_elements/recompute_ltr_release.py:273–317`; `analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json`; final caption Fig.4 |408 eligible;381 IQR/407 coverage;380 intersection/30species; zero-depth entire internal region;2000 species/branch/statistic seeded element bootstrap with SHA256 first16hex; central percentile intervals. |
| 4: LTR19 regression | `analyses/01_transposable_elements/ltr_genome_pgls.R:1–144`; `S/inputs/ltr_pgls/README.md`; `S/inputs/ltr_pgls/model_summary.csv` |19species/239elements; one unweighted median per species; loglog lambda[0,1], t_(17); specified BM/OLS/raw fits and19deletions; no measurement/tree propagation or calibrated deletion-rate interpretation. |
| 5: paired allometry and allocation | `S/scripts/nuclear_allometry.R:25–112`; `S/figures/figure_10/analysis/README.md`; `paired_animal_summaries.csv`, `paired_hierarchical_refits.csv` in same directory |Unstandardized log nucleus~log cell slope; slope1 t test; arithmetic means of animal medians; fraction within actual pairs. Four fits released/equal-animal ×24/22. Seed20260924;1000 animal→pair draws and1000refits/lambda re-estimation; conditional measurement ranges distinct from t confidence intervals;42animals/1152pairs/24species. |
| 6–7: body regressions and fixed scaling | `S/scripts/body_size_analysis.R:9,93–154,158–226`; `S/figures/figure_8/analysis/README.md`, `scaling.csv` |Species maximum SVL, not microscopy animal length. g~b,c~b,n~b;c~g,n~c,n~g;c~g+b,n~c+b,n~g+b;n~g+c;n~g+c+b. Full24 logmeans/sampleSDs fixed for subset22. Excludes exactly aeneus/wrighti; body/development/taxon changes confounded. Four upstream-body classes, original ten testable classes separately rerun. |
| 7: animal sensitivity | `S/scripts/animal_robustness.R:5,21–37,54–90`; `S/figures/figure_9/analysis/README.md`; `verification.json` |Morphology median-of-animal-medians, separate from allometry/fraction mean-of-animal-medians. IOD mean images within animal then mean animals. Seed2026092409,250 animal→image→row draws, independently resampled assays, paired morphology, all24, fixed scaling,24released-summary deletions. Four-model set only; cross-assay covariance and calibration/SVL/tree uncertainty omitted. |
| Original three-trait suite supporting6 | `analyses/03_phylogenetic_path/build_path_release.py:391–496`; `analyses/03_phylogenetic_path/run_path_analysis.R:32–39,180–204,691–805`; canonical `candidate_dags.csv`, `path_simulation_calibration.csv` |25DAGs/11classes/10testable;250 pooled-pair + image→nucleus measurement draws with recomputed fuscus denominator;200published trees on18species;24deletions;4covariance models;800simulations(100×4×2) seed20260710; species offsets documented in builder. Original standardization is per fitted dataset, not fixed full24 for every original sensitivity. |
| 6: expanded body-role equivalence classes | `S/scripts/body_role_analysis.R:23–139`; `S/figures/figure_13/analysis/METHODS.md`; `equivalence_classes.csv`, `model_comparison.csv` |Eight distinct classes; exact skeleton/unshielded-collider enumeration; H5four equivalent orientations. Full24/subset22; primary b,g,c,n and reverse n,c,g,b testing order, ancestral direction preserved; within cohort/order/set weights; direct structural equations fixed scaling; no expanded measurement/tree bootstrap. |
| Saturated diagram boundary | `S/scripts/figure_6_network_candidate.R:88`; `S/scripts/figure_6.R:280`; final `scripts/figure_6.R` |Fully connected three-trait diagram imposes no independence claim, was not ranked or assigned a scored class weight. Distinct from current H5 representative and the scored models. |
| 8: exact21 TE integration | `S/scripts/te_integration.R:120–193,248–258`; `S/figures/figure_11/analysis/README.md`; `regression_coefficients.csv`, `te_scaling.csv` |WithinLTR and allhit endpoints each fit g~t;g~t+b;c~g+b;c~g+t+b;n~c+g+b;n~c+g+t+b. Full34 logit TE means/sampleSDs, full24 logcellular/SVL scaling. Compare24,matched21,TEadjusted21; no imputation or significance-based denominator/model selection. Nuclear4predictors+intercept df16; conditional diagnostic rather than selected causal network. |
| Interval distinction | `analyses/03_phylogenetic_path/run_path_analysis.R:384–385,429–430`; `S/scripts/figures_3_5.R:156–162`; final `scripts/figures_5_7.R:88`; extension regression helpers |Selected frozen Fig5 model_summary (S26) brackets are t_(22), halfwidth/SE2.0738731, and fitted bands use t. Supporting detailed S28/S31 intervals are±1.96SE; their coefficient tests use t. Later body,TE,allometry,expanded structural intervals use qt(.975,n−p). Initial conflation was corrected; no caption/art change required. |
| Implemented calibration | `analyses/02_morphology/validate_release.py:135–170`; `iod_source_manifest.json:73–103`; `outputs/submission_preparation_20261013/provenance_recovery.md` |Released fuscus equal-image species estimate from target314/315/320/321 is denominator, not designated standards413/414. Anchor16.1/0.978pg; draw-wise target/reference ratio. Unknown run matching, approximately41% standard median difference, anchor/linearity/acquisition and staining uncertainty persist. |

## Remaining ambiguities and retained limits

- No model refits or data changes were required to reconcile the reported procedures. Numerical equivalence-class rankings and model conclusions remain retained outputs, not newly established scientific conclusions.
- The focal tree is author-supplied and its exact computational bytes are retained. Publication/archive identity, calibration method, units, and original source provenance remain unqualified.
- Designated standard identity and same-experimental-run staining were previously author-confirmed. The implemented denominator uses finalized target-species IOD, not those standards or a recovered run-match. No standard QC or acquisition comparability claim was added.
- Initial source tracing conflated detailed legacy S28/S31 normal-approximation intervals with selected S26/Figure5 intervals. The parent resolved the selected source: all three frozen Figure5 model-summary rows have CI halfwidth/SE=2.0738731 with df22, confirming t_(22). Final Figure5 caption is correct; no caption/art/ZIP correction is needed. Narrative and scaffold now explicitly separate selected S26/Figure5 t intervals from supporting S28/S31 normal-approximation tables.
- User's newly pasted historical draft mentions fastp, DeepTE/TEClass2, 60× oil IOD imaging, and additional specimen/slide sources. These have not been matched to the retained releases or confirmed as applied methods. Narrative reflects subsequent author clarification: all retained morphometry and IOD images are brightfield40×; poor40× oil-immersion reimages were excluded. The prior60× claim was an error. Pyron smear preservation remains unconfirmed; the general all-animals ethanol assertion was removed. No historical pipeline or60× claim was promoted. Author confirmation can refine these facts later.
- Historical dnaPipeTE configuration and two RepeatMasker variants are known. Actual applied variant, versions, container/library bytes, job/output mapping, and remaining read preprocessing remain unresolved. The lab/ethics/segmentation/checkpoint gaps remain in the checklist.
- Expanded eight-class models do not inherit original three-trait250measurement/200tree/800simulation qualification. Animal robustness qualifies the original four-body hypotheses only.

## Validation record

Checks are against the final text artifacts. No scientific refits were invoked.

- Initial structural scan: **FAIL**, solely because this report had not yet been created and the new README link pointed to it. The same scan found no missing matrix paths, forbidden narrative tokens, or obsolete integration statements. Final repeat below resolves that missing-report dependency.
- An initial `uv run --frozen` attempted to hydrate missing local Python3.15 dependencies and was interrupted with exit130 before the edit script ran. Subsequent pure-stdlib edits used `uv run --no-project python`; final checks use the parent's hydrated Python3.13 environment. This interruption was an environment/setup failure, not a scientific check.
- Final structural/path/prose checks are recorded below after completion.

## Supporting PCA statistics correction after independent review

Independent review identified stale supporting scaffold statistics. The current canonical `analyses/01_transposable_elements/data/te_pca_trait_association_tests.csv` supplies the corrected values. No analysis was rerun, and the researcher narrative was unchanged.

| Test ID | Current canonical cohort and retained values |
|---|---|
| `reproductive_strategy_permanova` | n=34; R²=0.2601016805953517; P=1e-5; adjusted P=1e-5 |
| `reproductive_strategy_permdisp` | n=34; R²=not applicable; P=0.70663; adjusted P=0.70663 |
| `direct_development_vs_aquatic_larvae_permanova` | n=34; R²=0.12030625975602982; P=3.5e-4; adjusted P=3.5e-4 |
| `egg_environment_among_aquatic_larvae_permanova` | n=31; R²=0.16764327291205847; P=1e-5; adjusted P=2e-5 |
| `size_adjusted_model_reproductive_strategy` | n=34; R²=0.19904284290161434; P=1e-5; adjusted P=1e-5 |
| `size_adjusted_model_log10_max_svl_mm` | n=34; R²=0.020444419162724357; P=0.55115; adjusted P=0.55115 |

The scaffold now reports all34 for the omnibus/direct-development contrast/size model and31 for the aquatic-versus-terrestrial contrast, with canonical rounded eight-decimal R² values. M07 directly links the canonical supporting source. Its SHA-256 is `04c521b17df1bf39ba594cc63018b293daf4f79dae68724d186bd319163f9332`.

## Final figure-export scope correction

The scaffold now labels Nimbus Sans, older size targets, and 300-ppi sRGB TIFF as the older seven-figure controlled release. The final eight-figure settings are embedded Arial,180mm collage width, at most245mm height,7.5pt working minimum, vectorPDF,300dpiRGBPNG and600dpiCMYKTIFF. Sources: `outputs/manuscript_polished_20260924_01/README.md:9`, `scripts/common.R:12–13,49–55` and `scripts/package_figures.py:18` within that run. M20 links these final sources. No art or researcher narrative was changed.

Final command: `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python` with pure-stdlib checks for every matrix source path and scaffold/README Markdown link, canonical supporting statistics, absence of superseded statistics, explicit export scopes, and unchanged researcher narrative SHA-256. Exit0. Earlier final structural/prose and selected Figure5 t22 checks passed before these supporting-only corrections.

```json
{
  "status": "PASS",
  "evidence_rows": 26,
  "unique_existing_local_source_paths": 78,
  "canonical_supporting_statistics": "reconciled from CSV rows by test_id",
  "figure_export_scopes": "older seven-figure and final eight-figure explicitly distinguished",
  "narrative_unchanged": true,
  "file_pins": {
    "MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md": "1573733a1577413881e00d63ef5a7da0832673dd53766bd58a271cc88f5cce72",
    "METHODS_BULLET_SCAFFOLD.md": "6179647455d06aa52e12a0a333c8823e46bee65790900077881cfc85c508c479",
    "METHODS_EVIDENCE_MATRIX.csv": "92b7a856b8c6a27a8bbbb94f5801fc5e9fae3bf28dc2376784ac1667af57f46b",
    "METHODS_REQUIRED_AUTHOR_INPUT.csv": "2036411b11d812ea83d8d1815e05b1eb02b8db98402f05c0ee25acd07759f1c3",
    "README.md": "6e90716dbeba8e270c425e068ce0c41e9a1ee5045d1a15b207d995d242ae7fdd"
  }
}
```

## Coordinator integration after writer verification

Independent review prompted a four-word qualifier limiting heart sampling to specimens sampled directly. Word/Markdown were resynchronized and final focused tests passed. The evidence matrix was normalized from CRLF to LF only. Earlier writer pins above identify the bytes checked at those stages; FINAL_PINS.json records the integrated deliverables. No scientific data or figures changed.

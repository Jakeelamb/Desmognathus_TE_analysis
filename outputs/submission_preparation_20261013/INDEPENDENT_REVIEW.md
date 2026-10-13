# Independent Methods and contract review

Reviewed on 2026-10-13. The repository revision is `15dfb3c18755cc8659b626fa213087c62c6db02f` with preexisting author changes and the current reconciliation. File hashes below identify the reviewed working revision. This is a separate native reviewer, not a cross-model-family review.

**Verdict.** PASS for the bounded scientific reconciliation, final Word consistency/layout, and repaired contract safeguards. The final revision passes execution and rejection controls. Submission readiness remains unqualified because calibration, acquisition, specimen, tree, upstream-run, and deposition evidence is incomplete. The formatting finding was resolved and rechecked. The excluded older auxiliary artwork remains a documented mismatch.

I changed no implementation, acceptance check, scientific data, Word draft, or figure. This report is my only assigned output. I applied Prove It Works, whose description says to "Verify against the real artifact." That changed the review from reading test summaries to inspecting the actual Word pages and retained auxiliary artwork.

## Findings and disposition

1. **Corrected during review.** The factual scaffold initially retained stale supporting PCA strategy numbers for 33 species. Its lines 63–65 now agree with `analyses/01_transposable_elements/data/te_pca_trait_association_tests.csv`. The omnibus test uses 34 species, R² 0.26010168, and PERMDISP P 0.70663. The direct-development contrast uses 34 species, R² 0.12030626, and Holm-adjusted P 0.00035. The egg-placement contrast uses 31 species, R² 0.16764327, and adjusted P 0.00002. The adjusted size comparison uses 34 species and SVL P 0.55115. The narrative already described the current panel correctly. The writer corrected the scaffold and linked the actual table. Confidence is high.
2. **Corrected during review.** The scaffold initially described Nimbus Sans and 300-ppi sRGB TIFFs as the unqualified publication style. Lines 165–166 now separate that older seven-figure release from the selected eight-figure package. The latter uses Arial, 180-mm collage width, maximum 245-mm height, minimum working text of 7.5 pt, 300-dpi RGB review PNGs, vector PDFs, and 600-dpi CMYK TIFFs. I checked the final package README and renderer specification. I did not repeat the coordinator's entire artwork audit. Confidence is high for the scope correction.
3. **Retained legacy-artifact mismatch.** `Publication/figures/Path_DAG_equivalence_class_reference.pdf` and its PNG still show "No directed paths" in the Independent panel. The current runtime contract requires only the 30 node labels across ten panels and three nodes with no edges in the Independent panel. I extracted text from the actual PDF and viewed the actual PNG. Thus the passing runtime check does not prove that this retained older auxiliary artwork matches the current renderer. The scientific caption caveats are present in the PDF. The coordinator explicitly excludes this older auxiliary from the selected polished eight-figure submission inventory and instructed that it remain unchanged. Confidence is high. Do not call this old artifact refreshed.
4. **Formatting gate initially failed.** The scoped `git diff --check` exited 2 because `METHODS_EVIDENCE_MATRIX.csv` had CRLF endings on all 27 lines. Its preserved baseline had LF endings. The author-input CSV remained LF. This was introduced by reconciliation and was reported to the coordinator for correction. It is unrelated to scientific correctness. Final disposition is recorded in the addendum.
5. **Collection scope corrected before acceptance.** The original direct heart-sampling sentence could be read as applying to supplied dried smears. The coordinator asked for a scope review and adopted the minimal qualifier, "For specimens sampled directly, erythrocytes were collected..." in narrative line 57. The final Word paragraph uses the same qualifier. This preserves the confirmed protocol while leaving supplied-smear preparation and retained-cohort mapping unresolved. It does not assert that Pyron smears were retained. Confidence is high for the reporting correction; the missing preparation evidence remains missing.

## Consequential source checks

`S` below means `outputs/manuscript_polished_20260924_01/inputs/source_run/`. The source-run figure numbers differ from the final eight-figure numbering.

| Requirement | Evidence read or executed | Observed conclusion |
|---|---|---|
| Actual conditional genome-size denominator | `analyses/02_morphology/provenance/iod_source_manifest.json:73–103`; `analyses/02_morphology/validate_release.py:140–169`; actual `nuclear_iod_objects.csv` filenames; narrative line 75 | The denominator is the finalized target-species D. fuscus IOD estimate from Processes 314, 315, 320, and 321. Designated standards 413/414, specimens 32469/32470, are not its current denominator. The prose distinguishes these and does not claim independent C-values. Target/standard staining-run mapping, standard QC, the approximately 41% difference, acquisition linearity, and anchor uncertainty remain unresolved. |
| Retained objective and image branch | Author clarification in this session; narrative line 57; scaffold line 96; actual rendered Word page 4 | Retained morphometry and IOD use brightfield 40× images. Poor 40× oil reimages were excluded. The historical 60× claim and blanket ethanol-fixation claim are absent. Pyron smear preservation remains unresolved in W01/W03. This review establishes correct reporting of the author account, not an independent microscope-log reconstruction. |
| Measurement and biological replication | Actual `cell_nucleus_objects.csv` and `nuclear_iod_objects.csv`; `S/scripts/animal_robustness.R:17–30`; narrative lines 63–73 | 1,152 paired morphology objects from 42 animals; 805 nuclei from 51 images and 50 animals. These are distinct assay cohorts. Species are regression observations; objects and images are not independent biological animals. Selected large-cell sampling is explicit. |
| Selected versus detailed interval methods | `S/figures/figure_5/analysis/model_summary.csv`; final `scripts/figures_5_7.R:87–88`; `analyses/03_phylogenetic_path/run_path_analysis.R:384–385,429–430`; narrative line 83 | All three frozen selected pairwise rows have 22 residual df and CI halfwidth/SE approximately 2.0738731. Project R gives t(.975,22) = 2.073873067904026. The selected Figure 5 t label is correct. Supporting detailed regression and path-edge tables use ±1.96 SE. No caption or artwork change was needed. |
| Allometry and allocation aggregation | `S/scripts/nuclear_allometry.R:33–45,48–94`; actual `paired_hierarchical_refits.csv` | Fractions are computed within actual pairs, then animal medians and species arithmetic means. The allometry alternative also uses arithmetic means of animal medians. There are 1,000 completed animal-then-pair draws and refits. The raw log–log slope and slope-one test are distinct from standardized pairwise coefficients. |
| Animal-aware body sensitivity | `S/scripts/animal_robustness.R:21–39,54–74`; actual `bootstrap_model_comparison.csv` | Morphology uses the median of animal medians. IOD uses image means within animal, then means across animals. There are 250 retained draws and four model rows per draw. Assays are independently resampled; cross-assay covariance is unmodeled. Full24 scaling is retained. These runs qualify the original four-body set, not the expanded eight classes. |
| Body/SVL equations and reduced cohort | `S/scripts/body_size_analysis.R:93–125,137–173`; narrative body-size section | Full24 and fixed subset22 differ only by excluding aeneus/wrighti. Body size is species maximum SVL rather than measured animal length. Scaling is fixed to full24, and the exclusion changes size range, development and taxon membership together. The described eleven equation specifications match the retained implementation. |
| Expanded body-role classes and test direction | `S/scripts/body_role_analysis.R:23–85,88–130`; actual `equivalence_classes.csv`, `model_comparison.csv`; narrative path-hypotheses section | Eight distinct classes, 32 ranking rows across two cohorts and two test orders. The cell-linked-SVL class has four graphical equivalent orientations. Primary order is b,g,c,n; reverse order is n,c,g,b. No saturated model or expanded measurement/tree bootstrap is claimed. Weights have only within-set/cohort/order meaning. |
| Original path sensitivity scope | `analyses/03_phylogenetic_path/build_path_release.py:436–496`; `run_path_analysis.R:32–39,686–805`; original candidate inventory | The original three-trait comparison has 25 DAGs, 11 classes, ten testable classes; 250 measurement datasets; 24 species deletions; 200 published trees on the 18-species overlap; and 800 simulations from 100 replicates × four generating classes × two residual regimes. The narrative does not transfer these qualifications to the expanded eight-class set. I read the retained implementation, not a new simulation run. |
| Exact21 TE integration | `S/scripts/te_integration.R:60–108,120–146,193–231`; actual score and coefficient tables | 34-species logit TE scaling; fixed full24 cellular scaling; exact21 join without imputation; six equations per endpoint. The most complex nuclear model has four predictors plus intercept and 16 residual df. The low-divergence all-hit fraction equals within-LTR fraction × LTR share. The alternative denominator is an explicit sensitivity, not a selected significant result. |
| Diversity21 and LTR19 inference | `S/scripts/figures_3_5.R:21–49`; `analyses/01_transposable_elements/ltr_genome_pgls.R:48–110`; `recompute_ltr_release.py:273–317` | Diversity is untransformed with lambda [0,1] and t19 intervals; LTR is log–log with lambda [0,1] and t17 intervals. LTR uses one unweighted species median, not element-count weighting. The bootstrap is element resampling with SHA-256-derived seeds and percentile intervals. Neither endpoint is a measured DNA-loss rate. |
| Evidence and checklist integrity | Actual Methods CSVs and README; direct path scan | 26 contiguous matrix rows and 78 distinct existing local source paths. W06 is partial. W15 is completed exploratory analysis with author review pending. Historical configured commands do not establish runtime versions, library bytes or applied execution variants. |

## Contract and Word review

The exact baseline comparison was against `outputs/submission_preparation_20261013/baseline/tests/`, not the full dirty diff from HEAD. The repository repair allows isolated review outputs while preserving other layout bans. It checks canonical source mappings, additional manifest dependencies and catalog paths, resolves paths to catch links into `outputs`, and keeps current artwork outside frozen input directories. Its negative control rejects an existing output source promoted into a canonical mapping. It does not independently verify asset content, checksums or scientific provenance.

The DAG repair moves case-sensitive source-literal caveats into an executed ggplot-caption check and checks S04's 24 conditional IOD-derived trait definitions. Existing candidate counts, representative classes, auxiliary placement, format existence and Independent-panel geometry checks remain. The negative control removes the C-value caveat from the executed caption and the unchanged test fails. This remains a renderer contract; it does not establish existing legacy-export freshness.

The Methods test changes retain the forbidden historical-path, audit-tag, unsupported pipeline, unique-Word-draft and overwrite protections. They extend evidence coverage from the old fixed 20 rows to contiguous IDs with M01–M26 required. The new heading test reads the generated document as a user would and checks visible text and Heading 2. It is an appropriate small behavior test. The builder recognizes `###` as a subsection, emits Heading 2, and checks that kind when comparing the document with Markdown. Default overwrite refusal remains before document creation. I found no weakened safeguard in these scoped changes.

The final DOCX matches Markdown. Of 40 baseline paragraphs, 29 are preserved exactly in the 61-paragraph final draft. Only `word/document.xml` changes in the ZIP package; all other package members are byte-identical. I inspected pages 1–8 of the coordinator's bundled-LibreOffice render at original image resolution. All are readable, equations and superscripts/subscripts are present, page numbers are visible, and there is no observed clipping, overlap, blank page or literal `###` heading. My separate default overwrite attempt on the actual author draft exits 1 and preserves its hash. An in-memory altered final paragraph is rejected by the synchronization check.

## Commands and results

All project Python commands used the coordinator's hydrated UV environment. No additional dependency installation, expensive producer, model suite, artwork build or scientific refit was performed by this reviewer. The full-suite execution belongs to the coordinator and is not claimed here.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python -m pytest -q tests/test_repository.py tests/test_path_dag_reference.py
```

Exit 0. Complete stdout was `13 passed in 5.91s` after thirteen dots.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python -m pytest -q -p no:cacheprovider tests/test_methods_release.py tests/test_repository.py tests/test_path_dag_reference.py
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python scripts/publication/build_methods_document.py --check
scripts/run_r.sh Rscript -e 'cat(sprintf("%.15f\n", qt(.975,22)))'
pdftotext -layout Publication/figures/Path_DAG_equivalence_class_reference.pdf -
```

Observed exits were 0, 0, 0 and 0. The tests reported `25 passed in 6.24s`, with no skips. The document check reported that visible DOCX prose matches the companion Markdown. R printed `2.073873067904026`. Actual PDF text includes `No directed paths` and the two scientific caveat sentences. The PNG was also viewed directly with `view_image`.

The first metadata comparison attempted `from scipy.stats import t` and exited 1 with `ModuleNotFoundError: No module named 'scipy'`. I did not install SciPy. I obtained the quantile from the maintained project R environment and reran the metadata comparison successfully. An initial lookup of a guessed `analyses/02_morphology/build_morphology_release.py` path failed because that file does not exist; the actual active implementation is `validate_release.py`, which I read. These were observation/setup issues, not failed scientific fits.

The retained-table and source-path comparison used the following assertions. Exit 0. This reads results and checks arithmetic; it does not fit models.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python - <<'PY'
from pathlib import Path
import pandas as pd
root = Path.cwd()
source = root/'outputs/manuscript_polished_20260924_01/inputs/source_run'
rows = pd.read_csv(source/'figures/figure_5/analysis/model_summary.csv')
ratios = (rows.pgls_ci_high - rows.pgls_ci_low)/(2*rows.pgls_standard_error)
assert len(rows) == 3 and rows.pgls_degrees_freedom.eq(22).all()
assert (abs(ratios - 2.073873067904026) < 1e-6).all()
assert (abs(ratios - 1.96) > .1).all()
morph = pd.read_csv(root/'analyses/02_morphology/data/cell_nucleus_objects.csv')
iod = pd.read_csv(root/'analyses/02_morphology/data/nuclear_iod_objects.csv')
assert len(morph) == 1152 and morph.specimen_id.nunique() == 42
assert len(iod) == 805 and iod.specimen_id.nunique() == 50 and iod.filename.nunique() == 51
fuscus = iod[iod.species.eq('D. fuscus')]
assert {x.split('_')[1] for x in fuscus.filename} == {'314','315','320','321'}
allometry = pd.read_csv(source/'figures/figure_10/analysis/paired_hierarchical_refits.csv')
assert len(allometry) == 1000 and allometry.fit_status.eq('success').all()
animal = pd.read_csv(source/'figures/figure_9/analysis/bootstrap_model_comparison.csv')
assert animal.id.nunique() == 250 and animal.model.nunique() == 4 and len(animal) == 1000
expanded = pd.read_csv(source/'figures/figure_13/analysis/model_comparison.csv')
equivalence = pd.read_csv(source/'figures/figure_13/analysis/equivalence_classes.csv')
assert len(equivalence) == 8 and len(expanded) == 32
assert equivalence.set_index('model').loc['cell_body_branch','equivalent_orientations'] == 4
te = pd.read_csv(source/'figures/figure_11/analysis/regression_coefficients.csv')
assert te.n_species.eq(21).all()
assert te[te.fit_id.eq('nucleus_te_adjusted')].residual_df.eq(16).all()
score = pd.read_csv(source/'figures/figure_12/analysis/repeat_score_species34.csv')
assert len(score) == 34
assert (abs(score.low_ltr_fraction_all_hits-score.within_ltr_fraction*score.ltr_fraction_all_hits) < 1e-15).all()
assert score.within_ltr_fraction.between(0,1,inclusive='neither').all()
methods = root/'Publication/methods'
evidence = pd.read_csv(methods/'METHODS_EVIDENCE_MATRIX.csv',keep_default_na=False)
paths = set()
for value in evidence.primary_sources:
    for path in value.split(';'):
        path = path.strip()
        if path and not path.startswith(('doi:','http://','https://')):
            assert (root/path).is_file(), path
            paths.add(path)
assert len(evidence) == 26 and len(paths) == 78
scaffold = (methods/'METHODS_BULLET_SCAFFOLD.md').read_text()
assert '33 coded TE34' not in scaffold and '0.2624' not in scaffold
assert all(v in scaffold for v in ['0.26010168','0.12030626','0.16764327','0.55115'])
assert 'older seven-figure' in scaffold.lower()
assert all(v in scaffold for v in ['Arial','600','CMYK'])
print('PASS: retained tables, CI arithmetic, denominator IDs, cohorts, completed draws, source paths and corrected scaffold')
PY
```

The reproducible negative and artifact comparison commands are retained below. They modify only in-memory objects, not implementation files or the final artifacts.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python - <<'PY'
import importlib.util
import json

def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

repository = load('review_repository', 'tests/test_repository.py')
original = repository.SUPPLEMENT_SOURCES
repository.SUPPLEMENT_SOURCES = {'review_negative_control':
    'outputs/manuscript_polished_20260924_01/inputs/source_run/scripts/te_integration.R'}
try:
    repository.test_isolated_review_outputs_remain_separate_from_canonical_sources()
except AssertionError as error:
    print(json.dumps({'negative_case': 'canonical source redirected into outputs',
                      'detected': True, 'assertion': str(error)}))
else:
    raise RuntimeError('Invalid canonical source was accepted')
finally:
    repository.SUPPLEMENT_SOURCES = original

dag = load('review_dag', 'tests/test_path_dag_reference.py')
original_run = dag.subprocess.run
observed = {}
def altered_run(command, **kwargs):
    command = list(command)
    code = command[-1]
    marker = 'caption <- tolower(Path_DAG_reference$labels$caption)'
    assert marker in code
    code = code.replace(marker,
        'Path_DAG_reference$labels$caption <- gsub("not a direct species C-value", '
        '"a species C-value", Path_DAG_reference$labels$caption, fixed=TRUE)\n' + marker, 1)
    command[-1] = code
    result = original_run(command, **kwargs)
    observed.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    return result
dag.subprocess.run = altered_run
try:
    dag.test_independent_dag_keeps_three_nodes_without_arrows_or_extra_annotation()
except AssertionError:
    assert observed['returncode'] != 0
    assert 'not a direct species c-value' in observed['stderr']
    print(json.dumps({'negative_case': 'remove direct-C-value caveat from executed ggplot caption',
                      'detected': True, **observed}))
else:
    raise RuntimeError('Missing scientific caveat was accepted')
finally:
    dag.subprocess.run = original_run
PY
```

Outer exit 0. Both controls printed `detected: true`. The canonical source rejection identifies the redirected output path. The executed R negative exits 1 with `Error: grepl("not a direct species c-value", caption, fixed = TRUE) is not TRUE` and `Execution halted`.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python - <<'PY'
from pathlib import Path
from zipfile import ZipFile
import hashlib, json, subprocess
from lxml import etree
from docx import Document
from runpy import run_path
root = Path.cwd()
current = root/'Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx'
baseline = root/'outputs/submission_preparation_20261013/baseline/Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx'
with ZipFile(baseline) as old, ZipFile(current) as new:
    assert set(old.namelist()) == set(new.namelist())
    changed = [p for p in old.namelist() if old.read(p) != new.read(p)]
    assert changed == ['word/document.xml'], changed
    ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    op = etree.fromstring(old.read('word/document.xml')).findall('w:body/w:p', ns)
    np = etree.fromstring(new.read('word/document.xml')).findall('w:body/w:p', ns)
    unchanged = sum(etree.tostring(p) in {etree.tostring(q) for q in np} for p in op)
    assert unchanged == 29
before = hashlib.sha256(current.read_bytes()).hexdigest()
result = subprocess.run(['uv','run','--frozen','--no-sync','--python','3.13','python',
    'scripts/publication/build_methods_document.py'], capture_output=True, text=True)
assert result.returncode == 1 and 'Refusing to overwrite an author-edited Word document' in result.stderr
assert hashlib.sha256(current.read_bytes()).hexdigest() == before
builder = run_path(root/'scripts/publication/build_methods_document.py')
modified = Document(current)
modified.paragraphs[-1].text = 'Deliberately altered prose for an in-memory negative control.'
check = builder['check_document_matches_source']
original = check.__globals__['Document']
check.__globals__['Document'] = lambda path: modified
try:
    check(root/'Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md', current)
except RuntimeError as error:
    assert 'out of sync' in str(error)
    mismatch_detected = True
else:
    raise RuntimeError('In-memory altered Word prose was accepted')
finally:
    check.__globals__['Document'] = original
print(json.dumps({'changed_word_package_parts': changed,
    'baseline_and_final_paragraphs': [len(op), len(np)],
    'unchanged_baseline_paragraphs': unchanged,
    'actual_word_overwrite_rejected_exit': result.returncode,
    'actual_word_bytes_preserved': True, 'in_memory_word_mismatch_detected': mismatch_detected,
    'docx_sha256': before}, indent=2))
PY
```

Outer exit 0. Observed package difference was only `word/document.xml`; paragraph counts were 40 and 61; 29 old paragraphs were preserved; overwrite rejection exit was 1; both negative controls passed. The actual DOCX hash was `f3fc7fc69dc81ac187e92d74377fc1adbfe5dd65eaf419d13a5a8ce5a758db5c`.

## Evidence paths and reviewed pins

- This report contains the independent command results and rejection controls.
- `outputs/submission_preparation_20261013/baseline/` contains the pre-reconciliation Methods and checks used for comparison.
- `outputs/submission_preparation_20261013/docx_render/` contains the initial eight PNGs and rendered PDF I inspected. The coordinator's initial render receipt is `checks/word-render.json`; its recorded command exit is 0 and log SHA-256 is `88d5740a373cb8b27148354d14015e541f4d016e6769d16a7c751a008cc911cb`.
- `outputs/submission_preparation_20261013/docx_render_final/` contains the actual final eight-page render after the sampling qualifier. I independently compared all eight PNG byte streams. Only page 4 changed. I viewed final page 4 directly at original resolution; it is clean and contains the qualifier. The other seven pages are byte-identical to those I inspected. `checks/word-render-final.json` records the final rendering.
- `Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md`: `2711baa93e0a3535946a7519fd2ee0497211bff35e79c0b304a20d7a63e31173`.
- `Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx`: `a827f7709334121637d9035414c560ce41514d46b75c5dfedf525d3114e00286`.
- `Publication/methods/METHODS_BULLET_SCAFFOLD.md`: `6179647455d06aa52e12a0a333c8823e46bee65790900077881cfc85c508c479`.
- `Publication/methods/METHODS_EVIDENCE_MATRIX.csv`: `d6c18ce6da8326e80105a0928ed0d177fea8bc741c1bb91f42080b0e9357d8a0`.
- `Publication/methods/METHODS_REQUIRED_AUTHOR_INPUT.csv`: `2036411b11d812ea83d8d1815e05b1eb02b8db98402f05c0ee25acd07759f1c3`.
- `Publication/methods/README.md`: `6e90716dbeba8e270c425e068ce0c41e9a1ee5045d1a15b207d995d242ae7fdd`.
- `scripts/publication/build_methods_document.py`: `3766a66ae5e58a023598b2434b28be5bd7d893d2c1b5a61292719c6e9615febf`.
- `tests/test_methods_release.py`: `8a4342150430d6c17210cd3efaa32c411166ed9db59b726c4d3161b9f19815a5`.
- `tests/test_repository.py`: `2a2b9a926aff7ce62fa88ed75f0d492a68aebd136b3be6f8842173eb92c9a175`.
- `tests/test_path_dag_reference.py`: `75fe0ba5d7ead8f972da85fe63a628a0b941745bc66517482041f73e440e62c6`.

## What remains unproven

Correct prose does not validate the absolute genome-size scale or establish a batch-corrected relative scale. The final model coefficients remain conditional on the unresolved staining calibration. The standard-target run map and standard remeasurement/QC are the next consequential evidence requirement.

The author-input checklist also retains specimen/voucher and ethics details, acquisition and sampling settings, Pyron smear preservation, focal-tree publication/calibration/branch-length provenance, segmentation checkpoint-to-mask attribution, actual upstream execution/library provenance, final author review, and public persistent deposition. Missing local records do not prove that these facts or records do not exist elsewhere.

This review does not establish causal direction, DNA-loss or recombination rates, genome occupancy from aligned-repeat fractions, random population sampling of erythrocytes, complete errors-in-variables inference, or combined tree/measurement/calibration uncertainty. It does not validate an assembled full manuscript, final bibliography, cover letter, author disclosures, or journal submission bundle. The full suite and final eight-figure artwork package are separately owned checks. No release or submission authority follows from the bounded PASS.

## Final addendum

All final hashes above were read from the actual files after the collection-scope correction. The final independent rerun of the same 25-test command exited 0 and reported `25 passed in 5.65s`, with no skips. Final Word `--check` exited 0. The retained-table/source-path script and ZIP/Word negative-control script above were executed directly from this report against the final files and both exited 0. The final Word negative-control output reports hash `a827f7709334121637d9035414c560ce41514d46b75c5dfedf525d3114e00286`; the earlier `f3fc...` output above is the retained pre-qualifier observation.

The matrix normalization was independently verified, not just read from the coordinator's receipt. Replacing each LF byte sequence in the final matrix with CRLF reconstructs the exact previously observed hash `92b7a856b8c6a27a8bbbb94f5801fc5e9fae3bf28dc2376784ac1667af57f46b`. The final matrix contains no CR bytes and 27 LF-terminated lines. Thus the cleanup changed only line endings. `matrix_line_endings.json` retains the coordinator's receipt. The final scoped command below exited 0 with no output.

```bash
git diff --check -- tests/test_methods_release.py tests/test_repository.py tests/test_path_dag_reference.py scripts/publication/build_methods_document.py Publication/methods/
```

The coordinator's full-suite command was `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python -m pytest -q`. I independently read `checks/full-suite.json`, verified the raw log against SHA-256 `7189e209c87612e4b0997ecdcb15b8c8421f47eb023a0b7f4991381df9a0854d`, and read its final result, `232 passed in 251.73s (0:04:11)`, exit 0. This full run precedes the final prose-only sampling qualifier and two subsequent Methods-test edits. I compared `checks/full_suite_methods_test_source.py` with the final test. The only differences are the equivalent set-comprehension syntax and the added direct-sampling scope assertion. The final independent 25-test run covers those changes and the actual deliverable. I did not duplicate the full suite.

The final eight-page observation and byte-preservation check cover the actual last Word revision. The bounded PASS does not close the documented scientific provenance or full-submission gaps, and it does not qualify the retained older auxiliary export as current.

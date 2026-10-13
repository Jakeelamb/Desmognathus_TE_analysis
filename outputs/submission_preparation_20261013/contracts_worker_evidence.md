# Repository and DAG contract repair evidence

Owned changes: `tests/test_repository.py`, `tests/test_path_dag_reference.py`, and this evidence file. Other author and agent changes were preserved. Starting repository revision: `15dfb3c18755cc8659b626fa213087c62c6db02f`, with existing dirty/untracked work. Baseline copies already retained by the coordinator were not changed.

The compact layout still excludes all retired paths and the obsolete TE output directory. `outputs/` is permitted as an isolated analysis/review surface. The added check follows actual canonical source registrations, publication dependencies, and the catalog to ensure none resolve into that surface, then verifies that the real eight-figure release keeps its frozen inputs and output artifacts separate.

The obsolete literal caption assertion was removed from the source-string check. The existing genuine R figure fixture now checks the assembled plot caption for the D. fuscus assembly anchor, pg/1C scale, non-direct-C-value qualification, equivalent arrows, and unscored saturated class. Its actual S04 input must retain the conditional nuclear-IOD definition. The plot caption itself does not spell out IOD; the IOD basis is checked in its released input definition. Auxiliary registration, 25-DAG/11-class/10-scored-class mapping, node labels, no-arrow independent panel, and released formats remain covered.

## Targeted checks


Command: `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python -m pytest tests/test_repository.py tests/test_path_dag_reference.py -q`

```text
Started: 2026-10-13T15:57:55Z
.............                                                            [100%]
13 passed in 6.38s

Exit status: 0
Finished: 2026-10-13T15:58:02Z
```

## Rejection controls

Command: `UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env uv run --frozen --no-sync --python 3.13 python -` (read-only in-memory mutations below)

Exact stdin payload:

```python
import runpy
import subprocess
from pathlib import Path

root = Path.cwd()
repository = runpy.run_path('tests/test_repository.py')
check = repository['test_isolated_review_outputs_remain_separate_from_canonical_sources']
check.__globals__['SUPPLEMENT_SOURCES'] = {
    'rejection_control': 'outputs/manuscript_polished_20260924_01/qa/source_baseline.json'
}
try:
    check()
except AssertionError as error:
    print('PASS: canonical source redirected into outputs rejected:', error)
else:
    raise RuntimeError('Canonical-output redirect was accepted')

dag = runpy.run_path('tests/test_path_dag_reference.py')
fixture = dag['test_independent_dag_keeps_three_nodes_without_arrows_or_extra_annotation']
code = next(value for value in fixture.__code__.co_consts
            if isinstance(value, str) and 'suppressPackageStartupMessages' in value)
original = 'caption <- tolower(Path_DAG_reference$labels$caption)'
replacement = '''caption <- tolower(gsub("not a direct species C-value", "a direct species C-value",
                             Path_DAG_reference$labels$caption, fixed=TRUE))'''
assert original in code
result = subprocess.run([str(root / 'scripts/run_r.sh'), 'Rscript', '-e',
                         code.replace(original, replacement)],
                        cwd=root, capture_output=True, text=True, timeout=90)
print('Mutated genuine R fixture exit:', result.returncode)
print(result.stdout, end='')
print(result.stderr, end='')
assert result.returncode != 0
assert 'grepl("not a direct species c-value", caption, fixed = TRUE) is not TRUE' in result.stderr
print('PASS: removal of direct-C-value qualification rejected by the assembled-caption check')
```

```text
Started: 2026-10-13T15:58:44Z
PASS: canonical source redirected into outputs rejected: outputs/manuscript_polished_20260924_01/qa/source_baseline.json
Mutated genuine R fixture exit: 1
Error: grepl("not a direct species c-value", caption, fixed = TRUE) is not TRUE
Execution halted
PASS: removal of direct-C-value qualification rejected by the assembled-caption check

Exit status: 0
Finished: 2026-10-13T15:58:48Z
```

## Final revision and diff checks

Command: `git rev-parse HEAD`

```text
15dfb3c18755cc8659b626fa213087c62c6db02f
```

Command: `git diff --check -- tests/test_repository.py tests/test_path_dag_reference.py`

```text
Exit status: 0
```

## Worker diff from preserved dirty baseline


Command: `diff -u outputs/submission_preparation_20261013/baseline/tests/test_repository.py tests/test_repository.py`

```diff
--- outputs/submission_preparation_20261013/baseline/tests/test_repository.py	2026-08-28 23:03:42.207613618 -0700
+++ tests/test_repository.py	2026-10-13 08:57:46.063497333 -0700
@@ -9,6 +9,8 @@
 import numpy as np
 import pandas as pd
 
+from scripts.project import SUPPLEMENT_SOURCES, TREE_SOURCES
+
 ROOT = Path(__file__).resolve().parents[1]
 
 
@@ -30,7 +32,6 @@
         "results",
         "path_analysis",
         "notebooks",
-        "outputs",
         "plans",
         "Dusky.yml",
         "paths.yaml",
@@ -40,6 +41,36 @@
     assert not (ROOT / "analyses/01_transposable_elements/output").exists()
 
 
+def test_isolated_review_outputs_remain_separate_from_canonical_sources() -> None:
+    outputs = (ROOT / "outputs").resolve()
+    canonical_sources = [*SUPPLEMENT_SOURCES.values()]
+    canonical_sources.extend(record[0] for record in TREE_SOURCES.values())
+    manifest = pd.read_csv(ROOT / "Publication/datasets/DATASET_MANIFEST.csv")
+    canonical_sources.extend(manifest["source_path"])
+    for dependencies in manifest["additional_source_paths"]:
+        canonical_sources.extend(json.loads(dependencies))
+    canonical_sources.extend(pd.read_csv(ROOT / "data/DATA_CATALOG.csv")["path"])
+    for relative in canonical_sources:
+        source = (ROOT / relative).resolve()
+        assert source.is_file(), relative
+        assert not source.is_relative_to(outputs), relative
+
+    review = outputs / "manuscript_polished_20260924_01"
+    frozen_inputs = (review / "inputs/source_run").resolve()
+    assert frozen_inputs.is_dir()
+    assert frozen_inputs.is_relative_to(outputs)
+    figure_manifests = list((review / "figures").glob("*/FIGURE_MANIFEST.csv"))
+    assert len(figure_manifests) == 8
+    for figure_manifest in figure_manifests:
+        figures = pd.read_csv(figure_manifest)
+        assert len(figures) == 1
+        for field in ("pdf", "png", "tiff"):
+            artifact = (figure_manifest.parent / figures.iloc[0][field]).resolve()
+            assert artifact.is_file(), artifact
+            assert artifact.is_relative_to(outputs), artifact
+            assert not artifact.is_relative_to(frozen_inputs), artifact
+
+
 def test_repository_validator() -> None:
     result = run_python("-m", "scripts.project", "validate")
     assert result.returncode == 0, result.stdout + result.stderr

Exit status: 1 (1 is expected for the retained patch)
```

Command: `diff -u outputs/submission_preparation_20261013/baseline/tests/test_path_dag_reference.py tests/test_path_dag_reference.py`

```diff
--- outputs/submission_preparation_20261013/baseline/tests/test_path_dag_reference.py	2026-09-17 18:29:04.754436588 -0700
+++ tests/test_path_dag_reference.py	2026-10-13 08:57:46.063497333 -0700
@@ -62,10 +62,6 @@
     assert "record_auxiliary_figure(" in block
     assert 'node_label = c("Genome\\nsize", "Nucleus\\narea", "Cell\\narea")' in block
     assert '"Relative\\nIOD"' not in block
-    assert "conditional D. fuscus assembly-anchored pg/1C estimate" in block
-    assert "not a direct species C-value" in block
-    assert "Arrows are not unique within classes" in block
-    assert "saturated class (6 members) unscored" in block
     for model in REPRESENTATIVES:
         assert f'"{model}"' in block
 
@@ -130,8 +126,18 @@
             identical(sort(points$x), c(.15, .5, .85)),
             identical(sort(points$y), c(.23, .23, .78)))
 }
-stopifnot(grepl("Arrows are not unique within classes", Path_DAG_reference$labels$caption,
-               fixed=TRUE))
+caption <- tolower(Path_DAG_reference$labels$caption)
+stopifnot(grepl("d. fuscus", caption, fixed=TRUE),
+          grepl("assembly-anchored", caption, fixed=TRUE),
+          grepl("pg/1c", caption, fixed=TRUE),
+          grepl("not a direct species c-value", caption, fixed=TRUE),
+          grepl("arrows are not unique within classes", caption, fixed=TRUE),
+          grepl("saturated class (6 members) unscored", caption, fixed=TRUE))
+traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
+stopifnot(nrow(traits) == 24L,
+          all(grepl("derived from nuclear IOD", traits$genome_size_definition, fixed=TRUE)),
+          all(grepl("conditional D. fuscus assembly-anchored", traits$genome_size_definition,
+                    fixed=TRUE)))
 '''
     result = subprocess.run(
         [str(ROOT / "scripts/run_r.sh"), "Rscript", "-e", code],

Exit status: 1 (1 is expected for the retained patch)
```

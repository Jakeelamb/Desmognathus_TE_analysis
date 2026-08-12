# 03 — Phylogenetic path analysis

This is the current implemented Path24 comparative-analysis surface. The fitted
model contains three traits: relative nuclear IOD, nucleus area, and cell area.
SVL, life-history, and TE-overlap tables are retained for inspection and
descriptive sensitivity work; they are not predictors in this fitted model.

The intended final project step may combine TE composition and collaborator
maximum SVL with morphology/IOD traits on an exact complete-case species panel.
That integrated model has not yet been specified or fitted. The current
three-trait results must not be presented as the completed all-data path
analysis.

## Inputs and outputs

- `data/path24_traits.csv`: exact 24-species model traits.
- `data/measurement_bootstrap_traits.csv`: measurement-bootstrap traits using
  the relative-IOD scale.
- `data/candidate_dags.csv`: candidate three-trait graph definitions. This file
  enumerates 25 possible three-node DAG members that collapse to 10 testable
  Markov-equivalence classes plus a 6-member saturated class that is unscored
  because it has no d-separation claim.
- `trees/source_time_tree_46.tre`: preserved dated source tree provided by Alex Pyron; its publication and calibration citation remains unresolved.
- `trees/path24_time_tree.nwk`: exact pruned tree used for Path24.
- `trees/published_*`: compact published-tree topology sensitivity inputs.
- `data/path_model_*`, `path_*`, `pairwise_*`, and `evolutionary_*`: current
  frozen fit and sensitivity outputs.
- `data/collaborator_max_svl.{xlsx,csv}` and organismal-trait tables: retained
  source/curation lane, outside the current three-trait fit.

Run a fast refit into ignored disposable output:

```bash
make path-quick
```

Render the compact review surface with:

```bash
make report
```

## Interpretation and provenance

These models are exploratory phylogenetic association/sensitivity analyses.
Markov-equivalent or otherwise close graphs cannot be used to assert a unique
causal arrow direction.

The unnumbered auxiliary DAG image placed immediately after Figure 7 in the
researcher-facing study report shows one representative DAG per scored
equivalence class. That gallery is descriptive only: relative nuclear IOD is
not absolute genome size, omitted members are Markov-equivalent
reorientations, and the saturated 6-member class is excluded because it has no
d-separation claim.

Relative nuclear IOD is an image-derived phenotype, not an independently
validated absolute C-value. Historical model files used a positive rescaling
of the same signal; log transformation followed by standardization makes the
saved fits invariant to that constant, so the active data and labels use the
defensible relative-IOD estimand.

Alex Pyron directly provided the focal 46-tip collaborator tree. Its publication
or archive identifier and calibration provenance are still unresolved and must
be supplied before publication. The exact pruned Path24 tree is reproducible
from that preserved source tree.

The collaborator-SVL aliases give *D. orestes* a maximum of 52 mm, while the
tree retains exact unsuffixed representative tips and prunes alternates. The
source-specific distinction is defined and tested in the
[shared identity-resolution contract](../../docs/IDENTITY_RESOLUTION.md).

`data/organismal_traits.csv` is the canonical one-row-per-species table for
SVL and vetted organismal traits. Its `reproductive_strategy` field preserves
Alex Pyron's three expert-coded groups as distinct values:

- `direct_development`
- `aquatic_eggs`
- `terrestrial_eggs_aquatic_larvae`

The coding covers 33 of the 37 canonical species. *D. brimleyorum*, *D.
folkertsi*, *D. fuscus*, and *D. ochrophaeus* remain blank because they were
not explicitly labeled in the supplied figure. The figure's *D. planiceps*
label is retained in the source record and is not transferred to the
morphological *D. fuscus* row. This trait is separate from adult aquaticity and
from the broader binary `development_mode` field.

Pyron supplied the interpretation after seeing the initial PCA. The categories
are retained as expert-vetted biological traits, but their apparent alignment
with that PCA is post hoc and must not be presented as an a priori or
confirmatory test. The original email date still needs to be recorded before a
formal personal-communication citation is finalized.

The TE–IOD overlap has 21 species. *D. brimleyorum*, *D. folkertsi*, and *D.
ochrophaeus* have Path24 measurements but no active TE resource, and their TE
values are not imputed.

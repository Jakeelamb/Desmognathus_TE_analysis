# 03 — Phylogenetic path analysis

This is the current implemented Path24 comparative-analysis surface. The fitted
model contains three traits: conditional *D. fuscus* assembly-anchored genome
size in pg/1C, nucleus area, and cell area.
SVL, life-history, and TE-overlap tables are retained for inspection and
descriptive sensitivity work; they are not predictors in this fitted model.

The intended final project step may combine TE composition and collaborator
maximum SVL with morphology/genome-size traits on an exact complete-case species panel.
That integrated model has not yet been specified or fitted. The current
three-trait results must not be presented as the completed all-data path
analysis.

## Inputs and outputs

- `data/path24_traits.csv`: exact 24-species model traits.
- `data/measurement_bootstrap_traits.csv`: measurement-bootstrap traits using
  the conditional genome-size scale.
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

The two path-analysis auxiliary images placed after the numbered figures in the
researcher-facing report serve different review purposes. The first aligns the
exact 24-tip tree with development mode and conditional genome-size,
nucleus-area, and cell-area intervals. It contains the two sampled direct
developers, *D. aeneus* and *D. wrighti*, plus 22 aquatic-larval species. *D.
organi* is absent from Path24 and no development mode or phenotype is imputed.
This alignment is descriptive and is not a life-history hypothesis test.

The second auxiliary image shows one representative DAG per scored equivalence
class. Its first node is labeled **Genome size**, meaning the conditional *D.
fuscus* assembly-anchored pg/1C sensitivity derived from IOD. Omitted members
are Markov-equivalent reorientations, the saturated 6-member class is excluded
because it has no d-separation claim, and the representative arrows do not
identify a unique causal direction.

The genome-size trait is produced by multiplying the target-to-*D. fuscus*
equal-image IOD ratio by the conditional 16.462167689-pg/1C assembly anchor. It
is the default analysis and display scale, but it is not an independent
genome-size assay or an independently validated absolute C-value. On the log10
scale this positive conversion adds a constant, and sample standardization
removes that constant, so the saved Path24 fit and model ranking are invariant
to the conversion. The IOD measurement and calibration audit remains in
Publication S21-S25; no relative-IOD analytical field remains in this analysis
directory. The DAGs contain one genome-size node and no separate IOD node or
IOD-to-genome-size arrow, because those labels describe one derived trait rather
than two causal variables.

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

The TE–genome-size overlap table is `data/te_genome_size_overlap.csv` and has 21
species. *D. brimleyorum*, *D. folkertsi*, and *D. ochrophaeus* have Path24
measurements but no active TE resource, and their TE values are not imputed.

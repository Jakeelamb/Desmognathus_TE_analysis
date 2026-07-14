# Desmognathus collaborator-review notebooks

This is the canonical eight-notebook review sequence for sharing the project. The
denominators stay explicit throughout: **TE34 n = 34** for vetted genomic-resource
descriptions, **Cell21 n = 21** for the established linked-cell panel,
**Morphology24 n = 24** for finalized microscopy, **Genome24 n = 24** for the
quality-matched genome-size estimates and three-trait path analysis, and
**path18 n = 18** only for the exact TE/cell/tree intersection. Reduced
microscopy sample sizes remain visible as interpretation fields and do not
exclude species from Genome24.

The TE and LTR figures are regenerated from audited current data using historical
R code/style. They are not screenshots or extracted slide panels. No historical
PDF/PNG is presented as a current result; the supplied manuscripts and lab-meeting
decks remain visual references only.

Six notebooks are thin review surfaces over frozen tables, manifests, trees, and
figures. Notebook 06
also retains and executes the complete lightweight analysis of the finalized
805-nucleus image-IOD panel, and Notebook 07 executes the lightweight exhaustive
DAG/input-design audit. The bundle can open optional local microscopy HTML
review pages. It does not rerun RepeatModeler,
RepeatMasker, dnaPipeTE, TE-library
construction, read mapping, cell/nucleus model training, segmentation inference, or
phylogenetic path-model fitting.

The checked-in static audit bundle is self-contained: every small data, tree,
manifest, and reused-figure dependency is snapshotted below `results/data/research_review/`
or `results/figures/research_review/`. The source-to-snapshot mapping, copy time, and
both SHA256 values are recorded in
`results/data/research_review/frozen_input_registry.csv`.
Complete figure provenance is the union of the R-runner manifest for newly
generated historical-style figures and this registry for reused audit figures.

## Review order

| Order | Notebook | Panel and purpose |
|---:|---|---|
| 1 | [Phylogenetic tree trimming](01_phylogeny_tree_trimming.ipynb) | path18 tree provenance, exact pruning, topology, calibration, and 200-tree uncertainty |
| 2 | [Data tables and provenance](02_data_tables_and_provenance.ipynb) | TE34/Cell21/path18 membership, accession identity, retry handling, source registry, and release gates |
| 3 | [Repeat analysis](03_repeat_analysis_te34.ipynb) | TE34 sequencing/assembly QC, mass accounting, composition, restored cross-taxon Shannon context, restored genome-size-estimate scatter on the exact path18 overlap, denominator sensitivity, PCA, scree, elbow, silhouette, clusters, and cluster/tree view |
| 4 | [LTR deletion footprint](04_ltr_deletion_footprint.ipynb) | 1,086 usable elements across 30 species; zero-aware terminal:internal depth, two source-corrupt exclusions, one explicitly audited mapping-pileup screen for presentation, full raw/log companions, coverage QC, and robust tests |
| 5 | [Cell modeling and measurement](05_cell_modeling_and_measurement.ipynb) | Cell21 traits, mask/model viewers, quality-screened selected-50 upper-tail estimand, segmentation checks, relative nuclear-IOD, pairwise traits, and tree tracks |
| 6 | [Genome-size estimation](06_genome_size_estimation.ipynb) | complete frozen quality-matched image-IOD analysis, conditional *D. fuscus*-anchored pg estimates, raw nuclei, image balance, aggregation sensitivity, pairwise traits, and measured phylogeny |
| 7 | [Cell–nucleus–genome path analysis](07_cell_nucleus_genome_path_analysis.ipynb) ([HTML](cell_nucleus_genome_path_analysis/index.html)) | all 24 finalized species; all 25 labeled DAGs/11 equivalence classes; primary d-separation/CICc fits; 250 measurement bootstraps, 200 published trees, 24 leave-one-out fits, 4 evolutionary models, and 800 simulations; cell-bridge class favored, causal orientation unidentified |
| 8 | [Integrated phylogenetic path analysis](08_integrated_phylogenetic_path_analysis.ipynb) | exploratory path18 TE/IOD/morphology DAG comparison, global fit, edges, measurement/tree/species sensitivity, and simulation calibration |

There is one notebook per distinct review domain and no second presentation tree.
The empty source notebook and HTML export formerly paired with the frozen-IOD
notebook were removed because the executed canonical notebook already retains
both its code and rendered outputs.

## Figure provenance

The exact current R runner is
[`scripts/processing/build_audited_historical_style_figures.R`](../../scripts/processing/build_audited_historical_style_figures.R).
It writes the compact review products to:

- `results/figures/research_review/`
- `results/data/research_review/`

The runner uses the preserved historical R plotting logic with the snapshotted,
audited TE34 and LTR inputs. In particular, `SRX19952890` and `SRX19952890R2` are normalized as
separate completed runs and then combined by equal-weight within-species averaging;
they remain one *D. orestes* species observation.

Notebook 03 and notebook 04 deliberately raise a clear missing-output error when
these lightweight R outputs are absent. First build or verify the portable
dependency snapshots and notebook sources:

```bash
scripts/run_in_dusky.sh python scripts/processing/build_research_review_notebooks.py
```

Then generate the lightweight R outputs with:

```bash
scripts/run_in_dusky.sh Rscript scripts/processing/build_audited_historical_style_figures.R
```

This command rebuilds figures and compact statistics from frozen current tables. It
does not rerun any expensive upstream repeat-discovery or read-processing tool.

Notebook 07 has its own fitted all-finalized-species path branch. To refit all declared
release sensitivities and regenerate its eight PNG/PDF figures plus notebook
source:

```bash
scripts/run_in_dusky.sh Rscript \
  path_analysis/scripts/run_cell_nucleus_genome_phylogenetic_path_analysis.R
scripts/run_in_dusky.sh python \
  path_analysis/scripts/build_cell_nucleus_genome_path_presentation.py
```

The full branch runs 250 measurement bootstraps, all 200 published trees, 24
leave-one-species-out fits, four evolutionary models, and 800 simulations. It
uses frozen compact inputs and does not rerun imaging or bioinformatics pipelines.

## Build and open

Rebuild the eight notebook sources, the lightweight genome-size products, and
the frozen dependency registry:

```bash
scripts/run_in_dusky.sh python scripts/processing/build_research_review_notebooks.py
```

Execute the notebooks only after both R figure bundles exist:

```bash
for notebook in notebooks/research_review/[0-9][0-9]_*.ipynb; do
  scripts/run_in_dusky.sh jupyter nbconvert \
    --to notebook --execute "$notebook" --inplace \
    --ExecutePreprocessor.timeout=600
done
```

Then open the review sequence:

```bash
scripts/run_in_dusky.sh jupyter lab notebooks/research_review
```

## Optional local microscopy viewers

The embedded static mask examples and all other audit figures are portable. The
large interactive HTML galleries remain in the optional sibling workspace
`/home/jake/Projects/cellprofiler_test` and are intentionally not copied into this
bundle. Notebook 05 reports which local viewers are available without failing when
they are absent. During a local interactive session, run
`start_optional_local_viewers()` in that notebook to start the loopback-only server
and open the galleries. The executed notebook intentionally stores no ephemeral
`127.0.0.1` URL or live iframe.

Verify the complete executed bundle before sharing:

```bash
uv run pytest -q scripts/python/tests/test_research_review_bundle.py
```

## Interpretation guardrails

- TE34 is not filtered by microscopy availability; Cell21 is not filtered by WGS
  availability; path18 alone is the integrated comparative denominator.
- Genome resources and microscopy specimens are not assumed to be the same sample.
- Current *fuscus* TE/LTR results use the expert-reidentified short-read resource
  `SRX20497025 / GCA_032353935.1`; the separate 2025 chromosome-level
  `GCA_050004315.1` assembly is validation-only and was not silently substituted.
- dnaPipeTE retry labels are computational runs, not extra species or biological
  replicates.
- The LTR terminal:internal statistic is a deletion-footprint/mapping proxy, not a
  measured ectopic-recombination rate.
- The corrected LTR30 primary table retains explicit zero-depth positions and has
  1,086 usable elements. Two truncated source files from the 1,088 selected set are
  excluded without repair; the old nonzero-only table is comparison-only.
- The artifact-screened LTR presentation retains 1,085 elements and all 30 species.
  Its single flagged one-sided terminal pileup remains visible in the row-level
  audit and unfiltered raw/log figures; no blanket IQR deletion is applied.
- The cross-taxon Shannon boxplot is descriptive context built on the published
  top-10-superfamily metric contract. Cross-study pipeline differences prevent it
  from serving as a homogeneous comparative test.
- The historical Shannon-versus-genome-size scatter is restored on the exact
  18-species TE/genome overlap. Its old all-16 placeholder has been replaced by
  variable fuscus-anchored image-IOD genome-size estimates with conditional
  bootstrap intervals. These are not direct C-values, and the pooled comparison
  with ten historical salamanders is descriptive because assay and TE-annotation
  workflows differ.
- The separate quality-matched relative nuclear-IOD panel remains an unscaled
  sensitivity index; it is not absolute genome size and is not calibrated to pg.
- Notebook 06 reports a distinct conditional rescaling of the quality-matched
  image-IOD measurements to *D. fuscus* = 16.36 pg. Those values are useful as
  genome-size estimates, but the anchor's exact literature provenance remains
  unresolved, anchor uncertainty is not propagated, and the estimates are not
  independent direct C-values.
- Notebook 07 fits the ten testable genome/nucleus/cell equivalence classes. The
  finalized 24-species panel favors the cell-bridge class; the three proposed
  nucleus-bridge arrow orientations remain Markov-equivalent, and arrow direction
  requires external identifying information.
- The selected-50 cell statistic uses the highest composite-ranked eligible linked
  cells; it is an upper-tail estimand, not a literal largest-mask sort or a typical
  erythrocyte size.
- The path analysis is exploratory. Shared-measurement algebra, tree uncertainty,
  model-selection uncertainty, and actual-tree calibration must remain visible.

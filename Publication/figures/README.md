# Genome Biology and Evolution figure release

The target journal is Genome Biology and Evolution (GBE). Official requirements
are recorded in `../journal_requirements/GBE_REQUIREMENTS.md`; the following
house-style choices apply within those requirements:

- Arial text, with no silent font fallback;
- 180-mm double-column and 88-mm single-column widths as project choices, not
  GBE figure-size mandates;
- 8-pt base text and a 7-pt minimum target at final size, above GBE's 6-pt floor;
- 1-pt axis lines as a project convention;
- newly rendered vector PDF masters with embedded Arial;
- review PNGs and 300-ppi CMYK TIFF revision derivatives for color figures;
- accessible colors plus non-color encodings;
- lowercase bold panel tags without parentheses;
- no figure-internal title unless scientifically necessary;
- one legend and one alt-text record for every finalized figure.

The displayed trait name is **Estimated genome size** throughout the review
gallery, including collage axes, standalone panel titles, captions and alt
text. Captions retain the IOD-derived, D. fuscus assembly-anchor and uncertainty
qualifications. Existing file basenames and source-column identifiers remain
stable so figure links and saved review notes keep their identities.

The curated placement and scientific role of every current figure are in
`GBE_FIGURE_PLAN.md`; diagnostic
and audit figures belong in the supplementary-information PDF or DOI repository.

The current pass generates the complete draft set before author-led refinement.
The review pool contains seven numbered candidates and five auxiliary/reference
candidates, shown as provisional review Figures 8-12. Every collage has isolated
panel exports. An additional Figure 4 PGLS companion has its own overview and
one isolated scatterplot. Two filter-sensitivity comparisons add matched violin and
PGLS collages, each with two isolated panels. Scientific values, shared palette, species ordering, axes, Arial,
and export properties are checked before review; draft completion does not
mean final manuscript selection or approval.

## Rebuild

For the first Wang-inspired review candidate, render Figure 1 alone:

```bash
scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R \
  --figure Figure_1 --output-dir Publication/figures/candidates/figure_1_wang
```

This writes the full collage's PDF, review PNG, CMYK TIFF, manifest, legend,
and alt text, plus separate exports of panels a-e in `panels/`, then stops
before Figure 2. Each isolated panel is 180 x 190 mm and reuses the exact
collage plot data, scales, species order, and colours. Species labels and
each panel's own legend are included so it can be read independently.
`panels/PANEL_MANIFEST.csv` and `panels/PANEL_LEGENDS_AND_ALT_TEXT.md` describe
the standalone files. Open the candidate's `index.html` to inspect the collage
and all five panels, with links to their PDFs, full-size PNGs, and CMYK TIFFs.

Render the next candidate, Figure 2, separately:

```bash
scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R \
  --figure Figure_2 --output-dir Publication/figures/candidates/figure_2_wang
```

Figure 2 contains the 34-species Bonsai composition topology, feature-resampling
recovery, and post-hoc strategy-association summary. Its 180 x 190-mm collage
and three isolated panels use Arial and the same strategy colours and shapes.
The topology constrains its species order; it is not a species phylogeny.
Recovery uses a full 0-100% axis, and the two tree support labels use exact
released refit counts (27.5% and 35.0%) rather than rounded Newick labels.

Render Figure 3 separately:

```bash
scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R \
  --figure Figure_3 --output-dir Publication/figures/candidates/figure_3_wang
```

Figure 3 replaces the low-divergence heatmap with three matched species dot
plots. Its 180 x 220-mm collage stacks LTR/LINE/TIR divergence profiles above
the 34-species comparison. Profiles share 0-50% divergence and 0-3% of all
species RepeatMasker-aligned hit bp per bin; the dot plots share 0-6% for the
sum below 5% divergence. Species retain the canonical supplied-tree order.
Gray curves show individuals, and coloured curves with distinct line types
show equal-weight means for the 7/24/3 reproductive-strategy groups. These
are descriptive repeat-contig assay summaries, not insertion ages or genome
occupancy. Isolated panels a and b are 180 x 85 and 180 x 170 mm.

Render Figure 4 separately:

```bash
scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R \
  --figure Figure_4 --output-dir Publication/figures/candidates/figure_4_wang
```

Figure 4 aligns the topology-only 30-species tree with horizontal violin
outlines, all 380 primary LTR element ratios as hollow points, short blue
median bars, and thin 95% intervals from 2,000 within-species element bootstraps.
Violin densities use log2 ratios, one pooled bandwidth, equal maximum widths,
and tails trimmed to each species' observed range. Right-hand labels pair
estimated genome size (Gbp*) and retained element count, separated by a slash.
The 19 available genome sizes use the released IOD-derived estimates converted
at 0.978 Gbp per pg; 11 unavailable sizes are marked `n/a` without imputation.
There are 3-50 retained elements per species; violin widths do not encode counts,
and sparse-species shapes are descriptive smoothing. Blue identifies summaries
in this figure; it does not encode
reproductive strategy. The log2 ratio axis and canonical species rows are
identical in the collage and isolated depth panel, which restores its own
species labels. The collage and both isolated panels are 180 x 190 mm with
embedded Arial. Median bars preserve the released medians of untransformed
ratios, then display them on the log2 axis. Decimal axis labels and the dotted
line at one mark equal terminal and internal depth. This remains a
mapping/deletion-footprint proxy, not a rate.

An exploratory companion to Figure 4 compares species median LTR ratios with
estimated genome size. Build it independently with:

```bash
scripts/run_r.sh Rscript scripts/publication/build_ltr_genome_pgls.R
```

`candidates/figure_4_pgls/` contains a single-panel overview plus the isolated
scatterplot, both 180 x 130 mm. The author omitted the former panel b
(leave-one-species-out plot) from review; the numerical influence audit remains
in `analysis/leave_one_species_out.csv`, and existing notes keep their stable
identities. The primary log-log Pagel-lambda
PGLS uses the 19 species with both measurements and retains the shared tree's
original branch lengths. Eleven LTR species without genome estimates are
excluded without imputation. One regression observation represents one species;
the 239 elements supplying these medians are not regression replicates.
The `analysis/` directory contains the exact joined data, exclusions, fitted
tree, model summaries, all 19 influence fits, prediction ribbon, input hashes,
software versions, and interpretation limits. Genome and median uncertainty
bars are descriptive; the PGLS does not propagate measurement error. Tree
provenance remains unresolved, so this is a review analysis, not final inference.

Build the five-domain versus four-domain visual comparisons independently:

```bash
scripts/run_r.sh Rscript scripts/publication/build_ltr_filter_comparison.R
```

`candidates/figure_4_filter_comparison/` pairs the current 380-element branch
with a 989-element four-domain sensitivity across the same 30 species. Both
retain the length, LTR/Gypsy, exact-annotation and regional coverage gates.
Within-species IQR fences are recalculated for the expanded eligible cohort,
adding 611 elements and excluding two formerly retained elements. The collage
uses shared species rows, 0.025-2.5 log2 ratio axes and the same smoothing
bandwidth; standalone panels also show estimated genome sizes. All retained
points are displayed. Current medians and intervals are preserved exactly;
four-domain intervals use 2,000 element-bootstrap replicates with recorded R
seeds (the same percentile method, a different random stream from the primary
NumPy implementation). Collage size is 180 x 210 mm; standalones are 180 x 190 mm.

`candidates/figure_4_pgls_comparison/` stacks matched PGLS plots for the same
19 species and supplied tree: 239 versus 606 contributing elements. The slopes
are 0.0242 and 0.0202, with P = 0.931 and 0.925, respectively. Both intervals
cross zero. Shared axes retain all measurement intervals. The collage is
180 x 265 mm and each standalone is 180 x 130 mm. Exact membership, gate
decisions, paired species summaries, seeds, model fits, predictions and input
hashes are recorded under the violin comparison's `analysis/` directory.
These are review-only sensitivity candidates; neither broadening eligibility
nor an unchanged PGLS establishes structural validity or calibrates one as an
empirical no-deletion baseline. Canonical datasets and existing figure exports
remain unchanged.

Build the remaining drafts independently of the retained Figures 1-4:

```bash
for number in {5..12}; do
  scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R \
    --figure "Figure_${number}" \
    --output-dir "Publication/figures/candidates/figure_${number}_wang"
done
```

The complete review set is:

| Review figure | Subject | Isolated panels |
|---|---|---:|
| 1 | Integrated TE landscape | 5 |
| 2 | TE composition and reproductive strategy | 3 |
| 3 | Repeat divergence profiles and species fractions | 2 |
| 4 | LTR terminal-to-internal depth | 2 |
| 4 PGLS companion | LTR median ratio and estimated genome size | 1 |
| 4 filter comparison | Five-domain and four-domain LTR violins | 2 |
| 4 PGLS filter comparison | Matched PGLS fits for both filters | 2 |
| 5 | Phylogeny and cellular traits | 4 |
| 6 | Estimated genome-size and morphology regressions | 3 |
| 7 | Path-model comparison and simulation calibration | 2 |
| 8 | All-34-species Wang-style integrated TE reference | 5 |
| 9 | TE diversity and estimated genome size | 2 |
| 10 | Path24 development and phenotype reference | 4 |
| 11 | Ten path-model equivalence-class diagrams | 10 |
| 12 | Wang et al. cross-study diversity reference | 4 |

Figures 8-12 retain their auxiliary status. Their review numbers are navigation
labels, not manuscript figure assignments. Cross-study values remain explicitly
method-labelled, and path diagrams remain representative equivalence classes.

Build the continuous review page from all candidate manifests:

```bash
review_args=()
for number in {1..12}; do
  review_args+=(--candidate "Publication/figures/candidates/figure_${number}_wang")
  if [[ "$number" == 4 ]]; then
    review_args+=(--candidate Publication/figures/candidates/figure_4_pgls)
    review_args+=(--candidate Publication/figures/candidates/figure_4_filter_comparison)
    review_args+=(--candidate Publication/figures/candidates/figure_4_pgls_comparison)
  fi
done
uv run --frozen --no-sync python scripts/publication/build_candidate_review.py "${review_args[@]}" \
  --output Publication/figures/candidates/index.html
```

The combined page shows each collage followed by its isolated panels, with
navigation links, captions, PDF/PNG/TIFF downloads, and a separate notes box
for every collage and panel. Start the local review server with:

```bash
uv run --frozen --no-sync python scripts/publication/serve_candidate_review.py
```

Open `http://127.0.0.1:8766/`. Use **Add / edit note** to jump to a figure's
notes. Changes save automatically to `Publication/figures/review_notes.json`,
outside the generated gallery. **Saved to project** confirms the save;
**Save now** retries pending changes and **Download notes** exports the current
comments as Markdown, including saved comments on panels removed from the
current gallery. The browser also keeps drafts if the server is unavailable.

Note identities follow stable figure/panel anchors, so rebuilding the gallery
does not erase comments. Figure-image and download links include a content
revision so a refreshed gallery cannot reuse an older cached export. A changed
note in another tab prompts a choice between
the saved version and the current draft. Panel manifests may set an optional
`display_panel_id` when visible letters change; `panel_id` stays fixed for
existing review anchors and notes. Figure 10 uses this after removing the
development-mode track: visible panels a-d keep the original a/c/d/e note keys.
Use this server for note-taking; a plain static-file server or a `file://`
preview cannot save notes into the project. The server binds only to localhost.

The current note-driven refinements add numeric estimates and a common circle
symbol to Figure 1b, tighten and label Figure 10's three trait axes, move
Figure 10's development-mode information into the caption, and remove the
redundant independent-DAG annotation in Figure 11. All released measurement
intervals are retained. Main-release exports remain unchanged pending selection.

Single-figure builds require a separate
output directory and cannot replace the main release or another candidate's
manifest. Selectors are available for Figures 1-12. The candidate pool
may contain more than eight figures; selection for GBE's main-text display
limit happens later.

The complete release can still be rebuilt with:

```bash
make setup-r  # first run only
make figures
```

`FIGURE_MANIFEST.csv` maps each figure to its source supplemental datasets,
release status, dimensions, font, legend, alt text, and PDF/PNG/TIFF files.
New builds require Arial and must fail if it is unavailable. Verify Arial
embedding in each newly rendered PDF and CMYK color space in its revision TIFF;
existing PDFs and manifest records are not evidence of compliance with this
updated contract.

Figure 2 is `Figure_2_te_diversity_and_bonsai`. It reads S42, S43, and
`trees/te34_bonsai_composition_tree_with_feature_support.nwk`. The Bonsai tree
is a non-phylogenetic representation of the 24-superfamily CLR composition;
the display root has no phylogenetic meaning, and feature resampling perturbs
compositional coordinates rather than biological specimens. S09-S11 retain the
PCA scores, variance, and loadings as a supplementary geometry audit; the
Wang-style diversity-versus-genome-size display remains an auxiliary candidate.

`AUXILIARY_FIGURE_MANIFEST.csv` records five separate unnumbered review figures,
all intentionally outside `FIGURE_MANIFEST.csv`:

- `Wang_style_integrated_TE_history_candidate.{pdf,png,tif}` is the integrated
  TE-history candidate. Its Figure 8 review export in `candidates/figure_8_wang/`
  now includes all 34 TE species, 21 estimated genome-size estimates with
  measurement intervals, and 13 explicit `n/a` entries. Stacked bars retain
  all repeat categories; matched LTR/LINE/TIR dots and unboxed divergence
  profiles use common 0-60% fraction and 0-3% profile-height scales. The full
  species aligned-hit denominator is retained. Repeat-hit fractions are not
  genome occupancy and divergence bins are not insertion ages. The root-level
  release export remains unchanged pending selection of the review version.
- `Wang_style_TE_diversity_genome_size_candidate.{pdf,png,tif}` adapts the
  reference paper's two-row small-multiple grammar to the exact overlap21.
  Pagel-lambda PGLS fits are shown for all species and the two sufficiently
  sampled reproductive-strategy strata; the two direct developers are shown
  without a fitted regression.
- `Path24_development_and_trait_phylogeny_reference.{pdf,png,tif}` aligns a
  topology-only display of the exact 24-tip tree with development mode and estimated genome-size,
  nucleus-area, and cell-area estimates. Path24 contains *D. aeneus* and *D.
  wrighti* as its two direct developers and 22 aquatic-larval species; *D.
  organi* is absent, and no trait is imputed.
- `Path_DAG_equivalence_class_reference.{pdf,png,tif}` shows one representative
  DAG for each of the 10 testable Path24 equivalence classes while leaving the
  saturated class unscored.
- `Wang2025_cross_study_TE_diversity_reference.{pdf,png,tif}` places the exact
  CC0 Dryad Figure 2 values beside current top-ten-superfamily *Desmognathus*
  summaries. It is descriptive only: the RepeatModeler/RepeatMasker and
  dnaPipeTE pipelines differ, current genome sizes are IOD-derived
  sensitivities, and no pooled regression or hypothesis test is shown.

The DAG node labeled **Genome size** is the *D. fuscus*
assembly-anchored pg/1C sensitivity derived from IOD. Positive rescaling leaves
the standardized Path24 fit invariant, but this is not an independent
genome-size assay. Representative DAG arrows are not a unique causal
orientation.

Figure 4 uses a phylogenetically ordered tree-over-distribution layout: hollow
points show the 380 retained elements, and blue points and horizontal intervals
show species medians and element-bootstrap
uncertainty, and the 1:1 line marks equal terminal and internal depth. Values
above one point in the direction expected under internal-sequence depletion;
the statistic is not an ectopic-recombination or DNA-loss rate.

Figure 6 uses estimated genome size rather than relative IOD. Its basename is
`Figure_6_genome_size_pairwise_pgls`. IOD remains visible only as the
measurement, quality, and calibration provenance for the conditional trait.
The integrated S40 TE-by-genome-size table remains available as supplemental
data but is no longer plotted as a publication figure.

The visual grammar follows Wang et al. (2025; doi:10.1038/s42003-025-08127-3)
where it is compatible with the journal contract: dense aligned small multiples,
raw observations beneath summaries, restrained framing, no grids, lowercase
panel tags, compact gutters, italic scientific names, and a color-blind-safe
blue/orange/green/purple palette. Figure 1 integrates 21 species across a
topology-only tree, reproductive strategy, estimated genome size, dnaPipeTE
composition, and RepeatMasker low-divergence signal. Numbered tree panels omit
branch lengths while publication identifier, calibration method, and
branch-length provenance remain pending.
Categorical species axes use one bottom-to-top order derived from the supplied
46-tip dated tree and preserve that relative order when a figure contains only
a subset. Figure 2's Bonsai tips are the explicit exception because their order
is constrained by the TE-composition topology rather than the species tree.

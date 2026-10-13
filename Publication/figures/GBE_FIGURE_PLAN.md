# Curated Genome Biology and Evolution figure plan

**As of:** 17 September 2026

This plan separates technical readiness from scientific placement. Genome
Biology and Evolution (GBE) is the target journal; the official requirements
authority is `../journal_requirements/GBE_REQUIREMENTS.md`. The shared house
style below implements the agreed figure conventions within those requirements.
The candidate collection is rendered under the Arial and CMYK contract.
Earlier main-release exports remain preserved until author selection; their
compliance is not inferred from the new candidates.

The author has requested the complete draft pool before the refinement pass.
The original pool contains 12 review figures and 46 isolated panels: seven
numbered candidates plus five auxiliary references, provisionally labelled
review Figures 8-12. The exploratory Figure 4 PGLS companion adds one collage
and one retained isolated panel. Two additional filter-comparison collages,
each with two isolated panels, bring the current pool to 15 overview figures
and 51 isolated panels.
Build and review candidates before choosing favorites; the eight combined figures-and-tables limit applies
to the eventual main-text submission, not to this candidate collection. The
placement recommendations below remain provisional.

The first revised candidate is in `candidates/figure_1_wang/`: an Arial,
180 × 190-mm integrated overlap21 display inspired by Wang et al. Figure 1.
It replaces the earlier low-divergence heatmap with 63 full LTR/LINE/TIR
profiles. Every profile uses 0–50% divergence and 0–3% of all species
RepeatMasker-aligned hit bp per 1-percentage-point bin. The nine superfamily
bars share a fixed 0–70% classified dnaPipeTE base scale. Its own manifest,
legend, and alt text accompany the candidate; the existing numbered release
is preserved until figure-by-figure selection.

**Review decision:** retain the Wang-inspired overview. Its aligned layout was
well received, but the many miniature plots make species differences difficult
to compare. The heatmap alternative was rejected and removed.

The next review candidate is `candidates/figure_2_wang/`: a 180 x 190-mm
Arial composition-tree and reproductive-strategy figure for all 34 TE species.
It includes three separately exported panels. The exact-edge recovery display
uses 0-100%, with counts separated from intervals, and the tree shows the exact
27.5% and 35.0% feature-refit recovery values. The global Fitch statistic uses
neutral charcoal. All source values and Bonsai tip order are retained; the
ordinary permutation tests remain post hoc and not phylogenetically corrected.
This candidate awaits author review.

Figure 3 is now a review candidate in `candidates/figure_3_wang/`: an Arial,
180 x 220-mm collage and two isolated panels. Three fixed-scale LTR/LINE/TIR
facets show 34 individual gray curves and the three descriptive strategy
means. A species dot plot replaces the heatmap, preserving the supplied-tree
order and using a common 0-6% axis across orders. Profile axes are 0-50%
divergence and 0-3% of all species RepeatMasker-aligned hit bp per bin, matching
Figure 1's denominator and display range. All observations are retained in
the 50-bin display, with absent bins zero-filled and no subset renormalization.
Figure 3 awaits author review.

Figure 4 is now a review candidate in `candidates/figure_4_wang/`: the
topology-only LTR30 tree and all 380 primary element observations, with
species medians and 95% element-bootstrap intervals. Following the author's
reference image, panel b now uses horizontal violin outlines, hollow points,
short vertical blue median bars, and thin horizontal uncertainty intervals.
The log2 densities use one pooled bandwidth and equal maximum widths, trimmed
to each species' observed range; widths do not encode sample size. The released
raw-ratio medians are transformed for display rather than recomputed in log
space. Decimal ratio labels and the dotted equal-depth line follow the reference
style while all 30 species stay aligned with the phylogeny. Its collage and two
standalone panels are 180 x 190 mm in Arial; the standalone depth panel
restores species labels and preserves the collage's axes and exact row order.
Right-hand labels pair available estimated genome size (Gbp*) and element
count (3-50), separated by a slash: 19 sizes are available and 11 are marked
`n/a`. Sizes are the released IOD-derived assembly-anchored estimates converted
at 0.978 Gbp per pg, not independently measured species C-values.
Blue identifies medians and intervals in this figure, not reproductive strategy.
The figure remains a descriptive mapping/deletion-footprint proxy and awaits
author review.

An additional exploratory companion, `candidates/figure_4_pgls/`, follows
Figure 4 in the gallery. Its single-panel overview and standalone panel a show the
19-species relationship between released median LTR depth ratios and
estimated genome size. The author omitted the leave-one-species-out PGLS plot
(former panel b) from review; its numerical results remain in the analysis
tables. The primary
log-log Pagel-lambda PGLS preserves the supplied dated-tree branch lengths:
slope 0.024, 95% t interval -0.561 to 0.609, P = 0.931, lambda near zero.
All 19 species-deletion intervals cross zero. The result is no detectable
association in this overlap, with uncertainty in both directions. It remains
exploratory because the dated-tree provenance is unresolved and trait
measurement uncertainty is displayed but not propagated into the fit.
The original Figure 4a/b exports and review-note identities remain intact.

Two review-only comparisons follow the PGLS companion:
`candidates/figure_4_filter_comparison/` and
`candidates/figure_4_pgls_comparison/`. The first contrasts the current
five-domain branch (380 elements) with four-domain eligibility (989 elements)
in aligned 30-species violin panels. Length, LTR/Gypsy, exact-annotation and
coverage criteria stay fixed; IQR fences are recalculated before the coverage
intersection, adding 611 elements and excluding two formerly retained elements.
Shared axes and a common bandwidth make the distribution changes inspectable.
The second uses the same 19 species and dated-tree branch lengths in both
PGLS fits: slopes 0.0242 and 0.0202, P = 0.931 and 0.925. All plotted points
and measurement intervals are retained. Every collage and isolated panel has
its own stable review-note identity. Primary data and existing exports are
preserved; the expanded cohort awaits author review and structural assessment.

The remaining drafts are in `candidates/figure_5_wang/` through
`candidates/figure_12_wang/`. Figures 5-7 show the Path24 phenotype tracks,
estimated genome-size regressions, and exploratory path-model sensitivity.
Review Figures 8-12 contain the full TE34 integrated reference, stratified
TE-diversity comparison, Path24 development reference, ten representative DAGs,
and method-labelled Wang et al. cross-study comparison. Their auxiliary
scientific status is retained. `candidates/index.html` collects the full set,
each collage followed by its individual panels, on one scrollable page.

Figure 8 was refined in response to the six saved review notes: all 34 TE
species now appear in the collage and each isolated panel. The genome-size
track retains the 21 observed assembly-anchored estimates and explicitly marks the
other 13 as unavailable. Stacked composition bars retain all repeat categories;
matched LTR/LINE/TIR dot plots and unboxed divergence profiles replace the
sparse seven-column miniature grids. Fractions retain the full aligned-hit
denominator, with common 0-60% fraction and 0-3% profile-height scales. The
original note anchors and comments are preserved for further author review.

The next saved-note pass refines Figures 1, 10, and 11. Figure 1b uses one
circle symbol, numeric genome-size labels, and an 8-27 pg/1C axis retaining
all conditional intervals; support cautions remain in S24 and the caption.
Figure 10 moves development mode into its caption and retains four panels:
phylogeny, genome size, nucleus area, and cell area. Trait estimates receive
numeric labels and tighter axes that include every released interval. Visible
letters are a-d; persistent review anchors retain the original panel identities
so existing notes stay with their plots. Figure 11 removes the redundant
"No directed paths" annotation from the independent diagram. Its three nodes
and absence of arrows are unchanged.

## Recommended main-text set

| Figure | Evidence panel | Narrative role | Why it belongs in the main text |
|---|---:|---|---|
| Figure 1 — integrated TE landscape | overlap21 | Aligns a topology-only phylogeny and reproductive strategy with estimated genome size, broad-order and top-superfamily dnaPipeTE composition, and low-divergence RepeatMasker signal. | Opens with the biological hierarchy and keeps all 21 species in one fixed phylogenetic row order. Assay denominators and genome-size limitations are explicit. |
| Figure 2 — Bonsai composition representation | TE34 | Shows the 34-species TE-composition topology, exact aquatic/direct strategy edges, feature-refit recovery, Fitch localization, and the PERMANOVA/PERMDISP diagnostic. | Central genomic result. It presents strong strategy-associated compositional structure without calling the topology a species phylogeny or claiming validated discrete clusters. PCA remains a supplementary geometry audit. |
| Figure 3 — repeat divergence landscape | TE34 | Shows individual and strategy-mean LTR, LINE, and TIR divergence profiles above phylogenetically ordered low-divergence dot plots. | Shared profile and dot-plot axes expose among-species and among-order differences; means zero-fill absent category/bin combinations and weight species equally within each strategy. |
| Figure 5 — phylogeny and cellular traits | Path24 | Aligns a topology-only phylogeny and development mode with estimated genome size, nucleus area, and cell area. | Shows the phenotype result in biological context while point shapes expose limited image or animal support. |

This four-figure core keeps the main narrative compact; Figures 4 and 6 are
strong sensitivity displays that can move into the main text if the manuscript
emphasizes deletion footprints or the genome-size/cell-size association.

## Recommended supplementary or sensitivity set

| Figure | Evidence panel | Placement | Claim boundary |
|---|---:|---|---|
| Figure 4 — terminal:internal LTR depth | LTR30 | Supplementary sensitivity | A topology-only 30-tip prune aligns the exact 380-element primary branch. Horizontal violin outlines show descriptive log2 distributions; hollow points are elements, short blue bars are released species medians, and thin horizontal intervals retain element-bootstrap uncertainty. Counts expose sparse support. The 1:1 line marks equal terminal and internal depth. Values above one point in the direction expected under internal-sequence depletion, but this mapping/deletion-footprint proxy is not an ectopic-recombination, solo-LTR, or DNA-loss rate. |
| Figure 6 — estimated genome-size pairwise PGLS | Path24 | Supplementary sensitivity | Estimated *D. fuscus* assembly-anchored genome size versus nucleus and cell area; point shape exposes three one-image/one-animal estimates and four failed underlying IOD quality-balance gates. The scale is IOD-derived, not an independent genome-size assay. |
| Figure 7 — phylogenetic path sensitivity | Path24 | Supplementary audit/sensitivity | Exploratory Markov-equivalence comparison shown on the *D. fuscus* assembly-anchored genome-size scale derived from IOD. Standardized fits are invariant to the positive conversion; this is neither an independent genome-size assay nor a basis for orienting causal arrows. |

Figure 7 should move into the main text only if the final manuscript treats the
path-model comparison as a core result and keeps its exploratory interpretation.
If so, replace another main display item rather than exceeding the combined cap.

Five companion review assets are unnumbered and stay outside
`FIGURE_MANIFEST.csv`; their records live in `AUXILIARY_FIGURE_MANIFEST.csv`:

- `Wang_style_integrated_TE_history_candidate.{pdf,png,tif}` adapts the
  reference paper's aligned phylogeny/genome-size/composition/landscape grammar
  to six deterministically selected overlap21 representatives. The composition
  denominator is aligned RepeatMasker hit bp rather than genome occupancy, and
  divergence bins are not insertion ages.
- `Wang_style_TE_diversity_genome_size_candidate.{pdf,png,tif}` shows the two
  canonical top-ten-superfamily diversity indices against estimated genome
  size for all 21 species and the 6/13/2 reproductive-strategy strata. The
  strata are post-hoc; the two direct developers have no fitted regression.
- `Path24_development_and_trait_phylogeny_reference.{pdf,png,tif}` aligns the
  exact Path24 tree with one development-mode track and estimated genome-size,
  nucleus-area, and cell-area intervals. The complete displayed panel contains
  two direct developers (*D. aeneus* and *D. wrighti*) and 22 species with
  aquatic larvae. *D. organi* has no Path24 phenotype record, is absent, and is
  not imputed. This is descriptive review material, not a development-mode
  hypothesis test.
- `Path_DAG_equivalence_class_reference.{pdf,png,tif}` summarizes one
  representative DAG for each of the 10 testable Path24 equivalence classes and
  leaves the saturated class explicitly unscored. Its **Genome size** node is
  the *D. fuscus* assembly-anchored pg/1C sensitivity derived from
  IOD. The standardized fit is invariant to this positive conversion, but the
  estimate is not an independent genome-size assay and the representative
  arrows do not establish a unique causal direction.
- `Wang2025_cross_study_TE_diversity_reference.{pdf,png,tif}` compares exact
  Wang et al. 2025 Dryad values with current top-ten-superfamily
  *Desmognathus* summaries in scatter and group-distribution panels. It is
  auxiliary descriptive context only: source pipelines are not harmonized,
  current genome sizes are IOD-derived sensitivities, and no
  pooled regression or hypothesis test is reported.

## Matching visual system

- Arial throughout, embedded in every newly rendered vector PDF; rendering must
  fail if Arial is unavailable rather than silently substitute another font.
- Eight-point base text with a seven-point minimum target at final dimensions;
  GBE's hard minimum is six points.
- Single-column 88-mm and double-column 180-mm widths are project house-style
  choices, not prescribed GBE figure dimensions.
- White background, minimal framing, no grids, lowercase unparenthesized panel
  tags, compact gutters, and a color-blind-safe blue/orange/green/purple palette.
- Wang et al. (2025; doi:10.1038/s42003-025-08127-3) is the layout benchmark:
  compact small multiples and raw observations layered beneath summaries are
  used wherever our released data support them, without copying their artwork
  or importing claims from a different pipeline.
- Italic scientific names and a single shared typography/spacing system across
  every panel.
- Shape or line-type redundancy wherever color encodes groups.
- No figure-internal titles; the first legend sentence states the result.
- New exports comprise an editable vector PDF master with embedded Arial, a
  review PNG, and a 300-ppi CMYK TIFF revision derivative for every color figure.
- One legend and concise alt text per figure in
  `FIGURE_LEGENDS_AND_ALT_TEXT.md`.
- Every categorical species axis follows one canonical bottom-to-top order from
  the supplied 46-tip dated tree, filtered without value-based reranking for the
  panel's species subset. Figure 2's Bonsai tip layout is topology-constrained
  by the TE-composition representation and is therefore the explicit exception.

## Figure-by-figure review

After generating the complete draft set, work through refinement with the
author: confirm each figure's purpose and reference layout, refine it using
the shared ggplot2 theme, and review its rendered result at final dimensions.
Check the
palette, species order, axes and units, typography, and export properties for
that figure. Existing files and manifest records remain pending this review;
changing the shared theme does not finalize the exported figure set.

## Submission packaging

- Embed only the selected main figures in the initial manuscript PDF.
- Put supplementary figures and their legends in the separate initial
  supplementary-information PDF.
- At revision, submit each requested figure as its own publication-quality file.
- Follow `../journal_requirements/GBE_REQUIREMENTS.md` and recheck the live GBE
  instructions immediately before submission.

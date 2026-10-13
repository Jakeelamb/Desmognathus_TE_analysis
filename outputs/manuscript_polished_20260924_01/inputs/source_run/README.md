# Isolated manuscript candidate run

This run contains thirteen manuscript figure candidates, a preserved composition supporting candidate, an alternative path display, and five consolidated tables. The central question is comparison of hypotheses linking genome size, nuclear size, cell size, and organismal size. Figures 9–13 extend the original eight with animal-level sampling sensitivity, nuclear–cell scaling, TE integration, phylogenetic divergence tests, and alternative roles for body size. The run uses copied inputs and has its own figure numbers, exports, review gallery, and notes. The existing publication files and the gallery at port 8766 are preserved.

Open the review page at **http://127.0.0.1:8767/**. Each collage is followed by its isolated panels, with PDF, PNG, and CMYK TIFF downloads. Five tables follow the figures and are supplied as CSV and one five-sheet Excel workbook in `tables/integrated/`. The earlier two-table workbook remains in `tables/`. All 64 note fields save only to this run's `review_notes.json`.

The primary filtering and original model definitions are retained. Additional models are separate analyses using frozen inputs and explicit formulas. No new specimens, images, upstream assemblies, or read mapping were generated. The divergence extension is included in Figures 11–12. Its primary endpoint is the fraction of LTR aligned-hit bp in divergence bins 0–4, divided by all LTR aligned-hit bp.

## Contents

- `inputs/`: copied publication datasets, source plotting scripts, retained LTR results, and previous review notes for reference.
- `figures/figure_1` through `figure_13`: collage, isolated panels, legends, plotted data, and analysis exports.
- `figures/support_composition`: preserved Bonsai composition candidate and its three panels, displayed as supporting S1.
- `tables/integrated/`: five consolidated table CSVs and `Manuscript_tables_integrated.xlsx`; earlier tables remain in `tables/`.
- `scripts/`: isolated rendering, table, gallery, and verification code.
- `qa/`: original-file preservation hashes, source hashes, figure and workbook checks, and browser evidence.
- `logs/`: execution logs.

## Rebuild

Run the R scripts from this directory using `/home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ`, with that environment's bin directory on PATH. Script order: `scripts/figures_1_2.R`, `scripts/figures_3_5.R`, `scripts/figure_4.R`, and `scripts/figure_6.R`. They write only into this run.

Figure 6 reuses S28 and S30–S37 without refitting its original models. The favored cell-bridge class contains three compatible causal directions. Model rankings and conditional associations distinguish this class from the nucleus-bridge class on the primary 24-species panel, while the 18-species subset separates them less clearly. See `MANUSCRIPT_STORY_REVISION_PATH_CENTRAL.md` for the original emphasis. The subsequent Figures 7–8 extension fits new isolated models as described below.

The weighted-arrow revision of Figure 6 labels each orientation with its matching S28 standardized phylogenetic regression coefficients, scales arrow width by absolute coefficient, and includes normal 95% intervals in panel c. The shared class weight is displayed separately. The previous Figure 6, its rendering script, and a snapshot of its notes are preserved at `revisions/figure_6_before_weighted_arrows/`, with a comparison link in the gallery. The earlier marginal/conditional comparison remains available among Figure 6 supporting results.

An additional reference-style candidate is available at `#figure-6-network` with its own notes and two isolated panels. Rebuild with `scripts/figure_6_network_candidate.R`. It uses the specified genome-first triangle: `cell ~ genome` and `nucleus ~ genome + cell`. Coefficients are 0.716664, 0.573623 and 0.007923. This fully connected model has no testable d-separation claim and is not assigned the cell-bridge class weight. Its files live in `figures/figure_6/alternatives/direct_indirect/`; existing main figures are retained.

From the repository root, use `uv run --frozen --no-sync python outputs/manuscript_story_20260924_01/scripts/build_integrated_tables.py`, followed by the bundled Node runtime on `scripts/build_tables.mjs`. Then run `build_review.py` with the project UV runtime. The Excel build requires the bundled `@oai/artifact-tool` dependency path. The original table builder and payload are retained in `revisions/before_integrated_extension/`; the original table exports remain untouched.

The separate local server is `uv run --frozen --no-sync python outputs/manuscript_story_20260924_01/scripts/serve_review.py --port 8767`. Do not point this server at the original notes file.

Run `verify_run.py verify` with the project UV runtime to check source preservation, figure manifests, assets, and embedded Arial fonts. The original-file baseline was recorded before rendering.

These are review candidates, not a replacement of the current publication release. Table estimates and figure legends state the measurement and inference boundaries. Final manuscript placement remains subject to figure-by-figure review.

All current figure labels and titles use “Genome size” following the author's requested terminology update. The prior exports and scripts remain at `revisions/before_genome_size_labels/`. This display change preserves the data, model results, measurement descriptions, and saved notes.

## Body-size incorporation: Figures 7–8

Rebuild from this run directory with the same R environment: `scripts/body_size_analysis.R`, then `scripts/figure_7.R` and `scripts/figure_8.R`; rebuild the gallery from the repository root using project UV as above. Analysis inputs and SHA-256 provenance are in `inputs/body_size/`. Numerical results, exact model formulas, warnings, tests, fixed full-panel scaling and methods are in `figures/figure_8/analysis/`. Previous figures, tables and notes remain unchanged; `qa/body_size_before.json` records their preservation baseline.

Figure 7 shows SVL for all 24 species in the established phylogenetic display order and its associations with genome size, cell area and nucleus area. Maximum SVL spans 29–122 mm and represents species maxima, not the lengths of microscopy animals. The complete panel is compared with a 22-species sensitivity excluding only D. aeneus and D. wrighti. All species are labeled in isolated scatterplots.

Figure 8 compares four specified hypotheses sharing body-size predictors of all three other traits. Candidate weights therefore compare residual genome–cell–nucleus structures conditional on that assumption; they do not demonstrate that body size is causally upstream or measure the improvement from adding body size. The cell-bridge candidate has weights 0.944 (24 species) and 0.828 (22 species) within this limited candidate set. All models, global independence-test P values and missing-edge constraints are exported. Weights from different candidate sets or species panels should not be compared as though they were on one common scale.

The genome–cell association remains positive after SVL adjustment: full-panel standardized beta 0.545 (95% t interval 0.192–0.899) across 24 species and 0.503 (0.041–0.965) across 22 species. The cell–nucleus association also persists. The full direct genome–nucleus coefficient, conditional on both cell area and SVL, remains imprecise and near zero; this does not prove the true effect is zero. These results extend the cellular-association story without identifying causal arrow direction.

Both panels use log10 traits standardized to the original 24-species means and SDs. Confidence intervals use t quantiles with regression residual degrees of freedom. No measurement-error propagation or new tree sensitivity was added. The two omitted taxa are simultaneously the smallest species and the only direct developers, so this sensitivity cannot separate size, developmental mode and species identity. No body-size cutoff search, exclusion of visual outliers, or separate two-species path fit was performed.

## Integrated extension: Figures 9–13

Rebuild the extension from this run directory with the same R environment: `scripts/animal_robustness.R`, `scripts/nuclear_allometry.R`, `scripts/te_integration.R`, and `scripts/body_role_analysis.R`. The first two also support their documented render-only workflow. Rebuild the consolidated tables and gallery afterward. Every analysis has its own frozen inputs and source checksums. No earlier figure needs rerendering for this extension.

Figure 9 compares released species summaries with the median of per-animal morphology medians and equal-animal IOD summaries. It adds 250 hierarchical measurement resamples and all 24 single-species deletions for the original four SVL-adjusted hypotheses. Figure 10 uses the mean of per-animal morphology medians as a separate equal-animal sensitivity, with 1,000 paired hierarchical refits. These are deliberately named distinct summaries, not one shared alternate dataset. Both retain all eligible objects.

Figures 11–12 use the within-LTR low-divergence share as the focal TE endpoint; the all-hit-denominator share is retained as a denominator sensitivity. Cellular traits retain full24 scaling; the logit TE score uses full34 scaling in the exact21 overlap. The configured dnaPipeTE sample is 1.5 Gb, not 15 Gb, and its quality measures have separate denominators from the RepeatMasker aligned-hit quantities. Strategy coding was provided after reviewing the original ordination, so these are exploratory associations.

Figure 13 broadens the candidate structures beyond body size as a common upstream predictor. Eight distinct equivalence classes are compared using a fixed test order and its reverse sensitivity. In phylopath, order determines the regression response for node pairs without a directed path; ancestral pairs retain the graph-implied direction. The earlier four-model set used the package default ordering. Scores and weights should be compared only within the same species panel, set and testing convention. Graphical equivalence is independently checked by skeleton and unshielded colliders. This expanded set has not received the Figure 9 bootstrap or a new tree sensitivity.

Read `INTEGRATED_RESULTS.md` for the resulting manuscript interpretation. The gallery retains the original path displays and restored composition supporting candidate so figure selection can proceed by review. Several new collages are deliberately tall review layouts; final journal-sized arrangements remain a refinement step.

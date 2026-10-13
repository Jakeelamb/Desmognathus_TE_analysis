# Isolated manuscript candidate run

This run reorganizes the existing Desmognathus results into six manuscript figures and two editable tables. The central question is comparison of causal hypotheses linking estimated genome size, nuclear size, and cell size. Figure 6 restores the phylogenetic path comparison as a main result, following the pairwise associations in Figure 5. It uses copied inputs and has its own figure numbers, exports, review gallery, and notes. The existing publication files and the gallery at port 8766 are preserved.

Open the new review page at **http://127.0.0.1:8767/**. Each collage is followed by its isolated panels, with PDF, PNG, and CMYK TIFF downloads. The two tables follow the figures and are also supplied as CSV and one two-sheet Excel workbook. Notes save only to this run's `review_notes.json`.

The primary filtering and model definitions are retained. Diversity fits are reproduced under the existing model specification to export slopes and intervals; the cellular and LTR regression displays reuse stored results. No new specimens, images, upstream assemblies, or read mapping were generated. The proposed additional phylogenetic divergence analysis is not part of this figure-rendering run.

## Contents

- `inputs/`: copied publication datasets, source plotting scripts, retained LTR results, and previous review notes for reference.
- `figures/figure_1` through `figure_6`: collage, isolated panels, legends, plotted data, and analysis exports.
- `tables/`: manuscript table CSVs and `Manuscript_tables.xlsx`.
- `scripts/`: isolated rendering, table, gallery, and verification code.
- `qa/`: original-file preservation hashes, source hashes, figure and workbook checks, and browser evidence.
- `logs/`: execution logs.

## Rebuild

Run the R scripts from this directory using `/home/jake/miniconda3/envs/desmognathus/bin/Rscript --no-environ`, with that environment's bin directory on PATH. Script order: `scripts/figures_1_2.R`, `scripts/figures_3_5.R`, `scripts/figure_4.R`, and `scripts/figure_6.R`. They write only into this run.

The path figure reuses S28 and S30–S37; no path models are refitted. The favored cell-bridge class contains three compatible causal directions. Model rankings and conditional associations distinguish this class from the nucleus-bridge class on the primary 24-species panel, while the 18-species subset separates them less clearly. See `MANUSCRIPT_STORY_REVISION_PATH_CENTRAL.md` for the revised emphasis.

The weighted-arrow revision of Figure 6 labels each orientation with its matching S28 standardized phylogenetic regression coefficients, scales arrow width by absolute coefficient, and includes normal 95% intervals in panel c. The shared class weight is displayed separately. The previous Figure 6, its rendering script, and a snapshot of its notes are preserved at `revisions/figure_6_before_weighted_arrows/`, with a comparison link in the gallery. The earlier marginal/conditional comparison remains available among Figure 6 supporting results.

An additional reference-style candidate is available at `#figure-6-network` with its own notes and two isolated panels. Rebuild with `scripts/figure_6_network_candidate.R`. It uses the specified genome-first triangle: `cell ~ genome` and `nucleus ~ genome + cell`. Coefficients are 0.716664, 0.573623 and 0.007923. This fully connected model has no testable d-separation claim and is not assigned the cell-bridge class weight. Its files live in `figures/figure_6/alternatives/direct_indirect/`; existing main figures are retained.

From the repository root, use `uv run --frozen --no-sync python outputs/manuscript_story_20260924_01/scripts/build_tables_payload.py`, followed by the bundled Node runtime on `scripts/build_tables.mjs`. Then run `build_review.py` with the project UV runtime. The Excel build requires the bundled `@oai/artifact-tool` dependency path.

The separate local server is `uv run --frozen --no-sync python outputs/manuscript_story_20260924_01/scripts/serve_review.py --port 8767`. Do not point this server at the original notes file.

Run `verify_run.py verify` with the project UV runtime to check source preservation, figure manifests, assets, and embedded Arial fonts. The original-file baseline was recorded before rendering.

These are review candidates, not a replacement of the current publication release. Table estimates and figure legends state the measurement and inference boundaries. Final manuscript placement remains subject to figure-by-figure review.

All current figure labels and titles use “Genome size” following the author's requested terminology update. The prior exports and scripts remain at `revisions/before_genome_size_labels/`. This display change preserves the data, model results, measurement descriptions, and saved notes.

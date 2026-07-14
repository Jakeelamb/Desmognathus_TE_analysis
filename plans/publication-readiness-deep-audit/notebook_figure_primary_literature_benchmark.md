# Primary-literature benchmark for the audit notebooks and figures

## Purpose

This document translates primary biological studies, primary methods papers, and Nature Portfolio reporting standards into concrete presentation requirements for six standalone *Desmognathus* audit notebooks:

1. phylogenetic-tree trimming and data completeness;
2. repeat discovery, repeat landscapes, diversity indices, and PCA;
3. the LTR terminal:internal deletion proxy;
4. cell/nucleus segmentation and morphometry;
5. relative nuclear IOD, nucleus size, and cell-size correlations across the phylogeny; and
6. phylogenetic path analysis.

The goal is not to imitate a journal's surface styling or to force the results to resemble a published result. The goal is to make every figure independently auditable: the biological unit, denominator, uncertainty, exclusions, and evidence boundary should be recoverable from the figure, caption, displayed table, and linked source data.

The most useful Nature-family precedents combine a dated tree with aligned trait/QC tracks, expose individual observations rather than bar-chart summaries, provide figure-specific source data, and place sensitivity analyses beside the headline result. These conventions are especially appropriate for this project because genomic resources, microscopy specimens, and tree tips are related at the species level but are not interchangeable physical samples.

## Portfolio-wide standard

### What every notebook should show before its biological results

Each notebook should begin with a compact audit header containing:

- the notebook purpose and the precise estimands it can support;
- the frozen 18-species panel and explicit exclusions;
- input paths, immutable file hashes, creation dates, and schema versions;
- software and package versions used to create the cached products;
- a statement that expensive upstream programs are **not rerun** by the notebook;
- the exact cached artifacts loaded instead;
- a data-unit statement distinguishing species, genomic resource, SRX run, technical retry, assembly, microscopy specimen, slide, image, cell, nucleus, and TE family/element;
- a visible release-gate box listing approved, exploratory, and blocked interpretations.

This implements the FAIR principle that data, algorithms, tools, and workflows should be findable and reusable, not merely the final plots ([Wilkinson et al. 2016](https://doi.org/10.1038/sdata.2016.18)). Nature Portfolio policy likewise treats access to data, code, and protocols as a condition of reproducibility ([Nature Portfolio reporting standards](https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards)).

### Figure and caption contract

Every inferential figure should have, either in its caption or an immediately adjacent audit table:

- exact `n` and what `n` represents;
- the statistical model and whether the test is one- or two-sided;
- effect estimate and interval, not only a significance marker;
- exact *P* value when one is computed, including nonsignificant values;
- the definition of every line, ribbon, box, violin, error bar, and point;
- units and transformation on every axis;
- the missing-data convention;
- the data exclusions and whether they were prespecified;
- the biological versus technical replicate definition;
- a direct link to the figure's machine-readable source table.

These requirements match Nature's instructions to define error bars, give exact sample sizes, identify replicates, and report exact significant and nonsignificant *P* values and test statistics ([Nature initial-submission guidance](https://www.nature.com/nature/for-authors/initial-submission)). Nature also recommends figures that remain comprehensible outside the immediate specialty and discourages unnecessary complexity ([Nature formatting guide](https://www.nature.com/nature/for-authors/formatting-guide)).

### Shared visual language

Use one semantic design system across notebooks:

- the same tip order everywhere, fixed by the final time tree;
- the same species colors only when color is genuinely needed;
- a colorblind-safe palette with redundant shape or line-style encoding;
- grey for contextual or unmodeled values, not for missingness;
- white or a hatched cell for missing data, with failed QC shown by a distinct symbol;
- solid lines for estimates supported by the specified model, dashed lines for uncertainty/sensitivity, and dotted lines only for explicit references such as 1:1 or zero;
- points for biological observations and lighter subordinate points for nested cells, images, or TE elements;
- no rainbow scales and no red-versus-green-only encoding;
- vector PDF/SVG for statistical figures and lossless PNG/TIFF for microscopy panels;
- scale bars on all microscopy images and one consistent physical calibration.

Nature explicitly asks authors to avoid rainbow scales, provide editable vector artwork, label axes with units, and design for color-vision accessibility ([Nature initial-submission guidance](https://www.nature.com/nature/for-authors/initial-submission); [Nature Communications checklist](https://www.nature.com/documents/ncomms-manuscript-checklist.pdf)).

### Shared table contract

Each displayed table should be a view of a saved CSV/TSV, not manually transcribed notebook output. Tables should have:

- one row per declared observational unit;
- stable identifiers rather than row numbers;
- species names, resource accessions, and specimen identifiers in separate fields;
- explicit units in column names or a linked dictionary;
- `NA` plus a reason field, never zero as a missing-data code;
- raw value, transformed value, estimate, uncertainty, and QC flag in separate fields;
- method/version/provenance columns for derived values;
- retry aliases separated from canonical biological-resource identity;
- a primary-key uniqueness assertion and a printed duplicate check;
- a hash or manifest entry for every exported table.

Small summary tables belong in the notebook; complete row-level tables should be downloadable and previewed with schema and validation summaries. Nature's FAIR and source-data conventions favor machine-actionable data underlying each display rather than values embedded only in a figure ([Wilkinson et al. 2016](https://doi.org/10.1038/sdata.2016.18)).

---

## Notebook 1 — Phylogenetic-tree trimming and data completeness

### Scientific purpose

This notebook should demonstrate exactly how the published time tree becomes the final 18-tip comparative tree. It should not infer a new tree. Its core claim is that every retained data row maps to exactly one intended tree tip, all exclusions and relabelings are explicit, branch lengths retain their time meaning, and tree uncertainty is carried forward.

### Required figure system

#### Figure 1. Tree-trimming audit flow

Show the source tree, name normalization, expert exclusions/reidentifications, exact-set intersection, final pruning, and validation gates as a left-to-right flow. Print counts at every transition and link the row-level crosswalk.

The accompanying table must list every source tip with `source_name`, `canonical_species`, `action`, `reason`, `evidence`, and `final_tip`. The trustworthy chromosome-level *fuscus* resource must remain distinguishable from discarded or reidentified resources. An exact set comparison should fail visibly on unmatched tips; silent dropping is unacceptable. This mirrors the explicit intersection behavior of comparative-data tools such as `geiger::treedata` ([official documentation](https://search.r-project.org/CRAN/refmans/geiger/html/treedata.html); [Pennell et al. 2014](https://doi.org/10.1093/bioinformatics/btu181)).

#### Figure 2. Before-and-after time trees

Place the untrimmed source tree and final tree side by side with:

- the same root orientation;
- a time scale in Ma;
- retained tips emphasized and removed tips muted;
- relabeled tips marked with a symbol;
- exact tip counts;
- source citation and tree filename;
- ultrametricity, duplicated-tip, zero/negative-branch, and root-to-tip checks printed below.

The primary Desmognathus tree source should be identified separately from the broad published uncertainty set. Stewart and Wiens provide a time-calibrated salamander phylogeny, calibration description, and supplementary tree products suitable for this provenance layer ([Stewart & Wiens 2025](https://doi.org/10.1016/j.ympev.2024.108272)).

#### Figure 3. Final time tree with aligned evidence tracks

Use one tip-aligned panel with columns for:

- genome/SRX resource and trust status;
- assembly availability and level;
- RepeatModeler/custom-library provenance;
- RepeatMasker output;
- dnaPipeTE output;
- divergence landscape;
- terminal:internal proxy availability;
- microscopy specimen and image count;
- morphometry availability;
- relative-IOD QC subset;
- inclusion in each path-model family.

Encode `present`, `missing`, `failed QC`, and `not applicable` separately. Do not display missing as zero. A strong Nature Communications precedent aligns a dated species tree with assembly size, TE fraction, BUSCO completeness, and gene-annotation tracks ([Cicconardi et al. 2023, Fig. 1](https://doi.org/10.1038/s41467-023-41412-5)). Another large comparative repeat study aligns time trees with genome size and repeat-content heatmaps ([Pasquesi et al. 2018](https://doi.org/10.1038/s41467-018-05279-1)).

#### Figure 4. Tree-uncertainty dashboard

Show:

- focal versus published-main topology;
- Robinson–Foulds distance across the 200 trees;
- branch-age or root-to-tip distributions where those vary;
- clade frequency for any variable split;
- frequency with which each downstream qualitative conclusion persists.

Do not show 200 overplotted trees as an unreadable cloud. A traitgram or a small number of consensus alternatives can be paired with quantitative frequency/interval panels. Comparative Nature work commonly summarizes analyses over tree sets using coefficient intervals and persistence frequencies rather than conditioning silently on one tree ([Comte et al. 2014](https://doi.org/10.1038/ncomms6053)).

### Reporting minimum

The notebook should print the complete kept/removed/relabelled table, branch-length units, rooting status, calibration source, number of alternative trees, and every exact-match assertion. A visually attractive tree without this audit layer is insufficient.

---

## Notebook 2 — Repeat discovery, landscapes, diversity, and PCA

### Scientific purpose

This notebook should explain how cached RepeatModeler, RepeatMasker, dnaPipeTE, and curated-library products became comparable species-level TE measurements. It should load their products, validate mass balance and identity, and regenerate inexpensive summaries and figures. It should not rerun expensive repeat discovery or masking.

### Required methods/provenance display

#### Figure 1. Cached repeat-analysis workflow

Use a pipeline diagram that identifies, for each species:

1. genomic resource and input reads/assembly;
2. custom-library discovery and curation;
3. RepeatMasker or dnaPipeTE quantification;
4. classification harmonization;
5. mass-accounting table;
6. landscape/diversity/PCA outputs.

For every expensive step, show the historical command or configuration, software/database version, completion marker, output directory, and hash. Technical retry suffixes belong in a provenance column and must converge on one canonical SRX/species identity.

dnaPipeTE was explicitly designed to assemble and quantify repeats from low-coverage read samples and produces repeat proportions and relative-age distributions ([Goubert et al. 2015](https://doi.org/10.1093/gbe/evv050)). RepeatModeler2 combines de novo family discovery with structural LTR discovery, but its output remains a candidate family library requiring clear classification/curation provenance ([Flynn et al. 2020](https://doi.org/10.1073/pnas.1921046117)).

### Required figure system

#### Figure 2. Mass-balance and method-concordance panel

For each species, show a 100% mass-accounting bar with mutually exclusive components such as known TE, unknown repeat, non-TE repeat, unresolved/conflicting classification, and unassigned/unmasked sequence. Print the row-sum residual and fail if it exceeds tolerance.

Beside it, show assembly-based versus read-based estimates for the same biological resource where both exist. Plot paired points, the 1:1 line, difference-versus-mean residuals, and label large disagreements. This is a method-comparison diagnostic, not a search for the highest correlation.

The Communications Biology salamander study reports the contribution of each repeat-mining method, the unknown fraction, low-coverage depth, and alternate library denominators rather than hiding these sources of variation ([Wang et al. 2025](https://doi.org/10.1038/s42003-025-08127-3)).

#### Figure 3. Phylogeny-aligned TE composition

Align the final time tree with:

- total identifiable repeat abundance;
- major TE-order composition;
- unknown/unclassified fraction;
- optionally the most abundant superfamilies as a heatmap.

Use a fixed ordering and palette for TE classes across every figure. Show absolute sampled-base percentages separately from within-TE proportions; do not mix their denominators. Pasquesi et al. provide a useful time-tree plus repeat-content heatmap precedent, while Wang et al. pair phylogeny, genome size, TE content, rank abundance, and proliferation history in one coherent figure system ([Pasquesi et al. 2018](https://doi.org/10.1038/s41467-018-05279-1); [Wang et al. 2025, Fig. 1](https://doi.org/10.1038/s42003-025-08127-3)).

#### Figure 4. Repeat landscapes

Use small multiples with:

- Kimura/substitution divergence on the x-axis;
- percentage of sampled bases or assembly bases on the y-axis;
- identical bins and class colors;
- shared y-scales whenever feasible, with any exceptions conspicuously marked;
- one panel per species in tree order;
- an adjacent table of included bases/hits and unresolved mass.

Caption the x-axis as divergence from a family consensus and a proxy for relative proliferation history, not calendar age. The axolotl Nature paper pairs repeat composition, TE-family phylogeny, family-length distributions, and Kimura-distance landscapes; it supplies quantitative source data and sample definitions ([Nowoshilow et al. 2018](https://doi.org/10.1038/nature25458)). The lungfish Nature paper similarly shows repeat composition, landscapes, and PCA and explicitly describes its divergence calculation and filtering ([Meyer et al. 2021](https://doi.org/10.1038/s41586-021-03198-8)).

#### Figure 5. Diversity profile, not one-number diversity

For every species, display at minimum:

- observed richness;
- Shannon entropy;
- exp(Shannon), if Hill numbers are used;
- Gini–Simpson, explicitly named as `1 - sum(p_i^2)`;
- inverse Simpson, if used, on a separately labeled scale;
- Pielou evenness;
- dominance/rank-abundance curve;
- the number of TE categories and abundance denominator.

Show species points, uncertainty from the defined resampling procedure, and a sensitivity for unknown/unclassified categories. Avoid bars containing only means. Shannon and Gini–Simpson have been used together in recent salamander TE studies, with TE superfamilies as community members and genomic percentages as abundance ([Wang et al. 2025](https://doi.org/10.1038/s42003-025-08127-3); [Decena-Segarra & Rovito 2024](https://doi.org/10.1093/molbev/msae225)). Haley and Mueller demonstrate that a biologically plausible result may be high and nearly unchanged diversity across genome expansion, reinforcing that the notebook must report the result rather than force a directional expectation ([Haley & Mueller 2022](https://doi.org/10.1007/s00239-022-10063-3)). Hill-number reporting makes the effective-number interpretation explicit and avoids treating indices with different units as interchangeable ([Jost 2006](https://doi.org/10.1111/j.2006.0030-1299.14714.x)).

#### Figure 6. CLR PCA scores, loadings, and robustness

The primary PCA display should include:

- CLR-transformed order-level composition;
- the stated zero-replacement rule and a zero-sensitivity result;
- species score plot with PC1/PC2 variance percentages;
- a separate loading plot or biplot with readable TE-order labels;
- the same species colors/order as the tree;
- leave-one-species-out score/loadings stability;
- phylogenetic PCA and 200-tree sensitivity in an adjacent panel;
- a clear designation of superfamily PCA as supplementary if sparsity makes it unstable.

Relative-abundance vectors are compositional, so ordinary PCA of raw percentages can produce closure artifacts. Log-ratio methods move ratio information to Euclidean space; the denominator and zero treatment are part of the estimand ([Gloor et al. 2017](https://doi.org/10.3389/fmicb.2017.02224)). A large Nature Communications benchmark likewise explains that compositional methods operate on ratios and defines the CLR denominator ([Nearing et al. 2022](https://doi.org/10.1038/s41467-022-28034-z)).

Do not add cluster ellipses or PERMANOVA labels unless a replicated grouping hypothesis exists and dispersion/permutation assumptions are satisfied. With one observation per species and related tips, an unrestricted species-label permutation is not evidence of biological clusters.

### Reporting minimum

Every repeat figure must identify its denominator, category level, input resource, and treatment of unknown/unclassified mass. Every diversity/PCA figure must expose the transformed source matrix and loading table. No rerun directory may appear as another biological species.

---

## Notebook 3 — LTR terminal:internal deletion-footprint proxy

### Scientific purpose and claim boundary

The notebook should use the term **terminal:internal deletion-footprint proxy** or **accumulated solo-LTR-producing deletion proxy**. A cross-sectional ratio of terminal and internal abundance is not an observed ectopic-recombination rate. It combines insertion history, element age, deletion, decay, assembly recovery, mapping, and family composition.

### Required figure system

#### Figure 1. Structural and computational estimand

Draw intact LTR structure with two terminals and an internal region, the solo-LTR product expected after intra-element recombination, and the exact numerator/denominator used in this project. Next to the schematic, show the cached workflow from candidate element/domain filtering through terminal/internal depth to species summary.

List aligner, mapping-quality rule, multimapper policy, length normalization, zero-depth policy, element-completeness rule, and coverage threshold. LTR_retriever demonstrates why LTR boundaries, terminal motifs, target-site duplications, internal domains, and counterfeit/truncated structures matter for high-confidence discovery ([Ou & Jiang 2018](https://doi.org/10.1104/pp.17.01310)).

#### Figure 2. Candidate-to-analysis attrition

For each species, show counts at every gate:

- candidate LTR sequences;
- terminal regions found;
- internal domains found;
- complete structures passing the declared domain rule;
- mapped terminal and internal targets;
- targets passing coverage;
- elements entering the species estimate.

Display absent assembly evidence as missing, not a ratio of zero. The notebook should fail closed below the prespecified evidence minimum.

#### Figure 3. Every element-level ratio

Plot all retained element/family ratios on a log y-axis, grouped by species in genome-size or phylogenetic order. Include:

- the individual points;
- median and bootstrap interval;
- exact element count per species;
- a dotted 1:1 reference;
- coverage or element-quality encoded redundantly;
- labels for high-influence elements.

Wang et al.'s direct salamander precedent displays each LTR/Gypsy family, species medians, a log scale, a 1:1 reference, and exclusion of a species with insufficient recovered consensus sequences ([Wang et al. 2025, Fig. 3](https://doi.org/10.1038/s42003-025-08127-3)). The wide within-species spread is scientifically important and should not be hidden behind one species bar.

#### Figure 4. Mapping and coverage sensitivity

Use paired/forest panels to show how each species estimate changes under:

- alternative defensible multimapper policies;
- minimum coverage thresholds;
- terminal definition choices;
- inclusion of five-domain versus stricter structural candidates;
- median versus family-weighted summaries;
- read/target resampling.

Also plot ratio versus total coverage and terminal versus internal depth to expose denominator instability. Frahry et al. used repeated low-coverage subsamples and family-level evidence in salamanders, providing the closest precedent for separating read-sampling from among-family uncertainty ([Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7)).

#### Figure 5. Structural validation where available

For chromosome/assembly-supported resources, show a small set of full-length and solo-LTR examples with boundaries, domains, target-site duplications, and read support. Summarize solo:intact or full-length:solo counts by scaffold/window only when those calls have passed structural validation. Whole-genome studies have related solo/full-length ratios to independently estimated recombination environments, a stronger structural estimand than pooled depth alone ([Ji & DeWoody 2016](https://doi.org/10.1007/s00239-016-9741-0)).

### Statistical reporting minimum

Families/elements within a species are nested evidence, not additional species. Any species-level predictor test must either use species summaries with uncertainty or a hierarchical model. Report intervals and influence, not only an omnibus ANOVA *P*. The conclusion should be phrased as more/fewer accumulated terminal-rich deletion footprints, never a measured contemporary rate.

---

## Notebook 4 — Cell/nucleus segmentation and morphometry

### Scientific purpose

This notebook must separate four layers that are easy to conflate:

1. image acquisition and staining;
2. cell and nucleus instance segmentation;
3. manual QC and cell–nucleus pairing;
4. biological morphometry and upper-tail selection.

A polished gallery is not validation by itself. Model identity, held-out truth, object-level error, accepted/rejected-object provenance, and nested sampling all need their own displays.

### Interactive species-mask viewers

The notebook should expose every species `index.html` through:

- a dropdown or link table keyed by canonical species;
- an inline `IFrame` when the notebook is served locally;
- an explicit relative path and “open in new tab” link;
- a static thumbnail/contact-sheet fallback for rendered notebook exports;
- a startup check that reports missing or broken viewer assets;
- a note explaining that browsers may block `file://` iframe content and that a local Jupyter/HTTP server is the supported route.

The viewer should retain image, cell-mask, nucleus-mask, overlay, object ID, paired-object ID, QC decision, and exclusion reason. It should not silently substitute a display-normalized image for the raw measurement image.

### Required figure system

#### Figure 1. Model and pipeline provenance

Show the exact production sequence from raw image to final linked measurements. Distinguish:

- pretrained/default model;
- any fine-tuned/custom model;
- the model actually used to generate the archived production masks;
- post-processing and pairing rules;
- manual decisions;
- the accepted-object table used downstream.

Cellpose 2.0 demonstrates the expected presentation: predictions from pretrained and successively trained models are shown on the same held-out image, with performance as a function of training amount rather than only a chosen final mask ([Pachitariu et al. 2022](https://doi.org/10.1038/s41592-022-01663-4)).

#### Figure 2. Segmentation QC atlas

For both cells and nuclei, show matched columns of:

- minimally processed raw image;
- independent truth mask;
- production prediction;
- color-coded overlay;
- accept/reject label and reason.

Include normal, overlapping, touching, broken, blurry, border-truncated, faint, unusually large, and irregular examples. Use the same LUT/contrast rule throughout, disclose any nonlinear processing, and include scale bars. Nature's microscopy policy requires whole-image adjustments, disclosure of thresholds/gamma/pseudocolor, acquisition settings, raw-image retention, and explicit treatment of images combined from different fields ([Nature image-integrity guidance](https://www.nature.com/nature/editorial-policies/image-integrity)).

#### Figure 3. Held-out performance dashboard

Report cells and nuclei separately, with metrics by tile/image, species, quality stratum, and size/shape bin:

- object precision, recall, and F1 at IoU 0.50 and 0.75;
- F1/AP as a curve across IoU thresholds;
- foreground IoU/Dice;
- object-count bias;
- foreground-area bias;
- split, merge, miss, and false-object rates;
- one-cell/one-nucleus pairing failures;
- bootstrap intervals over independent images/specimens.

Do not pool all pixels before scoring without also showing image-level distributions. Metric choice and aggregation can reverse model rankings; fixed IoU thresholds and per-image versus pooled aggregation must be named ([Hirling et al. 2024](https://doi.org/10.1038/s41592-023-01942-8)). The broader Metrics Reloaded framework recommends matching metrics to the exact segmentation task and failure mode rather than reporting one familiar overlap score ([Reinke et al. 2024](https://doi.org/10.1038/s41592-023-02150-0)).

Foundation/generalist microscopy models are evaluated across distinct image domains and with interactive correction/retraining workflows; their existence does not remove the need for local held-out validation ([Stringer et al. 2021](https://doi.org/10.1038/s41592-020-01018-x); [Archit et al. 2025](https://doi.org/10.1038/s41592-024-02580-4)).

#### Figure 4. Nested measurement distributions

Plot cells nested within image, slide/specimen, and species. A useful design is a muted cell-level raincloud or beeswarm with image/specimen estimates overlaid as larger points and species estimates with intervals. Print all four sample sizes. Do not let 50 cells visually imply 50 independent biological replicates.

Nature Methods emphasizes that population inference requires sampling the relevant level of variation and explicitly illustrates nested designs ([Altman & Krzywinski 2015](https://doi.org/10.1038/nmeth.3224)). CellProfiler's primary papers likewise make saved pipelines, traceable measurements, and segmentation validity central to reproducible downstream phenotyping ([Carpenter et al. 2006](https://doi.org/10.1186/gb-2006-7-10-r100); [McQuin et al. 2018](https://doi.org/10.1371/journal.pbio.2005970)).

#### Figure 5. Cell-selection sensitivity

For every species, show estimates and ranks under:

- all valid linked cells;
- image-balanced or specimen-balanced sampling;
- median and trimmed mean;
- upper quantiles;
- largest 50;
- largest 100.

Use a slopegraph/heatmap for rank changes and a forest plot for values/intervals. Also show eligible object count and the selected cells' percentile range. The largest-50 statistic must be called an upper-tail estimand, not typical cell size or a bias-free mean.

#### Figure 6. Pairing and review audit

Show one-cell/one-nucleus linkage rate, unlabeled/model-ranked versus manually reviewed fraction, reasons for removal, and sensitivity to inclusion of “maybe” objects. Preserve object IDs in the source table so every plotted point can be opened in the species viewer.

### Reporting minimum

The notebook must state whether the validation test species and slides overlap training data and whether any final-panel species are represented. A high-quality gallery cannot compensate for a non-independent test set. Any morphology conclusion should identify whether it concerns typical, balanced, or upper-tail erythrocytes.

---

## Notebook 5 — Relative IOD, nucleus size, and cell size across the phylogeny

### Scientific purpose and label boundary

This notebook should present the relationships among relative nuclear IOD, nuclear area, and cell area without relabeling relative IOD as absolute genome size. Unless DNA-stoichiometric staining, same-batch standards, optical-density conversion, and reference calibration are documented, axes and captions should use **relative nuclear IOD proxy** or **image-intensity phenotype**, not pg or C-value.

### Required figure system

#### Figure 1. Measurement anatomy and IOD decomposition

Show the equation used for each nucleus and a three-panel diagnostic:

- IOD versus nuclear area;
- IOD versus mean optical density;
- nuclear area versus mean optical density.

Color/shape by image or slide and print within-image and between-image variation. If the implementation is algebraically `area_px × mean_OD`, state that dependence directly. A path from IOD to nuclear area cannot then be interpreted as wholly independent biological evidence.

Feulgen image-analysis densitometry converts transmitted intensity to per-pixel optical density and sums it over a nucleus, with reference-standard calibration needed to obtain picograms ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)). The classic *Desmognathus* cytophotometric study used an internal DNA reference rather than arbitrary nuclear darkness ([Hally et al. 1986](https://doi.org/10.1007/BF00494802)).

#### Figure 2. Phylogeny plus three aligned traits

Pair the final time tree with columns for:

- relative nuclear IOD;
- nucleus area;
- linked cell area;
- measurement interval;
- cell/image/specimen counts;
- IOD QC status.

Use continuous, labeled scales and distinguish observed values from any reconstructed ancestral values. Do not color branches by reconstructed trait values unless the evolutionary model and reconstruction uncertainty are shown. The tree-plus-trait architecture follows strong comparative-genomics precedents in Nature Communications ([Cicconardi et al. 2023](https://doi.org/10.1038/s41467-023-41412-5); [Pasquesi et al. 2018](https://doi.org/10.1038/s41467-018-05279-1)).

#### Figure 3. Pairwise relationship matrix

Make three coordinated species-level panels:

1. relative IOD versus nucleus area;
2. relative IOD versus cell area;
3. nucleus area versus linked cell area.

Each panel should contain:

- all 18 species points, labeled or hoverable;
- measurement intervals on both axes where available;
- ordinary descriptive fit in neutral grey;
- PGLS fit and 95% interval in the primary color;
- effect estimate, 95% interval, exact *P*, Pagel's lambda, and exact species `n`;
- a leave-one-species-out coefficient inset;
- explicit identification of the microscopy estimator and IOD-QC subset.

The relevant salamander precedent reports log-transformed genome, nucleus, and cell-size relationships both before and after phylogenetic correction and shows that conclusions can change after correction ([Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)). That paper is more biologically relevant here than an unrelated high-profile visual template.

#### Figure 4. Raw versus phylogenetic effect forest

For each pairwise relationship and each defensible microscopy estimator, show the raw correlation/regression and PGLS coefficient side by side with 95% intervals. This makes the phylogenetic contribution visible rather than burying it in a caption.

#### Figure 5. Estimator and QC sensitivity

Use a heatmap or forest grid covering all morphology estimators and IOD subsets. Show coefficient sign, interval, sample size, and leave-one-out stability. Do not display only the selected-largest-50/image-QC-pass cell.

#### Figure 6. Calibration/agreement placeholder

If external calibrated C-values or future flow-cytometry values become available for the same specimens/batches, show:

- calibration curve with standards and unknowns;
- residuals and linearity/saturation diagnostics;
- Bland–Altman difference-versus-mean plot with bias and limits of agreement;
- prespecified acceptable agreement limits.

Correlation alone does not establish that two measurement methods are interchangeable ([Bland & Altman 1986](https://doi.org/10.1016/S0140-6736%2886%2990837-8)). Until such data exist, the panel should instead be a visible blocked-calibration gate, not an empty fitted genome-size graph.

### Reporting minimum

All inferential points are species-level, not cell-level. State the physical specimen mismatch between genomic and microscopy resources. Report the dependence of IOD on nuclear area and do not present a causal genome-size interpretation from relative intensity alone.

---

## Notebook 6 — Phylogenetic path analysis

### Scientific purpose

This notebook should compare a small, biologically prespecified set of directed acyclic graphs under phylogenetic correction. It should show whether the observed species-level covariance is consistent with these candidate graphs. It cannot prove causation, and it cannot promote relative IOD into absolute genome size.

### Required figure system

#### Figure 1. Candidate-DAG gallery

Plot every candidate family, including a full null model, with identical node coordinates and semantic arrow styling. Label each graph by its biological hypothesis, not only `model_1`. Print the testable d-separation basis set and flag saturated/no-basis-set models.

The `phylopath` primary paper identifies prespecification of the model set as the crucial confirmatory step and shows a model-set gallery before presenting rankings ([van der Bijl 2018](https://doi.org/10.7717/peerj.4718)). The original method frames phylogenetic path analysis as comparison among causal hypotheses using d-separation and PGLS, not data-driven arrow selection ([von Hardenberg & Gonzalez-Voyer 2013](https://doi.org/10.1111/j.1558-5646.2012.01790.x)).

#### Figure 2. Complete model-comparison display

Show every model's:

- number of independence claims;
- number of fitted parameters;
- Fisher's *C* and degrees of freedom;
- global-fit *P*;
- CICc;
- delta CICc;
- relative likelihood;
- model weight;
- fit/failure reason.

Use a table plus a weight/delta plot. Mark globally rejected models before interpreting their arrows. Treat delta CICc below two as competitive rather than proclaiming a unique winner.

#### Figure 3. Best-supported model set and coefficient forest

If one or more models remain competitive, show:

- standardized coefficients printed on arrows;
- edge width proportional to absolute coefficient;
- redundant color and line style for sign/support;
- a separate forest plot of standardized estimates and 95% intervals;
- direct, indirect, and total effects where identified;
- explicit model-averaging rule.

A recent Nature Communications phylogenetic-path study uses an a priori hypothesis diagram, coefficient labels, edge width, redundant sign encoding, a separate 95% interval effect panel, analysis over 100 trees, source data, and code ([Odom et al. 2025](https://doi.org/10.1038/s41467-025-60810-5)). Binary significant/nonsignificant arrows alone are not an adequate effect display.

#### Figure 4. Measurement-specification heatmap

For each candidate family, show top model, top weight, global-fit status, and key coefficient across every morphology estimator and IOD-QC subset. This is the clearest display of whether model ranking depends on one measurement decision.

#### Figure 5. Tree-uncertainty distributions

Across the main and 200 alternative trees, show:

- top-model frequency;
- model-weight distributions;
- standardized-coefficient distributions;
- fraction of globally rejected fits;
- topology/branch-length association with unstable outcomes, if present.

Do not summarize tree uncertainty as “results were similar” without showing the distributions. A 2026 Nature Communications analysis provides a strong current precedent: it reports raw/imputed sensitivity, two phylogenetic path frameworks, a 100-tree analysis, Fisher's *C*, model rankings, intervals, VIF checks, source data, and code ([Song et al. 2026](https://doi.org/10.1038/s41467-026-71230-4)).

#### Figure 6. Leave-one-species-out influence

Use a species-by-effect heatmap and/or coefficient forest showing:

- top-model changes;
- key coefficient changes;
- global-fit changes;
- delta-CICc changes;
- influential species labels.

The goal is not to delete influential species automatically. It is to expose whether a claim depends on one biological resource.

#### Figure 7. Exact-tree simulation calibration

For null and biologically relevant generating models, show:

- effect sizes spanning the observed interval;
- true-model recovery;
- false-selection rate under the null;
- coefficient bias and interval coverage;
- rejection/power rates;
- number of simulations and Monte Carlo uncertainty.

With only 18 tips, the structure of the actual tree may contain far less information than the nominal species count suggests. Phylogenetic Monte Carlo was developed specifically to assess comparative-model discrimination and false confidence on the actual phylogeny ([Boettiger et al. 2012](https://doi.org/10.1111/j.1558-5646.2011.01574.x)).

#### Figure 8. Release-gate summary

End with a machine-derived panel stating separately:

- computational fit status;
- stability across measurement specifications;
- stability across trees;
- leave-one-out stability;
- simulation calibration;
- upstream measurement validity;
- allowed claim language.

A robustly selected path involving relative IOD can still fail the biological genome-size gate. This separation must remain visible.

### Reporting minimum

Print the full basis sets and component PGLS equations, residual evolutionary model, phylogenetic parameter, transformations, sample size, model failures, and all candidate rankings. A nonsignificant global-fit test means a model was not rejected; it does not prove the graph true.

---

## Recommended figure order for tomorrow's scientific audit

For each notebook, use the same narrative rhythm:

1. **Question and evidence boundary** — what can this notebook establish?
2. **Provenance** — which cached artifacts and biological resources are loaded?
3. **Data-quality gate** — what passed, failed, or is missing?
4. **Primary descriptive figure** — all species or biological observations visible.
5. **Inferential figure** — effects and intervals with phylogenetic/statistical model stated.
6. **Sensitivity figure** — measurement, tree, threshold, or leave-one-out robustness.
7. **Source-data tables** — exact rows behind every figure.
8. **Verdict** — approved, exploratory, or blocked claim language.

This ordering resembles strong Nature-family articles in which the main integrated figure establishes sampling and genomic context before increasingly specific inferential panels. It is also much easier to audit live than a notebook organized by the chronology of scripts run months earlier.

## Figure practices to avoid

- Pie charts for 18-species TE comparisons; they make cross-species angle comparisons difficult. Use aligned bars or heatmaps.
- Bars with means and error bars but no points. Nature-family editorials explicitly encourage showing the observations ([“Show the dots in plots,” Nature Biomedical Engineering 2017](https://doi.org/10.1038/s41551-017-0079); [“Kick the bar chart habit,” Nature Methods 2014](https://doi.org/10.1038/nmeth.2837)).
- Significance stars without effect estimates, intervals, exact `n`, and exact *P* values.
- Regression at the cell or TE-element level when the stated claim is among species.
- Connecting genomic and microscopy records as if they were the same physical specimen.
- Calling a technical retry a replicate.
- Treating missing or failed-QC values as zero.
- Calling Kimura divergence an absolute insertion date.
- Calling terminal:internal depth an ectopic-recombination rate.
- Calling largest-50 morphometry an unbiased species mean.
- Calling relative nuclear IOD genome size without DNA-stoichiometric calibration.
- Coloring reconstructed branches without intervals and a declared evolutionary model.
- Showing one winning DAG without the candidate set, basis sets, global fit, and model weights.
- Hiding a failed validation gate because the resulting panel is visually unattractive.

## Compact citation map

| Notebook need | Primary or authoritative source | Reusable convention |
|---|---|---|
| Tree + aligned genomic/QC tracks | [Cicconardi et al. 2023](https://doi.org/10.1038/s41467-023-41412-5) | Dated tree aligned to assembly size, TE fraction, BUSCO, and annotation tracks |
| Tree + repeat heatmaps | [Pasquesi et al. 2018](https://doi.org/10.1038/s41467-018-05279-1) | Time tree, genome-size distributions, and tip-aligned repeat-content heatmaps |
| Salamander time-tree provenance | [Stewart & Wiens 2025](https://doi.org/10.1016/j.ympev.2024.108272) | Citable time tree, calibration description, supplementary tree set |
| Cached low-coverage repeat workflow | [Goubert et al. 2015](https://doi.org/10.1093/gbe/evv050) | Low-coverage repeat assembly, quantification, and divergence products |
| De novo TE-library workflow | [Flynn et al. 2020](https://doi.org/10.1073/pnas.1921046117) | Repeat-family discovery and structural LTR module provenance |
| Salamander TE integrated figure | [Wang et al. 2025](https://doi.org/10.1038/s42003-025-08127-3) | Phylogeny, genome size, composition, rank abundance, landscapes, diversity, and family-level T:I |
| Giant-genome TE landscape | [Nowoshilow et al. 2018](https://doi.org/10.1038/nature25458) | Repeat composition, TE phylogeny, length distributions, Kimura landscape, source data |
| Repeat landscape + PCA | [Meyer et al. 2021](https://doi.org/10.1038/s41586-021-03198-8) | Repeat composition, divergence landscape, PCA, explicit filtering |
| Salamander diversity alternatives | [Decena-Segarra & Rovito 2024](https://doi.org/10.1093/molbev/msae225); [Haley & Mueller 2022](https://doi.org/10.1007/s00239-022-10063-3) | Report diversity definitions and phylogeny without forcing direction |
| Compositional PCA | [Gloor et al. 2017](https://doi.org/10.3389/fmicb.2017.02224); [Nearing et al. 2022](https://doi.org/10.1038/s41467-022-28034-z) | Ratio-based analysis, explicit CLR denominator and zero handling |
| LTR structural qualification | [Ou & Jiang 2018](https://doi.org/10.1104/pp.17.01310) | Boundaries, motifs, TSDs, domains, intact/solo qualification |
| Salamander T:I sampling | [Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7) | Family distributions and low-coverage resampling |
| Structural solo/full-length analysis | [Ji & DeWoody 2016](https://doi.org/10.1007/s00239-016-9741-0) | Whole-genome structural counts and recombination comparison |
| Segmentation examples/training curve | [Pachitariu et al. 2022](https://doi.org/10.1038/s41592-022-01663-4) | Same-image masks across model versions and performance versus training amount |
| Segmentation metrics | [Hirling et al. 2024](https://doi.org/10.1038/s41592-023-01942-8); [Reinke et al. 2024](https://doi.org/10.1038/s41592-023-02150-0) | Multiple IoU thresholds, per-image aggregation, task-specific failure metrics |
| Microscopy integrity | [Nature image-integrity guidance](https://www.nature.com/nature/editorial-policies/image-integrity) | Raw-image retention, whole-image adjustments, LUT/settings/manipulation disclosure |
| Nested sampling | [Altman & Krzywinski 2015](https://doi.org/10.1038/nmeth.3224) | Show biological and nested technical variation separately |
| Salamander genome/nucleus/cell relationships | [Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010) | Raw and phylogenetically corrected relationships may differ |
| Image genome-size calibration | [Hardie et al. 2002](https://doi.org/10.1177/002215540205000601); [Hally et al. 1986](https://doi.org/10.1007/BF00494802) | Optical-density equation, reference calibration, same-assay standards |
| Measurement agreement | [Bland & Altman 1986](https://doi.org/10.1016/S0140-6736%2886%2990837-8) | Difference-versus-mean bias and limits of agreement, not correlation alone |
| Phylogenetic path implementation | [van der Bijl 2018](https://doi.org/10.7717/peerj.4718) | Candidate DAGs, full ranking, standardized paths, model averaging |
| Nature-family path presentation | [Odom et al. 2025](https://doi.org/10.1038/s41467-025-60810-5); [Song et al. 2026](https://doi.org/10.1038/s41467-026-71230-4) | A priori graphs, coefficients/intervals, model fit, tree and data sensitivity, source code/data |
| Comparative-model calibration | [Boettiger et al. 2012](https://doi.org/10.1111/j.1558-5646.2011.01574.x) | Simulation on the actual tree to measure model-selection power/error |

## Bottom line

The project does not need to rerun RepeatModeler, RepeatMasker, or dnaPipeTE during the scientific meeting. It does need to make their frozen inputs, commands, versions, outputs, identity mappings, and inexpensive downstream transformations inspectable. The strongest figure package is therefore not a set of polished conclusions alone. It is a linked system in which the tree fixes row order, every point opens into its source table or image viewer, every headline effect has an adjacent sensitivity display, and every blocked biological interpretation remains visibly blocked.

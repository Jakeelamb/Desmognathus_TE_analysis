# Upstream-methods primary-literature benchmark

## Scope and decision standard

This document benchmarks four upstream analysis domains needed before a phylogenetic path analysis of *Desmognathus* genome size, transposable elements (TEs), erythrocyte size, and nucleus size can be considered publication-ready. It is a literature benchmark, not a certification of the repository's present implementation.

The central distinction is:

- **Minimum release gate:** evidence needed to make a narrowly worded, reproducible claim without a known fatal validity problem.
- **Gold standard:** additional work that materially strengthens inference, robustness, and comparability with leading studies.

The benchmark does **not** require results to have the same direction as published results. Recent salamander studies have reported increasing, decreasing, and nearly unchanged TE diversity across different taxonomic scales, despite using related diversity metrics ([Decena-Segarra and Rovito 2024](https://doi.org/10.1093/molbev/msae225); [Haley and Mueller 2022](https://doi.org/10.1007/s00239-022-10063-3); [Wang et al. 2025](https://doi.org/10.1038/s42003-025-08127-3)). Matching quality means matching the clarity of the estimand, controls, phylogenetic treatment, uncertainty, and reproducibility—not forcing agreement with a favored biological result.

## Executive release gates

The path analysis should remain downstream-blocked until all of the following are true:

1. The ectopic-recombination analysis is explicitly labeled as a **cross-sectional proxy for accumulated LTR deletion footprints**, not a directly observed recombination or deletion rate.
2. Every LTR terminal:internal estimate has a documented sequence definition, denominator, mapping policy, family-level sample size, uncertainty interval, and coverage/quality sensitivity analysis.
3. TE diversity is calculated from one harmonized taxonomic level and denominator across species; Shannon, Gini-Simpson, richness, and dominance/evenness are reported together or their omission is justified.
4. PCA/ordination uses a transformation and distance appropriate for compositional abundance data. Raw percentage PCA is not the primary ordination.
5. Any PERMANOVA has a genuine replicated grouping or predictor design, reports effect size and dispersion diagnostics, and does not use unrestricted species-label permutations as though related species were exchangeable.
6. Every biological species is represented once in inferential TE analyses. Retried dnaPipeTE directories or run labels with appended suffixes such as `2` or `3` are recorded as **technical execution aliases**, not biological replicates.
7. The time tree has a citable source, branch lengths in time units, an exact tip-to-species crosswalk, and no silent taxon substitution or arbitrary polytomy resolution.
8. Phylogenetic path models are a small, a priori set of directed acyclic graphs; each has a nonempty testable basis set or is explicitly identified as saturated; global fit, CICc/model weights, component-model diagnostics, standardized coefficients, and uncertainty are reported.
9. Tree uncertainty is propagated over a distribution of plausible time trees, not represented solely by one consensus tree.
10. Cell and nucleus segmentation is validated on held-out images/slides/species using both pixel overlap and object-level detection/split/merge metrics, with error stratified by image quality.
11. Cells are not treated as independent species replicates. Cell-level measurements remain nested within image, slide, specimen, and species, and species-level uncertainty is passed downstream.
12. A “largest 50 cells” summary is labeled as an **upper-tail estimand**. It is not described as an unbiased estimate of typical cell or nucleus size, and it is accompanied by balanced random/all-valid-cell sensitivity analyses.
13. Nuclear image intensity is called genome size only when it comes from a DNA-stoichiometric assay with a same-batch or same-slide reference standard and an optical-density calculation. Arbitrary grayscale darkness, fluorescence intensity without calibration, or IOD from non-Feulgen staining is not a genome-size measurement.
14. Image-derived genome-size estimates are validated against an established method or defensible reference specimens; agreement is evaluated as agreement, not merely correlation.

---

## 1. LTR deletion, solo-LTR formation, and ectopic-recombination inference

### 1.1 What the analysis can estimate

Homologous recombination between the two LTRs of an element can delete the internal region and one LTR, leaving a solo LTR. Consequently, an excess of terminal sequence relative to internal sequence can record accumulated deletion footprints ([Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7); [Ji and DeWoody 2016](https://doi.org/10.1007/s00239-016-9741-0)). A snapshot of terminal and internal abundance does not observe events through time, however. It combines insertion history, family age, deletion, subsequent decay, detection, and sampling. The primary defensible term is therefore **terminal:internal deletion proxy** or **accumulated ectopic-recombination-mediated deletion footprint**, not “ectopic recombination rate.”

The closest direct salamander benchmark is Frahry et al. They estimated terminal:internal (T:I) ratios for abundant LTR families in low-coverage salamander data, tested 1% subsamples from complete genomes five times, and complemented shotgun estimates with manual inspection of salamander BACs ([Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7)). Their analysis also showed large among-family variation within a genome, making family identity and family-level replication essential rather than incidental.

The strongest recent salamander implementation selected genomic contigs longer than 3 kb annotated as LTR/Gypsy, identified terminal and internal regions with LTRpred, manually confirmed internal protein domains and terminal repeats, removed incomplete structures, and retained 12–25 contigs per analyzed species. Reads were mapped in two configurations, terminal and internal depth was length-normalized, and a T:I value was calculated per contig ([Wang et al. 2025, Methods: “Estimation of ectopic recombination-mediated deletion of LTR elements”](https://www.nature.com/articles/s42003-025-08127-3#Sec19)). That study excluded a species with only one recovered LTR/Gypsy contig rather than presenting a nominal estimate from inadequate structural replication ([Wang et al. 2025](https://doi.org/10.1038/s42003-025-08127-3)).

Whole-assembly structural counting provides a complementary estimand. Ji and DeWoody compared full-length:solo-LTR ratios among chromosomes and 1-Mb windows and tested their relationship with independent recombination maps ([Ji and DeWoody 2016](https://doi.org/10.1007/s00239-016-9741-0)). High-confidence intact-element discovery can use LTR boundaries, target-site duplications, terminal motifs, and internal protein domains to reject false structures; LTR_retriever was explicitly benchmarked on these features and documents solo, truncated, nested, and counterfeit structures ([Ou and Jiang 2018](https://doi.org/10.1104/pp.17.01310)).

### 1.2 Minimum release gate

The repository must add or expose all of the following in the shareable notebook and machine-readable outputs:

1. **A formal estimand and equation.** Define whether the quantity is terminal depth divided by internal depth, summed terminal bases divided by internal bases, solo:intact counts, or full-length:solo counts. State whether “terminal” combines both LTRs or uses one consensus terminal. Ratios from different definitions must not be placed on one scale without a conversion argument.
2. **Per-family/per-contig evidence.** Retain the individual family/contig ratios and identifiers, not only a pooled species value. Show the number recovered, number passing structure filters, and number used for each species. Wang et al. observed roughly 10- to 30-fold within-species variation among LTR/Gypsy families, so a species mean without the distribution hides biologically important uncertainty ([Wang et al. 2025](https://www.nature.com/articles/s42003-025-08127-3#Sec6)).
3. **Structural qualification.** Record minimum contig length, LTR boundary method, internal-domain evidence, completeness filter, nested-insertion treatment, and manual-review outcome. The rationale is directly supported by the false-positive and truncated-element benchmarks in [Ou and Jiang (2018)](https://doi.org/10.1104/pp.17.01310).
4. **Mapping qualification.** Report aligner/version, preset, minimum aligned length, mapping-quality threshold, paired/unpaired treatment, multimapper policy, and depth-normalization formula. Repeat-rich targets make these choices part of the estimand, not implementation trivia. At minimum, repeat the calculation under two defensible mapping policies, as Wang et al. did with unpaired versus merged-plus-unmerged reads ([Wang et al. 2025](https://www.nature.com/articles/s42003-025-08127-3#Sec19)).
5. **Coverage and recovery gate.** Establish an a priori minimum amount of usable family/contig evidence and fail closed below it. Do not substitute zero, `NA`, or a one-family point estimate for inadequate recovery. The exact threshold should be justified from this dataset's resampling behavior rather than copied mechanically from another organism.
6. **Technical-retry identity.** Map every suffixed or restarted dnaPipeTE execution to one canonical SRX resource and one biological species. Preserve the retry label in provenance, but collapse it before biological statistics. The dnaPipeTE method was designed to quantify repeats from a low-coverage genomic sample; repeated executions of the same input are computational repeatability checks, not new organisms ([Goubert et al. 2015](https://doi.org/10.1093/gbe/evv050)).
7. **Uncertainty from reads and families.** Report a species median and a bootstrap interval that resamples families/contigs; separately show read/subsampling sensitivity. Frahry et al.'s repeated 1% subsamples are the closest published precedent for evaluating low-coverage sampling behavior ([Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7)).
8. **Nested inference.** Do not turn multiple families from one species into independent species replicates. If testing a species-level predictor, use a hierarchical model with family/contig nested within species or reduce to a species estimate with uncertainty. Repeated observations from one unit violate ordinary independence and can strongly inflate false-positive rates ([Aarts et al. 2014](https://doi.org/10.1038/nn.3648)).
9. **Claim restraint.** Conclusions must say “consistent with more/fewer accumulated solo-LTR-producing deletions” or equivalent. A higher T:I ratio is not direct measurement of the contemporary molecular rate, and a nonsignificant relationship is not proof of equal rates.

### 1.3 Gold standard

- Run a second, assembly-structural analysis on every suitable long-read/chromosome-level assembly: identify intact, truncated, and solo LTRs; verify LTR boundaries and target-site duplications; and report solo:intact ratios by superfamily. Use high-confidence structural screening such as the criteria benchmarked by [Ou and Jiang (2018)](https://doi.org/10.1104/pp.17.01310).
- Validate the read-depth T:I metric against direct solo:intact counts in at least the reliable chromosome-level *D. fuscus* assembly and any other high-contiguity focal assembly. Agreement between the two approaches is much stronger evidence than either alone.
- Calibrate mapping and detection bias with in silico intact and solo elements spanning the observed lengths, divergence, copy numbers, and read depths. Recover the known ratios through the complete pipeline.
- Screen paired LTRs for gene conversion before interpreting LTR-LTR divergence as insertion age. Gene conversion can homogenize terminal sequences and bias ages downward ([Kijima and Innan 2010](https://doi.org/10.1093/molbev/msp295)).
- Where an independent recombination map or chromosome-scale proxy exists, test whether solo/intact structure varies with recombination environment, paralleling the direct map-based design of [Ji and DeWoody (2016)](https://doi.org/10.1007/s00239-016-9741-0).
- Propagate family and measurement uncertainty into downstream phylogenetic models rather than treating a species ratio as measured without error; comparative methods can be biased by ignored within-species or measurement variance ([Ives, Midford, and Garland 2007](https://doi.org/10.1080/10635150701313830)).

### 1.4 Required figures and tables

1. **Family-level deletion-proxy plot:** one point per family/contig; species median and bootstrap interval; log-scaled ratio; sample size printed for every species; a reference line whose biological meaning matches the exact ratio definition.
2. **Recovery/QC panel:** input candidates, complete candidates, mapped candidates, and final candidates per species, plus usable read depth.
3. **Sensitivity forest plot:** change in species estimate under mapping mode, minimum-alignment threshold, multimapper policy, family filter, and bootstrap scheme.
4. **Structural validation panel:** representative intact, truncated, and solo calls with LTRs, internal domains, target-site duplications, and read-depth tracks.
5. **Assembly-versus-read comparison:** paired T:I and solo:intact evidence for species with suitable assemblies, with disagreement discussed rather than hidden.

---

## 2. TE diversity indices, PCA/ordination, and PERMANOVA

### 2.1 Leading salamander coverage

Recent comparative salamander papers treat TE superfamilies as community categories and their genomic proportions as abundances. Decena-Segarra and Rovito calculated Shannon and Gini-Simpson diversity for 16 bolitoglossine salamanders and tested relationships with genome size using ordinary and phylogenetic regression ([Decena-Segarra and Rovito 2024](https://doi.org/10.1093/molbev/msae225)). Haley and Mueller used Shannon and Simpson-family metrics in four *Plethodon* species and emphasized that the indices weight rare and dominant groups differently ([Haley and Mueller 2022](https://doi.org/10.1007/s00239-022-10063-3)). Wang et al. analyzed 84 vertebrates, used the top ten superfamilies because rare-category abundance was uncertain, calculated both Shannon and Gini-Simpson indices, and used within-clade PGLS against genome size ([Wang et al. 2025, Methods: “Comparison of genomic TE community diversity across genome sizes”](https://www.nature.com/articles/s42003-025-08127-3#Sec18)).

These leading salamander analyses did **not** need PCA or PERMANOVA to support the primary diversity-versus-genome-size claim. Their core inferential structure was diversity indices plus phylogenetically controlled regression. PCA is valuable for visualizing which TE groups drive compositional differences; PERMANOVA is appropriate only for a prespecified multivariate hypothesis with valid replication and permutation units.

### 2.2 Diversity-index minimum gate

1. **One community definition.** Fix the category level (preferably superfamily for cross-species low-coverage comparisons), the accepted TE classes, the treatment of “unknown/unclassified,” and the abundance denominator. Report both classified-repeat composition and absolute-genome-load views when possible; they answer different questions.
2. **Mass accounting.** For every species, show total analyzed bases/reads, classified TE bases, unresolved repeat bases, non-TE bases, and the sum of all diversity proportions. A diversity index calculated after silently discarding a different unresolved fraction in each species is not directly comparable.
3. **Multiple diversity aspects.** Report richness, Shannon entropy, Gini-Simpson, and a dominance/evenness measure. State the exact Simpson convention; `sum(p_i^2)`, `1-sum(p_i^2)`, and its reciprocal have opposite or differently scaled interpretations. Convert to Hill effective numbers when a common “number of equally abundant categories” scale improves interpretation ([Jost 2006](https://doi.org/10.1111/j.2006.0030-1299.14714.x)).
4. **Rare-category sensitivity.** Recalculate after a documented prevalence/abundance filter and under top-*k* versus all-classified categories. Wang et al. limited the primary calculation to the top ten superfamilies because rare groups were uncertain ([Wang et al. 2025](https://www.nature.com/articles/s42003-025-08127-3#Sec18)); Haley and Mueller explain why Shannon is more sensitive to rare categories than Gini-Simpson ([Haley and Mueller 2022](https://doi.org/10.1007/s00239-022-10063-3)).
5. **Sampling-depth sensitivity.** Use repeated equal-depth subsampling or another coverage-standardization analysis and show stability of ranks and regression slopes. dnaPipeTE was developed for low-coverage data but loses sensitivity for old/divergent families, so equal nominal coverage does not eliminate detection bias ([Goubert et al. 2015](https://doi.org/10.1093/gbe/evv050)).
6. **Phylogenetic regression for claims.** Test diversity-versus-genome-size relationships with PGLS or another explicitly phylogenetic model, reporting effect estimate, confidence interval, residual diagnostics, and the fitted phylogenetic correlation parameter. Do not rely on ordinary Pearson correlation as the publication claim.

### 2.3 Compositional PCA/ordination minimum gate

TE relative-abundance vectors are compositions constrained to a constant sum. Standard Euclidean PCA on raw percentages can create covariance induced by closure rather than biological co-variation. Aitchison's log-ratio PCA was developed specifically for compositional data ([Aitchison 1983](https://doi.org/10.1093/biomet/70.1.57)).

The primary ordination should therefore be one of:

- centered-log-ratio (CLR) PCA with an explicit zero-replacement rule and a sensitivity analysis across defensible replacement values; or
- Hellinger-transformed abundance PCA, which was developed to make community-composition data suitable for Euclidean ordination ([Legendre and Gallagher 2001](https://doi.org/10.1007/s004420100716)); or
- a PCoA on a prespecified ecological distance, with any negative-eigenvalue correction documented.

Publication output must include the complete loading matrix, proportion of variance for every retained axis, a scree plot, species scores, the transformation and centering/scaling choices, and the list of categories included. Axis direction is arbitrary, so claims must be based on loadings and relative positions rather than the sign of PC1.

For interspecific data, phylogeny also affects covariance among traits. If the PC scores become inferential variables, include a phylogenetic PCA or a direct multivariate phylogenetic model as a sensitivity analysis; ignoring phylogeny during interspecific PCA or preprocessing can alter comparative inference ([Revell 2009](https://doi.org/10.1111/j.1558-5646.2009.00804.x)). Keep ordinary/compositional PCA as descriptive if its purpose is visualization only.

### 2.4 PERMANOVA minimum gate

PERMANOVA partitions variation in a chosen dissimilarity space and obtains significance by permutation; its inference depends on the distance and exchangeability of the permuted units ([Anderson 2001](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x)). Therefore:

1. Specify the biological hypothesis, response matrix, distance, formula, term order/type of sums of squares, permutation scheme, permutation count, and seed.
2. Report pseudo-*F*, partial or sequential *R*² as applicable, permutation *P*, and confidence/stability under leave-one-species-out analysis. A small *P* without effect size is incomplete.
3. Test multivariate dispersion with a distance-matched procedure and show a dispersion plot. PERMANOVA location effects can be difficult to distinguish from dispersion effects, for which Anderson developed a direct test ([Anderson 2006](https://doi.org/10.1111/j.1541-0420.2005.00440.x)).
4. Do not run a group PERMANOVA when a group has no genuine replication or when group membership is nearly identical to one clade. Related species are not freely exchangeable permutation units.
5. For the primary interspecific test, prefer a multivariate phylogenetic linear model or residual-randomization procedure conditioned on the phylogeny. Distance-based PGLS supports high-dimensional multivariate responses ([Adams 2014](https://doi.org/10.1111/evo.12463)), and RRPP has appropriate behavior for multivariate phylogenetic regression and ANOVA ([Adams and Collyer 2018](https://doi.org/10.1111/evo.13492)). A conventional PERMANOVA may remain a clearly labeled nonphylogenetic sensitivity analysis.

### 2.5 Gold standard

- Fit the diversity and compositional models across the posterior/time-tree set and summarize between-tree variation.
- Refit after equal read-depth subsampling, omission of each species, alternative unclassified-repeat handling, alternative taxonomic aggregation, and alternative zero replacement.
- Use multivariate phylogenetic regression directly on transformed TE composition, so the primary inference is not filtered through axes selected solely to maximize total variance.
- Propagate uncertainty in TE proportions from read resampling or dnaPipeTE iterations to diversity and ordination scores.
- Provide the exact community matrix used in every figure and test. Each technical retry alias must resolve to a single canonical resource row before matrix construction.

### 2.6 Required figures and tables

1. **Mass-accounted composition bars:** identical category order and denominator across species; unresolved and unclassified mass visible.
2. **Diversity-versus-genome-size panels:** Shannon, Gini-Simpson, richness, and dominance/evenness; PGLS estimate and 95% interval; species labels; raw points; leave-one-out slope range.
3. **Compositional biplot:** species scores, loadings for influential TE groups, variance labels, and phylogeny/taxonomic coloring that does not imply an untested group effect.
4. **Scree and loading table:** all axes, not only the visually convenient first two.
5. **Distance/dispersion panel:** PCoA plus distance-to-centroid distributions for every PERMANOVA factor.
6. **Sensitivity heat map:** changes in diversity rank, PC scores, and regression slope across depth, taxonomy, unresolved-mass, and zero-replacement choices.

---

## 3. Time-calibrated phylogeny and phylogenetic path analysis

### 3.1 Tree standard

The strongest currently available broad salamander resource contains 765 species and 503 genes and supplies **200 time-calibrated trees for comparative analyses**, specifically enabling propagation of phylogenetic uncertainty ([Stewart and Wiens 2025](https://doi.org/10.1016/j.ympev.2024.108272)). The *Desmognathus* nuclear phylogenomic literature also documents gene flow, reticulation, and mitochondrial/nuclear conflict, so taxonomy matching cannot be reduced to string cleanup ([Pyron et al. 2020](https://doi.org/10.1016/j.ympev.2020.106751); [Pyron et al. 2022](https://doi.org/10.1002/ece3.8574)).

Minimum tree provenance is: source publication and archive, exact tree filename/hash, whether it is a maximum-clade-credibility, consensus, or posterior draw, branch-length units, rooting, pruning commands, dropped tips, substituted tips, resolved polytomies, and the accession/specimen-to-tip crosswalk. If a time tree is newly estimated, calibration identity and placement must be justified and reproducible; fossil-calibration best practices require explicit specimen, apomorphy, phylogenetic placement, locality/stratigraphy, and age evidence ([Parham et al. 2012](https://doi.org/10.1093/sysbio/syr107)).

No species should be replaced by a congener simply to complete the panel without a labeled sensitivity analysis. No zero-length branch, duplicate tip, non-ultrametric tree, or arbitrary random polytomy resolution should pass silently.

### 3.2 Path-analysis minimum gate

Phylogenetic confirmatory path analysis tests whether observations are consistent with prespecified causal graphs while accounting for phylogenetic non-independence; ordinary path analysis can be spurious when shared ancestry is ignored ([von Hardenberg and Gonzalez-Voyer 2013](https://doi.org/10.1111/j.1558-5646.2012.01790.x)). The method remains confirmatory evidence among observational causal hypotheses, not experimental proof of causation.

Required components are:

1. **A priori model set.** Define a small set of biologically justified directed acyclic graphs before inspecting their path significance. Include a defensible null/alternative structure. `phylopath` identifies model definition as the crucial confirmatory step and implements the original PPA procedure reproducibly ([van der Bijl 2018](https://doi.org/10.7717/peerj.4718)).
2. **Testable basis sets.** Print every conditional-independence claim for every graph. A saturated model with no missing paths has no independence claims and cannot receive an evidential global-fit pass; it may estimate paths but cannot test whether the graph omits important relationships.
3. **Global fit.** Report the number of independence claims, every component *P*, Fisher's *C*, degrees of freedom, and global *P*. The d-separation approach combines the basis-set probabilities in Fisher's *C*, evaluated against a chi-square distribution ([Lefcheck 2016](https://doi.org/10.1111/2041-210X.12512)). A rejected global model cannot be rescued by highlighting significant individual arrows.
4. **Model comparison.** Report CICc, delta CICc, relative likelihood, and model weight for all candidate graphs. The PPA implementation uses small-sample-corrected C-statistic information criterion and exposes the number of claims and parameters ([van der Bijl 2018](https://doi.org/10.7717/peerj.4718)). When several graphs remain plausible, report model-averaged coefficients or preserve model uncertainty rather than declaring one graph uniquely true.
5. **Component-model evidence.** For every path and independence claim, report transformation, coefficient, standard error/interval, test statistic, *P*, sample size, estimated phylogenetic parameter such as Pagel's lambda, residual diagnostics, and influential species.
6. **Scale and collinearity.** Log-transform size variables when biologically/statistically justified, standardize coefficients on a declared scale, inspect nonlinearities, and report collinearity. Cell size and nucleus size are mechanically and biologically coupled; TE subclasses can also sum to a whole. A path coefficient does not remove those constraints.
7. **Measurement error.** Supply species-level sampling uncertainty for image-derived traits and TE proxies. Phylogenetic comparative methods that ignore within-species and measurement variance can yield biased parameters and overstated precision ([Ives, Midford, and Garland 2007](https://doi.org/10.1080/10635150701313830); [Silvestro et al. 2015](https://doi.org/10.1111/2041-210X.12337)).
8. **Small-sample power.** With a small species panel, use simulation on the actual tree and missing-data pattern to estimate power, false-selection rate, coefficient bias, and interval coverage for the candidate graphs. Information-criterion selection can have high error when a phylogeny contains insufficient information, and phylogenetic Monte Carlo was developed to diagnose this ([Boettiger, Coop, and Ralph 2012](https://doi.org/10.1111/j.1558-5646.2011.01574.x)). Do not rely on a universal observations-per-parameter rule.
9. **Influence analysis.** Refit after omitting each species, and separately after removing any taxonomically or measurement-problematic resource. Show whether path sign, support, and model rank are stable.

### 3.3 Tree-uncertainty minimum and gold standard

At minimum, repeat the complete model comparison and coefficient estimation over a representative set of time trees and report how often each model ranks first, the distribution of delta CICc/weights, and the distribution of each standardized path. Methods that condition on one fixed tree can understate uncertainty; comparative models can integrate a posterior tree set directly ([de Villemereuil et al. 2012](https://doi.org/10.1186/1471-2148-12-102)). A general frequentist combination using Rubin's rules is also available for phylogenetic and species-sampling uncertainty ([Nakagawa and de Villemereuil 2019](https://doi.org/10.1093/sysbio/syy089)).

The gold standard is a hierarchical Bayesian phylogenetic structural-equation model that jointly propagates trait error and tree uncertainty. PhyBaSE was designed to include trait and phylogenetic uncertainty and showed greater model-discrimination power than classic maximum-likelihood PPA in its simulations ([von Hardenberg and Gonzalez-Voyer 2025](https://doi.org/10.1111/2041-210X.70044)). This is an enhancement, not a substitute for a transparent classical PPA benchmark.

### 3.4 Required figures and tables

1. **Time tree with data completeness:** tips, time scale, support/uncertainty, and a matrix showing which specimen/resource supplies each trait.
2. **All candidate DAGs:** identical node placement, clear biological hypotheses, and no arrow added after viewing significance without being labeled exploratory.
3. **Model comparison table/plot:** global *C*, global *P*, CICc, delta, weight, number of claims, and number of parameters for every graph.
4. **Final path diagram:** standardized estimates and 95% intervals; line width by magnitude, not binary significance alone; direct, indirect, and total effects reported.
5. **Tree-uncertainty forest/violin:** coefficient and model-weight distributions across trees.
6. **Influence dashboard:** leave-one-species-out coefficients, global fit, and model rank.
7. **Simulation calibration:** power and false-selection rate under effect sizes spanning the observed intervals.

---

## 4. Erythrocyte cell/nucleus segmentation and image-derived genome size

### 4.1 Separate the two measurement problems

Morphometry and DNA-content densitometry use the same biological objects but are not the same assay.

- **Cell/nucleus morphometry** requires calibrated spatial scale and valid instance boundaries.
- **Genome-size densitometry** additionally requires DNA-stoichiometric staining, controlled illumination/exposure, optical-density conversion, and a reference standard of known DNA content.

CellProfiler's original methods paper states that segmentation accuracy determines downstream cell-measurement accuracy ([Carpenter et al. 2006](https://doi.org/10.1186/gb-2006-7-10-r100)); its modern implementation emphasizes saved, reusable pipelines and public ground-truth images for reproducibility ([McQuin et al. 2018](https://doi.org/10.1371/journal.pbio.2005970)). A custom model is acceptable only if its performance is demonstrated on this image domain.

### 4.2 Segmentation minimum gate

1. **Immutable image provenance.** Record specimen, species, slide, field, channel, objective, camera, pixel size, exposure, illumination, stain batch, acquisition date, and raw-file hash. Preserve raw images separately from display-normalized images.
2. **Leakage-proof splits.** Training, validation, and held-out test sets must split at least by slide/specimen, never by random cells or tiles from the same image. A stronger test holds out entire species or image-quality strata to measure domain shift.
3. **Independent ground truth.** Curate a blinded test set with cell and nucleus instance masks. Include overlapping, border-touching, broken, blurry, atypically shaped, faint, and densely packed examples in proportions that can be reported.
4. **Pixel and object metrics.** Report Dice or IoU for boundary overlap **and** object precision, recall, F1/average precision, count error, split rate, merge rate, missed-object rate, and false-object rate. Pixel overlap alone can rank an instance-segmentation result incorrectly when objects are missed, merged, or split; metric selection must match the biological task ([Reinke et al. 2024](https://doi.org/10.1038/s41592-023-02150-0)).
5. **Stratified performance.** Show metrics by species, slide, image quality, overlap status, size bin, and shape bin, with uncertainty. A single pooled Dice score can conceal failure in the exact large or irregular cells selected for biological inference.
6. **Pairing constraints.** Quantify one-cell/one-nucleus pairing failures, nuclei outside cells, multiple nuclei assigned to one cell, and cells lacking a nucleus. Make exclusion reasons machine-readable.
7. **Human reproducibility.** Have a second annotator review a stratified subset and report agreement or adjudication. Segmentation is the most consequential image-analysis step, and analyst choices can materially vary even on identical images ([Carpenter et al. 2006](https://doi.org/10.1186/gb-2006-7-10-r100)).
8. **Versioned model and pipeline.** Archive code, model weights, environment, inference parameters, post-processing, thresholds, and the exact accepted/rejected-object table. Cellpose itself was validated on a diverse held-out image collection and makes training data and code available, illustrating the expected reproducibility standard ([Stringer et al. 2021](https://doi.org/10.1038/s41592-020-01018-x)).

### 4.3 Cell-selection and statistical-unit gate

Selecting the largest 50 accepted cells per species is **not** a bias-free estimate of typical erythrocyte size. By definition it estimates an upper tail. If species have different total detected-cell counts, the largest 50 also represent different quantiles; and because merged or overlapping masks tend to be large, the selection can amplify residual segmentation errors.

Therefore:

1. Name the estimand explicitly, for example “mean area of the 50 largest QC-passing cells” or “upper-tail cell area.”
2. Make a central-tendency estimand—median, trimmed mean, or model-estimated specimen mean from all/balanced valid cells—the primary morphometric trait unless a specific biological hypothesis predicts the upper tail.
3. Show prespecified sensitivity for all valid cells, a balanced random sample, median, trimmed mean, 75th/90th percentile, largest 50, and largest 100.
4. Hold the number of eligible fields/cells balanced across species or use a quantile-based rule. Report how many objects were available before selection.
5. Preserve cell-nucleus pairing through every selection; do not independently select large nuclei after selecting cells.
6. Treat cells as nested observations. The independent biological sampling levels are specimens and species, not thousands of segmented cells. Multilevel methods are designed for this dependency, whereas treating nested observations as independent inflates type-I error ([Aarts et al. 2014](https://doi.org/10.1038/nn.3648)).
7. Produce specimen-level estimates and standard errors before species aggregation. If only one specimen represents a species, state that intraspecific variance is not estimable and treat that species' trait uncertainty conservatively.

This differs from standard erythrocyte morphometry, which commonly measures a random set of cells per individual; for example, a recent salamander hematology study measured 50 randomly extracted erythrocytes per individual ([Liu et al. 2023](https://doi.org/10.7717/peerj.15446)). The closest genome/cell/nucleus study analyzed 50 Feulgen-stained nuclei per individual, sampled 1–4 individuals per species, used species medians because of intraspecific variation and outliers, and repeated phylogenetic analyses with means ([Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)).

### 4.4 Genome-size image-density minimum gate

Feulgen image analysis densitometry is a quantitative genome-size method because the Feulgen reaction is DNA-specific and the bound stain is proportional to DNA; image pixels are converted from transmittance to absorbance and summed as integrated optical density (IOD) ([Hardie, Gregory, and Hebert 2002](https://doi.org/10.1177/002215540205000601)). Their protocol requires controlled Köhler illumination, white balance/background measurement, nonsaturated intensity, consistent exposure, undamaged nonoverlapping nuclei, and reference-standard calibration.

The exact salamander precedent used air-dried blood smears, overnight fixation, acid hydrolysis, freshly prepared Schiff reagent, dark storage, IOD from 50 nuclei per individual, and reference standards ([Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)). An earlier *Desmognathus* study used Feulgen cytophotometry and *Xenopus laevis* erythrocyte nuclei as an internal DNA reference ([Sessions and Kezer 1986](https://doi.org/10.1007/BF00494802)). These studies establish the minimum comparable assay.

Accordingly:

- If the slides were not stained with a validated DNA-stoichiometric stain, or no reference standard was co-processed under comparable conditions, nuclear darkness must **not** be reported as genome size. It can be retained as an image-intensity phenotype or QC covariate.
- IOD must be computed as a sum of per-pixel optical densities using local/slide background intensity, not as raw darkness, mean grayscale, area times an untransformed intensity, or display-adjusted pixel values. The equation and channel must be printed in the notebook ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)).
- Saturated, clipped, unevenly illuminated, blurred, overlapping, folded, broken, or out-of-focus nuclei require an objective fail rule. Show their frequency by species and slide.
- Reference and unknown nuclei must be stained, imaged, and analyzed together or under a demonstrated batch calibration. Report the reference's assumed 2C DNA amount and propagate its uncertainty.
- Genome-size conversion must use the ratio of unknown/reference IOD on the correct C-value scale. Do not calibrate to the same morphology correlation that will later be tested.
- Use specimen-level replicate slides/batches when available and quantify within-slide, among-slide, among-batch, and among-specimen variance.

### 4.5 Gold-standard validation

- Measure a representative subset with flow cytometry using a DNA-selective fluorochrome and internal reference. Simultaneous measurement of unknown and reference nuclei is the core precision principle of cytometric genome-size estimation ([Doležel, Greilhuber, and Suda 2007](https://doi.org/10.1038/nprot.2007.310)).
- Compare image-derived and cytometric values using a calibration plot **and** a Bland-Altman agreement plot. Correlation alone cannot establish interchangeability between measurement methods ([Bland and Altman 1986](https://doi.org/10.1016/S0140-6736%2886%2990837-8)). Predefine acceptable bias/limits from biological needs rather than after seeing the result.
- Include known salamander reference specimens spanning the observed genome-size range, not only one calibration point, to test linearity and saturation.
- Repeat Feulgen hydrolysis/staining for a subset to quantify batch repeatability. Mueller et al. reported that historical genome-size values could shift systematically with assumed standard values while relative species ratios remained stable, demonstrating why reference provenance matters ([Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)).
- Feed specimen/species measurement error into the comparative/path model using an error-aware phylogenetic method ([Ives et al. 2007](https://doi.org/10.1080/10635150701313830)).

### 4.6 Required figures and tables

1. **Segmentation QC atlas:** raw image, truth, prediction, overlay, and accept/reject reason for representative normal, overlapping, broken, blurred, border, small, large, and irregular cells/nuclei.
2. **Held-out performance dashboard:** object precision/recall/F1, merge/split rates, cell and nucleus IoU/Dice, area bias, and count bias by species/slide/quality.
3. **Selection-sensitivity plot:** species ranks and estimates under all-valid, balanced random, median, trimmed mean, upper quantiles, top-50, and top-100 summaries.
4. **Nested distribution plot:** cells within images within specimens within species, with biological replicate counts visibly distinct from cell counts.
5. **Feulgen QC plot:** background, exposure/saturation, nuclear area, mean OD, and IOD distributions by slide/batch; standards shown alongside unknowns.
6. **Calibration and agreement:** known/reference genome size versus IOD with residuals, plus Bland-Altman image-versus-cytometry agreement.
7. **Phylogenetic morphology plots:** genome size versus nucleus and cell size, raw species means with measurement intervals, PGLS fit/interval, and leave-one-species-out sensitivity. Mueller et al. found that ordinary and phylogenetically controlled conclusions can differ, which makes both raw and phylogenetic views necessary ([Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)).

---

## 5. Minimum publication package before path-analysis promotion

The shareable notebook should execute these upstream products in order and stop on failed gates:

1. **Canonical resource/specimen registry**
   - separate genome-resource accessions, dnaPipeTE technical retry aliases, microscopy specimens, slides, and phylogeny tips;
   - final species-panel membership and exclusion reason;
   - immutable input hashes and software versions.
2. **TE composition and diversity**
   - mass-accounted community matrix;
   - depth and rare-category sensitivities;
   - diversity indices and PGLS;
   - compositional ordination and, only if justified, phylogenetic multivariate inference.
3. **LTR deletion proxy**
   - structural candidates and manual checks;
   - family-level T:I estimates;
   - mapping/coverage sensitivities;
   - hierarchical/species estimates with intervals;
   - assembly structural validation where possible.
4. **Microscopy**
   - held-out segmentation validation;
   - object-level QC/exclusion ledger;
   - balanced cell-summary sensitivities;
   - Feulgen/reference validation or explicit withdrawal of the genome-size interpretation.
5. **Phylogeny**
   - exact crosswalk and pruned time-tree set;
   - topology/branch-length diagnostics;
   - tree-uncertainty representation.
6. **Only then: path analysis**
   - a priori DAGs and basis sets;
   - global fit and CICc comparison;
   - measurement-error, tree, influence, and simulation sensitivities;
   - claim-to-figure/table map.

## Bottom-line benchmark

A publication-quality version of this project can make unusually strong claims because it combines curated TE libraries, TE landscapes, LTR deletion proxies, time-calibrated phylogeny, and specimen-derived cell biology. The integration is also where the main validity risks lie. The defensible standard is not “all pipelines completed”; it is that each upstream variable has a named biological estimand, a canonical independent sampling unit, quantified measurement uncertainty, a method-specific negative/sensitivity control, and figures that expose rather than average away failure modes.

The most consequential hard boundary is genome-size imaging: **uncalibrated nuclear darkness is not a C-value**. The most consequential statistical boundary is replication: **reads, TE families, cells, and rerun folders are not additional species**. The most consequential phylogenetic boundary is model fit: **a path diagram without a testable basis set, global fit, tree uncertainty, and power/influence analysis is exploratory artwork, not a confirmed causal model**.

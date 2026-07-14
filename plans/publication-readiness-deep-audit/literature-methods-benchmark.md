# Literature and Methods Benchmark for Publication Readiness

Date: 2026-07-09

## Scope and evidence boundary

This is a read-only benchmark of the current Desmognathus_TE analysis against primary methods papers and closely relevant salamander studies. It does not remove, replace, or rerun an analysis. It asks which present results are already defensible, which claims need narrower wording, and which validation analyses should be completed before a public manuscript and shareable notebook are frozen.

Three evidence labels are used throughout:

- **Published evidence** means a claim supported directly by a linked primary paper or primary methods source.
- **Repository observation** means something directly visible in the current repository, current derived tables, or the linked CellProfiler workspace.
- **Audit interpretation** means an inference from comparing those two bodies of evidence. It is a recommendation, not a result already demonstrated by the data.

An item described as absent means it was not found in the repository surfaces inspected for this audit. That is a documentation gap, not proof that the laboratory step was never performed or recorded elsewhere.

## Executive verdict

The project has an unusually strong and potentially publishable biological foundation: custom repeat libraries, comparative repeatome summaries, direct paired-LTR sequence work, a time-calibrated phylogeny, curated erythrocyte masks, linked cell-nucleus measurements, explicit source hashes, and small prespecified path-model families. The repository also already states several important limitations honestly.

The current evidence does **not yet support treating every layer as one equally mature causal chain**. The most defensible publication structure is:

1. a primary comparative TE/genome-size analysis, after genome-size calibration and TE cross-species comparability are locked down;
2. a separately reported cell-size and nucleus-size analysis with the sampling estimand made explicit;
3. a genome-nucleus-cell path analysis only as an association/sensitivity analysis unless genome size comes from an assay independent of the nucleus images used to measure nucleus size;
4. paired-LTR divergence as an insertion-history proxy with conditional rate-calibrated ages;
5. terminal:internal LTR depth as an ectopic-deletion proxy, not a direct recombination rate, until it is validated against direct solo-LTR calls and mapping simulations.

Five issues should be treated as publication gates rather than polishing tasks:

- The “largest 50” design is selected predominantly on cell area and partly on nucleus darkness. It estimates the upper tail of accepted detected cells, not the average erythrocyte for a species.
- Genome size and nucleus size are derived from the same linked nuclear images. Their association is not independent evidence for a causal genome-to-nucleus edge.
- Fifteen of 21 current genome species have one image and one specimen; object-level resampling cannot recover missing specimen-level replication.
- The current phylogenetic path scaffold reports rankings and coefficient standard errors, but not the full basis-set tests, bootstrap confidence intervals, regression diagnostics, influence analysis, tree uncertainty, or model-recovery simulations needed for a small comparative panel.
- The TE, LTR-age, and ectopic layers need explicit estimator validation and cross-species technical comparability tables before mechanistic language is used.

These are tractable issues. None requires deleting the existing analyses. Most can be handled by declaring a primary estimand, adding alternative specifications, strengthening metadata and validation, and preserving the current outputs as sensitivity branches.

## Current repository snapshot

The following are current repository facts, not literature-derived claims.

| Surface | Current state | Publication implication |
|---|---:|---|
| TE species | 34 | Broadest comparative layer |
| Linked genome estimates | 21 | Smaller than the TE layer |
| TE + genome overlap | 18 | Effective sample for the core path family |
| TE + genome + ectopic overlap | 16 | Effective sample for the mechanism family |
| TE + genome + morphology overlap | 18 | Same linked-IOD dependence as genome layer |
| Stable genome estimates in the 18-species primary panel | 0 | Four are minor-caution; 14 are sensitivity-limited |
| Current genome support among all 21 species | 6 medium, 13 limited, 2 low | Support status must travel into every result |
| Species with one image and one specimen | 15 of 21 | Cell count is not biological replication |
| Species with two images and two specimens | 6 of 21 | Still shallow replication |
| Nuclei contributing to current primary genome estimates | 12 to 50 per species | Several “top-50” sets become smaller after OD QC |
| Current path runner | phylopath with default phylo_path call | Pagel-lambda PGLS is the package default |

Relevant repository anchors include:

- path_analysis/README.md:84-97, which already limits causal interpretation of the morphology chain;
- path_analysis/CELLPROFILER_PROVENANCE_AUDIT.md:81-96, which describes weighted medians, OD-QC filtering, and scaling to D. fuscus = 16.36 pg;
- path_analysis/data/derived/analysis_panel_summary.csv, which records 0 stable and 14 sensitivity-limited genome estimates in the 18-species panel;
- path_analysis/scripts/path_model_scaffold.R:494-551, which writes model ranking, best edges, averaged edges, and plots;
- ../cellprofiler_test/docs/MIXED_CELLPOSE_YOLO_WORKFLOW.md:127-132 and ../cellprofiler_test/scripts/build_top50_linked_pair_worklist.py:145-169,341-350, which define the current top-50 eligibility and ranking procedure;
- scripts/processing/ec.py:72-170,200-216,341-377, which define the terminal:internal depth statistic;
- results/reports/LTR_AGE_AUDIT.md:14-32, which records the paired-LTR substrate and successful sequence estimates.

## 1. Biological benchmark: genome size, nucleus size, and cell size

### Published evidence

[Hally et al. (1986)](https://doi.org/10.1007/BF00494802) measured DNA-Feulgen cytophotometry in four Desmognathus species using Xenopus laevis erythrocyte nuclei as an internal DNA standard. They found statistically detectable but small interspecific differences and did not find those differences sufficient to support the proposed life-history relationships. This is a particularly close benchmark because it is Desmognathus-specific and uses erythrocyte nuclei.

[Mueller et al. (2008)](https://doi.org/10.1016/j.zool.2007.07.010) tested genome, nucleus, and erythrocyte size in Batrachoseps. Raw species analyses showed positive relationships among genome, nucleus, and cell size, but after phylogenetic independent contrasts the genome-nucleus relationship remained while genome-cell relationships did not. The authors used species medians because within-species variation was substantial. Their result is a direct warning that a plausible raw three-variable chain may become a narrower relationship after phylogenetic correction.

[Sessions and Larson (1987)](https://doi.org/10.1111/j.1558-5646.1987.tb02463.x) established the broader plethodontid developmental-genome-size framework. [Itgen et al. (2022)](https://doi.org/10.1111/evo.14519) further showed that genome-size effects on morphology can be organ-specific rather than a single uniform organism-wide scaling relationship.

[The 2025 D. fuscus assembly paper](https://doi.org/10.1093/g3journal/jkaf157) notes that historical C-values assigned to Desmognathus and D. fuscus vary with method and taxonomic interpretation. That makes specimen identity, locality, ploidy convention, and calibration source material parts of the biological evidence, not clerical metadata.

### Accession-specific *D. fuscus* resource benchmark

The current repository does not use the 2025 `aDesFus1-2.1 / GCA_050004315.1` resource. Its *D. fuscus* TE, LTR-age, and ectopic rows trace to the independently deposited short-read resource `SRX20497025 / GCA_032353935.1`, publicly labeled *D. planiceps* but identified by Alex Pyron as actual *D. fuscus*. The separate older resource publicly labeled *D. fuscus*, `SRX19953421 / GCA_030265095.1`, is excluded as untrusted.

The published assembly is 16.118 Gb, whereas the current short-read assembly spans 0.653 Gb, only 4.05% as much assembled sequence. The current dnaPipeTE table reports 64.696% LTR and 16.360% LINE among classified order composition; Myers et al. report approximately 36% LTR and 15% LINE of the whole 2025 assembly, with approximately 75% total TE content. These values diagnose potentially important method and assembly sensitivity but are not interchangeable because the current table is closed relative composition while the paper reports genome occupancy. The fully accession-labeled calculation is preserved in [`fuscus_resource_benchmark.md`](fuscus_resource_benchmark.md) and [`fuscus_resource_benchmark.csv`](fuscus_resource_benchmark.csv).

The paper states that EarlGrey annotation outputs are deposited, but the live [Zenodo record](https://zenodo.org/records/15255946), checked 2026-07-09, exposes only the compressed assembly FASTA and AGP. The promised filtered-repeat GFF is absent. The most efficient next comparison therefore requires the EarlGrey summary directory, filtered-repeat GFF, repeat library, divergence table, raw LTR_FINDER GFF, commands, and checksums from the authors; it does not require downloading the full assembly first.

### Repository observation

The candidate-model document treats the classical genome-to-nucleus-to-cell sequence as a prespecified hypothesis and includes plausible alternatives. The repository also explicitly warns that genome size is calibrated from nucleus IOD and is not an independent assay for formal genome-to-nucleus causality.

The current full morphology DAG nevertheless uses genome size and nucleus area from the same linked objects. Current species values are then reduced to point estimates before phylopath is run.

### Audit interpretation

The biological hypothesis is well motivated, but the published record does not justify treating the complete chain as an expected fact. The analysis should be confirmatory about a small set of alternatives, and the manuscript should be prepared for a result in which genome size predicts nucleus size more clearly than cell size.

The following claim hierarchy is defensible:

- **Descriptive:** species differ in measured erythrocyte cell area, nucleus area, and calibrated nuclear IOD.
- **Associational:** phylogenetically corrected species values covary, with effect sizes and uncertainty reported.
- **Mechanistic:** TE state, DNA-loss proxies, genome size, nucleus size, and cell size form the proposed causal sequence.

The first level can be supported by the current measurements after sampling/QC validation. The second requires uncertainty-aware comparative models. The third requires stronger identification assumptions and, for the genome-nucleus edge, an independent genome-size assay.

### Publication gate

Before final interpretation:

1. State the genome-nucleus-cell DAGs as competing hypotheses, not as an assumed biological ordering.
2. Preserve the null, direct, and mediated alternatives.
3. Report raw and phylogenetically corrected relationships side by side.
4. If independent genome-size measurements cannot be obtained, label the genome-nucleus-cell family “linked-IOD association/sensitivity analysis” and avoid causal verbs such as drives, determines, or mediates.
5. State whether genome values are 1C, 2C, or a 1C estimate calibrated from 2C erythrocyte nuclei; use one convention consistently.

## 2. Feulgen/image cytometry and genome-size calibration

### Published evidence

[Hardie, Gregory, and Hebert (2002)](https://doi.org/10.1177/002215540205000601) is the central primary methods benchmark for Feulgen image analysis densitometry. Its relevant requirements include:

- Feulgen staining is DNA-specific and can be stoichiometric only when fixation, hydrolysis, staining, imaging, and analysis conditions are controlled.
- Camera response must be tested for linearity and saturation; exposure and illumination must remain in the valid range.
- Integrated optical density is the sum of per-pixel optical density across the entire nucleus, so segmentation, background choice, clipping, and image dynamic range all affect the measurement.
- A hydrolysis curve is needed to establish the staining optimum for the material and protocol.
- Standards and unknowns should be stained and imaged together, preferably with multiple standards spanning the unknown range.
- Field density, slide age, illumination, background, white balance/channel choice, and microscope settings can materially shift IOD.
- A sample of 50 **random** erythrocyte nuclei can estimate a slide-level mean efficiently under a controlled protocol; this is not support for selecting the 50 darkest or largest objects.

[Hally et al. (1986)](https://doi.org/10.1007/BF00494802) used Xenopus erythrocyte nuclei as an internal reference standard of 6.35 pg DNA per cell in their Desmognathus work. The core design feature is that standards experience the same staining and imaging conditions as unknowns.

### Repository observation

The current bridge:

- scales verified linked nucleus IOD to D. fuscus = 16.36 pg;
- selects OD-QC-pass linked nuclei for the primary genome estimate;
- retains all-selected and high-OD-QC sensitivity states;
- reports bootstrap intervals and support tiers;
- has only one or two images/specimens per species in the current genome table;
- explicitly acknowledges that the resulting genome estimate is not independent of nucleus imaging.

The current repository contains extensive computational provenance for selected masks and image traces. In the inspected project surfaces, a complete publication-grade wet-lab and imaging protocol was not found that ties each slide to fixation, stain batch, hydrolysis duration/temperature, staining duration, reference nuclei, microscope/camera settings, acquisition channel, slide age, field density, and calibrator provenance.

The exact source and specimen/taxonomic match for the constant D. fuscus = 16.36 pg is not documented next to the current calibration rule.

### Audit interpretation

The computational trace is stronger than the measurement-model trace. A narrow object-bootstrap interval around a species median quantifies variation among selected nuclei under the current image, but it does not include uncertainty from stain batch, calibrator value, camera linearity, slide, image, specimen, locality, taxonomic assignment, or the selection mechanism.

Scaling all species to one D. fuscus point reference also induces shared calibration uncertainty. Treating that scale as exact makes every downstream genome-size standard error too small. If the D. fuscus reference specimen is also the calibrator species in the comparative dataset, its fixed value is not an independently estimated observation.

### Required measurement record

Create one row per slide or acquisition batch with at least:

| Domain | Required fields |
|---|---|
| Biological sample | accepted species name, original field ID, specimen/voucher ID, locality, date, sex/stage if known, individual ID |
| Material | blood preparation, anticoagulant if any, fixation method and duration, storage/preservation |
| Feulgen reaction | reagent recipe, hydrolysis acid/concentration, temperature, duration, staining time, wash procedure, batch ID |
| Standard | standard species/tissue, specimen ID, accepted 1C or 2C value, source DOI/database record, standard and unknown co-location |
| Microscope | instrument, objective, numerical aperture, illumination/Köhler settings, filters/channel |
| Camera | model, bit depth, exposure, gain, gamma, white balance, saturation threshold, linearity-test result |
| Image | file ID, pixel size calibration, field density, background definition, acquisition date, slide age |
| Analysis | segmentation model/version/hash, nuclear mask rule, exclusion rules, IOD equation, channel, background correction |

### Publication gate

At minimum:

1. Document the full staining and imaging protocol and link it to each image.
2. Provide the exact provenance and 1C/2C interpretation of D. fuscus = 16.36 pg.
3. Demonstrate camera/channel linearity, absence of saturation, background stability, and batch effects.
4. Show whether standard and unknown nuclei were stained on the same slide/run. If they were not, state the limitation and quantify run-to-run uncertainty.
5. Propagate reference-value and reference-IOD uncertainty into every species estimate.
6. Validate a subset against an independent method or external C-value source matched by taxon and, ideally, population. Flow cytometry or an independently calibrated Feulgen run would break the current genome-nucleus measurement dependence.
7. Run leave-one-image and leave-one-specimen analyses where replication permits.
8. Show genome results under primary QC, all quality-eligible objects, and a calibration/batch sensitivity grid.

## 3. The largest-50 design, selection bias, and hierarchical replication

### Published evidence

[Mueller et al. (2008)](https://doi.org/10.1016/j.zool.2007.07.010) used species medians to reduce sensitivity to outliers because within-species variation was substantial; the paper did not define the largest cells as the species trait.

[Hardie et al. (2002)](https://doi.org/10.1177/002215540205000601) discusses efficient measurement of random erythrocyte nuclei after quality control, not outcome-selected nuclei.

[Ives, Midford, and Garland (2007)](https://doi.org/10.1080/10635150701313830) and [Felsenstein (2008)](https://doi.org/10.1086/587525) show that phylogenetic comparative analyses can be biased when within-species variation and measurement error are ignored. Species values are estimates, not error-free true means.

[Hurlbert (1984)](https://doi.org/10.2307/1942661) provides the classic pseudoreplication framework: many subsamples from one experimental or biological unit do not create many independent units.

### Repository observation

The linked worklist:

- requires one-to-one valid cell-nucleus pairs and several defensible mask/file checks;
- drops the bottom within-species quartile for nucleus darkness and clarity by default;
- then ranks candidates with a score that is 80% cell-size signal, 10% nucleus-darkness signal, 7% clarity, and 3% other terms;
- freezes up to 50 manually reviewed pairs per species;
- uses weighted species medians downstream;
- often has only one image and one specimen for a species.

This design intentionally targets large cells and dark nuclei. It is not simply “take 50 clean cells.”

### Audit interpretation

The current estimator is approximately the median of the largest accepted, well-focused, relatively dark-nucleus detected erythrocytes available to the ranking workflow. It is a legitimate biological quantity if that upper-tail phenotype is the declared target, but it is not an unbiased estimator of mean or median erythrocyte size.

Selection on cell area can exaggerate or otherwise alter between-species cell-size differences when species have different distribution shapes, different numbers of eligible cells, different segmentation failure rates, or different image areas. Selection on nucleus darkness is especially consequential because nuclear IOD is then reused for genome-size calibration. Conditioning on quality variables that also correlate with the predictor and outcome can change their observed association. The direction and magnitude cannot be assumed; it must be measured with alternative specifications.

Object-level bootstrap intervals are also conditional on the observed specimen/image clusters. With one specimen and one image, resampling 50 objects can estimate object-distribution uncertainty but cannot estimate between-individual, between-slide, or between-image variation. A hierarchical bootstrap becomes equivalent to object resampling when there is only one higher-level cluster.

### Required estimand comparison

Keep the current top-50 analysis, but compare it with:

1. **All quality-eligible linked pairs:** exclusion rules based only on prespecified segmentation/image quality, not cell size or IOD rank.
2. **Repeated random 50:** many seeded samples of 50 from the quality-eligible pool per species.
3. **Stratified random 50:** balanced across specimen and image, when possible.
4. **Top 50 by cell area only:** isolates the effect of tail selection.
5. **Top 50 under current score:** preserves the existing design.
6. **Median, trimmed mean, and distributional quantiles:** shows whether conclusions depend on the location statistic.
7. **Specimen-first species estimate:** compute image summaries, then specimen summaries, then species summaries so an image with many cells cannot dominate.

For every specification, report the exact estimand in the figure and methods. If top 50 remains primary, call it “median area of the 50 highest-ranked eligible linked erythrocytes” rather than “species cell size.”

### Required hierarchy and segmentation checks

- Report numbers of specimens, slides, images, eligible cells, accepted cells, and selected cells separately.
- Use specimen as the highest biological sampling unit and image/field as a nested technical unit.
- Propagate the hierarchy by bootstrap or a multilevel measurement model.
- Do not present narrow object-bootstrap intervals as total biological uncertainty for one-specimen species.
- Validate segmentation on held-out images stratified by species, staining/acquisition batch, and image quality.
- Report cell and nucleus Dice/IoU or boundary error, area bias, detection precision/recall, and linked-pair error.
- Include blinded manual repeat measurements for a subset and inter-reviewer or repeat-review agreement.
- Test whether segmentation error depends on cell area, nucleus area, darkness, focus, species, or batch.
- Preserve the rejected-object reason table and show species-specific rejection rates.

### Publication gate

The largest-50 analysis should not be removed. It should be demoted to sensitivity unless:

- the upper-tail erythrocyte phenotype is explicitly the biological estimand;
- the same count is available from a sufficiently large, similarly sampled eligible pool in every species;
- selection does not use nuclear darkness for any analysis in which IOD is a predictor;
- conclusions are stable across all-eligible and repeated-random-50 specifications.

## 4. Phylogenetic path analysis, d-separation, CICc, and small comparative samples

### Published evidence

[von Hardenberg and Gonzalez-Voyer (2013)](https://doi.org/10.1111/j.1558-5646.2012.01790.x) introduced phylogenetic confirmatory path analysis as a way to compare prespecified causal hypotheses while fitting phylogenetically corrected local regressions. Its evidential force comes from testing the conditional independencies implied by each DAG, not from drawing arrows after observing correlations.

[Shipley (2013)](https://doi.org/10.1890/12-0976.1) explains the d-separation basis set, Fisher's C statistic, and information-criterion comparison for path models. A non-rejected C test means the data do not reject the model's implied independencies; it does not prove the DAG or rule out unmeasured common causes.

[van der Bijl (2018)](https://doi.org/10.7717/peerj.4718) documents phylopath, including standardized coefficients, model ranking, averaging, and bootstrap confidence intervals. The package's official [introductory vignette](https://ax3man.github.io/phylopath/articles/intro_to_phylopath.html) notes that coefficient confidence intervals require an explicit boot argument; they are off by default.

[Boettiger, Coop, and Ralph (2012)](https://doi.org/10.1111/j.1558-5646.2011.01574.x) show that information-criterion model selection in phylogenetic comparative analyses can have high error rates depending on the exact tree, sample size, effect sizes, and model structure. They recommend simulation on the study's own design rather than a universal sample-size rule.

### Repository observation

Positive design choices already present:

- candidate families are intentionally small;
- one observed proxy is used per conceptual block in first-pass models;
- null, direct, additive, and mediated alternatives are included;
- morphology is explicitly labeled as staged/sensitivity work;
- data are pruned to a finite, common family-specific analysis table before the model set is fit.

Current limitations:

- effective sample sizes are about 18 for TE/genome, 16 for ectopic, 15 for LTR history, and 8-9 in strict-body panels;
- some local equations contain three predictors at those sample sizes;
- phylo_path is called without an explicit evolutionary model or bootstrap setting;
- output includes ranking, warnings, best/averaged coefficients, and standard errors, but not a complete exported basis-set table;
- the current runner does not export local-regression diagnostics, residual plots, leverage/influence, lambda estimates/boundary behavior, coefficient confidence intervals, or leave-one-tip results;
- a single dated tree file is used by the scaffold;
- species-level measurement uncertainty is not propagated into model ranking or coefficients.

### Audit interpretation

There is no defensible universal claim that n = 18 makes phylopath invalid. The correct conclusion is narrower: model discrimination and individual path coefficients may be fragile, and the exact fragility must be measured using the project's tree, models, trait covariance, missingness, and measurement errors.

CICc can rank weakly distinguishable models even when all coefficients are poorly identified. A best-model label without model weight, evidence ratio, bootstrap interval, influence analysis, and recovery simulation will look more decisive than the data warrant.

The strict-body panels with 8-9 tips are best treated as qualitative sensitivity analyses, especially for local regressions with several predictors.

### Required model-reporting bundle

For every model family, save:

1. exact species list and tree hash;
2. exact data table and transformations;
3. DAG edge list and biological rationale written before seeing its result;
4. basis-set conditional-independence claims;
5. every local PGLS equation, coefficient, standard error, confidence interval, lambda, sample size, and warning;
6. Fisher C, model-fit p value, parameter count, CICc, delta CICc, and model weight;
7. best-model and model-averaged coefficients with bootstrap confidence intervals;
8. residual diagnostics, influential species, collinearity checks, and convergence/boundary warnings;
9. leave-one-tip and leave-one-major-clade results;
10. alternative branch-length/time-tree results;
11. a primary/sensitivity status field.

### Required small-n validation

Run a model-recovery experiment on the exact observed tree:

1. Simulate data under each candidate DAG using plausible weak, medium, and observed effect sizes.
2. Add species-level measurement error drawn from morphology/genome bootstrap distributions and shared calibration uncertainty.
3. Reproduce the observed missingness and panel sizes.
4. Refit the complete candidate set.
5. Report the probability that the generating model is ranked first, false selection rates, coefficient bias, interval coverage, and the frequency of singular or boundary fits.
6. Repeat across plausible trees or posterior tree draws if those are available.

If the best model cannot be recovered reliably under observed-size effects, the appropriate manuscript conclusion is that the candidate models are not distinguishable with this panel—not that the top-ranked model is demonstrated.

### Additional safeguards

- Use the same species/data rows for all candidate DAGs within a ranked family.
- Keep exploratory model development separate from the frozen confirmatory set.
- Do not average across models that encode biologically contradictory directions unless the averaged coefficient has a coherent interpretation.
- Show raw scatterplots and phylogenetic residual relationships; a path diagram alone is insufficient.
- Propagate genome and morphology uncertainty through repeated-draw path analyses, not only through a standard error column that phylopath does not consume.
- Treat tree uncertainty and taxonomic crosswalk uncertainty as analytical uncertainty.
- Avoid language that a large d-separation p value “supports causality.” It only indicates that the implied independencies were not rejected.

## 5. TE annotation, landscapes, diversity, and compositional predictors

### Published evidence

[Goubert et al. (2015)](https://doi.org/10.1093/gbe/evv050) designed dnaPipeTE for repeat assembly and quantification from low-coverage raw reads. The source method requires samples below 1x coverage to avoid assembling nonrepetitive genome content, often uses 0.1-0.25x, and uses independent samples for assembly and quantification. Its benchmarks recover broad recent TE composition well but can underestimate total repeat content and are most informative for relatively low-divergence copies.

[Flynn et al. (2020)](https://doi.org/10.1073/pnas.1921046117) shows that de novo family discovery and accurate consensus compilation are critical and difficult, even with RepeatModeler2. Library fragmentation, family merging, contamination, and inconsistent classification can therefore change apparent family diversity.

[Ou and Jiang (2018)](https://doi.org/10.1104/pp.17.01310) validates structure-aware LTR_retriever calls and highlights the value of structural evidence for intact LTR retrotransposons. [Ou et al. (2018)](https://doi.org/10.1093/nar/gky730) introduced the LTR Assembly Index, which is useful for documenting whether an assembly is sufficiently continuous for LTR-based inference.

[Aitchison (1982)](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x) establishes that proportions constrained to sum to a constant require compositional treatment. Raw percentages cannot be interpreted as independent predictors.

[Sun et al. (2012)](https://doi.org/10.1093/gbe/evr139) provides the close salamander-specific benchmark that LTR retrotransposons contribute strongly to plethodontid genome gigantism, while also illustrating that repeat abundance and repeat history are distinct questions.

### Repository observation

The repository has:

- order- and superfamily-level dnaPipeTE summaries;
- custom repeat-library work;
- RepeatMasker-derived divergence/deletion landscapes;
- Pielou evenness and other diversity summaries;
- a targeted log LTR:LINE balance predictor;
- centered-log-ratio machinery for broader order compositions;
- source hashes and a path-analysis TE feature table;
- an explicit warning that the LTR-history axis is Gypsy-heavy rather than a global TE-age axis.

The path model uses compact TE predictors, which is appropriate for n = 18. The current CLR code replaces zeros using half the minimum observed positive value, a practical but assumption-sensitive choice.

The repository documents local derived tables more completely than the full upstream per-species raw-read/assembly execution context. A publication reader will need a single manifest proving that every species was processed comparably.

### Audit interpretation

The biggest TE threat is technical confounding across species, not the choice between one more diversity index. Read technology, read length, library preparation, raw depth, subsampling coverage, genome-size value used to calculate coverage, assembly span/continuity, library version, and unclassified fraction can all create species differences that resemble biology.

dnaPipeTE read-based abundance and RepeatMasker assembly annotation are different estimands. They can validate each other, but their percentages should not be silently mixed in one comparative column.

TE “landscape age” based on divergence to a family consensus is a relative copy-divergence distribution, not a direct insertion-time distribution. Old/divergent elements are harder to assemble, classify, and map, and the consensus itself is an estimated ancestral proxy.

Pielou evenness is interpretable only relative to a stable family/order definition and denominator. If custom libraries split one lineage's families more finely, apparent diversity can rise without any biological increase.

### Required per-species TE manifest

For every species, record:

- BioProject/BioSample/SRA or local source accession;
- voucher/individual and taxonomic crosswalk;
- sequencing platform, layout, read length, library type, and raw read count;
- trimming, contaminant removal, and mitochondrial filtering;
- genome-size value used to calculate target coverage, with source and uncertainty;
- dnaPipeTE target coverage, sampled read count, random seed, iterations, and independent assembly/quantification samples;
- dnaPipeTE, Trinity, RepeatMasker, Repbase/RepeatMasker library, RepeatModeler, TEsorter, and all custom-script versions;
- custom-library construction, classification, curation, redundancy reduction, contamination filtering, and checksums;
- classified, unclassified, simple-repeat, and total repeat fractions;
- if assembly-based: accession/version, assembly span, contig/scaffold N50, BUSCO, read coverage, and LAI where applicable;
- every denominator used in order/superfamily percentages;
- batch or source-study identifiers.

### Required TE validation

1. Rerun or verify all species at the same target coverage and with recorded seeds.
2. Demonstrate saturation/stability across at least two nearby coverage levels for representative small, medium, and large genomes.
3. Use the 16.1-Gb D. fuscus reference assembly as an anchor: compare dnaPipeTE abundance, assembly RepeatMasker abundance, order composition, unclassified fraction, and LTR history, while retaining its documented contiguity and repeat-resolution limitations.
4. Repeat the anchor comparison for any other sufficiently continuous assemblies.
5. Quantify batch/read-length/platform associations with each TE predictor.
6. Confirm that the genome-size values used to set dnaPipeTE coverage are not circularly derived from the same response being modeled, or show that modest coverage-value perturbations do not change the TE predictors.
7. Preserve unclassified repeats as an explicit component or exclusion category.
8. Report compositional sensitivity to zero replacement, denominator choice, CLR versus an interpretable ILR/pivot balance, and exclusion of rare orders.
9. Report diversity results under order, curated superfamily, and effective-number representations.
10. Treat consensus divergence as a relative landscape statistic and validate any age-language separately with paired LTRs.

### Statistical use

The targeted LTR:LINE log ratio is a defensible low-dimensional balance if it was biologically prespecified. The broader CLR is also directionally correct. Before publication:

- state the numerator and denominator in every balance;
- test at least two zero-replacement strategies;
- do not include several algebraically dependent composition coordinates in the same small model;
- show that the chosen balance is not driven by the unclassified component;
- report raw proportions for biological readability while fitting log-ratio coordinates;
- keep order-level primary models and superfamily-level analyses supplementary unless the effective sample increases.

## 6. Paired-LTR divergence and insertion-age calibration

### Published evidence

The paired-LTR clock rests on the expectation that the 5-prime and 3-prime LTRs are identical at insertion and then diverge. [SanMiguel et al. (1998)](https://doi.org/10.1038/1695) is an early primary demonstration of using paired LTR divergence to reconstruct retrotransposon history. The usual conversion is:

T = K / (2r)

where K is corrected divergence between the two LTRs and r is the substitution rate per site per year.

[Kimura (1980)](https://doi.org/10.1007/BF01731581) provides the two-parameter correction used for K2P distances. [Kijima and Innan (2010)](https://doi.org/10.1093/molbev/msp295) shows that gene conversion between LTRs can homogenize them and bias inferred insertion times toward younger values. Rate heterogeneity, alignment uncertainty, structural miscalls, and survivor bias add further uncertainty.

### Repository observation

The current paired-LTR branch is substantially stronger than a divergence-to-consensus heuristic:

- 1,371 valid paired elements in the master table;
- 1,291 successful paired-sequence divergence estimates;
- 238 successful high-confidence estimates;
- 1,363 of 1,371 candidate paired elements classified as Gypsy;
- median 363 comparable ungapped sites;
- K2P divergence stored directly;
- absolute ages shown under a recommended rate window from 1.0e-9 to 1.53e-9 substitutions/site/year, with 1.227e-9 as a central value;
- the rate anchors come from salamander/frog nuclear genes, not directly estimated LTR lineage rates.

### Audit interpretation

The direct paired-LTR divergence is publication-grade in concept. The absolute ages remain conditional calibration products. A single central age should not be reported without the rate window, sequence uncertainty, structural-confidence tier, and gene-conversion caveat.

Because the candidates are overwhelmingly Gypsy, the species summary is a retained intact Gypsy-history statistic. It cannot represent the insertion history of all LTR superfamilies or all TEs.

The intact elements that survive assembly and pass structural filters are not a random sample of all historical insertions. Older, fragmented, recombined, or assembly-collapsed elements are preferentially absent.

### Publication gate

1. Make K2P divergence the primary measured variable.
2. Report rate-calibrated ages as a range or sensitivity interval, not as exact dates.
3. Show the complete formula, time units, rate sources, and whether the rate is per lineage.
4. Define alignment software/parameters, minimum comparable sites, gap/ambiguous-base rules, and K2P failure handling.
5. Report all-pair and high-confidence-pair summaries.
6. Stratify by superfamily where counts permit; explicitly label the main result Gypsy-heavy.
7. Bootstrap elements within species and propagate substitution-rate uncertainty.
8. Test whether results change after requiring Complete = yes, 5+ domains, longer LTRs, or higher comparable-site thresholds.
9. Screen for or discuss gene conversion and show a sensitivity excluding suspiciously homogenized/discordant pairs if a defensible diagnostic can be specified.
10. Do not compare absolute ages across species as if assembly completeness and candidate recovery were identical without an assembly-quality sensitivity.

## 7. Terminal:internal LTR depth as an ectopic-deletion proxy

### Published evidence

[Frahry et al. (2015)](https://doi.org/10.1007/s00239-014-9663-7) is the closest salamander-specific benchmark. Ectopic recombination between the two LTRs of an element removes the internal region and one LTR, leaving a solo LTR. The study compared solo-LTR levels in four salamanders and five smaller-genome vertebrates and also used terminal:internal evidence at the LTR-family level. Equivalent read sample sizes and outlier handling were part of the comparison.

The observable is therefore a deletion signature or relative abundance pattern, not a directly measured per-generation recombination rate.

### Repository observation

The current implementation:

- starts from candidate LTR elements at least 3,000 bp long;
- calculates mean terminal and internal depth within each candidate element;
- discards every zero-depth position before calculating both means;
- defines the per-element ratio as mean terminal depth divided by mean internal depth;
- merges TEsorter annotations and retains elements with exactly five or six annotated domains;
- does not require Complete = yes in the final retained table;
- averages usable per-element ratios within species for the path-analysis predictor;
- has 16 species in the current TE + genome + ectopic overlap.

### Audit interpretation

The statistic is biologically motivated, but it is not yet equivalent to the Frahry et al. family-level estimator and is not a direct recombination rate.

Dropping zero-depth positions conditions both numerator and denominator on observed mapping. In repetitive sequence, zeros may arise from true absence, low coverage, masking, divergence, or multimapping policy. Removing them can change the ratio in a species- and region-dependent direction. Averaging per-element ratios also weights an element with low or patchy coverage the same as a well-covered element and differs from a ratio of pooled terminal/internal depth.

Requiring five or six domains selects protein-coding, structurally recognizable elements, while not requiring Complete = yes allows mixed structural states. The resulting species mean is therefore an estimator of the selected candidate set, not all LTR deletion events.

### Required validation ladder

1. **Reproduce the literature-like estimator.** Build family-level terminal:internal summaries after equalizing read count/coverage across species and applying prespecified outlier rules.
2. **Keep zero-depth positions.** Compare full-length denominators, nonzero-only denominators, and explicitly mappability-masked denominators.
3. **Compare aggregators.** Per-element mean ratio, median ratio, coverage-weighted mean ratio, pooled terminal depth / pooled internal depth, and family-balanced mean.
4. **Stratify structural confidence.** Complete = yes; 5+ domains; exact 5/6 domains; high-confidence paired-LTR set.
5. **Audit mapping.** Mapper/version, multimapping policy, MAPQ, duplicate handling, read length, GC, local mappability, library depth, and coverage normalization.
6. **Simulate.** Map reads from known mixtures of intact elements and solo LTRs over the observed read-length/depth range to quantify estimator bias and recovery.
7. **Validate directly.** On the 16.1-Gb D. fuscus reference assembly, call family-matched solo LTRs and intact LTRs with a structure-aware pipeline, then compare direct solo:intact abundance with every depth estimator while reporting assembly limitations.
8. **Cross-validate.** Repeat on any other continuous assemblies and test rank concordance across species.
9. **Quantify uncertainty.** Bootstrap reads/families/elements as appropriate and propagate it into comparative models.
10. **Rename the predictor.** Use “terminal:internal LTR depth ratio” or “ectopic-deletion proxy,” not “ectopic recombination rate.”

### Publication gate

Until the depth statistic is concordant with direct solo:intact calls and robust to mapping/denominator choices, it should remain a mechanistic sensitivity analysis. A non-significant result would not establish that ectopic deletion is unimportant; a significant result would not by itself measure its rate.

## 8. Reproducibility and the shareable Jupyter notebook

### Published evidence

[Sandve et al. (2013)](https://doi.org/10.1371/journal.pcbi.1003285) recommends recording every operation, avoiding manual data manipulation, versioning scripts and intermediate results, preserving raw data, recording seeds, and making the complete analysis executable.

[Wilkinson et al. (2016)](https://doi.org/10.1038/sdata.2016.18) defines FAIR research objects as findable, accessible, interoperable, and reusable, with persistent identifiers and rich metadata.

[Nosek et al. (2015)](https://doi.org/10.1126/science.aab2374) frames transparent data, code, materials, and analysis plans as part of research-evaluation standards.

### Repository observation

The project already has many building blocks:

- source IDs and hashes;
- derived-table builders;
- analysis panel and readiness tables;
- explicit primary/sensitivity naming;
- current methods notes;
- tests and guard scripts;
- a conda environment and project wrapper;
- provenance links to the sibling CellProfiler workspace.

Current reproducibility friction remains:

- several canonical inputs are reached through a sibling project and absolute local paths;
- upstream HPC/raw TE generation is not represented by one uniform per-species manifest and one command;
- the final notebook does not yet exist as a clean, top-to-bottom publication artifact;
- the current analysis state is distributed across many scripts, reports, and generated files;
- a clean-room rebuild from accession/raw input to every manuscript figure has not yet been demonstrated.

### Audit interpretation

The notebook should be the transparent narrative and inspection surface, not the only place where analysis logic lives. Large hidden notebook cells are difficult to test and review. Stable functions should remain in versioned scripts/modules; the notebook should call them, show their inputs/outputs, and make every table and figure traceable.

### Required release bundle

The publication release should contain:

- a frozen Git commit or tagged release;
- an environment lock file and, ideally, a container recipe;
- accession/download manifest for data that cannot be redistributed;
- checksums for every frozen input;
- a data dictionary with units, 1C/2C convention, estimands, missing-value semantics, and primary/sensitivity status;
- taxonomy/voucher/locality crosswalk;
- image/slide/specimen hierarchy table;
- TE per-species technical manifest;
- exclusion/QC logs and reason counts;
- frozen candidate DAG definitions;
- one command that executes the notebook from a clean kernel;
- one command or workflow that rebuilds its derived inputs;
- fixed random seeds or recorded seed streams;
- software and database versions in notebook output;
- machine-readable model tables and figure-source tables;
- a clean execution log and test report;
- an archived release with a persistent DOI.

### Recommended notebook structure

1. **Scope and claim tiers**
   - primary questions;
   - confirmatory versus exploratory models;
   - explicit limitations.

2. **Environment and provenance**
   - commit, environment, seeds, input checksums;
   - taxonomy and source manifests.

3. **Sample flow**
   - species/specimen/slide/image/cell counts;
   - exclusions and missingness.

4. **Segmentation validation**
   - held-out accuracy;
   - area bias and failure modes;
   - example masks chosen by a deterministic rule.

5. **Morphology estimands**
   - all eligible, random 50, and top 50;
   - hierarchical uncertainty;
   - within- and between-species distributions.

6. **Genome-size calibration**
   - wet-lab and image-QC metadata;
   - reference equation and uncertainty;
   - batch/quality/independent-validation results.

7. **TE provenance and validation**
   - raw-read/assembly manifests;
   - coverage saturation;
   - composition and diversity.

8. **LTR history**
   - structural filters;
   - K2P divergence;
   - conditional rate windows.

9. **Ectopic-deletion proxy**
   - estimator definitions;
   - zero/mapping/aggregation sensitivities;
   - direct solo:intact validation.

10. **Phylogeny**
    - dated tree provenance;
    - taxonomy reconciliation;
    - tree uncertainty.

11. **Primary comparative analyses**
    - raw relationships;
    - phylogenetic models;
    - full d-separation and CICc output.

12. **Robustness**
    - alternative cell selection;
    - genome calibration states;
    - TE denominators;
    - leave-one-tip/clade;
    - tree draws;
    - measurement-error draws;
    - model-recovery simulation.

13. **Final figures and tables**
    - generated only from frozen machine-readable results.

14. **Session and data availability**
    - sessionInfo/package versions;
    - accessions and DOI;
    - known limitations.

The notebook must restart and run all cells in order without relying on prior kernel state, untracked local files, or hard-coded /home/jake paths.

## 9. Gap matrix

### Required before publication

| Gap | Why it is required | Minimum closure evidence |
|---|---|---|
| Declare the morphology estimand | Largest-50 is an upper-tail estimator, not a neutral species mean | Methods text plus all-eligible/random-50/top-50 comparison |
| Separate quality exclusion from size/IOD selection | Current ranking uses both cell area and nucleus darkness | Frozen exclusion rules and specification table |
| Resolve pseudoreplication language | 15/21 species have one specimen/image | Hierarchy table; intervals labeled conditional where replication is absent |
| Document Feulgen/imaging protocol | IOD depends on stain, hydrolysis, optics, camera, and batch | Slide-level protocol/metadata table |
| Document D. fuscus calibrator | 16.36 pg is currently a fixed shared scale | Source, specimen/taxon match, 1C/2C convention, uncertainty |
| Validate genome-size scale | Current values are linked-IOD derived and often sensitivity-limited | Independent subset comparison or explicit sensitivity-only status |
| Keep genome-nucleus causality bounded | Genome IOD and nucleus area share the same nuclear images | Independent genome assay or associational wording |
| Propagate measurement uncertainty | Species points are estimated with unequal support | Repeated-draw comparative analysis with shared calibration error |
| Export complete path evidence | Ranking alone is insufficient | Basis sets, local PGLS tables, C/p/CICc/weights, bootstrap CIs |
| Demonstrate small-n resolvability | n is 8-18 depending on family | Exact-tree model-recovery and coefficient-coverage simulation |
| Test influence and tree uncertainty | One tip/clade/tree may determine the result | Leave-one-tip/clade and alternative-tree report |
| Lock TE technical comparability | Species technical differences can mimic biology | Per-species TE manifest and batch/coverage checks |
| Validate TE estimators | Read-based and assembly-based summaries differ | Assembly-anchor and coverage-saturation comparisons |
| Treat compositions correctly | Raw TE percentages are constrained | Prespecified balances/CLR plus zero/denominator sensitivity |
| Bound LTR-age claims | Rates are borrowed and gene conversion biases ages | K2P primary; rate window and structural sensitivity |
| Validate ectopic proxy | Current estimator differs from direct solo-LTR evidence | Direct solo:intact anchor plus simulations and estimator grid |
| Freeze a clean executable release | Results must be inspectable and regenerable | Clean-kernel notebook run, locked environment, manifests, checksums |

### Strongly recommended

| Gap | Benefit | Suggested evidence |
|---|---|---|
| Increase independent specimen replication | Estimates biological rather than object-only variation | At least several independently prepared individuals per key species where feasible |
| Use multiple DNA standards | Detects nonlinear calibration and brackets unknowns | Same-run standard curve |
| Blind segmentation/review validation | Reduces confirmation and curation bias | Held-out species/batch evaluation and repeat-review agreement |
| Preserve full trait distributions | Avoids reducing all biology to a single point | Long-form object/image/specimen table in supplement |
| Use posterior tree draws | Captures topology/branch-length uncertainty | Distribution of rankings/effects across trees |
| Use an ILR/pivot balance | Improves interpretable CoDA | Prespecified contrast and sensitivity to coordinate choice |
| Audit library fragmentation | Protects family-diversity claims | Curated merge/split sensitivity and unclassified fraction |
| Compare direct and proxy deletion metrics in more than one assembly | Establishes generality beyond the D. fuscus anchor | Rank/effect concordance across assemblies |
| Predeclare primary figures/models | Distinguishes confirmation from exploration | Timestamped analysis plan or frozen candidate-model file |
| Add negative controls | Reveals technical associations | Shuffled labels, technical covariates, and null simulations |

### Supplementary rather than blocking

| Extension | Why supplementary |
|---|---|
| Exhaustive superfamily-level path networks | Too many parameters for the current species sample |
| Interactive dashboards | Helpful for collaboration, not inferential evidence |
| Full ancestral-state visualizations | Contextual and model-dependent |
| Every substitution-rate value ever published | The defensible rate window is more useful |
| ML-based automatic curation beyond the frozen validation set | Useful operationally but not needed if the primary masks are auditable |
| Phylogenetic imputation to expand primary panels | Better retained as sensitivity because it does not add independent observations |

## 10. Recommended manuscript and figure hierarchy

### Primary result candidates

1. **TE composition/diversity versus independently defensible genome size**
   - exact 18-species panel;
   - one or two prespecified TE coordinates;
   - species-level genome uncertainty;
   - raw and phylogenetic relationships;
   - assembly/coverage validation.

2. **Genome size, nucleus size, and cell size**
   - only primary if genome size is independently validated;
   - otherwise title it as a linked-image association analysis;
   - show all-eligible/random-50/top-50 stability.

### Secondary mechanism results

- paired-LTR K2P history versus genome/TE state;
- terminal:internal depth versus genome/TE state after proxy validation;
- combined path families only when model recovery is adequate.

### Sensitivity/supplement

- current largest-50 specification;
- strict-body panels with 8-9 tips;
- phylogenetically imputed traits;
- absolute LTR ages under alternative rates;
- unvalidated ectopic depth formulations;
- superfamily-rich exploratory models.

### Figure set

1. **Study and sample-flow figure**
   - species through TE, genome, morphology, LTR, ectopic, and final panels;
   - specimens/images/cells made visible.

2. **Measurement-validation figure**
   - segmentation examples and held-out errors;
   - IOD linearity/batch/reference validation;
   - species support tiers.

3. **Morphology specification figure**
   - within-species distributions;
   - all eligible versus random 50 versus top 50;
   - specimen-level points, not only species summaries.

4. **Phylogenetic trait figure**
   - dated tree with genome, nucleus, cell, and selected TE coordinates;
   - uncertainty/support encoded without implying causality.

5. **Primary bivariate comparative figure**
   - raw species points with measurement intervals;
   - phylogenetic fit and uncertainty;
   - labels for influential tips.

6. **DAG/model-comparison figure**
   - complete candidate set;
   - C, CICc, delta, weight;
   - path coefficients with bootstrap intervals.

7. **TE validation/composition figure**
   - read-based versus assembly-based anchor;
   - selected log-ratio coordinate;
   - unclassified component.

8. **LTR history/deletion figure**
   - K2P distributions and conditional age scale;
   - direct solo:intact versus terminal:internal proxy validation.

9. **Robustness specification curve**
   - cell selection, genome state, TE denominator, tree, and tip-removal effects in one compact panel.

Every figure caption should state the analysis unit, number of species, number of specimens/images where relevant, transformations, phylogenetic model, interval type, and whether the panel is primary or sensitivity.

## 11. Concrete stop/go criteria

### Go for a primary claim when

- the estimator and biological target are the same;
- technical provenance is complete;
- species are the inferential units and lower-level hierarchy is propagated;
- effect direction and interpretation survive prespecified reasonable alternatives;
- the candidate model is recoverable on the exact tree at plausible effect sizes;
- no single species or clade determines the result;
- the coefficient interval is biologically informative;
- the claim language matches whether the evidence is descriptive, associational, or mechanistic.

### Stop or downgrade to sensitivity when

- a trait is calibrated from another variable in the same causal path;
- one image/specimen is presented as if 50 cells were 50 biological replicates;
- model rank changes with one tip, one QC state, or one reasonable denominator;
- CICc identifies a “best” model that simulations cannot reliably recover;
- an absolute LTR age is driven by a borrowed rate without an uncertainty window;
- a depth ratio is called a recombination rate without direct solo-LTR validation;
- the notebook cannot rebuild from a clean kernel and frozen inputs.

## 12. Source map

### Salamander genome/cell/nucleus biology

- [Hally et al. 1986. Cytophotometric evidence of variation in genome size of desmognathine salamanders.](https://doi.org/10.1007/BF00494802)
- [Sessions and Larson 1987. Developmental correlates of genome size in plethodontid salamanders.](https://doi.org/10.1111/j.1558-5646.1987.tb02463.x)
- [Mueller et al. 2008. Genome size, cell size, and the evolution of enucleated erythrocytes in attenuate salamanders.](https://doi.org/10.1016/j.zool.2007.07.010)
- [Sun et al. 2012. LTR retrotransposons contribute to genomic gigantism in plethodontid salamanders.](https://doi.org/10.1093/gbe/evr139)
- [Itgen et al. 2022. Genome size drives morphological evolution in organ-specific ways.](https://doi.org/10.1111/evo.14519)
- [Myers et al. 2025. First complete assembly for D. fuscus.](https://doi.org/10.1093/g3journal/jkaf157)

### Cytometry, hierarchy, and measurement error

- [Hardie, Gregory, and Hebert 2002. From Pixels to Picograms.](https://doi.org/10.1177/002215540205000601)
- [Hurlbert 1984. Pseudoreplication and the design of ecological field experiments.](https://doi.org/10.2307/1942661)
- [Ives, Midford, and Garland 2007. Within-species variation and measurement error in phylogenetic comparative methods.](https://doi.org/10.1080/10635150701313830)
- [Felsenstein 2008. Comparative methods with sampling error and within-species variation.](https://doi.org/10.1086/587525)

### Phylogenetic path analysis

- [von Hardenberg and Gonzalez-Voyer 2013. Phylogenetic confirmatory path analysis.](https://doi.org/10.1111/j.1558-5646.2012.01790.x)
- [Shipley 2013. AIC for path models compared using d-separation.](https://doi.org/10.1890/12-0976.1)
- [van der Bijl 2018. phylopath.](https://doi.org/10.7717/peerj.4718)
- [Boettiger, Coop, and Ralph 2012. Is your phylogeny informative?](https://doi.org/10.1111/j.1558-5646.2011.01574.x)

### TE annotation, composition, and LTR methods

- [Aitchison 1982. The statistical analysis of compositional data.](https://doi.org/10.1111/j.2517-6161.1982.tb01195.x)
- [Goubert et al. 2015. dnaPipeTE.](https://doi.org/10.1093/gbe/evv050)
- [Flynn et al. 2020. RepeatModeler2.](https://doi.org/10.1073/pnas.1921046117)
- [Ou and Jiang 2018. LTR_retriever.](https://doi.org/10.1104/pp.17.01310)
- [Ou et al. 2018. LTR Assembly Index.](https://doi.org/10.1093/nar/gky730)
- [Kimura 1980. A simple method for estimating evolutionary rates of base substitutions.](https://doi.org/10.1007/BF01731581)
- [SanMiguel et al. 1998. The paleontology of intergene retrotransposons of maize.](https://doi.org/10.1038/1695)
- [Kijima and Innan 2010. Estimation of LTR retrotransposon insertion time.](https://doi.org/10.1093/molbev/msp295)
- [Frahry et al. 2015. Low levels of LTR retrotransposon deletion by ectopic recombination in salamanders.](https://doi.org/10.1007/s00239-014-9663-7)

### Reproducibility

- [Sandve et al. 2013. Ten simple rules for reproducible computational research.](https://doi.org/10.1371/journal.pcbi.1003285)
- [Nosek et al. 2015. Promoting an open research culture.](https://doi.org/10.1126/science.aab2374)
- [Wilkinson et al. 2016. FAIR Guiding Principles.](https://doi.org/10.1038/sdata.2016.18)

## Bottom line

The project should not discard its present analyses. It should make their estimands and evidence tiers explicit.

The current repository is closest to publication when it treats:

- TE/genome relationships as the primary comparative result, conditional on genome and TE measurement validation;
- cell and nucleus measurements as hierarchical species estimates whose answer must be stable to selection design;
- genome-nucleus-cell paths as linked-measurement associations unless genome size is independently assayed;
- paired-LTR K2P as measured history and calibrated ages as conditional;
- terminal:internal depth as a proxy awaiting direct solo-LTR validation;
- CICc ranking as one component of evidence, completed by diagnostics, bootstrap intervals, influence tests, tree uncertainty, and exact-design simulations.

That framing is both more conservative and more scientifically informative: it will show readers exactly which conclusions arise from robust comparative signal, which depend on a particular measurement choice, and which are promising mechanistic hypotheses for follow-up.

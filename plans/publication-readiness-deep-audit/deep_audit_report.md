# Desmognathus TE Project: Publication-Readiness Deep Audit

**Audit date:** 2026-07-09
**Identity clarification:** Updated 2026-07-09 after direct guidance from Alex Pyron and an accession-level repository trace.
**Scope:** accession identity, TE processing, LTR analyses, erythrocyte segmentation and selection, genome-size image cytometry, phylogenetic comparative statistics, figures, documentation, and release-notebook design.
**Change boundary:** the initial pass was read-only. Subsequent remediation is
non-destructive: no historical analysis, data product, or figure was removed or
overwritten; corrected artifacts use versioned paths under `results/*/corrected/`,
`plans/publication-readiness-deep-audit/`, and `notebooks/`.

**Remediation update, 2026-07-09:** Historical products remain preserved, but
non-destructive corrected branches have now been built for accession identity,
dnaPipeTE mass accounting, RepeatMasker hit classification, divergence,
terminal:internal LTR depth, TE diversity/compositional PCA, focal and published
phylogenies, phylogenetic PCA, and microscopy. The microscopy release audit is
[`microscopy_release_audit_analysis18_v1.md`](microscopy_release_audit_analysis18_v1.md).
It scores the exact archived production cell masks, finds only diagnostic
two-tile segmentation support, defines the morphology trait as a composite
upper-tail sensitivity estimand, and withdraws the unsupported image-IOD
picogram conversion. Corrected products are under `results/data/corrected/`
and `results/figures/corrected/`. The corrected exploratory path audit and
executed release notebook are now complete:
[`corrected_path_analysis_audit_analysis18_v1.md`](corrected_path_analysis_audit_analysis18_v1.md)
and
[`../../notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb`](../../notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb).

The corrected path branch completed 958 fits with no fitting failures, including
all 200 published time trees and leave-one-species-out influence. Its actual-tree
simulation still blocks confirmatory promotion: strict unique false-selection
rates under independent traits are 19.0% for TE–IOD, 20.0% for IOD–morphology,
and 9.5% for the integrated family; observed-scale integrated-model recovery is
only 50%. These diagnostics reinforce rather than relax the genome-size and
causal-claim block.

## Executive verdict

The project has an unusually rich and potentially valuable evidence chain, but the current integrated causal analysis is **not yet publication-ready**. This is not mainly a plotting or notebook-polish problem. There is one confirmed upstream TE-classification error, one unresolved genomic-resource provenance/sensitivity gate, several measurement-validity gaps, and statistical outputs that are currently being interpreted more strongly than their diagnostics permit.

The two most urgent upstream gates are:

1. The accession decision behind the current `D. fuscus` row is not encoded in a machine-readable provenance record. The row uses `SRX20497025 / GCA_032353935.1`, which NCBI labels *D. planiceps* but Alex Pyron identifies as actual *D. fuscus*. A separate reliable 2025 *D. fuscus* assembly (`GCA_050004315.1`) exists but is not used by the current repository outputs. These resources must remain accession-explicit, and microscopy must remain an independent specimen stream.
2. The RepeatMasker merge assigns one dnaPipeTE classification to every RepeatMasker hit on a contig. In the current 9.53-million-hit table, 18.148% of hits with comparable known orders disagree with their own RepeatMasker order classification. Order-specific divergence and alignment-indel summaries therefore need to be rebuilt.

The central measurement issue is equally consequential. Current genome-size values are relative image-IOD estimates calibrated to a hard-coded `D. fuscus = 16.36 pg`; most species and the current reference are supported by only one slide/image/specimen, the exact 16.36-pg value is not sourced in the repository, and the wet-lab/acquisition metadata required to establish between-slide Feulgen comparability are absent. A current alternative-QC sensitivity changes all non-reference species upward by 19.8%–74.0% (median 49.7%). That alternative is itself not ready to replace the primary panel because it admits visibly broken cell masks; the result demonstrates estimator fragility, not which estimator is correct.

The safest path is to preserve every current output as a named pre-audit or sensitivity branch, encode genomic-resource identity and repair TE classification first, explicitly define the morphology and genome-size estimands, then rerun a frozen and simulation-validated phylogenetic analysis. Only after those gates should the project freeze a portable Jupyter release bundle.

## Release disposition by evidence layer

| Evidence layer | Current state | Publication disposition |
|---|---|---|
| Accession/taxon crosswalk | Public `planiceps` label versus expert *fuscus* identification is encoded as the current reidentified short-read resource; the separate chromosome-level 2025 *fuscus* resource benchmark is accession-explicit but lacks annotation assets for full parity | **Use the current accession-explicit short-read branch; retain the chromosome-level validation request** |
| dnaPipeTE relative composition | Useful processed signal, but closed to 100%; historical all-resource processing also sums one retry output | **Corrected final-panel mass ledger complete; absolute load still unresolved** |
| RepeatMasker divergence/indels | Corrected hit-native classification and conserved divergence branches are complete; historical tables remain preserved | **Use corrected final-18 branch; keep gap metrics out of DNA-loss language** |
| Paired-LTR K2P | Conceptually defensible for detectable paired LTRs, overwhelmingly Gypsy, with sparse support for some species | **Conditional/sensitivity until identity and ascertainment gates pass** |
| Terminal:internal LTR depth | Outlier-sensitive, zero-depth and mapping provenance gaps, not independently validated as ectopic recombination | **Exploratory proxy only** |
| Cell/nucleus morphology | Final-18 linkage, six estimands, rank sensitivity, hierarchical conditional intervals, and segmentation diagnostics are frozen; final-panel held-out validation remains absent | **Upper-tail sensitivity only; no general segmentation claim** |
| Image-IOD genome size | Calibration and between-slide comparability are not established; IOD is algebraically area × mean OD | **Relative nuclear-IOD proxy only; not pg or C-value** |
| Phylogenetic path analysis | Explicit-null families, global fit, CICc, 200-tree, species-influence, and simulation outputs are complete; measurement validity and model-recovery remain limiting | **Exploratory association sensitivity; no causal/genome-size claim** |
| Current notebook | Portable relative-path workbench executes from corrected artifacts and exposes a weird-result queue | **Approved for collaborator audit, not as evidence that blocked claims are valid** |

## Original pre-remediation numeric snapshot

These counts were recomputed from the live historical artifacts during the
initial audit. They are retained as the diagnosis that motivated remediation;
current corrected counts and verdicts are in the versioned audit reports.

- The dated tree has 46 tips and contains distinct `fuscus`, `fuscus_A`, `fuscus_E`, and `planiceps` tips. Its maximum root-to-tip difference is approximately `2e-6`, which is biologically negligible but causes a strict `ape::is.ultrametric()` check to return false unless tolerance or a documented correction is used.
- Current integrated panels contain 18 species for TE–genome, 15 for paired-LTR history, 16 for ectopic depth, 21 for genome–morphology, and 18 for TE–genome–morphology. Strict organismal panels contain only 8–9 species.
- The frozen morphology input contains 1,050 cell–nucleus pairs: exactly 50 for each of 21 species. Species candidate pools range from 76 to 701 rows.
- Only 440 frozen rows are `reviewed_keep`; 5 are `reviewed_maybe`; 605 are `model_ranked` and effectively unlabeled. The frozen output rejects explicit `discard` and `nucleus_only` rows but accepts `maybe` and unlabeled rows.
- The selected set overlaps the literal 50 largest cells by only 17–36 cells per species (mean 22.429). The median selected cell-area percentile ranges from 0.586 to 0.902, and candidate-pool size correlates with that percentile at Spearman `rho = 0.773`.
- Morphology summaries use a median of two images/specimens per species; three species have one. Current primary genome estimates are much thinner: most supported species use one image and one specimen.
- In the current 18-species TE–genome model family, `mediated_evenness` has model weight 0.999, but absolute fit is borderline (`Fisher C p = 0.0896`) and the `te_evenness -> genome size` edge is `-0.433 ± 0.225 SE`; an approximate 95% interval crosses zero.
- In the 16-species ectopic family, **all four candidate DAGs fail the d-separation goodness-of-fit test** (`p < 0.003`). Calling `ectopic_only` the winner because its relative model weight is 0.900 does not make the family adequate.
- In the 8-species strict ectopic panel, CICc is infinite or undefined for two models and all model weights are `NA`, yet a “best model” file is still exported.
- Current genome-QC labels contain no stable species: the live panel summary marks 14/18 TE–genome species as sensitivity-limited and four as minor caution.

The reproducible measurement calculations are in [`measurement_selection_metrics.py`](measurement_selection_metrics.py). The accession-labeled current-versus-2025 *fuscus* comparison is in [`fuscus_resource_benchmark.md`](fuscus_resource_benchmark.md), with its machine-readable values in [`fuscus_resource_benchmark.csv`](fuscus_resource_benchmark.csv). The fuller literature comparison is in [`literature-methods-benchmark.md`](literature-methods-benchmark.md). The prior focused genome-size audit remains relevant and is preserved at [`../size-estimation-deep-audit/deep_audit_report.md`](../size-estimation-deep-audit/deep_audit_report.md).

## P0: blockers before any integrated causal result is shared as final

### [P0-01] Encode the accession-specific *D. fuscus* decision and compare the current resource with the 2025 assembly

- **Status:** Biological analysis label resolved by expert guidance; provenance encoding and resource sensitivity remain open.
- **Evidence:** [`input_data/lookup_table.txt:18`](../../input_data/lookup_table.txt) maps `D.fuscus` to `SRX20497025` and `GCA_032353935.1`. Current NCBI records label both as *Desmognathus planiceps*: [SRA SRX20497025](https://www.ncbi.nlm.nih.gov/sra/SRX20497025), [assembly GCA_032353935.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032353935.1/). Alex Pyron identifies this resource as actual *D. fuscus*, so its use on the `fuscus` tree tip is the intended analysis decision rather than evidence that a *planiceps* row should replace it.
- **Evidence:** A distinct modern *D. fuscus* resource, `aDesFus1-2.1 / GCA_050004315.1`, was published in 2025 from PacBio data: [Myers et al. 2025](https://doi.org/10.1093/g3journal/jkaf157), [NCBI assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_050004315.1/). Repository-wide searches found no assembly, accession, read-run, EarlGrey annotation, or derived output from this resource. Current dnaPipeTE, RepeatMasker, LTR-age, and ectopic outputs all trace instead to `SRX20497025 / GCA_032353935.1` and its `JAUCNN` contigs.
- **Evidence:** The tree at [`input_data/phylogeny/desmo900dated_test.tre`](../../input_data/phylogeny/desmo900dated_test.tre) has separate `fuscus` and `planiceps` tips. [`scripts/processing/clean_tree_phylo.R:180`](../../scripts/processing/clean_tree_phylo.R) and [`path_analysis/scripts/path_model_scaffold.R:252`](../../path_analysis/scripts/path_model_scaffold.R) match lookup/data labels to exact tree-tip names; the current genomic values therefore land on the `fuscus` branch while `planiceps` is dropped.
- **Evidence:** The separate older resource labeled *D. fuscus* by NCBI, `SRX19953421 / GCA_030265095.1`, is one of the resources Alex Pyron does not trust. Its staged RepeatMasker file is explicitly excluded by [`scripts/processing/repeatmask.py`](../../scripts/processing/repeatmask.py), and it contributes no current dnaPipeTE, LTR-age, or filtered ectopic row.
- **Evidence:** The integrated species table combines species-level summaries across genomic and microscopy evidence. The erythrocyte images were collected from different biological samples than the genomic resources; any genomic accession shown with microscopy-derived species summaries is a species-level analysis reference, not a matched-specimen identifier.
- **Impact:** The current `fuscus` row is not shown to be biologically incoherent. The remaining risks are reproducibility of the expert reidentification, ambiguity about which resource generated each result, and assembly/method sensitivity when comparing the short-read resource with the 2025 PacBio assembly.
- **Effort:** M — encode the decision, keep specimen streams separate, and run an accession-labeled old-versus-2025 sensitivity comparison before choosing any replacement.
- **Risk:** Medium for current read-based TE composition; potentially high for assembly-dependent LTR and ectopic analyses because the available assemblies differ greatly in completeness and sequencing technology.
- **Confidence:** High on the repository lineage and current analysis decision; medium on how much assembly quality changes the LTR/ectopic results until the paired comparison is run.
- **Fix sketch:** Preserve the current `D. fuscus` outputs as the `SRX20497025 / GCA_032353935.1` analysis branch. Extend the existing source manifest and taxonomy crosswalk with public label, expert identification, voucher/BioSample, technology, analysis role, and inclusion decision. Run the 2025 assembly and published EarlGrey annotation as an explicitly separate validation/sensitivity branch. Do not infer genomic specimen identity for microscopy records, and do not create or analyze a *D. planiceps* TE row without independently selected evidence.

### [P0-02] Rebuild RepeatMasker summaries from each hit's own TE classification

- **Status:** Confirmed implementation error and current-data corruption.
- **Evidence:** [`scripts/processing/repeatmask.py:341`](../../scripts/processing/repeatmask.py) sorts dnaPipeTE annotations per contig and keeps only the longest classification; [`repeatmask.py:404`](../../scripts/processing/repeatmask.py) then joins that single class to every RepeatMasker hit on the contig. [`scripts/processing/divergence.py:48`](../../scripts/processing/divergence.py) groups by the inherited `Class`, `Order`, and `Superfamily`, not the hit's native `repeat_class`.
- **Evidence:** The first current merged record is an immediate example: RepeatMasker calls the hit `LTR/Gypsy`, while the inherited fields call it `LINE/L2`.
- **Evidence:** The reproducible shared-map audit in [`repeatmasker_hit_classification_audit.md`](repeatmasker_hit_classification_audit.md) scanned all 9,533,364 merged hits. It found 1,502,703 of 8,280,403 hits with comparable known orders mismatched (18.148%; 11.426% of comparable aligned bp). Known-superfamily mismatch was 20.381% of comparable hits and 12.762% of comparable aligned bp. Only 137 native RepeatMasker hits (0.0014%) used labels absent from the shared exact-label map; those 19 labels are retained in [`repeatmasker_unmapped_labels.csv`](repeatmasker_unmapped_labels.csv).
- **Impact:** Order- and superfamily-specific divergence, insertion, and deletion summaries are assigned to the wrong biological categories. Any path variable or figure derived from those summaries is currently unreliable.
- **Effort:** L — correct classification, rerun summaries, rebuild path inputs/models/figures, and re-audit conservation.
- **Risk:** High — core mechanistic predictors can change.
- **Confidence:** High.
- **Fix sketch:** Parse class/order/superfamily from each RepeatMasker hit's `repeat_class`. Keep the dnaPipeTE contig label only as a separate comparison field. Emit disagreement rates and unmatched mass. Add tests that preserve hit count and aligned bp through every split. Rename RepeatMasker `%del` and `%ins` as **alignment gap statistics relative to consensus**, not genomic DNA-loss or ectopic-recombination rates.

### [P0-03] Do not present current image-IOD values as publication-stable absolute genome sizes

- **Status:** Confirmed validation and provenance gap; the exact source of error among staining, acquisition, background, and mask selection is not yet separable.
- **Evidence:** The pipeline hard-codes `D. fuscus = 16.36 pg` (for example, [`path_analysis/scripts/build_balanced_genome_iod_sensitivity.py:39`](../../path_analysis/scripts/build_balanced_genome_iod_sensitivity.py)); repository-wide search found no DOI, database record, 1C/2C convention, or derivation for exactly 16.36 pg.
- **Evidence:** Published *D. fuscus* values are heterogeneous, while the classic Desmognathus cytometry study used an internal *Xenopus laevis* erythrocyte standard. Modern Feulgen protocols emphasize controlled fixation, hydrolysis, staining, camera linearity, and co-stained standards rather than assuming independent slides are directly comparable. See [Hally et al. 1986](https://doi.org/10.1007/BF00494802), [Hardie et al. 2002](https://doi.org/10.1177/002215540205000601), and [Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010).
- **Evidence:** The primary calibration is based on 23 selected nuclei from one *D. fuscus* image/specimen. Most species likewise have one image/specimen. The image metadata do not record fixation, hydrolysis duration, stain/reagent batch, co-stained standard, instrument, illumination, exposure/gain, scanner settings, saturation/linearity checks, or slide age. OME metadata retain pixel size and green-channel identity but not those acquisition variables.
- **Evidence:** The current balanced-QC sensitivity shifts every non-reference species upward by 19.821%–74.001% (median 49.743%; 16/20 at least 40%). It changes the reference from 23 nuclei/one image/one specimen to 50 nuclei/eight images/four specimens and changes the pg/IOD scale substantially. However, visual/metric audit shows that this alternative also admits broken cell masks and changes cell-area ranks sharply, so it is diagnostic rather than a replacement.
- **Impact:** Between-slide staining/acquisition variation and selection of the reference nuclei can change the species ranking and all genome-related path coefficients. Because genome IOD and nucleus area are derived from linked nuclear images, a genome-size-to-nucleus-size causal path is also partly circular unless genome size is independently assayed.
- **Effort:** L for independent validation; M to reconstruct methods and reframe the current values honestly.
- **Risk:** High — this affects the central genome–nucleus–cell claim.
- **Confidence:** High that absolute comparability is not established; medium on the contribution of each technical source.
- **Fix sketch:** Recover and cite the reference C-value and convention; reconstruct the complete wet-lab and imaging batch record; test camera linearity and saturation; model slide/batch explicitly; and validate a representative subset with flow cytometry, published C-values, or same-slide internal standards. Until then, call these **relative image-IOD genome-size sensitivity estimates** and keep `genome -> nucleus` causal language exploratory.

### [P0-04] Freeze an honest morphology estimand and finish the cell-level evidence state

- **Status:** Confirmed mismatch between method name and implemented estimator.
- **Evidence:** The candidate score is 80% cell-size percentile, 10% nucleus darkness, 7% clarity, and 3% other quality terms in `/home/jake/Projects/cellprofiler_test/scripts/build_top50_linked_pair_worklist.py:341-350`. The final grid prioritizes quality rank and pair-quality probability before the composite and cell area. The freeze script excludes only `discard` and `nucleus_only`, while its documentation calls the output an “exact reviewed top-50.”
- **Evidence:** Of 1,050 frozen pairs, 440 are reviewed keeps, five are reviewed maybe, and 605 are model-ranked/unlabeled. The literal-largest-50 overlap averages only 22.429 per species.
- **Evidence:** A fixed 50 samples different parts of each candidate distribution: candidate pools are 76–701, selected median percentiles are 0.586–0.902, and candidate-pool size is strongly associated with the selected percentile (`rho = 0.773`) and selected-versus-candidate size uplift (`rho = 0.751`).
- **Impact:** The current input estimates a quality-screened upper tail, not the literal maximum-cell distribution or population-typical erythrocyte size. Differential image yield and segmentation success can become a species-level covariate.
- **Effort:** M — no wholesale resegmentation is necessarily required, but the estimand, remaining review, and sensitivity set must be frozen.
- **Risk:** Medium-high — cell-size ranks and path coefficients can change under a different selection rule.
- **Confidence:** High.
- **Fix sketch:** Choose and name the biological estimand before looking at the path results. Compare at least: all eligible cells with hierarchical weighting; repeated random 50; a common upper quantile; literal largest 50; current quality-ranked upper-tail 50; and manual-only keeps. Resolve all `maybe` rows, distinguish human-accepted from model-accepted records, and validate the chosen estimator on a blinded manual sample. Preserve every current estimator as a labeled sensitivity panel.

### [P0-05] Treat model weights as relative comparison only; current path outputs do not establish the proposed causal mechanism

- **Status:** Confirmed statistical interpretation error/gap.
- **Evidence:** The 18-species TE model's top weight is 0.999, but its absolute d-separation fit is only `p = 0.0896`, and the focal `te_evenness -> gs` edge has an approximate interval crossing zero. Relative weight says which candidate is least unsupported within the supplied family; it does not validate an edge or the family.
- **Evidence:** Every model in the current 16-species ectopic family is rejected (`Fisher C p < 0.003`), yet a winner is still named. This indicates omitted dependencies or a misspecified family, not support for an ectopic-only mechanism.
- **Evidence:** The 8-species strict ectopic run contains infinite/undefined CICc values and all `NA` weights, but the scaffold still calls `best()` and writes best-model artifacts.
- **Evidence:** `te_evenness` and `ltr_balance` are derived from the same closed composition. A directed mediation relation between them can be a mathematical consequence of a shared denominator and is not causally identified by the current cross-species data. Several competing DAG directions are Markov-equivalent at this sample size.
- **Evidence:** Measurement uncertainty, tree uncertainty, exact-tree model-recovery simulations, and coefficient bootstrap intervals are absent. Lambda estimates frequently hit parameter boundaries.
- **Impact:** The current diagrams can look decisive while the focal edge is uncertain, a whole candidate family is inadequate, or the information criterion is undefined.
- **Effort:** L — refreeze hypotheses and rebuild the inferential layer after upstream corrections.
- **Risk:** High — causal conclusions may narrow or change.
- **Confidence:** High.
- **Fix sketch:** Define confirmatory DAGs before rerunning corrected data; reject a family when its global d-separation test fails; refuse winner export for non-finite CICc/weights; report all basis-set tests, standardized edges with bootstrap intervals, model fit, and influence. Use a single prespecified ILR/log-ratio TE contrast or latent/compositional block rather than causal arrows between mechanically related diversity summaries. Run simulation-based model recovery on the exact tree, sample sizes, missingness, effect ranges, and measurement error.

## P1: high-priority methods and reproducibility gaps

### [P1-01] Demonstrate segmentation and pair-quality generalization across species, slides, and acquisition conditions

- **Status:** Confirmed validation-design gap.
- **Evidence:** The Cellpose cell training set contains 12 annotated tiles: 10 train and two test, from only two species and two source whole-slide images (`D. ochrophaeus` Process_316 and `D. folkertsi` Process_366). Each test tile comes from one of those same source images.
- **Evidence:** The final nucleus YOLO dataset has 156 training and 32 validation patches with 536 instances, but all 12 source tiles trace to the same two images/species. Patch-level validation therefore does not test new taxa, specimens, slides, or acquisition batches.
- **Evidence:** The pair-quality classifier uses row-stratified five-fold CV, not grouped CV by species/image/specimen. Its balanced accuracy is 0.738; discard precision is 0.588 and discard recall 0.696. Those are useful triage metrics, not proof that 605 model-ranked final pairs are publication-quality.
- **Impact:** Cross-species differences in stain, morphology, focus, overlap, and background can produce species-specific segmentation error that mimics biological differences.
- **Effort:** L if annotation must expand; M for a carefully stratified audit sample.
- **Risk:** Medium-high.
- **Confidence:** High.
- **Fix sketch:** Annotate a taxonomically and technically stratified validation set across the 21 species, multiple slides, stain depths, cell sizes, and known failure modes. Hold out entire images/specimens/species. Report instance detection, mask IoU/boundary error, area bias, nucleus–cell pairing error, and downstream species-summary bias. Include failure-case panels, not only representative successes.

### [P1-02] Make specimen/slide the biological unit and distinguish technical from biological uncertainty

- **Status:** Confirmed inference gap.
- **Evidence:** Morphology has one to four images/specimens per species (median two); current genome estimates are usually one. The implemented bootstrap is hierarchical where replication exists, but collapses to resampling cells when a species has one specimen/image.
- **Impact:** Fifty cells from one animal do not estimate among-animal species variance. An object bootstrap interval is conditional technical/segmentation uncertainty, not a population-level confidence interval.
- **Effort:** M analytically; L if additional biological specimens must be collected.
- **Risk:** Medium — intervals may widen or become explicitly unavailable.
- **Confidence:** High.
- **Fix sketch:** Summarize or model cell measurements within image, image within specimen, and specimen within species. Where at least several specimens exist, use a hierarchical measurement model or specimen bootstrap. For one-specimen taxa, report the observed specimen estimate and conditional cell-level spread while marking biological uncertainty as unestimated. Propagate species-estimate draws into phylogenetic models.

### [P1-03] Preserve absolute TE mass and unresolved classification mass

- **Status:** Unresolved-mass accounting repaired; a configured-denominator absolute-load sensitivity is recovered for the final 18-species panel but is not confirmatory.
- **Evidence:** [`scripts/processing/dnaPipe.py:473`](../../scripts/processing/dnaPipe.py) aggregates dnaPipeTE aligned bases by retained category; [`dnaPipe.py:551`](../../scripts/processing/dnaPipe.py) divides by the retained row sum, forcing each species to 100%. The current predictors therefore encode relative TE composition, not genome-wide TE load.
- **Evidence:** Generic `LTR`, `LTR/ERV1`, `Unknown TIR`, and other unresolved classes are omitted before superfamily renormalization. Current order tables retain 91.64%–96.45% of aligned bases, but superfamily tables retain only 73.13%–82.49%. `LTR/ERV1` alone is 7.074% of all aligned bases and disappears from the superfamily denominator.
- **Implemented evidence:** [`dnapipete_mass_accounting_analysis18_v2.csv`](../../results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv) preserves total aligned-repeat bases and retained/unresolved mass for exactly the final 18 species. Its manifest matches every observed SRX to the active lookup and records the assembly accession separately. Order-level unresolved mass is 3.548%–6.179% (median 5.171%); superfamily unresolved mass is 17.507%–21.586% (median 19.995%). The retained-order percentages reproduce the historical table within 7.105e-15 percentage points. *D. orestes* and all other out-of-panel resources are absent. The earlier v1 snapshot remains preserved but is superseded for analysis.
- **Implemented evidence:** Git history recovers `dnaPipeTE.sh` at commit `cc64ffb134968250671ff0e3456255063a1aba97`, with a 15-Gb configured genome and 0.1× coverage. [`dnapipete_absolute_load_sensitivity_analysis18_v1.csv`](../../results/data/corrected/dnapipete/dnapipete_absolute_load_sensitivity_analysis18_v1.csv) therefore divides conserved aligned-repeat bases by the configured 1.5-Gb independent quantification sample. Repeat-aligned fractions are 60.619%–68.196% (median 65.744%); *D. fuscus* is 66.010% repeat-aligned and 62.536% TE-classified.
- **Impact:** Two genomes with very different total repeat fractions can have identical current predictors. Diversity/PCA/PERMANOVA differences can partly reflect unequal classification completeness.
- **Effort:** M to rerun/freeze the configured absolute layer reproducibly; explicit denominator and unresolved-bin reconstruction is complete.
- **Risk:** Medium.
- **Confidence:** High.
- **Remaining fix:** Recover the container digest/version, custom-library checksum, runtime logs, and replicate quantification samples, or rerun the final 18 accessions with those assets frozen. Treat the current absolute-load table as sensitivity-only. Analyze absolute repeat fraction separately from closed composition; for the corrected relative layer, carry unresolved mass and use CLR/ILR or prespecified log-ratios with transparent zero handling.

### [P1-04] Validate and robustify the terminal:internal LTR depth proxy

- **Status:** Corrected robustness and coverage audit implemented; biological interpretation remains blocked.
- **Evidence:** [`scripts/processing/ec.py:108`](../../scripts/processing/ec.py) removes every zero-depth position and imposes no covered-fraction or mapping-quality threshold. Mapping/reference, multimapper, secondary-alignment, duplicate, MAPQ, and zero-reporting settings are not preserved in the current provenance.
- **Evidence:** Species contribute 12–122 elements. The path input uses the species arithmetic mean ratio although a median is also available. For `D. intermedius`, one outlier makes the mean 3.185 and median 0.646, a 4.93-fold difference; the maximum element ratio is 57.773.
- **Implemented evidence:** [`ectopic_element_metrics_analysis18_v1.csv`](../../results/data/corrected/ectopic/ectopic_element_metrics_analysis18_v1.csv) contains 567 exact-coordinate, zero-aware 5/6-domain elements across 16 final-panel species. It reproduces the historical nonzero-only ratios within 4.441e-16, retains 2,872 explicit zeros, and shows that 565/567 elements pass 80% regional positive-position coverage. [`ectopic_species_robustness_analysis18_v1.csv`](../../results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv) adds median, geometric, trimmed and winsorized means, bootstrap intervals, complete-element/coverage branches, and leave-one-element-out influence.
- **Implemented evidence:** After zero-aware correction, *D. intermedius* still has arithmetic mean 3.114, median 0.646, and geometric mean 0.800; one element supplies 75.0% of its ratio sum and its removal shifts the arithmetic mean by 2.302. This is an estimand problem, not primarily a low-coverage problem.
- **Impact:** One element can dominate a causal predictor, unequal element counts are treated as equally precise, and repeat mappability/copy number can masquerade as solo-LTR formation.
- **Effort:** L for biological validation; robust statistical sensitivities and figures are complete.
- **Risk:** High for the ectopic mechanism, low for retaining it as an exploratory proxy.
- **Confidence:** High on instability; medium on exact biological bias.
- **Remaining fix:** Recover MAPQ/multimapper/reference/duplicate/secondary-alignment provenance and validate against direct intact-versus-solo LTR calls on assemblies or a manually checked subset. If retained as an exploratory mapping proxy, use the species median or mean element-level log-ratio with bootstrap and influence reporting—not the arithmetic mean of raw ratios. The published benchmark supports structural solo:intact evidence as the closer mechanistic quantity: [Frahry et al. 2015](https://doi.org/10.1007/s00239-014-9663-7).

### [P1-05] Replace unrestricted clade PERMANOVA and freeze the multiple-testing hierarchy

- **Status:** Historical test rejected for inference; descriptive CLR-PCA replacement and sensitivity products complete.

- **Status:** Confirmed invalid exchangeability and reporting gap.
- **Evidence:** Current PERMANOVA groups are defined from phylogenetic clades but species labels are permuted without respecting phylogenetic covariance. The permutation null therefore does not match the data-generating structure.
- **Evidence:** The only nominal superfamily Bray-Curtis result is `p = 0.037`; it would not survive a simple Benjamini-Hochberg correction across the four corresponding overall tests if those were one family and not prespecified.
- **Evidence:** Some candidate-model changes were made after seeing earlier results (for example, dropping weak aquaticity paths), which turns later comparisons into exploratory model development.
- **Impact:** Nominal community-level significance and selectively refined DAGs can overstate confirmatory evidence.
- **Effort:** M.
- **Risk:** Medium.
- **Confidence:** High.
- **Fix sketch:** Treat current results as discovery. Freeze one confirmatory family and endpoint hierarchy before corrected reruns. For composition, prefer phylogenetic regression/ordination on ILR coordinates or a justified phylogenetic permutation scheme; test multivariate dispersion separately. Report all tested families and adjust multiplicity within declared families.

### [P1-06] Quantify tree, leverage, and small-sample uncertainty

- **Status:** Tree structure, published 200-tree propagation, species influence,
  and model-recovery simulation are complete; collaborator focal-tree
  provenance and measurement-error propagation remain open.
- **Evidence:** The dated focal tree source/calibration record is still incomplete and its local source is ignored by Git. A versioned final-panel derivative now retains exactly 18 tips, uses no congener substitution, has positive bifurcating branches, and corrects at most two years of Newick rounding.
- **Evidence:** Stewart and Wiens (2025; DOI `10.1016/j.ympev.2024.108272`) supply an optimal time tree and 200 time-calibrated bootstrap trees. All 18 focal taxa are exact matches. Across the published set, final-panel crown age has median 13.859 Myr and 95% interval 12.166-16.403 Myr. The published topology differs from the collaborator focal tree by one rooted bipartition in the black-bellied complex; all 200 published trees share the published-main topology.
- **Implemented evidence:** Order-level phylogenetic PCA and all four corrected
  path families were refit across the published main tree plus all 200 time
  trees. All 804 tree-family fits completed; top-model identity was stable, and
  the terminal:internal family failed global fit on every tree.
- **Implemented evidence:** Seventy leave-one-species-out family fits export
  rankings and standardized edges. Actual-tree simulation includes 1,500
  selection replicates plus generating-edge bias, RMSE, and approximate interval
  coverage. It exposes non-negligible false selection and only 50% strict
  integrated-model recovery at observed coefficients.
- **Impact:** A few tips, alternative *fuscus* resource-derived measurements, or plausible tree/calibration variation can alter coefficients and model selection.
- **Effort:** L.
- **Risk:** High for direction-specific causal claims.
- **Confidence:** High.
- **Remaining fix:** Obtain the focal-tree citation/calibration record and
  propagate specimen-level microscopy uncertainty and TE uncertainty through an
  error-aware phylogenetic SEM. Treat the current candidate set as exploratory
  and freeze a prospective set only after new validation data exist. Boettiger
  et al. show why adequacy/power simulation is important for comparative model
  selection: [Boettiger et al. 2012](https://doi.org/10.1111/j.1558-5646.2011.01574.x).

### [P1-07] Bound paired-LTR claims to the detectable, mostly Gypsy subset

- **Status:** Confirmed ascertainment boundary; magnitude/direction of bias needs analysis.
- **Evidence:** Current successful estimates are 1,283/1,291 Gypsy. Species medians are unweighted, support can be sparse, and current summaries lack bootstrap intervals, minimum-pair gates, assembly-callable opportunity, or gene-conversion/boundary sensitivities.
- **Evidence:** Absolute ages use borrowed salamander/frog nuclear synonymous rates rather than a *Desmognathus* LTR-specific substitution rate.
- **Impact:** The result describes detectable paired, overwhelmingly Gypsy LTRs in assemblies, not the genome-wide insertion-age distribution. Fragmented assemblies preferentially recover a non-random subset.
- **Effort:** L.
- **Risk:** Medium — relative K2P may remain useful while absolute-age and burst language narrows.
- **Confidence:** Medium-high.
- **Fix sketch:** Use K2P divergence as the primary observed quantity; present calendar ages only across an explicit rate window. Bootstrap species medians, report pair counts and minimum-support gates, compare all versus structural/high-confidence pairs, quantify callable long-contig opportunity, and screen for gene conversion and boundary artifacts.

### [P1-08] Create a portable, immutable release bundle before creating the final notebook

- **Status:** Portable corrected collaborator notebook and manifest implemented;
  public archival packaging and environment locking remain open.
- **Evidence:** [`.gitignore:88`](../../.gitignore) excludes `input_data/`, `interim/`, and `results/`; a clean clone cannot inspect the canonical source layer. [`DATA_MANIFEST.yml`](../../DATA_MANIFEST.yml) records globs rather than per-file hashes, commands, parameters, tool/database versions, and parent-output relationships.
- **Evidence:** [`Dusky.yml`](../../Dusky.yml) leaves most dependencies unpinned and has no lock file. The live Dusky environment has `ape`, `phylopath`, and `phylolm`, but lacks `caper`, `vegan`, and `testthat`, so documented PGLS/PERMANOVA/test entrypoints do not all run now.
- **Evidence:** [`path_analysis/scripts/build_top50_size_notebook.py:73`](../../path_analysis/scripts/build_top50_size_notebook.py) hard-codes `~/Projects/cellprofiler_test`, mutable run directories, a local review URL, and local PDFs; it deletes old output files, records size/existence but not hashes, and globs model rankings rather than using a release allowlist.
- **Implemented evidence:**
  `notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb` uses relative
  repository paths, reads only versioned corrected tables for release results,
  executes without errors, embeds figures, and includes an automatic
  weird-result queue. Corrected manifests identify inputs/outputs and explicit
  claim blocks.
- **Remaining impact:** Large/raw inputs still depend on local acquisition and
  the environment is not fully locked for an archival clean-clone rerun.
- **Effort:** L for a complete evidence bundle; M after upstream outputs are frozen.
- **Risk:** Low — additive work, but essential for trust.
- **Confidence:** High.
- **Remaining fix:** Publish the canonical compact bundle and selected review
  exemplars in a versioned archive/DOI; manifest large inputs by accession and
  SHA-256; capture missing custom-library/container hashes; and lock the
  environment. Add strict read-time hash verification to the notebook before
  archival release.

## P2: cleanup and hardening findings

### [P2-01] Block TE-landscape regeneration until hit identity and mass are conserved

- **Evidence:** [`scripts/processing/parse_repeatmasker_landscape.py:76`](../../scripts/processing/parse_repeatmasker_landscape.py) uses `abs(end - start)` instead of inclusive length; [`:238`](../../scripts/processing/parse_repeatmasker_landscape.py) strips Trinity isoform suffixes before joining; [`:298`](../../scripts/processing/parse_repeatmasker_landscape.py) discards hit identity and deduplicates on coarse aggregate values.
- **Impact:** Distinct hits can be duplicated by many-to-many joins or collapsed by deduplication. Current `results/landscapes/` is empty, so this is a regeneration/figure blocker rather than evidence that current path inputs already contain this error.
- **Effort/Risk/Confidence:** M / medium / high.
- **Fix sketch:** Join on exact accession plus full contig, validate cardinality, retain unique hit ID and coordinates, use inclusive lengths, and assert source-versus-output hit-count and bp conservation.

### [P2-02] Keep dnaPipeTE computational retries from becoming biological replicates

- **Status:** Retry provenance clarified and excluded from the final 18-species corrected branch; archival selection between the two *D. orestes* realizations remains open.
- **Evidence:** The operator identifies suffixes such as `2` or `3` as computational retry labels added after failed/retried dnaPipeTE runs, not new specimens, SRA accessions, mates, or biological replicates. The staged repo contains one such case: `SRX19952890R2`.
- **Evidence:** `SRX19952890` and `SRX19952890R2` are both full-scale outputs (7,404,782 versus 7,404,964 summed read counts; 971,774,483 versus 971,590,301 aligned bases). Their classification profiles correlate at 0.99998, while stochastic contig identities differ. They are alternative computational realizations, not non-overlapping completion shards.
- **Impact:** Summing them double-weights one accession in the historical all-resource merge. This cannot affect the v2 analysis-scoped mass ledger because *D. orestes* is not processed, and the corrected writer rejects multiple source labels for any included species.
- **Effort/Risk/Confidence:** S / low / high; retry identity is operator-confirmed and file-scale evidence is consistent with full reruns.
- **Fix sketch:** Preserve retry suffixes as computational-run identifiers and select exactly one completed output per SRA using a run-selection manifest. Never sum or treat retries as biological replication. See `dnapipete_retry_output_audit.md`.

### [P2-03] Standardize the Simpson definition

- **Status:** Remediated and tested.
- **Evidence:** Both [`scripts/processing/diversity_stats.py`](../../scripts/processing/diversity_stats.py) and [`path_analysis/scripts/build_canonical_diversity_tables.py`](../../path_analysis/scripts/build_canonical_diversity_tables.py) now implement Gini-Simpson `1 - sum(p^2)` after normalization. Percent and proportion inputs are explicitly tested for scale invariance. Current stored values reproduce within `3.331e-16` at order and `4.441e-16` at superfamily level.
- **Impact:** The metric is now named and computed consistently; it remains conditional on the declared composition denominator.
- **Effort/Risk/Confidence:** S / low / high.
- **Fix sketch:** Complete. Preserve the explicit Gini-Simpson name and unresolved-mass sensitivity in all release outputs.

### [P2-04] Correct the source-material description

- **Evidence:** [`DATA_README.md`](../../DATA_README.md) says dnaPipeTE used RNA-seq and RepeatMasker used genome assemblies. Current NCBI metadata identify all 34 SRA experiments as paired genomic WGS, and the `SRX*_Trinity.align` files are RepeatMasker alignments on Trinity component contigs.
- **Impact:** The current methods describe the wrong molecule and substrate, preventing reviewers from assessing coverage and annotation scope.
- **Effort/Risk/Confidence:** M / low / high.
- **Fix sketch:** Document the actual shallow-WGS/dnaPipeTE/Trinity workflow, sampled bases/effective coverage, paired-read handling, custom library, classifier/database versions, and the distinction between read-supported relative composition and assembly-derived LTR analysis.

### [P2-05] Add tests at identity, denominator, and model-validity boundaries

- **Evidence:** The initial audit found only 16 Python tests and no R test
  directory. The expanded suite now covers accession-specific genomic-resource,
  hit/bp conservation, unresolved mass, retry-output guards, corrected
  microscopy, path-ranking completeness, 200-tree/leave-one-out coverage,
  simulation design, release manifests, and executed-notebook integrity. The
  corrected R scripts are exercised end to end; a dedicated R unit-test layer
  remains absent.
- **Impact:** The suite can remain green while the biological unit, denominator, and inferential result are wrong.
- **Effort/Risk/Confidence:** M / low / high.
- **Fix sketch:** Add fixtures/gates for accession-taxid-BioSample-voucher-tip concordance, exactly one declared input unit, retained-mass accounting, join cardinality, hit/bp conservation, robust ectopic summaries, finite CICc/weights, global-fit gates, and notebook release-manifest hashes. Add R tests and make `verify_setup.py` check every package used by documented entrypoints, not only `ape`.

## Published-method benchmark: what the repository still needs

The detailed evidence review and full source list are in [`literature-methods-benchmark.md`](literature-methods-benchmark.md). The most important comparisons are:

| Domain | Published-method expectation | Current repository gap | Required response |
|---|---|---|---|
| Salamander cell/nucleus/genome biology | Separate organismal replication from cell replication; make the C-value assay independently interpretable | Mostly 1–2 specimens; genome and nucleus derived from linked images | Hierarchical estimates and independent/validated genome assay |
| Feulgen image cytometry | Controlled fixation/hydrolysis/stain/acquisition, linear response, internal or co-stained standard, random/defined nucleus sampling | Missing batch/acquisition record; one-image reference; quality/size/darkness-selected nuclei | Reconstruct protocol, validate linearity and standards, declare selection estimand |
| Phylogenetic path analysis | Small a priori DAG family, d-separation/global fit, all basis tests and edges, adequacy/power analysis | Data-adaptive families, rejected family still ranked, undefined strict results, no model recovery | Freeze confirmatory DAGs, fail closed, simulate exact design |
| Low-coverage TE profiling | Preserve sampling depth, repeat fraction, classification completeness, library/tool provenance | Closed relative composition, missing absolute denominator/unresolved mass, missing upstream manifest | Rebuild absolute and relative layers with mass accounting |
| Repeat annotation/landscapes | Hit-level annotation identity and bp conservation; documented custom library | Contig class overwrites hits; landscape parser loses identity | Rebuild with conservation and concordance tests |
| Paired-LTR history | Report observed divergence and ascertainment; calibrate age cautiously | Mostly Gypsy, sparse pairs, borrowed rates, no callable-opportunity sensitivity | K2P primary, age windows and assembly/pair-support sensitivity |
| Ectopic deletion | Structural solo:intact evidence or a validated proxy | Unvalidated depth ratio, omitted zeros/mapping provenance, mean dominated by outliers | Robust estimator plus mapping QC and structural validation |
| Reproducible publication | Immutable inputs, environment, provenance, figures/tables from one release | Ignored data/results, mutable sibling paths, unpinned environment | Versioned evidence bundle and manifest-driven notebook |

## Recommended final analysis design

### 1. State the estimands before rebuilding

The manuscript should distinguish five quantities that the current language sometimes blends:

1. **Typical erythrocyte cell and nucleus size for a species** — ideally a specimen-level population estimand.
2. **Quality-screened upper-tail erythrocyte size** — a legitimate but different estimand if the biological question is maximum/near-maximum mature-cell size.
3. **Absolute haploid/diploid genome size** — requires an independently calibrated assay and explicit C-value convention.
4. **Absolute repeat abundance** — fraction of sampled/genomic bases attributable to repeats.
5. **Relative TE composition** — balance among TE categories conditional on repeat-assigned mass.

Do not ask one variable to stand in for another. In particular, relative LTR/LINE balance is not total TE load, RepeatMasker `%del` is not a deletion rate, terminal:internal depth is not yet an ectopic-recombination rate, and image IOD is not automatically an absolute C-value.

### 2. Use the actual sampling hierarchy

Represent `species -> specimen -> slide/image -> cell/nucleus` explicitly. Cell-level distributions are valuable for technical QC and within-specimen biology, but the species-level comparative unit is the specimen-derived species estimate. Record sex, life stage, collection/locality, fixation/stain batch, image/acquisition batch, and review state when available. Avoid generic “95% CI” labels when only conditional cell-level uncertainty was estimable.

### 3. Split confirmatory and exploratory analyses

- **Confirmatory:** a small frozen DAG family centered on independently supported genome size, nucleus size, and cell size, plus at most one prespecified TE absolute-load or ILR contrast.
- **Exploratory:** alternative TE summaries, LTR K2P, ectopic depth proxy, organismal covariates, landscapes, clade comparisons, and estimator sensitivities.
- Put all DAGs and endpoint families in a dated registry before the corrected run. Do not remove prior models; mark them `exploratory_pre_audit`.

### 4. Propagate uncertainty rather than substituting point estimates

Generate species-level draws from the measurement hierarchy, rerun phylogenetic coefficients/path models across those draws and a sample of plausible trees, and summarize the distribution of coefficients, signs, model fit, and rankings. Where the biological variance cannot be estimated, encode that limitation rather than treating 50 cells as independent animals.

### 5. Use compositional methods for TE balance

Analyze absolute repeat fraction separately. For relative composition, use prespecified log-ratios or ILR coordinates, retain unresolved mass, state the zero-replacement rule, and avoid causal arrows among metrics computed from the same closed vector unless a generative model justifies them.

### 6. Require model adequacy before model preference

For each DAG family:

- reject or revise the family if the global d-separation test fails;
- refuse CICc ranking when information criteria/weights are non-finite;
- show all basis-set claims and p-values;
- show standardized edge estimates with bootstrap/measurement/tree intervals;
- disclose lambda boundary warnings and residual diagnostics;
- show species and clade influence;
- demonstrate exact-design recovery/power in simulation;
- state what candidate DAGs are observationally equivalent.

## Release notebook blueprint

The final notebook should be a thin, auditable reader of frozen compact artifacts, not a hidden rerun of multi-gigabyte segmentation or TE pipelines. Heavy workflows should produce manifest-verified release tables first.

1. **Release identity and environment**
   - release ID, git commit, dirty-state declaration, environment lock, input/output SHA-256 table;
   - hard failure on a missing or mismatched file.
2. **Species/accession/voucher/tree crosswalk**
   - one row per biological unit;
   - automated taxid/BioSample/voucher/tip checks and explicit historical-name decisions.
3. **Sampling and missingness dashboard**
   - per species: specimens, images/slides, cells, review state, TE/LTR/ectopic support, genome assay, panel membership.
4. **Segmentation and curation QC**
   - stratified validation metrics, representative successes and failure cases, mask overlays, human/model decision counts.
5. **Morphology estimand comparison**
   - all eligible, random-50 repeats, common quantile, literal largest 50, current quality-ranked 50, manual-only;
   - species rank/effect sensitivity and candidate-pool diagnostics.
6. **Genome-size calibration and sensitivity**
   - reference derivation, 1C/2C convention, batch/slide metadata, standards, saturation/linearity, alternative reference/selection panels, independent validation.
7. **TE provenance and mass accounting**
   - per-species WGS depth, input hashes, custom library, absolute repeat fraction, unresolved mass, relative ILR/log-ratios, hit-class concordance.
8. **LTR-history and ectopic-proxy QC**
   - pair/element counts, K2P distributions, support gates, rate-window ages, assembly opportunity, mean/median/leave-one-out depth sensitivities.
9. **Tree and analysis panels**
   - tree provenance/calibrations, tip matching, ultrametric handling, panel sizes, complete-case consequences.
10. **Descriptive and phylogenetic associations side by side**
    - raw species scatterplots must be labeled descriptive;
    - paired panels show PGLS/path estimates, intervals, residual/influence diagnostics, and sample size.
11. **Confirmatory path analysis**
    - frozen DAG registry, basis sets, global fit, CICc validity, all edges, model recovery, tree/measurement sensitivity.
12. **Figures, tables, and anomaly drill-down**
    - every manuscript panel generated from notebook-visible data;
    - click/filter or static species dossiers for unusual results, with source images/masks and provenance links.

The notebook must use relative paths, never `127.0.0.1`, never glob an uncontrolled results directory, and never delete shared artifacts during execution. Write to a release-specific output directory.

## Recommended figure and statistics hierarchy

### Main figures

1. **Study design and phylogenetic coverage:** dated tree plus aligned heatmap of specimens/images/cells, genome-assay status, TE depth/classified mass, LTR pairs, ectopic elements, and missingness.
2. **Measurement validity:** stratified segmentation examples and errors; species-level morphology distributions with specimen/image layers; explicit comparison of selection estimands.
3. **Genome–nucleus–cell relationship:** species points with measurement intervals, raw descriptive relationships, and matched phylogenetic effect panels. If genome values remain image-IOD sensitivities, label them so in the axis and title.
4. **TE abundance and composition:** absolute repeat fraction if recovered; otherwise clearly label relative composition. Show unresolved mass and use interpretable log-ratios rather than stacked percentages alone.
5. **Confirmatory path result:** only after the gates pass; show a small DAG with standardized coefficients and intervals, sample size, global-fit statistic, and a nearby alternative-model table.
6. **Mechanism/sensitivity:** paired-LTR K2P and robust ectopic-proxy summaries with support counts and sensitivity intervals.

### Supplementary figures/tables

- Full accession/BioSample/voucher/taxid/tree-tip crosswalk.
- Full wet-lab and imaging batch table.
- All segmentation validation strata and failure cases.
- Candidate-yield and selected-quantile plots.
- Every basis-set test, model ranking, coefficient interval, lambda warning, and influence analysis.
- Tree and measurement-error sensitivity.
- TE classified/unclassified mass and RepeatMasker concordance by species/order.
- LTR pair-support, assembly-opportunity, alignment, and rate-window analyses.
- Ectopic element distributions, zero coverage, mapping QC, median/mean/trimmed/leave-one-out results.

### Minimum reporting for every statistical result

- biological and technical sample sizes separately;
- exact panel/species list and missing-data rule;
- estimand, transformation, units, and denominator;
- effect estimate and interval, not only p-value or model weight;
- phylogenetic covariance method and tree source;
- diagnostics/global-fit/basis-set results;
- multiplicity family and adjustment;
- influence and sensitivity result;
- whether the result is confirmatory, exploratory, or measurement sensitivity.

## Stop/go gates and execution order

No current analysis needs to be deleted. Freeze it with a manifest and use it as a comparison branch. Then proceed in this order:

| Gate | Required evidence | If it fails |
|---|---|---|
| 1. Biological identity | accession–BioSample–taxid–voucher–name–tree-tip decision recorded; genomic and microscopy specimen streams explicitly separate; current reidentified, excluded untrusted, and 2025 validation `fuscus` resources distinguished | Stop all integrated reruns |
| 2. TE integrity | hit-level RepeatMasker class; mass-conserving dnaPipeTE tables; one declared input unit/species | Keep current TE results only as pre-audit artifacts |
| 3. Morphology estimand | completed review; cross-species segmentation audit; frozen estimand and sensitivity comparison | Downgrade to exploratory morphology summaries |
| 4. Genome assay | sourced C-value convention, full protocol/batch record, validated comparability or independent assay | Call image IOD relative sensitivity only; no genome→nucleus causal claim |
| 5. Phylogenetic design | frozen DAGs, finite criteria, adequate global fit, exact-design model recovery, measurement/tree/influence sensitivity | Report associations only, not a causal path |
| 6. Mechanism layers | K2P support/ascertainment and validated robust ectopic proxy | Put mechanism results in exploratory/supplementary section |
| 7. Release | immutable manifest, hashes, locked environment, clean-clone notebook execution | Do not distribute the notebook as reproducible analysis |

The first implementation tranche should therefore be genomic-resource provenance and TE correctness, not figure redesign. The second should resolve measurement estimands/validation. The third should rebuild statistics and the release notebook. This sequencing avoids polishing results that are known to depend on unencoded resource provenance or incorrect classification.

## What remains scientifically usable now

- The raw microscopy images, linked masks, pair-level provenance, and extensive curation tooling are a strong foundation for a defensible morphology dataset.
- The current selected-cell table is useful as a **quality-ranked upper-tail sensitivity panel** once named accurately; it should not be thrown away.
- The read-based TE layer can support relative-composition questions after the taxon, input-unit, unresolved-mass, and denominator fixes.
- Paired-LTR K2P is a useful observed-history measure for detectable intact/paired, mostly Gypsy elements after species identity and support/ascertainment checks.
- The terminal:internal depth metric is useful as an exploratory mapping/structure proxy while validation is built.
- The dated tree is a useful working tree after source/calibration provenance, tolerance handling, and uncertainty are recorded.
- Existing path models are valuable as a discovery record and as simulation targets; they should be labeled pre-audit/exploratory rather than deleted.

## Verification performed

- Ran `scripts/run_tests.sh`: 16 Python tests passed; no R tests were discovered.
- Ran `scripts/run_in_dusky.sh python verify_setup.py`: the current check passed, but it verifies only the R package `ape`.
- Probed the live Dusky R environment: `ape`, `phylopath`, and `phylolm` are present; `caper`, `vegan`, and `testthat` are absent. Documented PGLS/PERMANOVA/test entrypoints therefore do not all execute in the current environment.
- Executed the current generated notebook end to end: all 20 code cells have outputs and no stored execution error. This establishes local execution only, not portability or scientific validity.
- Recomputed selection, review-state, candidate-pool, estimator-shift, and morphology-sensitivity metrics with [`measurement_selection_metrics.py`](measurement_selection_metrics.py).
- Inspected current model ranking, edge, warning, panel, support, and source-registry artifacts; recomputed the key model-fit and sample-size statements above.
- Inspected Cellpose, YOLO, and pair-classifier training/validation manifests in the sibling microscopy workspace.
- Reconciled all 34 accession pairs against current NCBI taxon/BioSample records: 33 matched the declared taxon; the current `fuscus` pair was the sole public-label mismatch. Direct expert guidance identifies that pair as actual *D. fuscus*. A separate accession trace established that the 2025 chromosome-level *D. fuscus* resource is not used in current repository outputs.
- Vetted the RepeatMasker merge on code, sample rows, and a full chunked current-table mismatch scan.
- Compared repository methods against the primary literature summarized in [`literature-methods-benchmark.md`](literature-methods-benchmark.md).

## Rejected or narrowed concerns

These were investigated and should not be repeated as stronger claims than the evidence supports:

- There is no general SRA-versus-assembly specimen mismatch: all 34 listed pairs share the same BioSample. For the current `fuscus` genomic pair, public labeling and expert identification differ, but the SRA and assembly still form one internally matched genomic resource. This does not make the independently collected microscopy specimens genomically matched samples.
- Current morphology inputs do not mix preservation states: all 43 contributing images are brightfield, dry, dried-blood preparations. The unresolved problem is between-slide/stain/acquisition comparability, not mixed preservation in the active panel.
- The bootstrap is not flat by design; it is hierarchical where replication exists. Its limitation is collapse to cell resampling when a species has only one specimen/image.
- OME metadata are not wholly absent: pixel size and green-channel identity are available. Instrument, exposure/illumination, staining, and batch metadata are the missing critical fields.
- Simple assembly length or scaffold N50 did not significantly predict current LTR-pair or ectopic-element counts in a screening analysis. Assembly ascertainment still needs callable-opportunity analysis, but the current evidence does not justify claiming those two coarse metrics already explain the results.
- The balanced genome-QC panel is not proven superior to the primary top-50 panel; it is currently valuable because its very large shifts expose sensitivity, while its broken-cell admissions show why it should not silently replace the primary panel.

## Decision point

The audit supports three self-contained implementation plans, in priority order:

1. **Provenance and TE correctness:** accession-specific *fuscus* provenance, current-versus-2025 sensitivity, RepeatMasker hit classification, analysis-panel dnaPipeTE mass accounting, retry-output isolation, and conservation tests are implemented. Annotation-level validation against the 2025 *fuscus* resource and recovery of absolute dnaPipeTE denominators remain open.
2. **Measurement validity:** microscopy validation set, completed row review, morphology estimand comparison, hierarchical uncertainty, genome-calibration protocol and independent validation plan.
3. **Confirmatory inference and release:** frozen DAGs, compositional features, robust ectopic estimator, tree/measurement uncertainty, exact-design model recovery, figures, manifest, and portable notebook.

Those should be planned and executed as separate gated tranches so that every downstream rerun consumes corrected, frozen inputs.

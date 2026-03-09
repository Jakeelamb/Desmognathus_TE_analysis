# TODO

This file records the highest-value analyses to run while the cell size,
nucleus size, and final genome size layers are still being generated.

The original audit queue is now complete.

Current state:

- No active TE-side audit items remain.
- The repo is in a waiting state for finalized cell size, nucleus size, and
  genome size estimates.
- The next major work block should begin only after those estimates are ready
  for integration into the full causal chain.

The goal is no longer more TE-side pipeline churn. The next goal is to merge
the finalized size estimates into the now-audited TE framework and rerun the
inside-out versus outside-in comparative tests.

## Wait State

- Status: waiting for finalized size estimates
- Blocking inputs:
  - final genome size estimates
  - final nucleus size estimates
  - final cell size estimates
- Next trigger:
  once those estimates are ready, reopen this file as a new integration queue
  for full path-model reruns, updated figure selection, and manuscript-facing
  synthesis.

## Working Rule

- Every new analysis must write:
  - a machine-readable result table under `results/data/` or `path_analysis/data/derived/`
  - a short tracked interpretation note in markdown
  - explicit source paths for every upstream dataset used
- Prefer analyses that strengthen trust in the current TE/genome story over
  broad exploratory side projects.
- Keep primary claims tied to observed/source-backed data. Sensitivity layers
  must stay clearly labeled as such.

## Completed Queue

### 1. Assembly-Bias Audit

- Status: completed 2026-03-08
- Question:
  Are the main TE signals tracking real biology, or are they mostly tracking
  assembly quality / fragmentation / low-coverage artifacts?
- Why this matters:
  This is the single most important remaining trust check because many genomes
  were assembled from low-coverage data.
- Existing inputs:
  - `results/data/ltr_age/genome_assembly_manifest.csv`
  - local FASTA assemblies under `input_data/genomes/`
  - `path_analysis/data/derived/te_path_features.csv`
  - `path_analysis/data/derived/path_input_master.csv`
  - `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
  - `results/data/ltr_age/ltr_age_species_summary.csv`
- Build:
  - per-assembly contig count, total assembled bp, N50, longest contig, gap fraction if applicable
  - correlations between assembly metrics and:
    - `ltr_line_logratio`
    - `order_pielou`
    - `weighted_te_divergence_p90`
    - `ectopic_log10_mean_ratio`
    - `ltr_history_median_k2p_distance`
    - `ltr_history_n_pairs_estimated`
- Outputs to create:
  - `results/data/assembly_bias/assembly_quality_metrics.csv`
  - `results/data/assembly_bias/assembly_bias_correlations.csv`
  - `results/data/assembly_bias/assembly_bias_models.csv`
  - `ASSEMBLY_BIAS_AUDIT.md`
- Pass condition:
  The main TE features and LTR-history variables are not dominated by assembly-quality covariates.
- Escalation condition:
  If a core TE variable is strongly explained by fragmentation/assembly size, it
  must be downgraded or residualized before manuscript use.
- Current result:
  Most core TE and paired-LTR history variables passed this screen. The main
  exception is `ltr_line_logratio`, which shows a moderate negative association
  with `assembly_top10_bp_fraction` in `ASSEMBLY_BIAS_AUDIT.md`.

### 2. Cross-Method Concordance Audit

- Status: completed 2026-03-08
- Question:
  Do the read-based `dnaPipeTE` summaries and assembly/RepeatMasker summaries
  tell the same species-level TE story?
- Why this matters:
  Agreement across methods would strongly support the biological validity of the
  main TE gradients.
- Existing inputs:
  - `results/data/dnaPipeTE_order_breakdown.csv`
  - `results/data/dnaPipeTE_superfamily_breakdown.csv`
  - `results/data/merged_repeatmasker_data.csv`
  - `results/data/divergence/divergence_summary_statistics_by_species.csv`
  - `path_analysis/data/derived/te_path_features.csv`
- Build:
  - matched species-by-order and species-by-superfamily summaries from both pipelines
  - concordance for:
    - order ranks
    - major balances (`LTR`, `LINE`, `TIR`)
    - evenness/richness summaries where comparable
    - divergence-linked summaries where method overlap exists
- Outputs to create:
  - `results/data/cross_method/cross_method_order_concordance.csv`
  - `results/data/cross_method/cross_method_superfamily_concordance.csv`
  - `results/data/cross_method/cross_method_species_disagreement.csv`
  - `CROSS_METHOD_CONCORDANCE_AUDIT.md`
- Pass condition:
  The same broad species-level TE axes recur across methods.
- Escalation condition:
  If concordance is poor, the main claims need to be narrowed to the more robust
  method-specific features.
- Current result:
  The main order-level TE features show good cross-method agreement in
  `CROSS_METHOD_CONCORDANCE_AUDIT.md`, including strong concordance for `LTR`,
  `DIRS`, `order_pielou`, and moderate concordance for `ltr_line_logratio`.

### 3. TE Age-Spectrum Shape Metrics

- Status: completed 2026-03-08
- Question:
  Are species differences driven by recent bursts, old retained TE tails, or a
  mix of both?
- Why this matters:
  This gets closer to mechanism than simple composition summaries.
- Existing inputs:
  - `results/landscapes/repeat_landscape_*.csv`
  - `results/landscapes/landscape_analysis_summary.csv`
  - `results/data/divergence/divergence_summary_statistics_by_species.csv`
  - `results/data/ltr_age/ltr_age_species_summary.csv`
- Build:
  - order-specific and all-TE shape metrics such as:
    - recent-mass fraction
    - old-tail fraction
    - peak divergence bin
    - skewness / tail weight
    - burst concentration or multimodality proxy
  - connect those summaries to paired-LTR history where possible
- Outputs to create:
  - `results/data/te_age_spectra/te_age_spectrum_metrics.csv`
  - `results/data/te_age_spectra/te_age_spectrum_correlations.csv`
  - `TE_AGE_SPECTRUM_AUDIT.md`
- Pass condition:
  At least one interpretable age-shape axis is stable enough to describe
  species-level TE tempo.
- Current result:
  The local landscapes support a stable recent-vs-old tempo axis in
  `TE_AGE_SPECTRUM_AUDIT.md`. This axis strongly tracks the existing divergence
  summaries, shows a weaker entropy/evenness relationship, and is largely
  independent of the paired-LTR age layer.

### 4. Path-Model Stability Diagnostics

- Status: completed 2026-03-08
- Question:
  Are the main winner families robust to species removal and clade subsampling?
- Why this matters:
  Current path-model sample sizes are informative but still modest.
- Existing inputs:
  - `path_analysis/data/derived/panels/*.csv`
  - `path_analysis/results/*model_ranking.csv`
  - `path_analysis/results/*best_model_edges.csv`
  - `input_data/phylogeny/desmo900dated_test.tre`
- Build:
  - leave-one-species-out reruns for:
    - `te_genome_primary_mediumplus`
    - `te_genome_organismal_primary_mediumplus`
    - `te_genome_ectopic_organismal_primary_mediumplus`
    - `te_genome_ltr_history_primary_mediumplus`
  - clade-jackknife sensitivity where group sizes allow
- Outputs to create:
  - `path_analysis/data/derived/model_stability_summary.csv`
  - `path_analysis/data/derived/model_stability_detail.csv`
  - `path_analysis/MODEL_STABILITY_AUDIT.md`
- Pass condition:
  Winner identity and main sign patterns are not driven by one or two species.
- Escalation condition:
  If a family is highly unstable, it should be downgraded to supplementary.
- Current result:
  `MODEL_STABILITY_AUDIT.md` shows that `te_genome` and
  `te_genome_ectopic_organismal` are highly stable, `te_genome_ltr_history` is
  moderately stable, and `te_genome_organismal` is near-tied and unstable
  enough to keep below the strongest result tier.

### 5. TE-History Mechanism Audit

- Status: completed 2026-03-08
- Question:
  Does the paired-LTR historical axis help explain current TE state, or is it
  largely independent?
- Why this matters:
  This is the cleanest next mechanistic test available before the full
  cell/nucleus/genome-size chain is finalized.
- Existing inputs:
  - `path_analysis/data/derived/ltr_history_features.csv`
  - `path_analysis/data/derived/te_path_features.csv`
  - `results/data/ltr_age/ltr_divergence_feature_correlations.csv`
  - `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
- Build:
  - compact comparative tests among:
    - `ltr_history_median_k2p_distance`
    - `ltr_line_logratio`
    - `order_pielou`
    - `ectopic_log10_mean_ratio`
    - selected divergence/deletion summaries
  - optional small path sensitivity family if the simpler tests justify it
- Outputs to create:
  - `results/data/ltr_history_mechanism/ltr_history_mechanism_tests.csv`
  - `LTR_HISTORY_MECHANISM_AUDIT.md`
- Pass condition:
  The historical axis either yields a coherent mechanistic link to current TE
  state or is clearly documented as mostly orthogonal.
- Current result:
  `LTR_HISTORY_MECHANISM_AUDIT.md` shows that paired-LTR age is mostly
  orthogonal to current TE state, while recoverable paired-LTR abundance/density
  tracks landscape recency and lower old-tail burden.

### 6. Phylogenetic Residual / Outlier Audit

- Status: completed 2026-03-08
- Question:
  Which species are unusually TE-heavy, ectopic-heavy, or old-LTR-rich relative
  to phylogenetic expectation?
- Why this matters:
  Outlier species often produce the sharpest biological insight and help avoid
  overgeneralized clade stories.
- Existing inputs:
  - `input_data/phylogeny/desmo900dated_test.tre`
  - `path_analysis/data/derived/path_input_master.csv`
  - `results/data/ltr_age/ltr_age_species_summary.csv`
  - `results/data/ectopic_recombination_species_summary.csv`
- Build:
  - phylogenetic residuals for the main TE and LTR-history variables
  - ranked outlier tables with source-linked notes on species coverage / support
- Outputs to create:
  - `results/data/phylo_residuals/phylo_residual_summary.csv`
  - `results/data/phylo_residuals/phylo_outlier_species.csv`
  - `PHYLOGENETIC_OUTLIER_AUDIT.md`
- Pass condition:
  The project gains a small set of biologically interpretable species-level
  exceptions worth discussing in the paper.
- Current result:
  `PHYLOGENETIC_OUTLIER_AUDIT.md` identifies species-level terminal deviations
  worth carrying into the paper, including recurrent outlier behavior for
  `aeneus`, strong recentness in `fuscus`, unusual ectopic behavior in
  `intermedius` and `mavrokoilius`, and unusual TE-balance/evenness in
  `catahoula`.

## Execution Plan

### Phase 1: Trust Layer

Run first:

1. Assembly-bias audit
2. Cross-method concordance audit

Reason:

- These are the highest-value validation steps for the main TE story.
- If either fails badly, later comparative interpretation has to be narrowed.

Deliverable for phase completion:

- A short decision note stating which TE variables remain primary-paper safe.

Current phase outcome:

- Broadly passed.
- Keep the TE feature set in play, but treat `ltr_line_logratio` as the most
  technically sensitive core variable and avoid overrelying on it alone.

### Phase 2: Temporal Mechanism Layer

Run second:

3. TE age-spectrum shape metrics
4. TE-history mechanism audit

Current phase outcome:

- Passed.
- The landscape branch supports a stable recent-vs-old TE tempo axis.
- Paired-LTR age itself is mostly independent of that axis.
- Recoverable paired-LTR abundance behaves more like a present-day recency /
  recoverability signal than a deep historical age signal.

Reason:

- These exploit the newly rebuilt landscape and paired-LTR layers.
- They sharpen whether TE history reflects recent bursts, old retention, or both.

Deliverable for phase completion:

- A source-linked statement about TE tempo across species that can sit beside
  the LTR calibration results.

### Phase 3: Model Robustness Layer

Run third:

5. Path-model stability diagnostics
6. Phylogenetic residual / outlier audit

Current phase outcome:

- Passed with useful qualification.
- The strongest primary-family support is now for `te_genome` and
  `te_genome_ectopic_organismal`.
- `te_genome_organismal` is too winner-fragile to overstate.
- `te_genome_ltr_history` is informative but still sensitivity-grade.
- The project now has concrete species-level phylogenetic exceptions to discuss
  instead of only clade-level summaries.

Reason:

- Once the TE variables themselves are validated, the next job is to assess
  robustness and identify biologically meaningful exceptions.

Deliverable for phase completion:

- A revised statement of which path families are strong enough for the main
  text, and which belong in supplement/sensitivity only.

## Implemented Scripts

Implemented in this queue:

1. `scripts/processing/audit_assembly_bias.py`
2. `scripts/processing/audit_cross_method_concordance.py`
3. `scripts/processing/build_te_age_spectrum_metrics.py`
4. `scripts/processing/audit_ltr_history_mechanism.py`
5. `path_analysis/scripts/audit_model_stability.R`
6. `scripts/processing/audit_phylogenetic_outliers.R`

## Next Queue

When the size estimates land, the next queue should include:

1. rebuild the merged size-aware analysis inputs
2. rerun the full TE -> genome size -> nucleus size -> cell size families
3. compare inside-out versus outside-in support with the finalized size layer
4. refresh manuscript tables, figure selection, and methods/results wording

## Stop Conditions

- Stop adding new TE-derived predictors if the trust-layer audits show strong
  assembly or method artifacts.
- Stop adding new path families if the main existing families fail basic
  stability checks.
- Do not promote any new historical TE variable to the primary claim tier
  unless it survives both the trust-layer and stability-layer audits.

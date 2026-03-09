# Candidate Models

These model families are the starting point for development. They are intentionally small because the current comparative sample sizes are modest.

## Modeling Rules

- Use one observed proxy per conceptual block in the first-pass DAGs.
- Start with order-level TE summaries, not full superfamily matrices.
- Keep candidate sets small enough to compare cleanly with CICc.
- Treat morphology families as scaffolded and provisional until independent final genome-size estimates are ready.

## Family 1: `te_genome`

Purpose: test whether broad TE characteristics are associated with genome-size variation across the TE/genome overlap set.

Nodes:

- `gs`: scaled log genome size
- `ltr_balance`: log LTR:LINE balance
- `te_evenness`: order-level Pielou evenness

Candidate models:

1. `genome_null`
   No directed paths among the observed variables.
2. `ltr_balance_only`
   `gs <- ltr_balance`
3. `evenness_only`
   `gs <- te_evenness`
4. `additive_load_evenness`
   `gs <- ltr_balance + te_evenness`
5. `mediated_evenness`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness`

Interpretation goal:

- Does genome size track broad TE-composition balance?
- Is evenness informative beyond a simple LTR-heavy vs LINE-heavy contrast?

## Family 2: `te_genome_ectopic`

Purpose: add the current ectopic-recombination proxy as a possible constraint or mediator.

Nodes:

- `gs`
- `ltr_balance`
- `te_evenness`
- `ectopic_index`

Candidate models:

1. `ectopic_only`
   `gs <- ectopic_index`
2. `ltr_to_ectopic`
   `ectopic_index <- ltr_balance`
   `gs <- ectopic_index`
3. `evenness_and_ectopic`
   `ectopic_index <- ltr_balance`
   `gs <- te_evenness + ectopic_index`
4. `full_mechanism`
   `ectopic_index <- ltr_balance + te_evenness`
   `gs <- ltr_balance + te_evenness + ectopic_index`

Interpretation goal:

- Does an LTR-rich TE landscape associate with the ectopic proxy?
- Does the ectopic proxy explain genome-size variation directly, or only alongside TE composition?

## Family 3: `genome_morphology`

Purpose: stage the classical nucleotypic hypothesis.

Nodes:

- `gs`
- `ns`: scaled log nucleus area
- `cs`: scaled log cell area

Candidate models:

1. `morphology_null`
   No directed paths.
2. `genome_to_nucleus`
   `ns <- gs`
3. `genome_to_cell_direct`
   `cs <- gs`
4. `mediated_cell_size`
   `ns <- gs`
   `cs <- ns`

Interpretation goal:

- Is cell size better explained as a direct correlate of genome size, or as a nucleus-mediated consequence?

Caveat:

- This family should remain a planning or sensitivity analysis until the final genome estimates are independent of nucleus-area scaling.

## Family 4: `te_genome_morphology`

Purpose: stage the full mechanistic chain suggested by the chapter framing.

Nodes:

- `ltr_balance`
- `te_evenness`
- `gs`
- `ns`
- `cs`

Candidate models:

1. `te_to_genome_to_nucleus_to_cell`
   `gs <- ltr_balance + te_evenness`
   `ns <- gs`
   `cs <- ns`
2. `te_to_genome_partial_cell`
   `gs <- ltr_balance + te_evenness`
   `ns <- gs`
   `cs <- gs + ns`
3. `te_direct_to_genome_nucleus`
   `gs <- ltr_balance`
   `ns <- gs`
   `cs <- ns`
4. `te_evenness_path`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness`
   `ns <- gs`
   `cs <- ns`

Interpretation goal:

- This is the eventual integrated chapter analysis, but not the first result to prioritize.

## Family 5: `te_genome_organismal`

Purpose: test whether the core TE-genome association remains after adding the most defensible organismal covariates from the curated trait layer.

Nodes:

- `ltr_balance`
- `te_evenness`
- `body_size`: scaled log adult-oriented body-size proxy
- `aquaticity`: scaled aquaticity index
- `gs`

Design note:

- `development_mode` is intentionally omitted from this first organismal family because the current TE+genome overlap contains only three direct developers and those are largely nested inside the lowest aquaticity class.

Candidate models:

1. `te_baseline`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness`
2. `body_size_additive`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + body_size`
3. `aquaticity_additive`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + aquaticity`
4. `organismal_additive`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + body_size + aquaticity`
5. `aquaticity_confounds_body_and_te`
   `body_size <- aquaticity`
   `te_evenness <- ltr_balance + aquaticity`
   `gs <- te_evenness + body_size`

Interpretation goal:

- Does body size or broad aquaticity explain genome-size variation better than the current TE summary?
- Does the TE signal survive after a compact organismal adjustment?

## Family 6: `te_genome_ectopic_organismal`

Purpose: test whether body size still matters after adding the current ectopic-recombination proxy, or whether the ectopic signal absorbs the same variance.

Nodes:

- `ltr_balance`
- `te_evenness`
- `ectopic_index`
- `body_size`
- `gs`

Design note:

- `aquaticity` is omitted because the first organismal family showed little support for it.
- The key comparison here is between `ectopic_index`, `body_size`, and the previous TE-only path.

Candidate models:

1. `ectopic_baseline`
   `gs <- ectopic_index`
2. `ectopic_body_size_additive`
   `gs <- ectopic_index + body_size`
3. `te_body_size_baseline`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + body_size`
4. `te_ectopic_body_size`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + ectopic_index + body_size`
5. `ltr_to_ectopic_body_size`
   `ectopic_index <- ltr_balance`
   `gs <- ectopic_index + body_size`

Interpretation goal:

- Does ectopic signal remain the best genome-size predictor once body size is included?
- Does the TE-evenness path still help once ectopic and body size are both available?

## Family 7: `te_genome_ltr_history`

Purpose: test whether the species-level historical LTR divergence layer explains
additional genome-size variation beyond the current TE-state variables.

Nodes:

- `ltr_history`: species-level median paired-LTR K2P divergence
- `ltr_balance`
- `te_evenness`
- `gs`

Design note:

- This family uses divergence rather than a hard-dated age because any fixed
  positive substitution-rate calibration is only a scalar rescaling of K2P and
  therefore does not change the standardized path-model term.
- The biological timescale interpretation is carried separately in
  `LTR_SUBSTITUTION_RATE_CALIBRATION.md`.

Candidate models:

1. `te_baseline`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness`
2. `history_direct`
   `gs <- ltr_history`
3. `history_additive`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness + ltr_history`
4. `history_to_balance`
   `ltr_balance <- ltr_history`
   `te_evenness <- ltr_balance`
   `gs <- te_evenness`
5. `history_to_evenness`
   `te_evenness <- ltr_balance + ltr_history`
   `gs <- te_evenness`

Interpretation goal:

- Does retained historical LTR divergence add explanatory value beyond the
  current TE composition summary?
- Is any history signal better understood as upstream of current TE state or as
  an additional direct correlate of genome size?

## What Not To Do Yet

- Do not put `morph_nc_ratio` in the same DAG as both `ns` and `cs`.
- Do not add many raw TE order percentages to the same small-sample DAG.
- Do not mix provisional area-derived genome size with nucleus size and then over-interpret directionality.
- Do not fold the TE dilution simulation into the same path-analysis family. That is a separate appendix/supplement question.

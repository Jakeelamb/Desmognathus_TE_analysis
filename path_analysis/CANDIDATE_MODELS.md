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

## What Not To Do Yet

- Do not put `morph_nc_ratio` in the same DAG as both `ns` and `cs`.
- Do not add many raw TE order percentages to the same small-sample DAG.
- Do not mix provisional area-derived genome size with nucleus size and then over-interpret directionality.
- Do not fold the TE dilution simulation into the same path-analysis family. That is a separate appendix/supplement question.

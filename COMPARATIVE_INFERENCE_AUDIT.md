# Comparative Inference Audit

This note records the current audited state of the legacy PGLS and PERMANOVA
branches used for TE composition comparisons outside the phylogenetic
path-analysis workflow.

## Bottom Line

- The previous PGLS pairwise screen was not methodologically safe as written,
  because it regressed raw TE compositions against each other even though the
  order and superfamily tables are perfectly closed at `100%` per species.
- The PGLS workflow now runs on CLR-transformed compositions with explicit zero
  replacement, which makes it a defensible compositional screening layer rather
  than a raw-proportion correlation scan.
- The previous PERMANOVA branch used an arbitrary `cutree(k = 5)` split that
  produced tiny groups (`38/4/2/1/1`) and overstated pairwise clade contrasts.
- The PERMANOVA workflow now uses curated named clades, collapses undersized
  named groups into `other`, and restricts pairwise testing to named clades
  with adequate sample size.

## PGLS Status

Canonical script:
- `scripts/processing/pgls_analysis.R`

Current output contract:
- `results/data/pgls/pgls_order_pairwise.csv`
- `results/data/pgls/pgls_superfamily_pairwise.csv`

Current method:
- CLR-transform order and superfamily composition matrices before model fitting
- Use half the smallest positive value as the deterministic zero replacement
- Fit pairwise PGLS models as a compositional screening layer

Important consequence of the audit:
- The old-looking `LINE ~ LTR` signal is no longer supported once the workflow is
  run in CLR space. In the current output it is weak and non-significant
  (`p = 0.295`, `p_adjusted = 0.442`, `R^2 = 0.034`).

Current CLR-screen summary:
- Order-level pairwise models passing BH at `0.05`: `13 / 36`
- Major-superfamily pairwise models passing BH at `0.05`: `12 / 15`

Interpretation standard:
- These PGLS results are now acceptable as an exploratory compositional
  screening layer.
- They should not be treated as direct evidence of raw abundance tradeoffs,
  because CLR coordinates are still compositional coordinates rather than
  independent unconstrained traits.

## PERMANOVA Status

Canonical script:
- `scripts/processing/permanova_analysis.R`

Current output contract:
- `results/data/permanova/permanova_summary.csv`
- `results/data/permanova/permanova_order_pairwise.csv`
- `results/data/permanova/permanova_superfamily_pairwise.csv`
- `results/data/permanova/species_clade_assignments.csv`

Current clade scheme:
- curated named clades from the Desmognathus tree/literature scaffold
- undersized named clades collapsed into `other`
- current group sizes in the TE dataset: `other = 25`, `mountain_clade = 7`,
  `coastal_clade = 5`, `southern_appalachian = 5`, `small_clade = 4`

Current PERMANOVA summary:
- Order Bray-Curtis: `R2 = 0.171`, `p = 0.201`
- Order CLR-Euclidean: `R2 = 0.193`, `p = 0.075`
- Superfamily Bray-Curtis: `R2 = 0.240`, `p = 0.037`
- Superfamily CLR-Euclidean: `R2 = 0.170`, `p = 0.154`

Pairwise status after restricting to named clades only:
- Order pairwise contrasts surviving BH at `0.05`: `0 / 6`
- Superfamily pairwise contrasts surviving BH at `0.05`: `0 / 6`

Interpretation standard:
- Order-level clade structure is not robust across distance spaces.
- Superfamily-level clade structure is present in Bray-Curtis space but does
  not replicate in CLR-Euclidean space, so it should be described as
  distance-sensitive rather than fully robust.
- Pairwise clade separation should not currently be claimed as significant.

## Recommended Paper Use

- Main inference: keep relying on the path-analysis layer and the audited TE
  feature engineering rather than the legacy pairwise PGLS screen.
- PGLS: use as exploratory/supplementary compositional screening only.
- PERMANOVA: use the overall superfamily Bray result cautiously and report its
  lack of replication in CLR space.
- Avoid strong statements that TE order composition is separated by clade.

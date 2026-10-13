# Bonsai publication representation and PCA comparison

This directory compares the canonical TE34 CLR-PCA with the Bonsai tree
representation described by de Groot et al. (2026). The comparison and full
validation remain here as the analysis authority; the compact tip annotations,
validation summary, and supported Newick are promoted to S42, S43, and Figure 2
by the controlled publication builder. Canonical PCA outputs are unchanged and
remain available as the supplementary S09-S11 geometry audit.

The comparison uses the same 34 species and 24 ubiquitous-superfamily CLR
features as the final PCA. No feature standard deviations are supplied because
the current data do not provide a defensible species-by-feature error matrix.
Bonsai therefore uses its documented near-zero-error fallback. Feature-variance
rescaling is disabled so that both methods start from the same unscaled
CLR/Aitchison geometry.

## Result

On this dataset, Bonsai is complementary rather than uniformly better. Its
global squared-distance correlation is 0.881, compared with 0.918 for PC1-PC2,
and the median species-specific correlation is 0.875 versus 0.913. Bonsai does
better on five-nearest-neighbor recall (0.724 versus 0.488) and on normalized
scale-adjusted stress (0.278 versus 0.397). PC1-PC6 is better than both on every
reported fidelity metric (global correlation 0.994; neighbor recall 0.824;
stress 0.099). All 23 PCs reproduce the full CLR geometry exactly.

These metrics answer different questions. Correlation measures distance
ordering, neighbor recall emphasizes local membership, and stress measures
absolute reconstruction after one global scale factor. There is no single
"Bonsai wins" result here.

## Reproductive-strategy validation

The strategy-colored pattern is much stronger than randomized labels on the
fixed fitted tree. Across 99,999 label permutations preserving the observed
group sizes, the tree requires only 2 Fitch parsimony steps versus a null mean
of 9.40 (`P = 0.00001`; Holm-adjusted `P = 0.00009` across all nine reported
tests). The exact aquatic-egg group (7 species) has adjusted `P = 0.00009`, and
the exact direct-development group (3 species) has adjusted `P = 0.00312`.

That localization is not the same as stable topology. In 200 perturbation runs
that resample the 24 CLR superfamily coordinates and refit a bounded Bonsai
search, the exact aquatic-egg edge recurs in 55 trees (27.5%) and the exact
direct-development edge in 70 trees (35.0%). The terrestrial-egg + aquatic-larva
strategy occupies the intervening backbone in the full tree rather than one
exact edge side and is exact in 2/200 resamples (1.0%). Median normalized
Robinson-Foulds distance from the full fit is 0.774, so the detailed topology is
feature-sensitive even though the full fitted tree localizes the labels very
strongly.

These are post hoc exploratory results. Fixed-tree label permutation does not
establish causation or phylogenetically independent association, and resampling
correlated compositional coordinates is a feature perturbation—not a
biological-specimen bootstrap.

## Reproduce the upstream run

The exact read-only publication replay is:

```bash
make te-bonsai
```

The wrapper fetches the exact upstream commit, prepares the general-data input,
runs the authors' implementation through `uv`, and writes scratch results below
the ignored `output/` directory:

```bash
analyses/01_transposable_elements/bonsai_comparison/run_upstream_bonsai.sh
```

The full run takes about three minutes on the current laptop. The tracked result
can be replayed without network access or rerunning Bonsai:

```bash
uv run --frozen --no-sync python \
  analyses/01_transposable_elements/bonsai_comparison/compare_bonsai_to_pca.py check

scripts/run_r.sh Rscript \
  analyses/01_transposable_elements/bonsai_comparison/plot_bonsai_comparison.R
```

The declared validation campaign is resumable and uses the same pinned upstream
commit in an isolated environment:

```bash
BONSAI_BOOTSTRAP_JOBS=8 \
  analyses/01_transposable_elements/bonsai_comparison/run_feature_bootstrap.sh

uv run --frozen --no-sync python \
  analyses/01_transposable_elements/bonsai_comparison/validate_strategy_localization.py check

scripts/run_r.sh Rscript \
  analyses/01_transposable_elements/bonsai_comparison/plot_bonsai_validation.R
```

## Interpretation boundary

The Bonsai tree is a visualization of TE-composition similarity, not a species
phylogeny. The direct visual benchmark is PC1-PC2; all 23 released PCs remain an
exact representation of the 24-feature CLR distances. The companion validation
quantifies label localization and bounded feature-resampling sensitivity, but
it does not supply biological replication or a feature-error model. Both
analyses remain exploratory only.

Primary sources:

- de Groot et al. (2026), *Nature Biotechnology*,
  <https://doi.org/10.1038/s41587-026-03220-2>
- Official implementation, pinned in the manifest:
  <https://github.com/dhdegroot/Bonsai-data-representation>

The upstream implementation is CC-BY-NC-4.0 and is not vendored here.

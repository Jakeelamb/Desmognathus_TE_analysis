# Path analysis as the central manuscript question

This revision applies only to this isolated candidate run. The earlier editorial draft, original publication figures, and all existing data remain preserved. Figures 1–5 and their review-note identities remain in place; Figure 6 adds the focal path-analysis result.

## Question

Which causal hypotheses can explain covariation among estimated genome size, nuclear size, and cell size in Desmognathus?

The repeat analyses establish the genomic context, and the pairwise regressions establish the observed associations. The path analysis then asks whether the pattern fits nuclear size as the sole connection between genome size and cell size, cell size as the connection between genome size and nuclear size, or other candidate structures. This comparison belongs in the main Results.

## Main inference from the existing analysis

The primary 24-species phylogenetic path comparison favors the cell-bridge class (CICc 13.383; model weight 0.917; global-fit P = 0.976). The nucleus-bridge class has delta CICc 13.772 and fails its global-fit test (P = 0.000997). Weights describe relative support within the ten scored classes; they are not probabilities that a causal explanation is true. The six fully connected DAGs form an unscored saturated class because they imply no conditional-independence constraint.

The biological discrimination is visible in the conditional associations. Estimated genome size remains associated with cell area when nuclear area is included (standardized coefficient 0.575, stored 95% normal interval 0.280 to 0.870). Its coefficient for nuclear area after cell area is included is near zero (0.008, interval -0.492 to 0.508). The broad latter interval does not establish an exactly zero effect. These coefficients are from S28, whose normal intervals differ from the t intervals used in Figure 5 and Table 2.

The strongest concise conclusion is:

> Phylogenetic path analysis favors models in which cell size connects genome-size and nuclear-size variation and disfavors complete mediation of genome–cell covariation through nuclear size in the primary 24-species panel. The directions of the supported connections remain unresolved.

Three alternative directions share the favored class's conditional-independence prediction: genome to cell to nucleus; nucleus to cell to genome; and cell to both genome and nucleus. They are displayed together as competing causal hypotheses. No single direction is selected or assigned a causal probability.

## How much of this result is stable?

The cell-bridge class ranks first in 246/250 measurement resamples, 24/24 species-removal analyses, and 4/4 evolutionary-model choices. It also ranks first on 200/200 published bootstrap trees on an 18-species subset, but nucleus bridge is competitive and globally admissible in all 200. On the reduced focal and published trees, weights are approximately 0.434 for cell bridge and 0.322 for nucleus bridge (delta CICc 0.593). The species subset therefore weakens discrimination. Near-zero fitted phylogenetic covariance also limits what the repeated tree result establishes.

Simulation recovery is a separate measure of model discrimination: true cell bridge ranks first in 78/100 Brownian-residual and 80/100 independent-residual simulations. Recovery is lower for other generating classes. These rates are available with the figure's supporting validation output; they do not translate into a probability that the empirical causal explanation is correct.

## Manuscript emphasis

Lead the biological motivation with the alternatives linking genome, nucleus, and cell dimensions. Present their measured associations and path comparison together in the main Results. Use the repeat architecture, divergence, diversity, and LTR-ratio results to describe the genomic setting and evaluate the particular repeat hypotheses tested here. The current three-trait path model does not include TE variables or reproductive strategy, so the narrative should not draw an estimated TE-to-cell causal path.

No additional specimens or images are required to display or interpret these existing comparisons. The IOD-derived genome-size estimates, selected-large-cell sampling, and phylogeny provenance belong in the Methods and figure legend. A common multiplicative genome-size anchor cannot change the standardized log-scale path ranking; species-specific measurement errors remain relevant.

Sources: snapshotted S28, S30, S32–S37, and `inputs/path/candidate_dags.csv`; the primary path source is `analyses/03_phylogenetic_path/README.md`. Method context: [phylopath model-comparison documentation](https://ax3man.github.io/phylopath/reference/plot.phylopath_summary.html).

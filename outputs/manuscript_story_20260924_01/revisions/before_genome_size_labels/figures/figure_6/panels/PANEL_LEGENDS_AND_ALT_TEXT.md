# a Weighted causal hypotheses linking genome, cell and nucleus

The cell-bridge class contains three distinct directed acyclic graphs, reproduced from the frozen candidate-DAG enumeration. Estimated genome size -> cell area -> nucleus area, its reverse chain, and cell area -> both traits all imply genome–nucleus conditional independence given cell area. Arrow labels are direction-specific standardized phylogenetic regression coefficients from S28; every child in these diagrams has one parent, so its displayed coefficient comes from the matching two-variable regression. Arrow width is proportional to the absolute coefficient. Traits are log10-transformed and standardized to unit sample standard deviation across 24 species. The genome-first orientation has genome -> cell beta = 0.71666 and cell -> nucleus beta = 0.57930. Reverse and fork coefficients are retrieved separately, rather than copied or averaged; they agree after rounding to two decimals. These orientations are statistically indistinguishable in this three-trait comparison and share one class weight of 0.91674; their common conditional structure is the supported result. The highlighted orientation is a hypothesis, not a uniquely identified causal direction.

Alt text: Three equally compatible causal hypotheses have shared class weight 0.917. Direction-specific genome–cell arrows have coefficient 0.72 and cell–nucleus arrows 0.58.

# b Primary comparison of ten path classes

Ten testable path-equivalence classes are ranked by CICc using 24 species, standardized log10 traits, the supplied dated tree and Pagel-lambda regressions. The cell-bridge class has CICc = 13.3829, Fisher C = 0.04959, one independence claim, global-fit P = 0.97551, and relative weight = 0.91674. Its tested claim is the residual genome–nucleus association after conditioning on cell area (P = 0.97551). All nine other tested classes fail the P >= 0.05 global-fit criterion; the nucleus-bridge class has delta CICc = 13.77171 and global-fit P = 0.0009971. Weights compare only the supplied ten classes; six saturated DAGs impose no testable d-separation claim and are not ranked. All primary component regressions have fitted lambda approximately 1e-7. The class weight measures relative model support and is distinct from the direction-specific path coefficients in panel a.

Alt text: Ten models are ranked by delta CICc; only cell bridge passes global fit, with weight 0.917. Nucleus bridge is 13.77 units behind.

# c Path coefficients and uncertainty intervals

Upper rows show the two component regressions in the genome-first hypothesis: genome -> cell beta = 0.71666 (95% normal interval 0.42523 to 1.00810), and cell -> nucleus beta = 0.57930 (0.23869 to 0.91992). Lower rows show separate conditional regressions: genome -> nucleus adjusted for cell area, beta = 0.00792 (-0.49195 to 0.50780); genome -> cell adjusted for nucleus area, beta = 0.57486 (0.27994 to 0.86977). Circles identify displayed single-parent paths and triangles identify conditional comparisons. All four fits use 24 species, standardized log10 traits and Pagel-lambda phylogenetic regression. Intervals are normal approximations from S28 (beta +/- 1.96 SE), not bootstrap intervals or estimates of indirect effects. The broad conditional genome–nucleus interval does not establish an exactly zero effect.

Alt text: Genome-to-cell and cell-to-nucleus path coefficients are 0.72 and 0.58. Conditional genome-to-nucleus is 0.008 with a wide interval; conditional genome-to-cell remains positive.

# d Path-class ranking across sensitivity analyses

Cell bridge ranks first in 246/250 measurement resamples, 24/24 leave-one-species-out fits, 4/4 evolutionary-model specifications and 200/200 published-tree fits. The first and third checks retain 24 species, species omission retains 23 per fit, and the published-tree analysis contains only 18 species. All 200 published-tree fits also retain nucleus bridge as globally admissible and within delta CICc <= 2, so top ranking does not indicate unique support. On the focal tree pruned to the same 18 species, cell-bridge weight is 0.43357 and nucleus-bridge weight is 0.32240 (delta CICc = 0.59252); the published main tree gives essentially the same values. These reduced-species results primarily expose dependence of discrimination on taxon coverage, rather than establishing strong topology robustness.

Alt text: Cell bridge ranks first in 246 of 250 resamples, all 24 omissions, all four evolutionary models, and all 200 trees. The 18-species tree fits also support nucleus bridge.


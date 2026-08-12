# Curated GBE figure plan

**As of:** 10 August 2026

This plan separates technical readiness from scientific placement. The shared
house style and exported files satisfy the checked GBE/OUP artwork constraints,
but no format can guarantee editorial acceptance.

## Recommended main-text set

| Figure | Evidence panel | Narrative role | Why it belongs in the main text |
|---|---:|---|---|
| Figure 1 — focal phylogeny, study design, and genomic QC | TE34 / Path24 / overlap21 | Shows the supplied 24-tip dated Path24 tree and defines the three real denominators and accession/QC boundary. | Establishes phylogenetic context while preventing the 21-species overlap from being mistaken for a universal study panel; tree calibration and publication provenance remain visibly pending. |
| Figure 2 — TE diversity and CLR ordination | TE34 | Establishes classified order-level Shannon and Gini-Simpson diversity plus major superfamily compositional structure. | Central genomic result; panel A shows the 34-species order-level Shannon and Gini-Simpson distributions without presenting observed richness as a diversity index. S08 contains the primary order stratum and a superfamily sensitivity/audit stratum after classified categories are reclosed to one; unresolved aligned-base mass is excluded from diversity and remains only in S06-S07 dnaPipeTE QC. S09-S11 come from one retry-averaged 24-feature CLR-PCA fit restricted to superfamilies detected in all 34 species. |
| Figure 3 — repeat divergence landscape | TE34 | Shows category-specific divergence profiles and among-species young-repeat fractions. | Central repeat-landscape result; means now zero-fill absent category/bin combinations and weight all 34 species equally. |
| Figure 5 — reviewed cell and nucleus morphology | Path24 | Shows the finalized upper-tail image phenotype and its conditional uncertainty. | Point shape exposes the seven species represented by one animal, so object counts are not mistaken for biological replication. |

This four-figure main set leaves four slots for main-text tables or a
scientifically essential replacement under GBE's eight-item combined cap.

## Recommended supplementary or sensitivity set

| Figure | Evidence panel | Placement | Claim boundary |
|---|---:|---|---|
| Figure 4 — terminal:internal LTR depth | LTR30 | Supplementary sensitivity | Mapping/deletion-footprint proxy only; displayed support counts use the exact 380-element primary branch from 408 LTR/Gypsy candidates with at least five correctly parsed TEsorter domains, retained by both the within-species two-sided 1.5-IQR rule and the ≥80% positive-depth requirement in each terminal and internal region. |
| Figure 6 — relative IOD pairwise PGLS | Path24 | Supplementary sensitivity | Relative image-derived nuclear IOD; point shape exposes three one-image/one-animal estimates and four failed quality-balance gates. |
| Figure 7 — phylogenetic path sensitivity | Path24 | Supplementary audit/sensitivity | Exploratory relative-IOD Markov-equivalence comparison; neither an absolute-genome-size analysis nor a basis for orienting causal arrows. |
| Figure S1 — TE diversity × relative IOD | overlap21 | Supplementary sensitivity | Exact TE34 × Path24 intersection; uses full classified-superfamily Shannon entropy and relative IOD, with no pg/Gb or C-value claim. |

Figure 7 should move into the main text only if the final manuscript treats the
path-model comparison as a core result and keeps its exploratory interpretation.
If so, replace another main display item rather than exceeding the combined cap.

`Path_DAG_equivalence_class_reference.{pdf,png,tif}` is a companion review
asset, not a numbered manuscript figure. It stays outside
`FIGURE_MANIFEST.csv`, summarizes representative DAGs for the 10 testable
Path24 equivalence classes, and leaves the saturated class explicitly
unscored.

## Matching visual system

- Nimbus Sans, a Helvetica-compatible sans serif, embedded in vector PDF.
- Seven-point minimum target text at final dimensions.
- White background, complete black panel boxes, no grids, parenthesized panel
  tags, compact gutters, and a reference-inspired coral/cyan/viridis palette.
- Italic scientific names and a single shared typography/spacing system across
  every panel.
- Shape or line-type redundancy wherever color encodes groups.
- No figure-internal titles; the first legend sentence states the result.
- PDF master, review PNG, and 300-ppi CMYK TIFF derivative for every figure.
- One legend and concise alt text per figure in
  `FIGURE_LEGENDS_AND_ALT_TEXT.md`.

## Submission packaging

- Embed only the selected main figures in the initial manuscript PDF.
- Put supplementary figures and their legends in the separate initial
  supplementary-information PDF.
- At revision, submit each requested figure as its own publication-quality file.
- Recheck the live journal instructions immediately before submission; the
  evidence snapshot is `../journal_requirements/GBE_REQUIREMENTS.md`.

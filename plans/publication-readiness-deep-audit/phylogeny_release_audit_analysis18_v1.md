# Dated-phylogeny release audit: final 18-species panel

This audit preserves both historical tree files. It derives a versioned analysis tree from the exact 46-tip local source and does not infer that genomic and microscopy resources are matched specimens.

## Structural verdict

- The local source has 46 tips, is rooted and fully bifurcating, has positive branch lengths, and has a root age of 30.574569 Myr.
- Its root-to-tip range is 0.000002000 Myr (2.000 years). It is ultrametric at `1e-5` Myr tolerance but fails a strict floating-point check because the Newick branch lengths are rounded.
- The final tree keeps exactly the declared 18 tips, including `fuscus` and excluding `planiceps`. No congener substitution, duplicate-tip collapse, or random polytomy resolution is used.
- The only correction pads terminal edges to the maximum source root-to-tip depth; the maximum change is 2.000 years. The serialized output is strictly ultrametric and retains all internal ages and topology.
- Corrected final-panel root age: 16.708337 Myr.

**Approved:** structural validation, exact final-panel pruning, the accession-specific `fuscus` tip decision, and the rounding-only correction.

## Provenance and uncertainty verdict

**Not approved for confirmatory publication inference yet.** The repository describes the source only as a published time-calibrated *Desmognathus* phylogeny. It does not contain the exact source citation/archive, original tree identifier, calibration record, support values, tree type (posterior draw/consensus/MCC), or posterior/bootstrap tree set. The local source is ignored by Git even though current path scripts read it directly.

The tracked `results/phylogeny/processed_phylogeny.nwk` is not an alternative topology or uncertainty draw. On its 34 common tips it has the same rooted clade signature as the local source, but every pairwise distance is scaled by 0.800000; maximum residual from exact scaling is 1.421e-14. No justification for this 0.8 scale was found, so it is not approved for analysis or release figures.

One corrected point tree is a clean computational input, not propagation of phylogenetic uncertainty. That omission is now addressed by the separately provenance-locked Stewart-Wiens (2025) main tree and 200 time-calibrated bootstrap-tree sensitivity in `phylogeny_published_uncertainty_analysis18_v1.md`. The focal tree still needs its own exact source/calibration record, and neither tree set represents reticulate-network uncertainty.

## Data completeness

- TE composition: 18/18 species.
- Terminal:internal LTR proxy: 16/18 species; unavailable for kanawha, valtos.
- Microscopy morphology: 18/18 species.
- Image-IOD sensitivity stream: 18/18 species. Availability does not make it an approved absolute C-value.

## Release products

- `results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk`
- `results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_patristic_distances_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_tree_input_comparison_v1.csv`
- `results/figures/corrected/phylogeny/phylogeny_data_completeness_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_patristic_heatmap_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_root_to_tip_rounding_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_processed_scale_comparison_v1.png`

## Required external provenance from the tree provider

Request the source paper/DOI, archive URL, exact original filename and checksum, branch-length units, tree-estimation method, calibration identities/priors, support annotations, whether this is a posterior draw/consensus/MCC tree, and the posterior/bootstrap tree sample. Until those arrive, keep all path-model results exploratory with respect to tree uncertainty.

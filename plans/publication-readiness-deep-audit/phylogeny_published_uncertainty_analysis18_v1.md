# Published time-tree uncertainty audit: final 18-species panel

Stewart and Wiens (2025; DOI `10.1016/j.ympev.2024.108272`) provide an optimal time-calibrated tree (Supplementary File S3) and 200 time-calibrated bootstrap trees (Supplementary File S4). The exact publisher archives, extracted files, and SHA-256 hashes are recorded in the source registry. All 18 focal taxa are present by exact binomial name.

## Structural and topology results

- All 201 published trees prune to exactly 18 rooted, bifurcating tips with positive branches. Newick rounding required at most 3.000 years of terminal padding; topology and internal ages were not altered.
- All 200 bootstrap trees have the same rooted final-panel topology as the published main tree (`rooted RF = 0`). They differ from the collaborator focal tree by one rooted bipartition (`rooted RF = 2`).
- Focal-only clade: `kanawha;mavrokoilius`. Published-main/all-bootstrap-only clade: `intermedius;marmoratus;mavrokoilius`. This is a real black-bellied-complex topology sensitivity, not a label substitution.
- Focal-versus-published-main patristic correlation is 0.996, indicating broad relative-distance agreement despite the one topology difference and their different crown ages.

## Time uncertainty

- Final-panel crown age across 200 trees: median 13.859 Myr; 95% interval 12.166-16.403 Myr; range 11.134-17.717 Myr.
- The collaborator focal crown age is 16.708 Myr; the published main crown age is 13.912 Myr.
- Every pairwise patristic distance and every published-main clade age is exported with bootstrap mean, SD, and 95% interval.

## Verdict

**Approved as a published tree-uncertainty sensitivity set.** The repository no longer needs to represent phylogenetic uncertainty with one point tree. The collaborator focal tree remains usable as the focal topology after its exact source/calibration record is supplied; until then, report Stewart-Wiens main plus 200-tree results alongside it. The 200-tree set represents bootstrap/dating sensitivity, not gene-tree conflict or reticulate-network uncertainty.

## Outputs

- `results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk`
- `results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex`
- `results/data/corrected/phylogeny/phylogeny_source_registry_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_pairwise_uncertainty_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_clade_uncertainty_analysis18_v1.csv`
- `results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv`
- `results/figures/corrected/phylogeny/phylogeny_published_bootstrap_density_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_root_age_uncertainty_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png`

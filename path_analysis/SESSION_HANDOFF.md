# Session Handoff

This document is the secondary restart log for the Desmognathus phylogenetic path-analysis work. It captures what this workspace is for, what has already been built, the current dataset state, and the exact next steps.

For the main evergreen analysis map, read `PATH_ANALYSIS_STATUS.md` before this file.

## What We Are Trying To Do

The chapter-level goal is to build a publication-grade comparative framework for genome evolution in `Desmognathus`, centered on:

- TE characteristics and genome size
- genome size and nucleus size
- nucleus size and cell size
- broader organismal/ecological correlates where coverage is good enough to defend them

The immediate path-analysis goal is to have analysis-ready, traceable input tables ready before the final independent genome-size estimates arrive from `~/Projects/cellprofiler_test`.

The main causal backbone remains:

- `TE features -> genome size`
- `ectopic recombination / DNA loss -> genome size`
- `genome size -> nucleus size -> cell size`

Important caveat:

- The current CellProfiler `genome_size_pg` values are still provisional and area-derived.
- That means morphology-linked models are useful for planning and sensitivity analysis, but not yet the final manuscript basis for a clean `genome -> nucleus` claim.

## What Has Been Built

The work is now consolidated in [`path_analysis/`](/home/jake/Projects/Desmognathus_TE/path_analysis).

Core planning and workflow files:

- [`PATH_ANALYSIS_STATUS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/PATH_ANALYSIS_STATUS.md)
- [`README.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/README.md)
- [`INPUT_PREPARATION_PLAN.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/INPUT_PREPARATION_PLAN.md)
- [`TE_MODEL_INPUTS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/TE_MODEL_INPUTS.md)
- [`ANALYSIS_DATASETS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/ANALYSIS_DATASETS.md)
- [`SOURCE_TRACKING.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/SOURCE_TRACKING.md)
- [`CANDIDATE_MODELS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/CANDIDATE_MODELS.md)

Core build scripts:

- [`build_master_dataset.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/build_master_dataset.py)
- [`prepare_te_features.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/prepare_te_features.py)
- [`prepare_nc_biodiversity_traits.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/prepare_nc_biodiversity_traits.py)
- [`prepare_curated_organismal_traits.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/prepare_curated_organismal_traits.py)
- [`build_path_input_master.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/build_path_input_master.py)
- [`build_analysis_panels.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/build_analysis_panels.py)
- [`audit_source_traceability.py`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/audit_source_traceability.py)
- [`path_model_scaffold.R`](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/path_model_scaffold.R)

Traceability and input products:

- [`organismal_traits_curated.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/organismal_traits_curated.csv)
- [`te_model_feature_panel.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/te_model_feature_panel.csv)
- [`path_input_master.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/path_input_master.csv)
- [`analysis_panel_summary.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_panel_summary.csv)
- [`analysis_species_readiness.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_species_readiness.csv)
- [`source_traceability_gaps.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/source_traceability_gaps.csv)

## What Has Been Done So Far

### 1. Workspace and model scaffold

- Built a consolidated `path_analysis/` workspace instead of scattering scripts.
- Added a `phylopath` scaffold with prespecified DAG families.
- Extended the `phylopath` scaffold so it can run directly on the panel CSVs via `--panel`, writing panel-specific result files instead of overwriting the legacy family outputs.
- Ran the initial path-model families and saved outputs to [`results/`](/home/jake/Projects/Desmognathus_TE/path_analysis/results).

Current first-pass ranking winners:

- `te_genome`: `mediated_evenness`
- `te_genome_ectopic`: `ectopic_only`
- `genome_morphology`: `mediated_cell_size`
- `te_genome_morphology`: `te_evenness_path`

Model ranking tables:

- [`te_genome_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_model_ranking.csv)
- [`te_genome_ectopic_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ectopic_model_ranking.csv)
- [`genome_morphology_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/genome_morphology_model_ranking.csv)
- [`te_genome_morphology_model_ranking.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_morphology_model_ranking.csv)

Panel-aware reruns now available in [`results/`](/home/jake/Projects/Desmognathus_TE/path_analysis/results):

- `te_genome_primary_mediumplus`: best model remains `mediated_evenness`
- `te_genome_primary_strict_body`: best model remains `mediated_evenness`
- `te_genome_ectopic_primary_mediumplus`: best model remains `ectopic_only`
- `te_genome_ectopic_primary_strict_body`: best model remains `ectopic_only`
- `te_genome_morphology_primary_mediumplus`: best model remains `te_evenness_path`
- `te_genome_morphology_primary_strict_body`: not interpretable as a ranking result because `n = 9` and CICc is unavailable for models where `q >= n`
- `te_genome_organismal_primary_mediumplus`: best model is now `body_size_additive`
- `te_genome_organismal_primary_strict_body`: best model is also `body_size_additive`
- `te_genome_ectopic_organismal_primary_mediumplus`: best model is `te_body_size_baseline`
- `te_genome_ectopic_organismal_primary_strict_body`: best model is also `te_body_size_baseline`

### 2. TE preparation

- Built a compact TE feature panel for comparative modeling.
- The main usable TE variables are already prepackaged rather than requiring raw compositional tables at model time.

Recommended TE inputs for the main DAGs:

- `ltr_line_logratio`
- `order_pielou`
- one turnover/loss axis at a time
- `ectopic_log10_mean_ratio` where ectopic support exists

Avoid as the main manuscript approach:

- dumping many raw TE-order fractions directly into the path model
- all-order PCA as the only composition summary
- using multiple nearly redundant turnover/deletion summaries in the same small-`n` DAG

### 3. Source tracking

- Built a source manifest and source-file inventory.
- Added source ids, file hashes, and source usage summaries through the derived tables.
- Added a traceability audit so every source id used in analysis inputs can be checked.

Current status:

- [`source_traceability_gaps.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/source_traceability_gaps.csv) is empty.

### 4. Organismal trait curation

- Built a curated organismal trait table that merges NC Biodiversity HTML snapshots, AmphiBIO, manual extractions, and other source-linked inputs.
- Improved the NC parser to catch more adult body-size phrasings and elevation wording.
- Added stronger body-size handling:
  - adult SVL ranges
  - sex-specific adult SVL ranges
  - adult TL range midpoints instead of TL maxima when only TL is available
- Added a species-level panel-readiness audit.

Notable recent upgrades:

- `marmoratus` is now a high-confidence adult SVL range from NC text.
- `folkertsi` is now a high-confidence sex-specific adult SVL range.
- TL-based proxies now use range midpoints where possible instead of maxima.
- `auriculatus` now uses a source-linked transformed-specimen mean SVL from Graham et al. 2010 plus explicit blackwater-swamp lifestyle coding.
- `apalachicolae` now uses a source-linked transformed-specimen mean SVL from Graham et al. 2010 plus explicit seepage/streamside lifestyle coding.
- `valtos` now parses NC `transformed specimens` SVL phrasing as a medium-confidence transformed-body proxy instead of collapsing it into the lowest mixed-stage bucket.
- `monticola` now has a source-linked mean SVL from the 2023 revision abstract, which is enough to move it into medium-plus panels while still keeping it out of strict-body sets.
- New online rescue pass added source-linked adult SVL evidence for `aeneus`, `ocoee`, `orestes`, `organi`, `welteri`, and `wrighti`.
- `welteri` now also has thesis-backed high-confidence lifestyle coding (`stream_aquatic`, aquaticity `2`) from the Felix 2001 Marshall thesis.
- `prepare_curated_organismal_traits.py` now supports `adult_svl_mean_mm` so adult-only morphometric means can be used as medium-confidence body-size inputs when ranges are unavailable.
- `pascagoula` is no longer missing core organismal traits: the Plazi treatment now provides an adult SVL range, aquatic-larval development, and swamp-coastal-plain lifestyle coding.
- `ochrophaeus` now has a source-linked adult SVL range from the Canadian status report, upgrading it from a low-confidence AmphiBIO body-size proxy.

## Current Data State

### Overlap counts

These are the main comparative overlap sizes currently available:

- TE species in repo: `34`
- genome-estimate species from CellProfiler: `31`
- TE + genome overlap: `27`
- TE + genome + ectopic overlap: `24`
- genome + linked morphology overlap: `21`
- TE + genome + linked morphology overlap: `18`
- TE + genome + ectopic + linked morphology overlap: `16`

### Organismal trait status

From [`organismal_traits_curated.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/organismal_traits_curated.csv):

- body-size confidence: `9 high / 26 medium / 2 low / 1 missing`
- missing body size: `1`
- missing development: `1`
- missing lifestyle: `1`
- non-missing elevation midpoints: `10`

Interpretation:

- body size is the main factor limiting stricter organismal panels
- elevation is still too sparse for a main-model covariate
- development and lifestyle are now mostly in good shape

### Current panel sizes

From [`analysis_panel_summary.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_panel_summary.csv):

- `te_genome_all`: `27`
- `te_genome_primary_mediumplus`: `27`
- `te_genome_primary_strict_body`: `16`
- `te_genome_ectopic_all`: `24`
- `te_genome_ectopic_primary_mediumplus`: `24`
- `te_genome_ectopic_primary_strict_body`: `15`
- `te_genome_morphology_all`: `18`
- `te_genome_morphology_primary_mediumplus`: `18`
- `te_genome_morphology_primary_strict_body`: `9`

The key consequence of the latest rescue pass is that there are no longer any low-body-confidence exclusions inside the current `TE + genome` overlap. The medium-plus panels now equal the full overlap panels for `te_genome` and `te_genome_ectopic`, and the strict-body sets gained six species because the new rescues are adult-only SVL metrics rather than TL or mixed-stage proxies.

### Species-level readiness

Use [`analysis_species_readiness.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_species_readiness.csv) to see:

- whether a species is eligible for each panel
- why it is excluded when it is not eligible
- whether exclusion is driven by low-confidence body size, mixed-stage proxies, missing morphology, etc.

This file should be treated as the main species-triage sheet.

## Current Constraints and Rules

### 1. Genome-size caution

Do not over-interpret morphology-linked models yet.

- `TE -> genome size` models are currently the cleanest mechanistic target.
- `genome size -> nucleus size -> cell size` remains scaffolded, but should wait for independent final genome estimates before being used as the formal manuscript result.

### 2. Traceability is mandatory

Every new trait source needs:

- a `source_id`
- a raw file saved locally when possible
- a manifest entry
- a path or locator
- enough metadata to survive manuscript revision and methods writing

If a source cannot be stored locally, it still needs to be entered with a precise URL, citation, and access date.

### 3. Small-`n` comparative discipline

The path models should stay small and prespecified.

- Prefer one observed variable per conceptual block.
- Keep TE predictors compact.
- Keep elevation out of the core model unless coverage improves a lot.
- Keep `N/C ratio` out of DAGs that already include cell area and nucleus area.

### 4. Phylogenetic backfilling is sensitivity-only

The repo now has a separate phylogenetic-imputation layer built from:

- [`results/phylogeny/processed_phylogeny.nwk`](/home/jake/Projects/Desmognathus_TE/results/phylogeny/processed_phylogeny.nwk)
- [`organismal_traits_curated.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/organismal_traits_curated.csv)

Outputs:

- [`phylogenetic_trait_imputation_long.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/phylogenetic_trait_imputation_long.csv)
- [`organismal_traits_phylo_inference.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/organismal_traits_phylo_inference.csv)
- [`phylogenetic_trait_imputation_summary.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/phylogenetic_trait_imputation_summary.csv)

Rules:

- do not overwrite `organismal_traits_curated.csv` with phylogenetic predictions
- keep inferred values in explicit `phylo_` columns only
- use phylogenetic fills for sensitivity analysis or manual-review triage, not as silent replacements for observed facts

Current practical limit:

- the processed tree currently excludes `bairdi`, `brimleyorum`, `folkertsi`, `gvnigeusgwotli`, `imitator`, and `ochrophaeus`
- that means the present phylogenetic backfill only touches `abditus` and `catahoula`
- `abditus` now has a sensitivity-only phylogenetic body-size fill plus phylogenetic support for its low-confidence lifestyle coding
- `catahoula` now has sensitivity-only phylogenetic fills for `development_mode`, `aquaticity_index`, and `microhabitat_class`
- [`phylofill_panel_comparison.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/phylofill_panel_comparison.csv) currently shows no species-set change for `te_genome`, `te_genome_ectopic`, or `te_genome_morphology`: all `primary_phylofill` panels are identical to the observed-only `primary_mediumplus` panels
- the same no-change result now also holds for `te_genome_organismal` and `te_genome_ectopic_organismal`

## Tomorrow: Best Next Steps

### Priority 1. Interpret the new panel-aware path-model reruns

The panel-aware reruns are now on disk, so the next work is interpretation rather than more infrastructure.

Immediate next analyses:

- compare coefficient tables and averaged-path outputs between `all`, `primary_mediumplus`, and `primary_strict_body`
- compare the new organismal-augmented family against the TE-only baseline, especially whether body size is a robust additive covariate or mostly a proxy for the same variance
- compare the new ectopic-plus-organismal family against both `te_genome_ectopic` and `te_genome_organismal`, because the current winner suggests ectopic loses priority once body size is allowed to compete
- decide whether the strict-body morphology family should be dropped from sensitivity reporting because `n = 9` is too small for stable CICc ranking
- decide whether the unchanged winning models across panel definitions are strong enough to summarize in the chapter draft

Reason:

- the highest-value question is now model stability across panel definitions, not further species rescue inside the current overlap
- the TE and ectopic families remained stable in winner identity after the organismal rescue pass
- the first organismal-augmented family is now in place, and it currently points to body size rather than aquaticity as the only organismal covariate that materially competes with the TE-only baseline
- the combined ectopic-plus-organismal family goes further: once body size is allowed in the candidate set, the winning models drop `ectopic_index` entirely and revert to the TE-evenness-plus-body-size structure
- the morphology strict-body subset is visibly underpowered and should not be treated like the larger sensitivity sets

### Priority 2. Treat remaining data gathering as coverage expansion, not rescue within current panels

The current bottlenecks are now mostly outside online organismal backfilling.

Targets:

- missing genome-size coverage for species such as `adatsihi`, `balsameus`, `carolinensis`, and `conanti`
- missing morphology coverage for `aeneus`, `orestes`, `organi`, and `wrighti`
- any future primary supplements that would become useful once new genome or morphology estimates arrive

Important note:

- direct command-line fetch attempts against Dryad file streams returned `403 Forbidden`
- this likely needs browser/manual retrieval or a different authenticated download method tomorrow

### Priority 3. Mine additional organismal covariates only after model reruns

Do this after the model reruns:

- elevation
- reproductive timing / clutch size
- fecundity
- additional habitat/lifestyle variables

These are valuable, but they are now lower leverage than exploiting the newly expanded panel sizes.

### Priority 4. Revisit TE presentation only as a modeling/input question

The TE data are already usable, but tomorrow’s TE-focused work should be about manuscript presentation and model defensibility:

- decide the final “core TE block” for the main model set
- decide which TE variables are main-text versus supplement/sensitivity only
- possibly create a short TE-input methods figure or note from [`TE_MODEL_INPUTS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/TE_MODEL_INPUTS.md)

## Exact Restart Workflow

If restarting fresh later, open these first:

1. [`PATH_ANALYSIS_STATUS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/PATH_ANALYSIS_STATUS.md)
2. [`SESSION_HANDOFF.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/SESSION_HANDOFF.md)
3. [`analysis_species_readiness.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_species_readiness.csv)
4. [`TE_MODEL_INPUTS.md`](/home/jake/Projects/Desmognathus_TE/path_analysis/TE_MODEL_INPUTS.md)

Then rebuild if needed:

```bash
python3 path_analysis/scripts/build_master_dataset.py
python3 path_analysis/scripts/prepare_te_features.py
python3 path_analysis/scripts/prepare_amphibio_traits.py
python3 path_analysis/scripts/prepare_external_morphometrics.py
python3 path_analysis/scripts/prepare_nc_biodiversity_traits.py
python3 path_analysis/scripts/prepare_curated_organismal_traits.py
Rscript path_analysis/scripts/build_phylogenetic_trait_imputation.R
python3 path_analysis/scripts/build_te_model_panel.py
python3 path_analysis/scripts/build_path_input_master.py
python3 path_analysis/scripts/build_analysis_panels.py
python3 path_analysis/scripts/refresh_source_file_inventory.py
python3 path_analysis/scripts/audit_source_traceability.py
```

Then check the two must-pass outputs:

- [`analysis_panel_summary.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_panel_summary.csv)
- [`source_traceability_gaps.csv`](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/source_traceability_gaps.csv)

If touching the path models again:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate Dusky
Rscript path_analysis/scripts/path_model_scaffold.R --dataset te_genome
Rscript path_analysis/scripts/path_model_scaffold.R --dataset te_genome_ectopic
```

## Useful Commits

The main recent checkpoints are:

- `28be836` `Add path analysis workspace and initial phylopath scaffold`
- `1eed61f` `Expand path analysis input pipeline and traceability`
- `a2e154a` `Strengthen organismal trait curation and panel readiness`

## Bottom Line

The project is no longer at the “loose notes and drifting scripts” stage.

What exists now is:

- a consolidated path-analysis workspace
- a runnable comparative model scaffold
- a compact TE feature panel
- curated organismal inputs
- panelized analysis datasets
- species-level readiness/exclusion tracking
- publication-grade source traceability

Tomorrow’s highest-value work is to interpret the new panel-aware `phylopath` outputs, especially whether the unchanged winners across medium-plus and strict-body subsets are robust enough to narrate. There are no obvious remaining online body-size rescues left that would enlarge the current `TE + genome` medium-plus panels; the next data push, if needed, should target missing genome/morphology coverage or richer ecological covariates.

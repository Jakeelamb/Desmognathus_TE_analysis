# Source Tracking

This workspace treats traceability as part of the dataset, not as an afterthought.

## Principle

Every trait value used in the chapter should be reconstructable back to:

- the exact source
- the exact name concept used in that source
- the exact measurement type
- the exact curation decision that made it usable

## Files

- `data/templates/source_manifest.csv`
  Master bibliography and source registry for datasets, papers, supplements, and web sources.
- `data/templates/literature_trait_extraction.csv`
  Long-form raw extraction table. One row per source-derived trait observation.
- `data/templates/species_taxonomy_crosswalk.csv`
  Maps current repo names to legacy names and revision-backed name concepts.
- `data/derived/source_file_registry.csv`
  Hash registry for repo-local, external raw, and CellProfiler input files currently used by the path-analysis pipeline.
- `data/derived/source_usage_summary.csv`
  Table-level audit of every source id currently used in derived path-analysis inputs.
- `data/derived/source_traceability_gaps.csv`
  Machine-generated list of missing manifest entries, missing local file records, or populated values lacking source ids.

## Minimum standard for a usable source

A usable source entry should have:

- `source_id`
- short citation key
- full citation or title
- DOI and/or URL
- source type
- taxonomic scope
- notes on what trait blocks it can inform
- if stored locally, a file hash in `source_file_inventory.csv`

## Minimum standard for a usable extraction row

A usable extraction row should have:

- `current_species`
- `trait_name`
- `value`
- `unit`
- `source_id`
- `original_taxon_name`
- enough notes to recover the context

## Curation guidance

- If a source reports values for a broader species complex, do not coerce it into a current split species without justification.
- If a web source is used for habitat or elevation, keep the access date.
- If multiple sources disagree, keep all raw rows and resolve the conflict in a separate curated summary step.
- If a source is only used to validate a coding decision, it still belongs in `source_manifest.csv`.
- If a raw external file is stored locally, record its filename and sha256 hash.
- Rebuild the audit outputs after major source or table changes with `python3 path_analysis/scripts/audit_source_traceability.py`.

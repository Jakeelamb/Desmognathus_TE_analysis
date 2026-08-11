# Canonical data

`identity/` is the shared join and provenance layer for all three analyses.
Analysis measurements live beside their component under `analyses/*/data`.
There is deliberately no second processed-data hierarchy.

`DATA_CATALOG.csv` records the component, relative path, publication supplement
ID where applicable, dimensions, byte size, and SHA-256 of every canonical
identity, analysis, provenance, and tree file. Refresh it with:

```bash
make catalog
make identity
make validate
```

Join species through the canonical `species` or `species_id` fields and consult
`identity/species_taxonomy_crosswalk.csv`. Exact source aliases live only in
`identity/source_species_aliases.csv`; trait aggregation lives only in
`identity/trait_registry.csv`. Never infer a mapping from a suffix or apply an
accession-level genomic correction to another evidence stream. The full rule
and current cases are in [docs/IDENTITY_RESOLUTION.md](../docs/IDENTITY_RESOLUTION.md).

Panel definitions:

- `identity/te34_panel.csv`: vetted genomic-resource panel.
- `identity/path24_panel.csv`: finalized microscopy/path panel.
- `identity/analysis_availability.csv`: explicit analysis availability across
  the union of study species.

`identity/source_manifest.csv`, `source_species_aliases.csv`,
`trait_registry.csv`, `genomic_accessions.tsv`, and the taxonomy crosswalk
preserve the source and naming decisions needed to interpret those panels.

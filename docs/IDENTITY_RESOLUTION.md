# Species identity and source-label resolution

**Decision status:** active, 10 August 2026

This is the single explanatory authority for species-name resolution in the
working repository. The machine-readable authorities are
`data/identity/species_taxonomy_crosswalk.csv`,
`data/identity/source_species_aliases.csv`, `data/identity/trait_registry.csv`,
`data/identity/source_manifest.csv`, and the exact `tree_tip` values in
`analyses/03_phylogenetic_path/data/path24_traits.csv`.

Inspect the complete contract without writing files:

```bash
make identity
```

`make validate` runs the same audit as part of the repository-wide release
check. A passing data catalog alone is not enough: it proves file identity, not
that aliases, trait aggregation, tree selection, and accessions still agree.
Conversely, `make identity` checks semantics rather than frozen-file hashes;
`make validate` supplies both halves of the contract.

## Author decision and operational rule

The project author clarified on 10 August 2026 that “tips that vary only by
suffix get collapsed.” The repository implements that decision at the
species-data boundary while preserving the raw evidence:

1. Raw source rows and the 46-tip source tree are never renamed or overwritten.
2. A source label maps to a canonical species only through an exact,
   source-scoped row in `data/identity/source_species_aliases.csv`. There is no
   automatic suffix stripping, regular expression, case folding, or inferred
   synonymy.
3. Multiple resolved source rows use the aggregation registered for that trait.
   Maximum SVL uses `max`, so the calculation is a maximum of reported maxima.
4. A phylogenetic tree cannot safely average or merge branches. Path24 therefore
   retains the exact unsuffixed representative tip and prunes alternate tips.
   It does not rename tips or merge their most recent common ancestor (MRCA).
5. Accession-specific expert reidentification remains separate from both trait
   aliases and tree-tip selection.

This is one species-level policy with three evidence-specific operations, not a
universal string-replacement rule.

## Current resolutions

| Evidence stream | Source labels | Active resolution |
| --- | --- | --- |
| Collaborator maximum SVL | `orestes-ac=52`, `orestes-b=50` | Exact aliases to `orestes`; registered `max` gives 52 mm and both labels remain in provenance. |
| Focal tree | `orestes`, `orestes_B` | Retain exact `orestes`; prune `orestes_B`. Their MRCA also contains `ochrophaeus`, so an MRCA merge would be biologically and topologically wrong. |
| Focal tree | `fuscus`, `fuscus_A`, `fuscus_E` | Retain exact `fuscus`; prune the two alternate tips. |
| Focal tree and SVL | `planiceps` | Preserve as a distinct source label; it is outside Path24 and its 48-mm SVL is not transferred to `fuscus`. |
| Genomic resource | public `Desmognathus planiceps` | Map only `SRX20497025 / GCA_032353935.1` to the `fuscus` analysis row, based on Alex Pyron's 9 July 2026 written decision. |
| Source-only labels | `imitator`, `gamma`, `marm_G`, `quad_C`, `quad_G` | Preserve and exclude where the active panel does not select them; infer no mapping from spelling. |

The planiceps-to-fuscus genomic decision is not general synonymy. It does not
rename the collaborator `planiceps` SVL row, the source-tree tip, a morphology
specimen, or the reproductive-strategy figure label.

## Panel consequence

`overlap21` is always the exact set intersection `TE34 ∩ Path24`; an identity
correction does not redefine that panel. Resolving the two *D. orestes* SVL
rows gives *D. orestes* a canonical maximum SVL of 52 mm, so all 21 overlap
species currently have TE composition, reviewed morphology/IOD, and maximum
SVL. A final integrated model is still unselected and unfitted.

## What the audit proves

`make identity` checks all of the following from the preserved sources:

- the 40-row collaborator CSV is the exact ordered `Genus == Desmognathus`
  projection of the 56-row workbook;
- every one of the 37 organismal species has the correct source labels,
  max-of-maxima SVL, source ID, and source path;
- the released Path24 tree is the exact 24-tip source-tree prune by membership
  and rooted topology, with only known Newick rounding tolerance in distances;
- every Path24 species is bound row-by-row to its exact tip, all TE34 tips exist,
  and an unregistered suffix-like active tip forces a new curator decision;
- the exact `orestes` and `fuscus` tips are retained and their alternates absent;
- `planiceps` never becomes a global alias for `fuscus`; and
- the active fuscus SRA/assembly pair and Alex Pyron decision agree across the
  source manifest, TE34, Path24, and the availability table.

The command emits no report file or derived dataset. The stdout report is a
human inspection surface; the CSVs and source files remain the authorities.

## Changing the current identity contract

1. Preserve the original source unchanged and add its provenance to
   `data/identity/source_manifest.csv` if it is new.
2. Record the curator or expert decision as a source-manifest row.
3. For another collaborator-SVL alias, add one exact
   `(collaborator_desmog_svl, source_label)` row to
   `data/identity/source_species_aliases.csv`; never add a suffix pattern.
4. Preserve every contributing SVL source label; `max_svl_mm` continues to use
   the registered `max` aggregation.
5. For a tree, select the exact representative tip explicitly; do not put tree
   labels in the alias table.
6. A new source type or a different trait aggregation requires a deliberate
   extension of `scripts/identity.py` and a failing regression test first; the
   current module does not claim to be a universal resolver.
7. Run `make identity`, `make catalog`, and `make test`. A failure requires a
   deliberate resolution decision, not a silent fallback.

The exact publication/archive identifier, calibration method, and branch-length
provenance for the collaborator tree remain unresolved publication inputs.

# Blueprint: split phylogeny and data-provenance review notebooks

> Planned against commit `6895c2d` on 2026-07-09. The working tree was already
> heavily modified and contains the corrected audit products described below.
> This blueprint is intentionally non-destructive: the notebooks must load
> frozen outputs and run only lightweight checks and presentation plots. They
> must not invoke dnaPipeTE, RepeatMasker, RepeatModeler, TEsorter, read mapping,
> segmentation inference, or the full path-model fit.

## Outcome

Create two independently executable, presentation-ready notebooks:

1. `notebooks/01_phylogeny_trimming_provenance_uncertainty_analysis18_v1.ipynb`
2. `notebooks/02_analysis_data_tables_provenance_analysis18_v1.ipynb`

Both notebooks should work from the repository root or from `notebooks/`, use
only relative paths after locating the root, render their tables and figures
inline, contain no hidden state from the existing omnibus notebook, and retain
all outputs after execution. Each notebook should begin with a prominent
statement that it reads frozen products rather than rerunning the expensive
upstream pipeline.

The two notebooks should share the same final-panel declaration:

```python
FINAL18 = [
    "amphileucus", "anicetus", "apalachicolae", "auriculatus", "bairdi",
    "campi", "fuscus", "gvnigeusgwotli", "intermedius", "kanawha",
    "marmoratus", "mavrokoilius", "monticola", "ocoee", "perlapsus",
    "tilleyi", "valtos", "welteri",
]
```

Do not use `str.removeprefix`: the project Dusky environment is Python 3.8.

## Shared notebook setup

Use the following root-location and display convention in both notebooks:

```python
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import Image, Markdown, display

ROOT = Path.cwd().resolve()
if ROOT.name == "notebooks":
    ROOT = ROOT.parent
assert (ROOT / "DATA_MANIFEST.yml").exists(), "Open from repo root or notebooks/"

pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 200)
sns.set_theme(style="white", context="notebook")

def load_csv(relative_path, **kwargs):
    path = ROOT / relative_path
    assert path.exists(), f"Missing frozen artifact: {path}"
    return pd.read_csv(path, **kwargs)

def show_png(relative_path, width=1100):
    path = ROOT / relative_path
    assert path.exists(), f"Missing frozen figure: {path}"
    display(Image(filename=str(path), width=width))

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
```

Presentation conventions:

- Italicize *Desmognathus* and use `D. species` for display labels, while
  preserving lowercase epithets as join keys.
- Use a restrained colorblind-safe palette. Use one color for the focal tree,
  one for the published main tree, and translucent teal for the 200-tree set.
- Put units on every time axis (`Myr`) and identify whether a value is a crown
  age, root-to-tip distance, or pairwise patristic distance.
- Use direct annotations for the one disputed clade and the *fuscus* resource
  exception. Do not encode claim status only by color.
- Every displayed frozen figure should be preceded by the source table path and
  followed by one sentence stating the permitted interpretation.
- End each notebook with a compact `PASS / LIMITATION / BLOCKED` table rather
  than allowing caveats to live only in prose.

---

# Notebook 01: phylogeny trimming, provenance, and uncertainty

## Scientific purpose

Show exactly how the 46-tip local tree becomes the frozen 18-tip analysis tree,
prove that the operation uses exact taxon identities, display the terminal-edge
rounding correction, and compare the focal tree with the publisher-archived
Stewart-Wiens main tree and 200 dated bootstrap trees. This notebook should
teach the audit trail; it should not merely display a final tree image.

## Exact inputs and frozen outputs

### Inputs used only for lightweight inspection

| Role | Path | Current status |
|---|---|---|
| Collaborator-supplied 46-tip local tree | `input_data/phylogeny/desmo900dated_test.tre` | Local/ignored; SHA-256 `1cf561526405c794b70d5767182eefe84dfa1552f9f9c0be6b36c465b19301aa`; exact citation and calibration provenance unresolved |
| Declared analysis panel | `path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv` | 18 rows; use only its `species` column as the pruning declaration because the table also retains historical image-IOD fields |
| Stewart-Wiens publisher main tree | `input_data/phylogeny/stewart_wiens_2025/Supplementary_File_S3.nex` | Local/ignored publisher supplement; 796 full-tree tips |
| Stewart-Wiens 200 time trees | `input_data/phylogeny/stewart_wiens_2025/Supplementary_File_S4.nex` | Local/ignored publisher supplement; exactly 200 trees |

The notebook should prefer the frozen derivatives below for ordinary display.
Parsing the 46-tip input is needed only for the in-memory pruning demonstration.
It should not rerun either audit writer.

### Frozen tree products

| Role | Path |
|---|---|
| Corrected focal 18-tip tree | `results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk` |
| Published main 18-tip tree | `results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk` |
| Published 200-tree final-panel set | `results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex` |
| Tree structure metrics | `results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv` |
| Tip/accession/data-stream crosswalk | `results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv` |
| Focal patristic-distance matrix | `results/data/corrected/phylogeny/phylogeny_patristic_distances_analysis18_v1.csv` |
| Historical tree comparison | `results/data/corrected/phylogeny/phylogeny_tree_input_comparison_v1.csv` |
| Published-tree source registry | `results/data/corrected/phylogeny/phylogeny_source_registry_v1.csv` |
| Focal/main/200-tree metrics | `results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv` |
| Pairwise distance uncertainty | `results/data/corrected/phylogeny/phylogeny_pairwise_uncertainty_analysis18_v1.csv` |
| Clade-age uncertainty | `results/data/corrected/phylogeny/phylogeny_clade_uncertainty_analysis18_v1.csv` |
| Focal versus published topology difference | `results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv` |
| Structural manifest | `results/data/corrected/phylogeny/phylogeny_release_audit_analysis18_v1.manifest.json` |
| Published uncertainty manifest | `results/data/corrected/phylogeny/phylogeny_published_uncertainty_analysis18_v1.manifest.json` |

### Reports to link in notebook markdown

- `plans/publication-readiness-deep-audit/phylogeny_release_audit_analysis18_v1.md`
- `plans/publication-readiness-deep-audit/phylogeny_published_uncertainty_analysis18_v1.md`
- `plans/publication-readiness-deep-audit/upstream_methods_primary_literature_benchmark.md`
- Stewart and Wiens (2025), DOI `10.1016/j.ympev.2024.108272`

### Frozen figures

Primary presentation sequence:

1. `results/figures/corrected/phylogeny/phylogeny_data_completeness_analysis18_v1.png`
2. `results/figures/corrected/phylogeny/phylogeny_root_to_tip_rounding_analysis18_v1.png`
3. `results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png`
4. `results/figures/corrected/phylogeny/phylogeny_published_bootstrap_density_analysis18_v1.png`
5. `results/figures/corrected/phylogeny/phylogeny_root_age_uncertainty_analysis18_v1.png`

Supplementary/QC sequence:

- `results/figures/corrected/phylogeny/phylogeny_patristic_heatmap_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png`
- `results/figures/corrected/phylogeny/phylogeny_processed_scale_comparison_v1.png`

## Cell-by-cell presentation flow

### 1. Scope and release boundary

State:

- exact final 18 taxa;
- `fuscus` retained and `planiceps` excluded;
- exact pruning, no congener substitution, no duplicate collapse, no random
  polytomy resolution;
- focal tree structurally approved after a rounding-only terminal correction;
- publisher tree set approved for sensitivity, not for reticulation or
  gene-tree-conflict uncertainty;
- focal citation/calibration provenance remains open.

### 2. Artifact integrity gate

Load both JSON manifests, recursively inspect their `outputs`, and assert that
all 22 listed artifacts exist and match their SHA-256 values. Display a compact
table with `artifact`, `exists`, and `sha256_match` rather than printing hashes
as unstructured text.

Expected current result: `22 / 22` exist and `22 / 22` hash-match.

### 3. Show the untrimmed and declared panels

Load the local tree with `Bio.Phylo`, show:

- source tip count: 46;
- frozen tip count: 18;
- panel count: 18;
- all panel names are present in the source;
- source contains a `planiceps` tree tip but the frozen tree does not;
- frozen tree contains `fuscus`.

Display the sorted 18-species list. Do not imply that the public label on the
genomic *fuscus* accession changes the tree taxonomy.

### 4. Reproduce pruning in memory and compare to the frozen output

Import only the pure helper functions from the guarded audit module:

```python
from Bio import Phylo
from scripts.processing.audit_phylogeny_release import (
    canonical_species,
    internal_clade_signature,
    pad_terminal_rounding,
    pairwise_distances,
    prune_to_species,
    root_to_tip_depths,
)

panel = load_csv(
    "path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv",
    usecols=["species"],
)
source = Phylo.read(ROOT / "input_data/phylogeny/desmo900dated_test.tre", "newick")
pruned = prune_to_species(source, panel["species"])
in_memory, padding = pad_terminal_rounding(pruned)
frozen = Phylo.read(
    ROOT / "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk",
    "newick",
)
for tip in frozen.get_terminals():
    tip.name = canonical_species(tip.name)

in_memory_distance = pairwise_distances(in_memory).sort_index().sort_index(axis=1)
frozen_distance = pairwise_distances(frozen).sort_index().sort_index(axis=1)

pruning_checks = {
    "exact_tip_set": set(in_memory_distance.index) == set(FINAL18),
    "frozen_tip_set": set(frozen_distance.index) == set(FINAL18),
    "same_rooted_clades": (
        internal_clade_signature(in_memory) == internal_clade_signature(frozen)
    ),
    "max_pairwise_abs_difference_myr": float(
        np.max(np.abs(in_memory_distance - frozen_distance))
    ),
    "max_terminal_padding_years": float(padding["terminal_padding_years"].max()),
    "frozen_root_to_tip_range_myr": float(
        np.ptp(list(root_to_tip_depths(frozen).values()))
    ),
}
display(pd.Series(pruning_checks, name="value").to_frame())
assert pruning_checks["exact_tip_set"]
assert pruning_checks["frozen_tip_set"]
assert pruning_checks["same_rooted_clades"]
assert pruning_checks["max_pairwise_abs_difference_myr"] < 1e-10
assert pruning_checks["max_terminal_padding_years"] <= 3
assert pruning_checks["frozen_root_to_tip_range_myr"] <= 1e-10
```

Expected current values: maximum pairwise difference
`5.329070518200751e-15 Myr`, maximum terminal padding approximately `2 years`,
and zero frozen root-to-tip range at serialized precision.

This cell is the live proof that pruning can be explained and verified without
regenerating any tree-analysis product.

### 5. Tree structure before and after pruning

Display selected columns from `phylogeny_tree_metrics_analysis18_v1.csv`:

- `tree_id`, `n_tips`, `rooted_by_root_degree`, `fully_bifurcating`;
- `all_branch_lengths_positive`, `root_age_myr`;
- `root_to_tip_range_myr`, `strict_ultrametric_at_1e_10_myr`;
- `maximum_terminal_padding_myr`.

Then show the completeness tree and the rounding figure. The narration must
distinguish the source root age (`30.574569 Myr`) from the pruned focal crown age
(`16.708337 Myr`).

### 6. Tip identity, accessions, and independent samples

Display these crosswalk columns:

```text
species, tree_tip, tree_match_type, te_sra_accession,
te_assembly_accession, microscopy_n_specimens, microscopy_n_images,
genome_iod_n_specimens, genome_iod_n_images,
genome_and_microscopy_are_independent_samples, taxon_decision_note
```

Assertions:

```python
crosswalk = load_csv(
    "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv"
)
assert crosswalk.shape[0] == 18
assert set(crosswalk["species"]) == set(FINAL18)
assert crosswalk["tree_match_type"].eq(
    "exact_after_genus_prefix_removal"
).all()
assert crosswalk["genome_and_microscopy_are_independent_samples"].all()
assert crosswalk["te_sra_accession"].is_unique
assert crosswalk["te_assembly_accession"].is_unique
assert "planiceps" not in set(crosswalk["tree_tip"])
assert "fuscus" in set(crosswalk["tree_tip"])
```

The *fuscus* callout must say:

> The current comparative TE resource is `SRX20497025 / GCA_032353935.1`.
> Its public taxon label is *D. planiceps*, but the accession-specific expert
> identity decision assigns this resource to the *fuscus* comparative row.
> This is not a general synonymy and does not substitute a `planiceps` tree tip.
> Microscopy specimens are independent samples.

Do **not** call `GCA_032353935.1` the published chromosome-level *fuscus*
assembly. The separately published PacBio/chromosome-level resource recorded in
the source manifest is `GCA_050004315.1`, and it is currently
`validation_only`, absent from the active comparative outputs.

### 7. Published tree provenance

Display all three rows of `phylogeny_source_registry_v1.csv`. Also display the
original focal-source path and SHA from
`phylogeny_release_audit_analysis18_v1.manifest.json`, because the focal row in
the CSV registry points to the derived 18-tip working tree rather than the
original 46-tip local file.

### 8. Topology comparison

Display the two rows in
`phylogeny_topology_difference_analysis18_v1.csv`, then the cophylogeny figure.
State the exact result:

- rooted RF focal versus published main: 2;
- focal-only clade: `kanawha;mavrokoilius`;
- published-main/all-200-only clade:
  `intermedius;marmoratus;mavrokoilius`;
- focal versus published-main patristic correlation: `0.9957`.

Do not describe RF=2 as two independent biological conflicts; it reflects one
rooted bipartition replacement.

### 9. Time-tree uncertainty

Use `phylogeny_tree_uncertainty_metrics_analysis18_v1.csv` to calculate and
display the bootstrap crown-age median, 2.5th/97.5th percentiles, and range.
Show the density and crown-age figures.

Expected values:

- 200 published time trees;
- median final-panel crown age `13.8587 Myr`;
- 95% interval `12.1661–16.4026 Myr`;
- published main crown age `13.912 Myr`;
- focal crown age `16.708 Myr`;
- all 200 final-panel bootstrap topologies have RF=0 to the published main and
  RF=2 to the focal tree.

### 10. Distance uncertainty and historical 0.8-scale tree

Show the pairwise uncertainty heatmap as supplementary evidence. Then display
the one-row historical input comparison and the processed-scale figure. State
that `results/phylogeny/processed_phylogeny.nwk` is not an uncertainty draw: on
34 common tips it has the same rooted clades and every pairwise distance is
scaled by `0.8`, with no recovered justification. It is not approved for
analysis or release figures.

### 11. Final verdict table

| Component | Verdict | Permitted statement |
|---|---|---|
| Exact final-18 pruning | PASS | Exact taxon pruning; `fuscus` retained, `planiceps` excluded |
| Rounding correction | PASS | Terminal padding of at most about 2 years; topology/internal ages unchanged |
| Published 200-tree sensitivity | PASS | Bootstrap/dating sensitivity across Stewart-Wiens time trees |
| Focal tree provenance | LIMITATION | Focal topology is structurally usable, but its exact citation/archive/calibrations remain unresolved |
| Reticulation/gene-tree conflict | NOT REPRESENTED | The tree set does not model network or locus-conflict uncertainty |
| Historical processed 0.8 tree | BLOCKED | Do not use |

---

# Notebook 02: analysis data tables and provenance

## Scientific purpose

Provide a navigable data contract for the final 18-species project: which
tables are raw/local, historical bridges, corrected analysis products, or
release gates; what one row means in each table; how the genomic accessions map
to a species-level comparative row; and why genomic and microscopy evidence
must not be described as matched specimens.

This notebook should not duplicate the domain-specific graphs from the repeat,
ectopic, microscopy, or path notebooks. Its figures are provenance and
completeness views.

## Exact provenance and contract inputs

| Role | Path | Unit / warning |
|---|---|---|
| Share-facing data contract | `DATA_MANIFEST.yml` | Human-readable source/derived boundary |
| Data acquisition guide | `DATA_README.md` | Explains local large inputs and retry aliases |
| Variable dictionary | `path_analysis/DATA_DICTIONARY.md` | Required semantic reference |
| Active genomic lookup | `input_data/lookup_table.txt` | 34 species/resource rows; local/ignored; not a microscopy map |
| Source/bibliography manifest | `path_analysis/data/templates/source_manifest.csv` | 40 source rows; includes explicit *fuscus* resource decisions |
| Taxonomy crosswalk | `path_analysis/data/templates/species_taxonomy_crosswalk.csv` | 38 taxon concepts; *fuscus* exception is resource-specific |
| Source file/hash registry | `path_analysis/data/derived/source_file_registry.csv` | 55 file rows; currently all `exists=True` |
| Declared final-18 panel | `path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv` | Use as panel declaration/historical bridge, not as the corrected measurement table |
| Final-18 resource crosswalk | `results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv` | Canonical presentation roster |
| Release gate matrix | `results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv` | Nine machine-readable claim boundaries |

## Corrected analysis tables to inventory

| Layer | Path | Current shape | Unit of observation |
|---|---|---:|---|
| dnaPipeTE mass ledger | `results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv` | 18 × 20 | One final-panel species/genomic resource |
| Mass-accounted order composition | `results/data/corrected/dnapipete/dnaPipeTE_order_breakdown_mass_accounted_analysis18_v2.csv` | 18 rows | One species; relative aligned-repeat composition plus unresolved mass |
| TE diversity sensitivity | `results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv` | 72 × 12 | Species × taxonomic level × mass treatment |
| TE order PCA scores | `results/data/corrected/diversity_pca/te_pca_scores_analysis18_v1.csv` | Load for navigation | Species × ordination specification |
| Terminal:internal proxy robustness | `results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv` | 62 × 23 | Species × filtering branch; 16 species in the all-elements branch |
| Microscopy support | `results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv` | 18 × 20 | One species; nested specimen/image/pair support |
| Microscopy estimator sensitivity | `results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv` | 324 rows | Species × six estimators × three morphology metrics |
| Relative IOD sensitivity | `results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv` | Load for navigation | Species-level relative proxy across declared image-QC subsets |
| Corrected path input cube | `results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv` | 324 × 42 | Species × six morphology estimators × three IOD subsets |
| Measurement specifications | `results/data/corrected/path_analysis/corrected_path_measurement_specifications_analysis18_v1.csv` | 18 rows | One morphology-estimator/IOD-subset combination |

## Cell-by-cell presentation flow

### 1. Scope, layer vocabulary, and no-expensive-rerun contract

Define the layers explicitly:

- **local upstream inputs**: large or ignored resources needed only to rebuild;
- **historical bridge**: preserved legacy tables, including historical
  `genome_size_pg`, that are not current claim-ready measurements;
- **corrected derived**: versioned audit products loaded by these notebooks;
- **release gate**: table that limits permitted manuscript language.

State that this notebook runs only file existence, hash, shape, uniqueness,
join, and missingness checks.

### 2. Artifact catalog and current availability

Create a small hard-coded catalog with columns:

```text
layer, artifact, path, expected_rows, unit_of_observation,
claim_status, expensive_to_rebuild
```

Add live columns `exists`, `size_mb`, and `sha256`. Sort missing artifacts to
the top and assert that every presentation-critical corrected product exists.
Do not assert that all local upstream inputs are portable; report that status.

### 3. Final-18 canonical roster

Use `phylogeny_tip_crosswalk_analysis18_v1.csv`, not the 133-column historical
master table, as the presentation roster. Display:

```text
species, tree_tip, te_sra_accession, te_assembly_accession,
has_te_composition, has_terminal_internal_proxy,
has_microscopy_morphology, has_image_iod_sensitivity,
microscopy_n_specimens, microscopy_n_images,
genome_and_microscopy_are_independent_samples, taxon_decision_note
```

Assertions:

- exactly 18 unique species/tips;
- 18 unique SRX accessions and 18 unique assembly accessions;
- exact species set equals `FINAL18`;
- `planiceps`, `folkertsi`, and `ochrophaeus` are absent;
- `fuscus` is present with `SRX20497025 / GCA_032353935.1`;
- every independence flag is true;
- terminal:internal support is present for 16 species and absent exactly for
  `kanawha` and `valtos`.

### 4. Visualize stream completeness

Build a lightweight heatmap directly from the four boolean stream columns in
the roster. Add a narrow right-hand annotation for `microscopy_n_specimens`.
Use direct labels (`TE composition`, `T:I proxy`, `Morphometry`, `Relative IOD
stream`) and a footnote that availability does not imply claim approval.

This is the only new figure required in the data-tables notebook. Keep it
in-memory during execution; the dedicated phylogeny notebook already displays
the frozen tree-plus-completeness figure.

### 5. Genomic resource identity and *fuscus* decision ledger

Filter `source_manifest.csv` to rows where `tree_tip == "fuscus"` or
`analysis_taxon_name == "Desmognathus fuscus"`. Display:

```text
source_id, sra_accessions, assembly_accession, public_taxon_name,
analysis_taxon_name, voucher_or_isolate, analysis_role,
inclusion_decision, decision_basis_source_id, notes
```

The notebook must distinguish three resources:

1. `SRX20497025 / GCA_032353935.1`: current comparative TE/LTR/terminal:internal
   resource, public label *planiceps*, accession-specific expert reassignment
   to the `fuscus` analysis row;
2. older `SRX19953421 / GCA_030265095.1`: excluded as untrusted;
3. published PacBio/chromosome-level `GCA_050004315.1`: reliable validation and
   sensitivity candidate, but absent from current inputs/outputs.

This cell must explicitly correct the wording in the existing omnibus notebook
builder (`scripts/processing/build_publication_audit_notebook.py`), which calls
`GCA_032353935.1` “chromosome-level.” The source manifest identifies
`GCA_050004315.1` as the separate published chromosome-level resource.

### 6. Genomic and microscopy samples are independent

Display a two-column conceptual mapping rather than inventing a specimen join:

| Genomic namespace | Microscopy namespace |
|---|---|
| `te_sra_accession` | independently collected `specimen_id` |
| `te_assembly_accession` | source image / slide / mask IDs |
| one selected genomic resource per species | one or more specimens/images per species |
| species-level comparative join only | never populated from genomic BioSample/voucher |

Use the all-true independence flag as executable evidence. Do not join
`voucher_or_isolate` into a microscopy record and do not claim that a cell came
from the sequenced individual.

### 7. dnaPipeTE retry aliases are computational runs

Link and summarize:

- `plans/publication-readiness-deep-audit/dnapipete_retry_output_audit.md`;
- `path_analysis/TE_PROVENANCE_AUDIT.md`;
- `results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.manifest.json`.

If the local `input_data/dnaPipeTE/` directory exists, build a display-only run
inventory. Resolve filename prefixes against the known SRX values in
`input_data/lookup_table.txt`; do **not** try to strip appended digits with a
generic regular expression, because an appended `2` or `3` is syntactically
indistinguishable from the accession's digits without the lookup boundary.

Current staged example:

- `SRX19952890_reads_per_component_and_annotation`;
- `SRX19952890R2_reads_per_component_and_annotation`;
- both map to *D. orestes* / `SRX19952890`;
- neither belongs to `FINAL18`;
- neither enters the corrected 18-species mass ledger.

Show that the mass ledger has exactly 18 unique `te_sra_accession` values and
no *orestes* row. State the release rule: choose one completed computational
output per canonical SRX; never sum retries or treat them as biological
replicates.

### 8. Table shape, key, and semantic audit

Create a compact live audit table with one row per corrected table and columns:

```text
artifact, rows, columns, species_n, duplicate_key_rows,
missing_primary_measurements, interpretation
```

Use the following key contracts:

- mass ledger: `species` unique, 18 species;
- diversity sensitivity: unique
  `species + te_level + composition_mode`, 72 rows;
- ectopic robustness: unique `species + analysis_branch`, 62 rows; all-elements
  branch has 16 species;
- microscopy support: `species` unique, 18 species;
- corrected path cube: unique
  `species + morphology_estimator + iod_subset`, 324 rows;
- path specifications: unique `morphology_estimator + iod_subset`, 18 rows.

Expected path cube layout is six morphology estimators × three IOD subsets ×
18 species. Assert:

```python
assert not path_input["iod_proxy_is_absolute_genome_size"].any()
assert not path_input["terminal_internal_proxy_is_ectopic_rate"].any()
assert set(path_input.loc[
    path_input["terminal_internal_ratio_median"].isna(), "species"
]) == {"kanawha", "valtos"}
```

### 9. Historical versus corrected measurement table

Display only column names and a five-row semantic comparison, not the full
historical 69-column panel:

| Historical bridge field | Corrected release field | Current interpretation |
|---|---|---|
| `genome_size_pg` | `relative_nuclear_iod_proxy` | Historical pg conversion withdrawn; corrected field is relative sensitivity only |
| `order_simpson` | `gini_simpson` plus `simpson_dominance` | Name ambiguity resolved |
| `ectopic_log10_mean_ratio` | robust median/geometric terminal:internal fields | Deletion-footprint/mapping proxy, not a rate |
| single morphology summary | six declared estimators | Upper-tail estimand sensitivity |
| one point tree | focal + published main + 200 trees | Tree sensitivity propagated |

Warn prominently that
`path_analysis/data/derived/panels/te_genome_primary_mediumplus.csv` and
`master_species_table.csv` retain historical `genome_size_pg` fields. They are
useful for provenance reconstruction, but must not be presented as the current
analysis-ready genome-size measurements.

### 10. Release gates

Display all nine rows of
`publication_release_gate_matrix_analysis18_v1.csv`, including
`permitted_language`. End with assertions that:

- TE diversity and PCA are descriptive;
- time trees are sensitivity inputs;
- cell/nucleus linkage is provenance-approved;
- absolute genome size is blocked;
- terminal:internal values are blocked as a recombination rate;
- path analysis is exploratory only.

### 11. Download/navigation links

Use notebook-relative Markdown links to the exact CSVs, manifests, data
dictionary, and audit reports. For large tables, display a sortable subset and
also provide the path so the reviewer can open the complete CSV. Do not embed
the multi-gigabyte raw inputs in the notebook.

### 12. Final verdict table

| Component | Verdict | Presenter wording |
|---|---|---|
| Final-18 roster and genomic accessions | PASS | Exact, one row/resource per species |
| Species-level genomic/microscopy join | PASS WITH BOUNDARY | Joined by taxon only; samples are independent |
| dnaPipeTE retry identity | PASS | Technical aliases; one completed output selected per SRX |
| Corrected table keys/shapes | PASS | Versioned, explicit units of observation and claim flags |
| Historical `genome_size_pg` | BLOCKED | Retained for reconstruction, not current absolute genome size |
| Portable sharing through Git alone | BLOCKED UNTIL BUNDLED | Critical local inputs and all `results/` products are ignored |

---

# Lightweight verification gates

## Existing evidence checked while writing this blueprint

- `scripts/run_in_dusky.sh python -m unittest scripts.python.tests.test_phylogeny_release_audit`
  completed with `4` tests passing.
- Both phylogeny manifests currently enumerate `22` outputs; all `22` exist and
  all `22` match their stored SHA-256 values.
- Live in-memory prune check currently returns:
  - exact final-18 tips: `True`;
  - same rooted clade signature as frozen output: `True`;
  - maximum pairwise-distance difference:
    `5.329070518200751e-15 Myr`;
  - maximum terminal padding: approximately `2 years`;
  - serialized final-tree root-to-tip range: `0`.

## Required verification after notebook implementation

Execute each notebook independently from a clean kernel:

```bash
scripts/run_in_dusky.sh jupyter nbconvert \
  --to notebook --execute \
  notebooks/01_phylogeny_trimming_provenance_uncertainty_analysis18_v1.ipynb \
  --inplace --ExecutePreprocessor.timeout=180

scripts/run_in_dusky.sh jupyter nbconvert \
  --to notebook --execute \
  notebooks/02_analysis_data_tables_provenance_analysis18_v1.ipynb \
  --inplace --ExecutePreprocessor.timeout=180
```

Then programmatically verify for each notebook:

- no code-cell output has `output_type == "error"`;
- every code cell has an execution count;
- both notebooks contain the exact `FINAL18` set and current analysis version;
- the phylogeny notebook records `22/22` manifest hash matches;
- the data notebook records 18 unique species/SRX/assembly rows, the 324-row
  path cube, and the all-false absolute-genome/rate flags;
- output cells are retained for tomorrow's presentation.

Run the existing focused test and full repository test command after builders
or notebook-audit code is added:

```bash
scripts/run_in_dusky.sh python -m unittest \
  scripts.python.tests.test_phylogeny_release_audit
scripts/run_tests.sh
```

## STOP conditions

Stop and report rather than silently adapting if:

- any manifest-listed artifact is missing or hash-mismatched;
- final-panel membership differs from the declared 18 species;
- an active genomic resource does not match the crosswalk's SRX/GCA pair;
- a notebook requires regenerating an expensive upstream result to render;
- a notebook is about to call the historical `genome_size_pg` an absolute
  genome-size measurement;
- the *fuscus* current comparative resource is described as the separately
  published chromosome-level assembly;
- sharing requires paths outside the repository but no frozen snapshot is
  included in the bundle.

---

# Key blockers before external sharing

1. **The repository is not a portable share bundle.** `input_data/` and
   `results/` are ignored by Git. Both notebooks can execute on the current
   machine, but a collaborator receiving only the Git checkout will not receive
   the local trees, corrected CSVs, or figures. Package the exact required
   notebook inputs/figures/reports with checksums, or publish a versioned data
   archive, before sending it away.

2. **Focal-tree provenance is unresolved.** The 46-tip source has a locked hash
   and a transparent transformation, but the exact source paper/archive/tree
   identifier, calibration record, support values, and tree type remain
   missing. Keep the Stewart-Wiens set as the provenance-complete sensitivity
   anchor.

3. **Current versus chromosome-level *fuscus* resources are easy to conflate.**
   Active comparative outputs use `SRX20497025 / GCA_032353935.1`. The published
   PacBio/chromosome-level `GCA_050004315.1` is registered as a validation-only
   candidate and was not used to generate the current comparative tables.

4. **The focal source-registry row names a derived file.** In
   `phylogeny_source_registry_v1.csv`, the collaborator focal row's
   `local_source_file` is the derived 18-tip tree. The original 46-tip path/hash
   lives in the structural JSON manifest. The notebook must show both so a
   reviewer does not mistake the derivative for the original source.

5. **Historical panel tables still expose withdrawn fields.** The final-panel
   declaration CSV contains `genome_size_pg`, and the historical master table
   contains the same legacy field. Use those files only for panel/provenance
   reconstruction; use the corrected path cube for current measurements and
   enforce its `iod_proxy_is_absolute_genome_size == False` flag.

6. **Source registries include machine-local paths.** The microscopy source
   registry deliberately records absolute sibling-workspace paths. A portable
   presentation bundle needs frozen copies or a documented same-machine demo;
   links alone will not work on the scientist's computer.

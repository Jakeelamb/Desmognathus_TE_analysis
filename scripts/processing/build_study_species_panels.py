#!/usr/bin/env python3
"""Build explicit TE34, cell21, and integrated-path18 panel contracts.

The TE resource panel is not inferred from microscopy completeness. The linked
cell panel is independent of WGS availability. The path panel is their exact
evidence-complete intersection. No trait is imputed to move a species between
panels.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DECLARATION = ROOT / "path_analysis/data/templates/analysis_species_panels.csv"
LOOKUP = ROOT / "input_data/lookup_table.txt"
READINESS = ROOT / "path_analysis/data/derived/analysis_species_readiness.csv"
OUTPUT_DIR = ROOT / "path_analysis/data/derived/panels"
TE_OUTPUT = OUTPUT_DIR / "study_te_resource_panel34_v1.csv"
CELL_OUTPUT = OUTPUT_DIR / "study_cell_linked_panel21_v1.csv"
PATH_OUTPUT = OUTPUT_DIR / "study_integrated_path_panel18_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "study_species_panels_v1.manifest.json"

EXPECTED_CELL_ONLY = {"brimleyorum", "folkertsi", "ochrophaeus"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(values: pd.Series) -> pd.Series:
    return (
        values.astype(str)
        .str.strip()
        .str.replace("Desmognathus ", "", regex=False)
        .str.replace("D. ", "", regex=False)
        .str.replace("D.", "", regex=False)
        .str.lower()
    )


def build_panels() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    declared = pd.read_csv(DECLARATION)
    required = {
        "species",
        "te_resource_panel34_v1",
        "cell_linked_panel21_v1",
        "integrated_path_panel18_v1",
        "decision_basis",
    }
    missing = sorted(required.difference(declared.columns))
    if missing:
        raise ValueError(f"Panel declaration is missing columns: {missing}")
    declared["species"] = canonical_species(declared["species"])
    if declared["species"].duplicated().any():
        raise ValueError("Panel declaration contains duplicate species")
    for column in [
        "te_resource_panel34_v1",
        "cell_linked_panel21_v1",
        "integrated_path_panel18_v1",
    ]:
        if declared[column].dtype != bool:
            declared[column] = declared[column].map({"True": True, "False": False})
        if declared[column].isna().any():
            raise ValueError(f"Panel declaration contains invalid booleans in {column}")

    te_species = set(declared.loc[declared["te_resource_panel34_v1"], "species"])
    cell_species = set(declared.loc[declared["cell_linked_panel21_v1"], "species"])
    path_species = set(declared.loc[declared["integrated_path_panel18_v1"], "species"])
    if len(te_species) != 34 or len(cell_species) != 21 or len(path_species) != 18:
        raise ValueError(
            "Declared panel sizes must be TE34/cell21/path18, observed "
            f"{len(te_species)}/{len(cell_species)}/{len(path_species)}"
        )
    if path_species != te_species.intersection(cell_species):
        raise ValueError("Integrated path18 must equal the TE34/cell21 intersection")
    if cell_species - te_species != EXPECTED_CELL_ONLY:
        raise ValueError(
            f"Unexpected cell-only species: {sorted(cell_species - te_species)}"
        )

    lookup = pd.read_csv(LOOKUP, sep="\t")
    lookup["species"] = canonical_species(lookup["Species"])
    lookup = lookup.loc[lookup["species"].isin(te_species)].copy()
    if len(lookup) != 34 or lookup["species"].duplicated().any():
        raise ValueError("Active lookup does not map TE34 one-to-one")
    lookup = lookup.rename(
        columns={
            "SRA_Accension": "te_sra_accession",
            "Genome_Accension": "te_assembly_accession",
        }
    )[["species", "te_sra_accession", "te_assembly_accession"]]

    readiness = pd.read_csv(READINESS)
    readiness["species"] = canonical_species(readiness["species"])
    readiness_columns = [
        "species",
        "has_tree_tip",
        "has_te",
        "has_genome",
        "has_ectopic",
        "has_morphology",
        "has_ltr_history",
        "has_ltr_high_confidence",
        "ltr_history_n_pairs_high_confidence",
    ]
    readiness = readiness[readiness_columns]
    joined = (
        declared.merge(lookup, on="species", how="left", validate="one_to_one")
        .merge(readiness, on="species", how="left", validate="one_to_one")
        .sort_values("species")
        .reset_index(drop=True)
    )
    te = joined.loc[joined["te_resource_panel34_v1"]].copy()
    cell = joined.loc[joined["cell_linked_panel21_v1"]].copy()
    path = joined.loc[joined["integrated_path_panel18_v1"]].copy()
    if not te[["has_tree_tip", "has_te"]].all().all():
        bad = te.loc[~(te["has_tree_tip"] & te["has_te"]), "species"].tolist()
        raise ValueError(f"TE34 species lack tree or TE evidence: {bad}")
    if not cell[["has_genome", "has_morphology"]].all().all():
        bad = cell.loc[~(cell["has_genome"] & cell["has_morphology"]), "species"].tolist()
        raise ValueError(f"Cell21 species lack current linked microscopy evidence: {bad}")
    if not path[["has_genome", "has_morphology"]].all().all():
        bad = path.loc[~(path["has_genome"] & path["has_morphology"]), "species"].tolist()
        raise ValueError(f"Path18 species lack current microscopy evidence: {bad}")

    te["analysis_role"] = "te_resource_descriptive"
    cell["analysis_role"] = "cell_linked_descriptive"
    path["analysis_role"] = "microscopy_integrated_path"
    manifest = {
        "analysis_id": "study_species_panels_v1",
        "declaration": str(DECLARATION.relative_to(ROOT)),
        "declaration_sha256": sha256(DECLARATION),
        "lookup": str(LOOKUP.relative_to(ROOT)),
        "lookup_sha256": sha256(LOOKUP),
        "readiness": str(READINESS.relative_to(ROOT)),
        "readiness_sha256": sha256(READINESS),
        "te_resource_panel": {
            "path": str(TE_OUTPUT.relative_to(ROOT)),
            "n_species": len(te),
            "species": te["species"].tolist(),
            "n_with_terminal_internal_proxy": int(te["has_ectopic"].sum()),
            "n_with_ltr_history": int(te["has_ltr_history"].sum()),
        },
        "cell_linked_panel": {
            "path": str(CELL_OUTPUT.relative_to(ROOT)),
            "n_species": len(cell),
            "species": cell["species"].tolist(),
        },
        "integrated_path_panel": {
            "path": str(PATH_OUTPUT.relative_to(ROOT)),
            "n_species": len(path),
            "species": path["species"].tolist(),
        },
        "te_only_species_without_current_linked_microscopy": sorted(
            te_species - path_species
        ),
        "cell_only_species_without_active_te_resource": sorted(cell_species - te_species),
        "microscopy_traits_imputed_for_te_only_species": False,
        "analysis18_outputs_preserved": True,
    }
    return te, cell, path, manifest


def main() -> None:
    te, cell, path, manifest = build_panels()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    te.to_csv(TE_OUTPUT, index=False)
    cell.to_csv(CELL_OUTPUT, index=False)
    path.to_csv(PATH_OUTPUT, index=False)
    manifest["te_resource_panel"]["sha256"] = sha256(TE_OUTPUT)
    manifest["cell_linked_panel"]["sha256"] = sha256(CELL_OUTPUT)
    manifest["integrated_path_panel"]["sha256"] = sha256(PATH_OUTPUT)
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {TE_OUTPUT.relative_to(ROOT)} ({len(te)} species)")
    print(f"Wrote {CELL_OUTPUT.relative_to(ROOT)} ({len(cell)} species)")
    print(f"Wrote {PATH_OUTPUT.relative_to(ROOT)} ({len(path)} species)")
    print(f"Wrote {MANIFEST_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

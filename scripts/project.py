#!/usr/bin/env python3
"""Small command surface for the cleaned Desmognathus paper repository."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd

from .identity import audit_identity
from .publication_metadata import (
    SUPPLEMENT_DEPENDENCIES,
    SUPPLEMENT_METADATA,
    column_definition,
)

ROOT = Path(__file__).resolve().parents[1]
PUBLICATION_DIR = ROOT / "Publication" / "datasets"
CATALOG_PATH = ROOT / "data" / "DATA_CATALOG.csv"
ARCHIVE_PATH = ROOT.with_name(f"{ROOT.name}_archive")

# One canonical analysis-ready source for each publication supplement.
SUPPLEMENT_SOURCES = {
    "S01": "data/identity/analysis_availability.csv",
    "S02": "data/identity/te34_panel.csv",
    "S03": "data/identity/path24_panel.csv",
    "S04": "analyses/03_phylogenetic_path/data/path24_traits.csv",
    "S05": "analyses/01_transposable_elements/data/assembly_quality.csv",
    "S06": "analyses/01_transposable_elements/data/dnapipete_input_quality.csv",
    "S07": "analyses/01_transposable_elements/data/dnapipete_mass_accounting.csv",
    "S08": "analyses/01_transposable_elements/data/te_diversity.csv",
    "S09": "analyses/01_transposable_elements/data/te_pca_scores.csv",
    "S10": "analyses/01_transposable_elements/data/te_pca_variance.csv",
    "S11": "analyses/01_transposable_elements/data/te_pca_loadings.csv",
    "S12": "analyses/01_transposable_elements/data/repeatmasker_divergence_landscape.csv",
    "S13": "analyses/01_transposable_elements/data/repeatmasker_species_inventory.csv",
    "S14": "analyses/01_transposable_elements/data/ltr_element_metrics.csv",
    "S15": "analyses/01_transposable_elements/data/ltr_species_robustness.csv",
    "S16": "analyses/01_transposable_elements/data/ltr_resource_coverage.csv",
    "S17": "analyses/01_transposable_elements/data/ltr_excluded_elements.csv",
    "S18": "analyses/02_morphology/data/cell_nucleus_objects.csv",
    "S19": "analyses/02_morphology/data/cell_nucleus_species_estimates.csv",
    "S20": "analyses/02_morphology/data/cell_mask_review_audit.csv",
    "S21": "analyses/02_morphology/data/nuclear_iod_objects.csv",
    "S22": "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
    "S23": "analyses/02_morphology/data/nuclear_iod_by_image.csv",
    "S24": "analyses/02_morphology/data/nuclear_iod_quality_balance.csv",
    "S25": "analyses/02_morphology/data/nuclear_iod_quality_diagnostics.csv",
    "S26": "analyses/03_phylogenetic_path/data/pairwise_pgls.csv",
    "S27": "analyses/03_phylogenetic_path/data/path_standardized_traits.csv",
    "S28": "analyses/03_phylogenetic_path/data/pairwise_phylogenetic_regressions.csv",
    "S29": "analyses/03_phylogenetic_path/data/evolutionary_model_sensitivity.csv",
    "S30": "analyses/03_phylogenetic_path/data/path_model_comparison.csv",
    "S31": "analyses/03_phylogenetic_path/data/path_edge_estimates.csv",
    "S32": "analyses/03_phylogenetic_path/data/path_model_stability.csv",
    "S33": "analyses/03_phylogenetic_path/data/path_tree_sensitivity.csv",
    "S34": "analyses/03_phylogenetic_path/data/path_measurement_sensitivity.csv",
    "S35": "analyses/03_phylogenetic_path/data/path_leave_one_species_out.csv",
    "S36": "analyses/03_phylogenetic_path/data/path_simulation_calibration.csv",
    "S37": "analyses/03_phylogenetic_path/data/path_dsep_tests.csv",
    "S38": "analyses/03_phylogenetic_path/data/phylogeny_edges.csv",
    "S39": "analyses/03_phylogenetic_path/data/phylogeny_nodes.csv",
    "S40": "analyses/03_phylogenetic_path/data/te_relative_iod_overlap.csv",
}

TREE_SOURCES = {
    "trees/desmo900dated_source46.tre": (
        "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre",
        "TRE",
    ),
    "trees/desmognathus_time_tree_path24_v1.nwk": (
        "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk",
        "NWK",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def write_if_changed(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256(source) == sha256(destination):
        return
    shutil.copy2(source, destination)


def build_publication() -> None:
    manifest_path = PUBLICATION_DIR / "DATASET_MANIFEST.csv"
    manifest = read_csv(manifest_path)
    expected_ids = [f"S{number:02d}" for number in range(1, 41)]
    if manifest["supplement_id"].tolist() != expected_ids:
        raise ValueError("Publication manifest must contain ordered S01-S40 rows")

    inventory_rows: list[dict[str, object]] = []
    updated_rows: list[dict[str, object]] = []
    for row in manifest.to_dict(orient="records"):
        supplement_id = str(row["supplement_id"])
        source_relative = SUPPLEMENT_SOURCES[supplement_id]
        source = ROOT / source_relative
        destination = PUBLICATION_DIR / str(row["filename"])
        if not source.exists():
            raise FileNotFoundError(source)
        write_if_changed(source, destination)
        frame = read_csv(source)
        dependencies = list(SUPPLEMENT_DEPENDENCIES.get(supplement_id, ()))
        dependency_hashes = {}
        for dependency_relative in dependencies:
            dependency = ROOT / dependency_relative
            if not dependency.is_file():
                raise FileNotFoundError(dependency)
            dependency_hashes[dependency_relative] = sha256(dependency)
        row.update(SUPPLEMENT_METADATA.get(supplement_id, {}))
        row.update(
            source_path=source_relative,
            source_sha256=sha256(source),
            additional_source_paths=json.dumps(dependencies, separators=(",", ":")),
            additional_source_sha256=json.dumps(
                dependency_hashes,
                sort_keys=True,
                separators=(",", ":"),
            ),
            output_sha256=sha256(destination),
            rows=len(frame),
            columns=len(frame.columns),
            bytes=destination.stat().st_size,
        )
        updated_rows.append(row)
        for column in sorted(frame.columns):
            inventory_rows.append(
                {
                    "supplement_id": supplement_id,
                    "filename": row["filename"],
                    "column": column,
                    "pandas_dtype": str(frame[column].dtype),
                    "source_column": column,
                    "definition_status": column_definition(supplement_id, column),
                }
            )

    pd.DataFrame(updated_rows, columns=manifest.columns).to_csv(manifest_path, index=False)
    inventory_path = PUBLICATION_DIR / "COLUMN_INVENTORY.csv"
    pd.DataFrame(inventory_rows).to_csv(inventory_path, index=False)

    non_csv_rows = []
    for destination_relative, (source_relative, file_format) in TREE_SOURCES.items():
        source = ROOT / source_relative
        destination = PUBLICATION_DIR / destination_relative
        write_if_changed(source, destination)
        non_csv_rows.append(
            {
                "filename": destination_relative,
                "format": file_format,
                "reason_not_csv": (
                    "Newick/tree text preserves the source or exact Path24 phylogenetic structure."
                ),
                "sha256": sha256(destination),
                "bytes": destination.stat().st_size,
            }
        )
    (PUBLICATION_DIR / "NON_CSV_FILES.json").write_text(json.dumps(non_csv_rows, indent=2) + "\n")

    manifest = read_csv(manifest_path)
    summary = {
        "publication_dataset_count": len(manifest),
        "csv_files": len(manifest),
        "non_csv_tree_files": len(non_csv_rows),
        "total_csv_rows": int(manifest["rows"].sum()),
        "total_csv_bytes": int(manifest["bytes"].sum()),
        "release_status_counts": dict(Counter(manifest["release_status"]).most_common()),
        "manifest_sha256": sha256(manifest_path),
        "column_inventory_sha256": sha256(inventory_path),
    }
    (PUBLICATION_DIR / "release_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Rebuilt {len(manifest)} publication CSVs and {len(non_csv_rows)} trees")


def component_for(path: Path) -> str:
    relative = path.relative_to(ROOT)
    text = str(relative)
    if text.startswith("data/identity/"):
        return "shared_identity_and_provenance"
    if text.startswith("analyses/01_transposable_elements/"):
        return "transposable_elements"
    if text.startswith("analyses/02_morphology/"):
        return "morphology"
    if text.startswith("analyses/03_phylogenetic_path/"):
        return "phylogenetic_path"
    raise ValueError(relative)


def catalog_rows() -> list[dict[str, object]]:
    reverse_supplement = {path: key for key, path in SUPPLEMENT_SOURCES.items()}
    paths = sorted((ROOT / "data" / "identity").glob("*"))
    paths += sorted((ROOT / "analyses").glob("*/data/*"))
    paths += sorted((ROOT / "analyses").glob("*/provenance/*"))
    paths += sorted((ROOT / "analyses").glob("*/trees/*"))
    rows = []
    for path in paths:
        if not path.is_file():
            continue
        relative = str(path.relative_to(ROOT))
        n_rows: int | str = ""
        n_columns: int | str = ""
        if path.name.endswith(".csv") or path.name.endswith(".csv.gz"):
            frame = read_csv(path)
            n_rows, n_columns = len(frame), len(frame.columns)
        elif path.suffix == ".tsv":
            frame = pd.read_csv(path, sep="\t")
            n_rows, n_columns = len(frame), len(frame.columns)
        rows.append(
            {
                "component": component_for(path),
                "path": relative,
                "publication_id": reverse_supplement.get(relative, ""),
                "rows": n_rows,
                "columns": n_columns,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return rows


def build_catalog() -> None:
    rows = catalog_rows()
    pd.DataFrame(rows).to_csv(CATALOG_PATH, index=False)
    counts = Counter(row["component"] for row in rows)
    print(f"Cataloged {len(rows)} files: {dict(sorted(counts.items()))}")


def show_identity() -> None:
    audit = audit_identity(ROOT)
    print(audit.render())
    if audit.errors:
        raise SystemExit(1)


def validate() -> None:
    errors: list[str] = []
    manifest = read_csv(PUBLICATION_DIR / "DATASET_MANIFEST.csv")
    if manifest["supplement_id"].tolist() != list(SUPPLEMENT_SOURCES):
        errors.append("publication manifest IDs are not exact ordered S01-S40")

    for row in manifest.to_dict(orient="records"):
        supplement_id = str(row["supplement_id"])
        source_relative = SUPPLEMENT_SOURCES[supplement_id]
        source = ROOT / source_relative
        publication = PUBLICATION_DIR / str(row["filename"])
        if not source.exists() or not publication.exists():
            errors.append(f"{supplement_id}: missing source or publication file")
            continue
        source_hash = sha256(source)
        publication_hash = sha256(publication)
        if source_hash != publication_hash:
            errors.append(f"{supplement_id}: canonical and publication files differ")
        if source_hash != str(row["source_sha256"]):
            errors.append(f"{supplement_id}: manifest source hash is stale")
        if publication_hash != str(row["output_sha256"]):
            errors.append(f"{supplement_id}: manifest output hash is stale")
        for field, expected_value in SUPPLEMENT_METADATA.get(supplement_id, {}).items():
            if str(row[field]) != expected_value:
                errors.append(f"{supplement_id}: curated {field} is stale")
        frame = read_csv(source)
        if (len(frame), len(frame.columns)) != (int(row["rows"]), int(row["columns"])):
            errors.append(f"{supplement_id}: manifest dimensions are stale")
        try:
            dependency_paths = json.loads(str(row["additional_source_paths"]))
            dependency_hashes = json.loads(str(row["additional_source_sha256"]))
        except json.JSONDecodeError:
            errors.append(f"{supplement_id}: additional-source metadata is not valid JSON")
            continue
        if not isinstance(dependency_paths, list) or not isinstance(dependency_hashes, dict):
            errors.append(f"{supplement_id}: additional-source metadata has the wrong shape")
            continue
        expected_dependencies = list(SUPPLEMENT_DEPENDENCIES.get(supplement_id, ()))
        if dependency_paths != expected_dependencies or set(dependency_hashes) != set(
            expected_dependencies
        ):
            errors.append(f"{supplement_id}: additional-source dependency set is stale")
            continue
        for dependency_relative in expected_dependencies:
            dependency = ROOT / dependency_relative
            if not dependency.is_file():
                errors.append(f"{supplement_id}: missing dependency {dependency_relative}")
            elif dependency_hashes[dependency_relative] != sha256(dependency):
                errors.append(f"{supplement_id}: dependency hash is stale for {dependency_relative}")

    inventory_path = PUBLICATION_DIR / "COLUMN_INVENTORY.csv"
    if not inventory_path.is_file():
        errors.append("Publication/datasets/COLUMN_INVENTORY.csv is missing")
    else:
        inventory = read_csv(inventory_path)
        required_inventory_columns = {
            "supplement_id",
            "filename",
            "column",
            "pandas_dtype",
            "source_column",
            "definition_status",
        }
        if not required_inventory_columns.issubset(inventory.columns):
            errors.append("Publication/datasets/COLUMN_INVENTORY.csv has an invalid schema")
        elif inventory.duplicated(["supplement_id", "column"]).any():
            errors.append("Publication/datasets/COLUMN_INVENTORY.csv has duplicate keys")
        else:
            expected_inventory: dict[tuple[str, str], dict[str, str]] = {}
            manifest_filenames = manifest.set_index("supplement_id")["filename"].astype(str)
            for supplement_id, source_relative in SUPPLEMENT_SOURCES.items():
                frame = read_csv(ROOT / source_relative)
                for column in frame.columns:
                    expected_inventory[(supplement_id, column)] = {
                        "filename": manifest_filenames.loc[supplement_id],
                        "pandas_dtype": str(frame[column].dtype),
                        "source_column": column,
                        "definition_status": column_definition(supplement_id, column),
                    }
            observed_keys = set(
                inventory[["supplement_id", "column"]].itertuples(index=False, name=None)
            )
            if observed_keys != set(expected_inventory):
                errors.append("Publication/datasets/COLUMN_INVENTORY.csv keys are stale")
            else:
                inventory_index = inventory.set_index(["supplement_id", "column"])
                for key, expected_fields in expected_inventory.items():
                    for field, expected_value in expected_fields.items():
                        if str(inventory_index.loc[key, field]) != expected_value:
                            errors.append(
                                f"{key[0]}: column inventory {field} is stale for {key[1]}"
                            )

    te_panel = read_csv(ROOT / "data/identity/te34_panel.csv")
    path_panel = read_csv(ROOT / "data/identity/path24_panel.csv")
    te_species = set(te_panel["species"].astype(str))
    path_species = set(path_panel["species"].astype(str))
    if (len(te_species), len(path_species), len(te_species & path_species)) != (34, 24, 21):
        errors.append("panel contract must be TE34, Path24, overlap21")

    identity_audit = audit_identity(ROOT)
    errors.extend(f"identity: {error}" for error in identity_audit.errors)

    if not CATALOG_PATH.exists():
        errors.append("data/DATA_CATALOG.csv is missing")
    else:
        actual = {str(row["path"]): row for row in catalog_rows()}
        recorded = {
            str(row["path"]): row
            for row in read_csv(CATALOG_PATH).fillna("").to_dict(orient="records")
        }
        if set(actual) != set(recorded):
            errors.append("data catalog paths are stale; run `make catalog`")
        else:
            for path, expected in actual.items():
                observed = recorded[path]
                if str(observed["sha256"]) != str(expected["sha256"]) or int(
                    observed["bytes"]
                ) != int(expected["bytes"]):
                    errors.append(f"data catalog metadata is stale for {path}; run `make catalog`")

    if errors:
        raise SystemExit("Validation failed:\n- " + "\n- ".join(errors))
    print(
        f"PASS: 40 publication tables; TE34={len(te_species)}; "
        f"Path24={len(path_species)}; overlap={len(te_species & path_species)}; identity=PASS"
    )


def status() -> None:
    rows = catalog_rows()
    for component in sorted({str(row["component"]) for row in rows}):
        selected = [row for row in rows if row["component"] == component]
        print(
            f"{component}: {len(selected)} files, "
            f"{sum(int(row['bytes']) for row in selected) / 1024 / 1024:.1f} MiB"
        )
    availability = "available" if ARCHIVE_PATH.is_dir() else "not mounted"
    print(f"archive: {ARCHIVE_PATH} ({availability})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = {
        "catalog": build_catalog,
        "identity": show_identity,
        "publication": build_publication,
        "status": status,
        "validate": validate,
    }
    parser.add_argument("command", choices=tuple(commands))
    args = parser.parse_args()
    commands[args.command]()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge reviewed new-species microscopy objects with immutable legacy freezes."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd


KEY_COLUMNS = ["species", "review_key"]
RELEASE_METADATA_COLUMNS = [
    "decision",
    "decision_source",
    "review_status",
    "review_decision",
    "review_decision_source",
    "review_decision_file",
    "source_panel",
    "finalization_cohort",
]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DERIVED_ROOT = PROJECT_ROOT / "path_analysis/data/external/derived"
MORPH_ROOT = DERIVED_ROOT / "largest_cell_mask_review"
IOD_ROOT = DERIVED_ROOT / "image_quality_matched_genome_iod"
NEW_SPECIES = ("D. aeneus", "D. orestes", "D. wrighti")
IOD_REVIEW_SPECIES = (*NEW_SPECIES, "D. ochrophaeus")

DEFAULT_MORPH_FROZEN = MORPH_ROOT / "frozen_largest_cell_mask_top50.csv.gz"
DEFAULT_MORPH_EXTENSION = MORPH_ROOT / "largest_cell_mask_review_extension_20260714.csv.gz"
DEFAULT_MORPH_DECISIONS = (
    PROJECT_ROOT
    / "input_data/specimen_slides/literal_largest_cell_mask_review_extension_20260714_decisions_aenues_orestes_wrighti.csv"
)
DEFAULT_MORPH_OUTPUT = MORPH_ROOT / "finalized_largest_cell_masks_all_reviewed_species.csv.gz"
DEFAULT_MORPH_AUDIT = MORPH_ROOT / "finalized_largest_cell_masks_new_species_review_audit.csv"

DEFAULT_IOD_FROZEN = IOD_ROOT / "image_quality_matched_nuclei_frozen_reviewed.csv.gz"
DEFAULT_IOD_EXTENSION = IOD_ROOT / "image_quality_matched_review_extension_20260714.csv.gz"
DEFAULT_IOD_DECISIONS = (
    PROJECT_ROOT
    / "input_data/specimen_slides/image_quality_matched_review_extension_20260714_decisions_aeneus_wrighti_orestes_ochrophaeus.csv"
)
DEFAULT_IOD_OUTPUT = IOD_ROOT / "finalized_quality_matched_nuclei_all_reviewed_species.csv.gz"
DEFAULT_IOD_AUDIT = IOD_ROOT / "finalized_quality_matched_nuclei_new_species_review_audit.csv"
DEFAULT_MANIFEST = DERIVED_ROOT / "microscopy_review_finalization_20260714.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        frame.to_csv(
            path,
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
    else:
        frame.to_csv(path, index=False)


def write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def finalize_review_table(
    *,
    frozen: pd.DataFrame,
    extension: pd.DataFrame,
    decisions: pd.DataFrame,
    new_species: Sequence[str],
    panel: str,
    decision_file: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return immutable legacy rows plus reviewed-keep rows for new species."""
    species_set = set(new_species)
    if set(frozen["species"]) & species_set:
        raise ValueError("Legacy frozen table already contains a requested new species.")

    new_decisions = decisions.loc[decisions["species"].isin(species_set)].copy()
    new_decisions["decision"] = new_decisions["decision"].astype(str).str.strip().str.lower()
    unsupported = set(new_decisions["decision"]) - {"keep", "problem", "unsure"}
    if unsupported:
        raise ValueError(
            "Unsupported review decisions: " + ", ".join(sorted(unsupported))
        )
    missing_species = species_set - set(new_decisions["species"])
    if missing_species:
        raise ValueError(
            "Requested new species are missing exported decisions: "
            + ", ".join(sorted(missing_species))
        )
    if new_decisions.duplicated(KEY_COLUMNS).any():
        raise ValueError("New-species decisions contain duplicate review keys.")

    new_extension = extension.loc[extension["species"].isin(species_set)].copy()
    new_extension = new_extension.rename(
        columns={
            column: f"upstream_{column}"
            for column in RELEASE_METADATA_COLUMNS
            if column in new_extension.columns
        }
    )
    if new_extension.duplicated(KEY_COLUMNS).any():
        raise ValueError("New-species extension contains duplicate review keys.")

    source_index = pd.MultiIndex.from_frame(new_extension[KEY_COLUMNS])
    decision_index = pd.MultiIndex.from_frame(new_decisions[KEY_COLUMNS])
    missing = decision_index.difference(source_index)
    if len(missing):
        raise ValueError(f"Decision keys missing from extension table: {len(missing)}")

    audit = new_decisions.merge(
        new_extension[KEY_COLUMNS],
        on=KEY_COLUMNS,
        how="left",
        validate="one_to_one",
        indicator="source_match",
    )
    audit["included_in_finalized"] = audit["decision"].eq("keep")
    audit["decision_file"] = decision_file

    keep_decisions = new_decisions.loc[new_decisions["decision"].eq("keep")].copy()
    selected = new_extension.merge(
        keep_decisions[KEY_COLUMNS + ["decision", "decision_source"]],
        on=KEY_COLUMNS,
        how="inner",
        validate="one_to_one",
    )

    legacy = frozen.copy()
    legacy = legacy.rename(
        columns={
            column: f"upstream_{column}"
            for column in RELEASE_METADATA_COLUMNS
            if column in legacy.columns
        }
    )
    legacy["source_panel"] = legacy.get("panel", "")
    legacy["review_status"] = "reviewed_keep"
    legacy["review_decision"] = "keep"
    legacy["review_decision_source"] = "legacy_frozen"
    legacy["review_decision_file"] = ""
    legacy["finalization_cohort"] = "previously_frozen"

    selected["source_panel"] = selected.get("panel", "")
    selected["review_status"] = "reviewed_keep"
    selected["review_decision"] = selected.pop("decision")
    selected["review_decision_source"] = selected.pop("decision_source")
    selected["review_decision_file"] = decision_file
    selected["finalization_cohort"] = "new_species_reviewed"

    merged = pd.concat([legacy, selected], ignore_index=True, sort=False)
    merged["panel"] = panel
    merged = merged.sort_values(KEY_COLUMNS, kind="mergesort").reset_index(drop=True)
    if merged.duplicated(KEY_COLUMNS).any():
        raise ValueError("Finalized table contains duplicate review keys.")
    return merged, audit.sort_values(KEY_COLUMNS, kind="mergesort").reset_index(drop=True)


def write_finalized_review_release(
    *,
    frozen_path: Path,
    extension_path: Path,
    decisions_path: Path,
    output_path: Path,
    audit_path: Path,
    new_species: Sequence[str],
    panel: str,
) -> dict[str, Any]:
    # Capture source hashes before writing. Incremental finalization may use the
    # previous finalized output as the next frozen input and write back to the
    # same path.
    frozen_sha256 = sha256_file(frozen_path)
    extension_sha256 = sha256_file(extension_path)
    decisions_sha256 = sha256_file(decisions_path)
    frozen = pd.read_csv(frozen_path, low_memory=False)
    extension = pd.read_csv(extension_path, low_memory=False)
    decisions = pd.read_csv(decisions_path, low_memory=False)
    merged, audit = finalize_review_table(
        frozen=frozen,
        extension=extension,
        decisions=decisions,
        new_species=new_species,
        panel=panel,
        decision_file=str(decisions_path.resolve()),
    )
    write_csv(merged, output_path)
    write_csv(audit, audit_path)
    selected = merged.loc[merged["finalization_cohort"].eq("new_species_reviewed")]
    raw_decision_counts = audit.groupby(["species", "decision"]).size().unstack(fill_value=0)
    decision_counts = {
        str(species): {str(decision): int(count) for decision, count in row.items()}
        for species, row in raw_decision_counts.iterrows()
    }
    selected_counts = selected.groupby("species").size().to_dict()
    quality_status = {}
    if "quality_match_status" in selected:
        quality_status = {
            str(species): sorted(set(group["quality_match_status"].dropna().astype(str)))
            for species, group in selected.groupby("species")
        }
    return {
        "panel": panel,
        "new_species": list(new_species),
        "legacy_frozen_rows": int(len(frozen)),
        "legacy_frozen_species": int(frozen["species"].nunique()),
        "new_species_selected_counts": {
            str(species): int(selected_counts.get(species, 0)) for species in new_species
        },
        "new_species_decision_counts": decision_counts,
        "new_species_quality_match_status": quality_status,
        "ignored_legacy_decision_rows": int(len(decisions) - len(audit)),
        "finalized_rows": int(len(merged)),
        "finalized_species": int(merged["species"].nunique()),
        "frozen_source": str(frozen_path.resolve()),
        "frozen_source_sha256": frozen_sha256,
        "extension_source": str(extension_path.resolve()),
        "extension_source_sha256": extension_sha256,
        "decisions_source": str(decisions_path.resolve()),
        "decisions_source_sha256": decisions_sha256,
        "output_csv": str(output_path.resolve()),
        "output_sha256": sha256_file(output_path),
        "audit_csv": str(audit_path.resolve()),
        "audit_sha256": sha256_file(audit_path),
    }


def add_target_summary(summary: dict[str, Any], target: int) -> None:
    selected = summary["new_species_selected_counts"]
    summary["review_target_per_species"] = int(target)
    summary["new_species_shortfall_to_target"] = {
        species: max(0, int(target) - int(count)) for species, count in selected.items()
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--morph-frozen", type=Path, default=DEFAULT_MORPH_FROZEN)
    parser.add_argument("--morph-extension", type=Path, default=DEFAULT_MORPH_EXTENSION)
    parser.add_argument("--morph-decisions", type=Path, default=DEFAULT_MORPH_DECISIONS)
    parser.add_argument("--morph-output", type=Path, default=DEFAULT_MORPH_OUTPUT)
    parser.add_argument("--morph-audit", type=Path, default=DEFAULT_MORPH_AUDIT)
    parser.add_argument("--iod-frozen", type=Path, default=DEFAULT_IOD_FROZEN)
    parser.add_argument("--iod-extension", type=Path, default=DEFAULT_IOD_EXTENSION)
    parser.add_argument("--iod-decisions", type=Path, default=DEFAULT_IOD_DECISIONS)
    parser.add_argument("--iod-output", type=Path, default=DEFAULT_IOD_OUTPUT)
    parser.add_argument("--iod-audit", type=Path, default=DEFAULT_IOD_AUDIT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--morph-species",
        action="append",
        default=None,
        help="Species to append to morphology. Repeat once per species.",
    )
    parser.add_argument(
        "--iod-species",
        action="append",
        default=None,
        help="Species to append to IOD. Repeat once per species.",
    )
    parser.add_argument(
        "--skip-morphology",
        action="store_true",
        help="Leave the existing morphology release and manifest entry unchanged.",
    )
    parser.add_argument(
        "--skip-iod",
        action="store_true",
        help="Leave the existing IOD release and manifest entry unchanged.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.skip_morphology and args.skip_iod:
        raise ValueError("At least one release layer must be finalized.")

    required: list[Path] = []
    if not args.skip_morphology:
        required.extend(
            [args.morph_frozen, args.morph_extension, args.morph_decisions]
        )
    if not args.skip_iod:
        required.extend([args.iod_frozen, args.iod_extension, args.iod_decisions])
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing finalization inputs: " + ", ".join(missing))

    manifest: dict[str, Any] = {
        "release_id": "microscopy_review_finalization_20260714",
        "status": "finalized_reviewed_keeps_only",
        "legacy_frozen_tables_modified": False,
        "new_species": [],
    }
    if args.manifest.is_file():
        manifest.update(json.loads(args.manifest.read_text(encoding="utf-8")))

    released_species = set(manifest.get("new_species", []))
    if not args.skip_morphology:
        morph_species = tuple(args.morph_species or NEW_SPECIES)
        morphology = write_finalized_review_release(
            frozen_path=args.morph_frozen,
            extension_path=args.morph_extension,
            decisions_path=args.morph_decisions,
            output_path=args.morph_output,
            audit_path=args.morph_audit,
            new_species=morph_species,
            panel="finalized_largest_cell_masks_all_reviewed_species",
        )
        add_target_summary(morphology, 50)
        manifest["morphology"] = morphology
        released_species.update(morph_species)
    if not args.skip_iod:
        iod_species = tuple(args.iod_species or IOD_REVIEW_SPECIES)
        iod = write_finalized_review_release(
            frozen_path=args.iod_frozen,
            extension_path=args.iod_extension,
            decisions_path=args.iod_decisions,
            output_path=args.iod_output,
            audit_path=args.iod_audit,
            new_species=iod_species,
            panel="finalized_quality_matched_nuclei_all_reviewed_species",
        )
        add_target_summary(iod, 40)
        manifest["iod"] = iod
        released_species.update(iod_species)

    manifest["new_species"] = sorted(released_species)
    write_json(manifest, args.manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

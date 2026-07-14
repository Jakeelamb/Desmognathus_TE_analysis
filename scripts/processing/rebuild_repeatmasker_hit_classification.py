#!/usr/bin/env python3
"""Build a corrected RepeatMasker table without replacing the pre-audit output."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.processing.te_classification import (  # noqa: E402
    annotate_repeatmasker_hits,
    prepare_dnapipete_context,
)


DEFAULT_REPEATMASKER_INPUT = (
    PROJECT_ROOT / "results" / "data" / "merged_repeatmasker_data.csv"
)
DEFAULT_DNAPIPETE_INPUT = (
    PROJECT_ROOT / "results" / "data" / "dnaPipeTE_merged_classifications.csv"
)
DEFAULT_SPECIES_FILE = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "derived"
    / "panels"
    / "te_genome_primary_mediumplus.csv"
)
DEFAULT_LOOKUP_FILE = PROJECT_ROOT / "input_data" / "lookup_table.txt"
PRE_AUDIT_CANONICAL_OUTPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "repeatmasker_detailed_classification_combined.csv"
)
DNAPIPETE_COLUMNS = [
    "dnaPipeTE_contig_name",
    "Source",
    "Class",
    "Order",
    "Superfamily",
    "#reads",
    "aligned_bases",
    "RM_hit_length_bp",
    "RM_annotation",
    "RM_classification",
    "hitlength_contiglength",
]


def _portable(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def _load_analysis_species(path: Path) -> list:
    table = pd.read_csv(path, dtype=str)
    species_column = "species" if "species" in table.columns else table.columns[0]
    species = sorted(
        {
            _canonical_species(value)
            for value in table[species_column].dropna()
            if _canonical_species(value)
        }
    )
    if not species:
        raise ValueError(f"Analysis species file contains no species: {path}")
    return species


def filter_hits_to_analysis_species(
    hits: pd.DataFrame, analysis_species: set
) -> pd.DataFrame:
    """Retain only hits belonging to the declared final analysis panel."""

    if "Desmognathus_Species" not in hits.columns:
        raise ValueError("RepeatMasker hits lack Desmognathus_Species")
    hit_species = hits["Desmognathus_Species"].map(_canonical_species)
    return hits.loc[hit_species.isin(analysis_species)].copy()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--repeatmasker-input", type=Path, default=DEFAULT_REPEATMASKER_INPUT
    )
    parser.add_argument("--dnapipete-input", type=Path, default=DEFAULT_DNAPIPETE_INPUT)
    parser.add_argument(
        "--species-file",
        type=Path,
        default=DEFAULT_SPECIES_FILE,
        help="Final analysis panel whose species are retained in the corrected table.",
    )
    parser.add_argument("--lookup-file", type=Path, default=DEFAULT_LOOKUP_FILE)
    parser.add_argument("--chunksize", type=int, default=250_000)
    parser.add_argument(
        "--limit-rows",
        type=int,
        default=None,
        help="Smoke-test limit; omit for the full table.",
    )
    return parser.parse_args()


def verify_corrected_file(path: Path, chunksize: int = 250_000) -> dict:
    """Read back a corrected table and fail if aliases or conserved totals drift."""

    usecols = [
        "Class",
        "Order",
        "Superfamily",
        "repeatmasker_te_class",
        "repeatmasker_order",
        "repeatmasker_superfamily",
        "repeatmasker_classification_status",
        "query_start",
        "query_end",
        "SRX_ID",
        "Desmognathus_Species",
    ]
    n_hits = 0
    aligned_bp = 0
    statuses: Counter = Counter()
    srx_ids = set()
    species = set()
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=chunksize, low_memory=False):
        for canonical, explicit in (
            ("Class", "repeatmasker_te_class"),
            ("Order", "repeatmasker_order"),
            ("Superfamily", "repeatmasker_superfamily"),
        ):
            left = chunk[canonical].astype("string").fillna("<missing>")
            right = chunk[explicit].astype("string").fillna("<missing>")
            if not left.equals(right):
                raise RuntimeError(
                    f"Corrected RepeatMasker alias {canonical} differs from {explicit}"
                )
        start = pd.to_numeric(chunk["query_start"], errors="coerce")
        end = pd.to_numeric(chunk["query_end"], errors="coerce")
        bp = (end - start).abs() + 1
        if bp.isna().any() or (bp <= 0).any():
            raise RuntimeError("Corrected RepeatMasker readback has invalid coordinates")
        n_hits += len(chunk)
        aligned_bp += int(bp.sum())
        statuses.update(
            chunk["repeatmasker_classification_status"].value_counts().to_dict()
        )
        srx_ids.update(chunk["SRX_ID"].dropna().astype(str).unique())
        species.update(chunk["Desmognathus_Species"].dropna().astype(str).unique())
    return {
        "n_hits": n_hits,
        "inclusive_aligned_bp": aligned_bp,
        "classification_status_counts": dict(sorted(statuses.items())),
        "n_srx_ids": len(srx_ids),
        "n_species": len(species),
        "species_labels": sorted(species),
        "canonical_aliases_match": True,
    }


def main() -> None:
    args = parse_args()
    analysis_species = _load_analysis_species(args.species_file)
    analysis_species_set = set(analysis_species)
    lookup = pd.read_csv(args.lookup_file, sep="\t", dtype=str)
    lookup["canonical_species"] = lookup["Species"].map(_canonical_species)
    selected_lookup = lookup.loc[lookup["canonical_species"].isin(analysis_species_set)]
    missing_species = sorted(analysis_species_set.difference(selected_lookup["canonical_species"]))
    if missing_species:
        raise ValueError(f"Analysis species are absent from the TE lookup: {missing_species}")
    selected_srx = set(selected_lookup["SRA_Accension"].dropna())
    output = args.output.resolve()
    if output == PRE_AUDIT_CANONICAL_OUTPUT.resolve():
        raise ValueError(
            "Refusing to replace the pre-audit canonical output. Choose a new "
            "accession- or correction-labeled --output path."
        )
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    partial_output = output.with_suffix(output.suffix + ".partial")
    occupied = [path for path in (output, manifest_path, partial_output) if path.exists()]
    if occupied:
        raise FileExistsError(
            "Non-destructive rebuild requires unused output paths: "
            + ", ".join(str(path) for path in occupied)
        )

    dnapipete = pd.read_csv(args.dnapipete_input, usecols=DNAPIPETE_COLUMNS)
    dnapipete_rows_read = len(dnapipete)
    dnapipete = dnapipete.loc[dnapipete["Source"].isin(selected_srx)].copy()
    context_input_rows = len(dnapipete)
    context = prepare_dnapipete_context(dnapipete)

    output.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    aligned_bp = 0
    statuses: Counter = Counter()
    first_chunk = True
    reader = pd.read_csv(
        args.repeatmasker_input,
        chunksize=args.chunksize,
        nrows=args.limit_rows,
        low_memory=False,
    )
    for hits in reader:
        hits = filter_hits_to_analysis_species(hits, analysis_species_set)
        if hits.empty:
            continue
        annotated = annotate_repeatmasker_hits(hits, context)
        start = pd.to_numeric(annotated["query_start"], errors="coerce")
        end = pd.to_numeric(annotated["query_end"], errors="coerce")
        bp = (end - start).abs() + 1
        if bp.isna().any() or (bp <= 0).any():
            raise ValueError("Corrected RepeatMasker output contains invalid coordinates")

        annotated.to_csv(
            partial_output,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            index=False,
        )
        first_chunk = False
        n_rows += len(annotated)
        aligned_bp += int(bp.sum())
        statuses.update(annotated["repeatmasker_classification_status"].value_counts().to_dict())

    if n_rows == 0:
        raise RuntimeError("Corrected RepeatMasker rebuild produced zero rows")
    readback = verify_corrected_file(partial_output, chunksize=args.chunksize)
    if readback["n_hits"] != n_rows or readback["inclusive_aligned_bp"] != aligned_bp:
        raise RuntimeError("Corrected RepeatMasker readback failed row/bp conservation")
    if readback["classification_status_counts"] != dict(sorted(statuses.items())):
        raise RuntimeError("Corrected RepeatMasker readback status counts changed")
    readback_species = sorted(_canonical_species(value) for value in readback["species_labels"])
    if readback_species != analysis_species:
        raise RuntimeError(
            "Corrected RepeatMasker output species differ from the declared analysis panel"
        )
    partial_output.replace(output)

    manifest = {
        "output_path": _portable(output),
        "output_sha256": _sha256(output),
        "source_state": "corrected_hit_level_repeatmasker_classification",
        "repeatmasker_input": _portable(args.repeatmasker_input),
        "dnapipete_input": _portable(args.dnapipete_input),
        "species_file": _portable(args.species_file),
        "species_file_sha256": _sha256(args.species_file),
        "analysis_species": analysis_species,
        "analysis_species_count": len(analysis_species),
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "eligible_for_path_analysis": True,
        "pre_audit_output_preserved": _portable(PRE_AUDIT_CANONICAL_OUTPUT),
        "n_hits": n_rows,
        "inclusive_aligned_bp": aligned_bp,
        "dnapipete_rows_read": dnapipete_rows_read,
        "dnapipete_context_input_rows": context_input_rows,
        "dnapipete_context_unique_rows": len(context),
        "classification_status_counts": dict(sorted(statuses.items())),
        "readback_verification": readback,
        "limit_rows": args.limit_rows,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()

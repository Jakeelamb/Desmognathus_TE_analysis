#!/usr/bin/env python3
"""Build a compact divergence landscape from frozen corrected RepeatMasker hits.

This is a downstream cache builder. It does not invoke RepeatMasker or any
repeat-discovery software. Inclusive RepeatMasker query-coordinate bp are
binned by native RepeatMasker order and percent divergence, then normalized
within species to the total corrected aligned hit bp represented in the table.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv"
INPUT_MANIFEST = Path(str(INPUT) + ".manifest.json")
OUTPUT_DIR = ROOT / "results/data/corrected/repeat_landscape"
OUTPUT = OUTPUT_DIR / "repeatmasker_divergence_landscape_analysis18_v1.csv"
MANIFEST = OUTPUT_DIR / "repeatmasker_divergence_landscape_analysis18_v1.manifest.json"

FINAL18 = {
    "amphileucus", "anicetus", "apalachicolae", "auriculatus", "bairdi", "campi",
    "fuscus", "gvnigeusgwotli", "intermedius", "kanawha", "marmoratus", "mavrokoilius",
    "monticola", "ocoee", "perlapsus", "tilleyi", "valtos", "welteri",
}


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
        .str.replace("D.", "", regex=False)
        .str.lower()
    )


def summarize_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    frame = chunk.copy()
    frame["species"] = canonical_species(frame["Desmognathus_Species"])
    frame["order"] = frame["Order"].fillna("Unclassified").astype(str).str.strip()
    frame.loc[frame["order"].eq(""), "order"] = "Unclassified"
    frame["percent_divergence"] = pd.to_numeric(frame["percent_divergence"], errors="coerce")
    start = pd.to_numeric(frame["query_start"], errors="coerce")
    end = pd.to_numeric(frame["query_end"], errors="coerce")
    # RM_hit_length_bp belongs to the dnaPipeTE context row and is not the
    # length of this RepeatMasker hit. The corrected-table manifest conserves
    # inclusive RepeatMasker query coordinates, so use the same quantity here.
    frame["hit_bp"] = (end - start).abs() + 1
    frame = frame[
        frame["species"].isin(FINAL18)
        & frame["percent_divergence"].ge(0)
        & frame["hit_bp"].gt(0)
    ].copy()
    frame["divergence_bin_start_pct"] = np.floor(
        frame["percent_divergence"].clip(upper=50)
    ).astype(int)
    return (
        frame.groupby(
            ["species", "order", "divergence_bin_start_pct"], as_index=False, observed=True
        )
        .agg(hit_count=("hit_bp", "size"), hit_bp=("hit_bp", "sum"))
    )


def build_landscape(chunksize: int = 400_000) -> pd.DataFrame:
    summaries: List[pd.DataFrame] = []
    columns = [
        "Desmognathus_Species",
        "percent_divergence",
        "query_start",
        "query_end",
        "Order",
    ]
    for chunk in pd.read_csv(INPUT, usecols=columns, chunksize=chunksize):
        summaries.append(summarize_chunk(chunk))
    combined = pd.concat(summaries, ignore_index=True)
    combined = (
        combined.groupby(
            ["species", "order", "divergence_bin_start_pct"], as_index=False, observed=True
        )[["hit_count", "hit_bp"]]
        .sum()
    )
    species_total = combined.groupby("species")["hit_bp"].transform("sum")
    combined["species_total_hit_bp"] = species_total
    combined["fraction_species_hit_bp"] = combined["hit_bp"] / species_total
    combined["percent_species_hit_bp"] = 100 * combined["fraction_species_hit_bp"]
    combined["divergence_bin_label"] = combined["divergence_bin_start_pct"].map(
        lambda value: "50+" if value == 50 else f"{value}-{value + 1}"
    )
    return combined.sort_values(
        ["species", "order", "divergence_bin_start_pct"]
    ).reset_index(drop=True)


def validate(frame: pd.DataFrame, expected_bp: int) -> Dict[str, object]:
    observed_species = set(frame["species"])
    if observed_species != FINAL18:
        raise ValueError(f"Landscape species mismatch: {sorted(observed_species ^ FINAL18)}")
    observed_bp = int(round(frame["hit_bp"].sum()))
    if observed_bp != expected_bp:
        raise ValueError(f"Hit-bp conservation failed: {observed_bp} != {expected_bp}")
    fractions = frame.groupby("species")["fraction_species_hit_bp"].sum()
    if not np.allclose(fractions.values, 1.0, atol=1e-12):
        raise ValueError("Within-species landscape fractions do not sum to one")
    return {
        "n_rows": int(len(frame)),
        "n_species": int(frame["species"].nunique()),
        "n_orders": int(frame["order"].nunique()),
        "n_divergence_bins": int(frame["divergence_bin_start_pct"].nunique()),
        "hit_count": int(frame["hit_count"].sum()),
        "hit_bp": observed_bp,
        "maximum_fraction_sum_error": float((fractions - 1).abs().max()),
    }


def main() -> None:
    if not INPUT.exists() or not INPUT_MANIFEST.exists():
        raise FileNotFoundError("Corrected RepeatMasker hit table or manifest is missing")
    source_manifest = json.loads(INPUT_MANIFEST.read_text())
    expected_bp = int(source_manifest["inclusive_aligned_bp"])
    landscape = build_landscape()
    validation = validate(landscape, expected_bp)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    landscape.to_csv(OUTPUT, index=False)
    manifest = {
        "analysis_id": "repeatmasker_divergence_landscape_analysis18_v1",
        "input_path": str(INPUT.relative_to(ROOT)),
        "input_sha256_from_source_manifest": source_manifest["output_sha256"],
        "input_rows_from_source_manifest": source_manifest["n_hits"],
        "input_hit_bp_from_source_manifest": expected_bp,
        "classification_basis": "corrected native RepeatMasker hit order",
        "abundance_denominator": (
            "total inclusive corrected RepeatMasker query-coordinate aligned bp "
            "within species; overlapping hits are not interval-deduplicated"
        ),
        "divergence_bin_width_percent": 1,
        "overflow_bin": "50+ percent",
        "expensive_upstream_tools_executed": False,
        "validation": validation,
        "output_path": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": sha256(OUTPUT),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()

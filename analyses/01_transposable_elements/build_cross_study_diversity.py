"""Build the descriptive cross-study top-10 TE-diversity table.

The Wang et al. rows are copied from the lossless normalized Dryad snapshot.
Current Desmognathus indices are copied from the canonical top-10-superfamily
diversity release. No regression or hypothesis test is performed here.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(__file__).resolve().parent / "data"
PUBLISHED_INPUT = DATA_DIR / "wang2025_figure2_te_diversity.csv"
DIVERSITY_INPUT = DATA_DIR / "te_diversity.csv"
IOD_INPUT = ROOT / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
OUTPUT = DATA_DIR / "te_diversity_cross_study.csv"

DATASET_DOI = "10.5061/dryad.zpc866tkv"
PUBLISHED_FILE = "Fig_2/Data_TE_diversity.xlsx"
PUBLISHED_FILE_SHA256 = "367a8102fd18f9c17251c442b0a0112af6afc008b51bf4b110f72a98fcc05c70"
CURRENT_FILES = (
    "analyses/01_transposable_elements/data/te_diversity.csv;"
    "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
)
CURRENT_METRIC_SCOPE = "current_top10_superfamilies_dnapipete_retry_averaged"
PUBLISHED_METRIC_SCOPE = "published_top10_superfamilies_repeatmodeler_repeatmasker"
PG_TO_GB = 0.978

OUTPUT_COLUMNS = [
    "data_origin",
    "source_record_id",
    "source_doi",
    "source_file",
    "source_file_sha256",
    "metric_scope",
    "scientific_name",
    "common_name",
    "family",
    "comparison_group",
    "genome_size_gb",
    "genome_size_ci_low_gb",
    "genome_size_ci_high_gb",
    "shannon_top10_superfamilies",
    "gini_simpson_top10_superfamilies",
    "figure_inclusion",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def published_comparison_group(row: pd.Series) -> str:
    if row["figure_group"] == "Fishes":
        shape = row["fish_point_shape"]
        if pd.isna(shape):
            return "Bony fishes"
        return "Lobe-finned fishes" if int(shape) == 8 else "Cartilaginous fishes"
    return {
        "Anura": "Frogs",
        "Caudata": "Published salamanders",
        "Non-avian Reptillia": "Non-avian reptiles",
        "Aves": "Birds",
        "Mammalia": "Mammals",
    }[row["figure_group"]]


def build_published_rows() -> pd.DataFrame:
    source = pd.read_csv(PUBLISHED_INPUT, float_precision="round_trip")
    result = pd.DataFrame(
        {
            "data_origin": "wang2025_published",
            "source_record_id": source["source_row"].map(lambda value: f"wang2025_{value:03d}"),
            "source_doi": DATASET_DOI,
            "source_file": PUBLISHED_FILE,
            "source_file_sha256": PUBLISHED_FILE_SHA256,
            "metric_scope": PUBLISHED_METRIC_SCOPE,
            "scientific_name": source["species"],
            "common_name": source["english_name"],
            "family": source["family"],
            "comparison_group": source.apply(published_comparison_group, axis=1),
            "genome_size_gb": source["genome_size_gb"],
            "genome_size_ci_low_gb": np.nan,
            "genome_size_ci_high_gb": np.nan,
            "shannon_top10_superfamilies": source["shannon_top10_superfamilies"],
            "gini_simpson_top10_superfamilies": source["gini_simpson_top10_superfamilies"],
            "figure_inclusion": np.where(
                source["figure_group"].eq("Caudata"),
                "published_scatter_and_group",
                "published_group_only",
            ),
        }
    )
    return result[OUTPUT_COLUMNS]


def build_current_rows() -> pd.DataFrame:
    metrics = pd.read_csv(
        DIVERSITY_INPUT,
        float_precision="round_trip",
    )
    expected_columns = [
        "species",
        "n_superfamilies_included",
        "shannon_entropy",
        "gini_simpson",
    ]
    if list(metrics.columns) != expected_columns or metrics.shape != (34, 4):
        raise ValueError("canonical diversity must have the exact 34 x 4 S08 schema")
    if metrics["species"].duplicated().any():
        raise ValueError("canonical diversity species must be unique")
    if not metrics["n_superfamilies_included"].eq(10).all():
        raise ValueError("canonical diversity must include ten superfamilies per species")
    metric_values = metrics[["shannon_entropy", "gini_simpson"]].to_numpy(dtype=float)
    if not np.isfinite(metric_values).all():
        raise ValueError("canonical diversity metrics must be finite")
    metrics = metrics.rename(
        columns={
            "shannon_entropy": "shannon_top10_superfamilies",
            "gini_simpson": "gini_simpson_top10_superfamilies",
        }
    )
    metrics = metrics.drop(columns="n_superfamilies_included")
    iod = pd.read_csv(IOD_INPUT, float_precision="round_trip").copy()
    iod["species"] = iod["species"].str.removeprefix("D. ")
    iod = iod[
        [
            "species",
            "genome_size_pg_fuscus_anchored",
            "genome_size_pg_ci_low",
            "genome_size_pg_ci_high",
        ]
    ]
    joined = metrics.merge(iod, on="species", how="left", validate="one_to_one")
    diversity_hash = sha256(DIVERSITY_INPUT)
    iod_hash = sha256(IOD_INPUT)
    source_hashes = f"{diversity_hash};{iod_hash}"

    has_genome_size = joined["genome_size_pg_fuscus_anchored"].notna()
    result = pd.DataFrame(
        {
            "data_origin": "desmognathus_current",
            "source_record_id": joined["species"].map(lambda value: f"desmo_{value}"),
            "source_doi": "this_study",
            "source_file": CURRENT_FILES,
            "source_file_sha256": source_hashes,
            "metric_scope": CURRENT_METRIC_SCOPE,
            "scientific_name": "Desmognathus " + joined["species"],
            "common_name": np.nan,
            "family": "Plethodontidae",
            "comparison_group": "Current Desmognathus",
            "genome_size_gb": joined["genome_size_pg_fuscus_anchored"] * PG_TO_GB,
            "genome_size_ci_low_gb": joined["genome_size_pg_ci_low"] * PG_TO_GB,
            "genome_size_ci_high_gb": joined["genome_size_pg_ci_high"] * PG_TO_GB,
            "shannon_top10_superfamilies": joined["shannon_top10_superfamilies"],
            "gini_simpson_top10_superfamilies": joined["gini_simpson_top10_superfamilies"],
            "figure_inclusion": np.where(
                has_genome_size,
                "current_scatter_and_group",
                "current_group_only",
            ),
        }
    )
    return result[OUTPUT_COLUMNS]


def build() -> pd.DataFrame:
    result = pd.concat([build_published_rows(), build_current_rows()], ignore_index=True)
    if result.shape != (118, 16):
        raise ValueError(f"expected a 118 x 16 table, observed {result.shape}")
    if result["source_record_id"].duplicated().any():
        raise ValueError("source_record_id values must be unique")
    return result


def csv_bytes(frame: pd.DataFrame) -> bytes:
    stream = io.StringIO(newline="")
    frame.to_csv(stream, index=False, lineterminator="\n", float_format="%.17g")
    return stream.getvalue().encode("utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="replace the canonical output instead of verifying it read-only",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rendered = csv_bytes(build())
    if args.write:
        OUTPUT.write_bytes(rendered)
        print(f"wrote {OUTPUT.relative_to(ROOT)}")
        return 0
    if not OUTPUT.exists():
        print(f"missing canonical output: {OUTPUT.relative_to(ROOT)}", file=sys.stderr)
        return 1
    if OUTPUT.read_bytes() != rendered:
        print(
            f"canonical output differs: rerun with --write to replace {OUTPUT.relative_to(ROOT)}",
            file=sys.stderr,
        )
        return 1
    print(f"verified {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

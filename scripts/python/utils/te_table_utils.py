#!/usr/bin/env python3
"""
Helpers for normalizing TE feature tables to a species-by-feature matrix.
"""

from __future__ import annotations

import re
from typing import Iterable, Tuple

import pandas as pd


SPECIES_FIRST_COLUMN_HINTS = {
    "species",
    "taxon",
    "taxa",
    "sample",
    "sampleid",
    "sample_id",
    "unnamed0",
}


def _normalize_label(label: object) -> str:
    """Normalize a column label for orientation checks."""
    return re.sub(r"[^a-z0-9]+", "", str(label).strip().lower())


def _looks_like_species(value: object) -> bool:
    """Heuristic for Desmognathus species/sample identifiers."""
    return bool(re.match(r"^(D[._]|Desmognathus[._ ])", str(value).strip()))


def coerce_feature_table_to_species_matrix(
    df: pd.DataFrame,
    feature_axis_labels: Iterable[str],
) -> Tuple[pd.DataFrame, str]:
    """
    Convert a TE table into a canonical species-by-feature numeric matrix.

    Returns the normalized matrix and a short orientation label describing
    the input layout that was detected.
    """
    if df.empty:
        raise ValueError("Input table is empty.")

    feature_hints = {_normalize_label(label) for label in feature_axis_labels}
    first_col = df.columns[0]
    first_col_norm = _normalize_label(first_col)

    if first_col_norm in SPECIES_FIRST_COLUMN_HINTS:
        matrix = df.rename(columns={first_col: "Species"}).set_index("Species")
        orientation = "species_rows"
    elif first_col_norm in feature_hints:
        matrix = df.rename(columns={first_col: "feature"}).set_index("feature").transpose()
        matrix.index.name = "Species"
        orientation = "feature_rows"
    else:
        first_values = df.iloc[: min(len(df), 5), 0]
        species_like_rows = sum(_looks_like_species(value) for value in first_values)
        species_like_cols = sum(_looks_like_species(label) for label in df.columns[1:])
        if species_like_rows >= max(1, len(first_values) // 2):
            matrix = df.rename(columns={first_col: "Species"}).set_index("Species")
            orientation = "species_rows_inferred"
        elif species_like_cols >= max(1, len(df.columns[1:]) // 2):
            matrix = df.rename(columns={first_col: "feature"}).set_index("feature").transpose()
            matrix.index.name = "Species"
            orientation = "feature_rows_inferred"
        else:
            raise ValueError(
                "Could not infer TE table orientation from the first column or labels."
            )

    matrix.index = matrix.index.map(str)
    matrix.columns = [str(column) for column in matrix.columns]

    drop_columns = [
        column
        for column in matrix.columns
        if _normalize_label(column) in {"total", "totals", "sum"}
    ]
    if drop_columns:
        matrix = matrix.drop(columns=drop_columns)

    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    return matrix, orientation

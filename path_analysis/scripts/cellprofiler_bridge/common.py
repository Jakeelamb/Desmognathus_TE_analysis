from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


REFERENCE_GENOME_PG = 16.36
GENOME_GB_PER_PG = 0.978


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(value: object) -> str | pd.NA:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return pd.NA
    text = str(value).strip()
    if not text:
        return pd.NA
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip() or pd.NA


def nonempty_pct(series: pd.Series) -> float:
    values = series.fillna("").astype(str).str.strip()
    return float(values.ne("").mean() * 100.0) if len(values) else 0.0


def existing_pct(series: pd.Series) -> float:
    values = series.fillna("").astype(str).str.strip()
    if len(values) == 0:
        return 0.0
    exists = values.map(lambda value: Path(value).exists() if value else False)
    return float(exists.mean() * 100.0)


def file_record(label: str, path: Path) -> dict[str, object]:
    return {
        "label": label,
        "path": str(path.resolve()),
        "exists": path.exists(),
        "sha256": sha256_for_file(path) if path.exists() else None,
    }


def _json_safe(value: object) -> object:
    if value is pd.NA:
        return None
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and pd.isna(value):
        return None
    return value


def write_json(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_safe(data), indent=2, sort_keys=True) + "\n")

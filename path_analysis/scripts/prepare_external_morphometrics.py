#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "raw"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

CONANTI_FUSCUS_FILE = RAW_DIR / "Pyron_Beamer_conanti_fuscus_Appendix_S1.csv"
OUTPUT_FILE = DERIVED_DIR / "conanti_fuscus_morphometrics_summary.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_conanti_fuscus(df: pd.DataFrame, source_hash: str) -> pd.DataFrame:
    numeric_cols = [
        "SVL",
        "TL",
        "AG",
        "CW",
        "FL",
        "HL",
        "SG",
        "TW",
        "TO",
        "FI",
        "HW",
        "ED",
        "IN",
        "ES",
        "ON",
        "IO",
        "IC",
    ]

    summary = (
        df.groupby("Species")
        .agg(
            source_id=("Species", lambda _: "pyron_beamer_2023_conanti_fuscus_appendix_s1"),
            source_measurement_context=("Species", lambda _: "all_specimens_mixed_life_stages"),
            source_filename=("Species", lambda _: CONANTI_FUSCUS_FILE.name),
            source_file_sha256=("Species", lambda _: source_hash),
            morph_n_specimens=("Specimen", "count"),
            morph_n_counties=("County", "nunique"),
            morph_n_states=("State", "nunique"),
            svl_mean_all_mm=("SVL", "mean"),
            svl_sd_all_mm=("SVL", "std"),
            svl_q90_all_mm=("SVL", lambda s: float(np.quantile(s, 0.9))),
            svl_max_mm=("SVL", "max"),
            tl_mean_all_mm=("TL", "mean"),
            tl_q90_all_mm=("TL", lambda s: float(np.quantile(s, 0.9))),
            tl_max_mm=("TL", "max"),
        )
        .reset_index()
        .rename(columns={"Species": "species"})
    )

    measurement_means = (
        df.groupby("Species")[numeric_cols]
        .mean()
        .add_prefix("mean_")
        .reset_index()
        .rename(columns={"Species": "species"})
    )
    summary = summary.merge(measurement_means, on="species", how="left")
    return summary


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)

    if not CONANTI_FUSCUS_FILE.exists():
        raise FileNotFoundError(
            f"Missing raw source file: {CONANTI_FUSCUS_FILE}\n"
            "Download the appendix and place it in path_analysis/data/external/raw/."
        )

    df = pd.read_csv(CONANTI_FUSCUS_FILE)
    source_hash = sha256_for_file(CONANTI_FUSCUS_FILE)
    summary = summarize_conanti_fuscus(df, source_hash)
    summary.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Species: {len(summary)}")


if __name__ == "__main__":
    main()

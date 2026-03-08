#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "raw"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

AMPHIBIO_ZIP = RAW_DIR / "AmphiBIO_v1.zip"
OUTPUT_FILE = DERIVED_DIR / "amphibio_desmognathus_traits.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def standardize_species(species: pd.Series) -> pd.Series:
    return species.astype(str).str.strip().str.split().str[-1].str.lower()


def derive_development_mode(row: pd.Series) -> str | None:
    direct = row.get("Dir")
    larval = row.get("Lar")
    viv = row.get("Viv")
    if pd.isna(direct) and pd.isna(larval) and pd.isna(viv):
        return None
    if direct == 1 and larval != 1:
        return "direct_development"
    if larval == 1 and direct != 1:
        return "aquatic_larva"
    if viv == 1:
        return "viviparous"
    return "mixed_or_unclear"


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)

    if not AMPHIBIO_ZIP.exists():
        raise FileNotFoundError(
            f"Missing raw source file: {AMPHIBIO_ZIP}\n"
            "Download the AmphiBIO archive and place it in path_analysis/data/external/raw/."
        )

    source_hash = sha256_for_file(AMPHIBIO_ZIP)
    with zipfile.ZipFile(AMPHIBIO_ZIP) as archive:
        with archive.open("AmphiBIO_v1.csv") as handle:
            df = pd.read_csv(handle, encoding="latin1")

    df = df[df["Genus"].astype(str).str.lower() == "desmognathus"].copy()
    df["species"] = standardize_species(df["Species"])
    df["development_mode"] = df.apply(derive_development_mode, axis=1)

    keep_cols = [
        "species",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Body_mass_g",
        "Age_at_maturity_min_y",
        "Age_at_maturity_max_y",
        "Body_size_mm",
        "Size_at_maturity_min_mm",
        "Size_at_maturity_max_mm",
        "Longevity_max_y",
        "Litter_size_min_n",
        "Litter_size_max_n",
        "Reproductive_output_y",
        "Offspring_size_min_mm",
        "Offspring_size_max_mm",
        "Ter",
        "Aqu",
        "Arb",
        "Diu",
        "Noc",
        "Crepu",
        "Wet_warm",
        "Wet_cold",
        "Dry_warm",
        "Dry_cold",
        "Dir",
        "Lar",
        "Viv",
        "OBS",
        "development_mode",
    ]
    out = df[keep_cols].copy()
    out = out.rename(columns={c: f"amphibio_{c.lower()}" for c in keep_cols if c != "species"})
    out["amphibio_source_id"] = "amphibio_2017"
    out["amphibio_source_filename"] = AMPHIBIO_ZIP.name
    out["amphibio_source_file_sha256"] = source_hash
    out.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Species: {len(out)}")


if __name__ == "__main__":
    main()

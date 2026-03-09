#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results" / "data" / "ltr_age"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "derived"

SPECIES_SUMMARY_FILE = RESULTS_DIR / "ltr_age_species_summary.csv"
CALIBRATED_FILE = RESULTS_DIR / "ltr_age_species_summary_calibrated.csv"
CANDIDATES_FILE = RESULTS_DIR / "ltr_substitution_rate_candidates.csv"
OUTPUT_FILE = DERIVED_DIR / "ltr_history_features.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(value: object) -> str | pd.NA:
    if value is None or pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if not text:
        return pd.NA
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def load_species_summary() -> pd.DataFrame:
    df = pd.read_csv(SPECIES_SUMMARY_FILE)
    df["species"] = df["species"].map(canonical_species)
    rename_map = {
        "n_pairs_attempted": "ltr_history_n_pairs_attempted",
        "n_pairs_estimated": "ltr_history_n_pairs_estimated",
        "n_high_confidence_pairs_estimated": "ltr_history_n_pairs_high_confidence",
        "median_p_distance": "ltr_history_median_p_distance",
        "median_k2p_distance": "ltr_history_median_k2p_distance",
        "median_comparable_sites": "ltr_history_median_comparable_sites",
    }
    out = df.rename(columns=rename_map)[["species"] + list(rename_map.values())].copy()
    return out.dropna(subset=["species"]).drop_duplicates(subset=["species"]).reset_index(drop=True)


def load_calibrated_summary() -> pd.DataFrame:
    df = pd.read_csv(CALIBRATED_FILE)
    df["species"] = df["species"].map(canonical_species)
    keep = df[
        df["rate_label"].isin(
            [
                "salamander_rag1_floor",
                "frog_cmyc_midpoint",
                "frog_cmyc_high",
            ]
        )
    ].copy()
    keep["age_col"] = keep["rate_label"].map(
        {
            "salamander_rag1_floor": "ltr_history_age_low_mya",
            "frog_cmyc_midpoint": "ltr_history_age_central_mya",
            "frog_cmyc_high": "ltr_history_age_high_mya",
        }
    )
    keep["age_hc_col"] = keep["rate_label"].map(
        {
            "salamander_rag1_floor": "ltr_history_age_low_high_conf_mya",
            "frog_cmyc_midpoint": "ltr_history_age_central_high_conf_mya",
            "frog_cmyc_high": "ltr_history_age_high_high_conf_mya",
        }
    )

    age_wide = keep.pivot(index="species", columns="age_col", values="median_age_mya")
    hc_wide = keep.pivot(index="species", columns="age_hc_col", values="median_age_high_confidence_mya")
    out = pd.concat([age_wide, hc_wide], axis=1).reset_index()
    return out.dropna(subset=["species"]).drop_duplicates(subset=["species"]).reset_index(drop=True)


def load_candidate_metadata() -> dict[str, object]:
    df = pd.read_csv(CANDIDATES_FILE)
    recommended = df[df["included_in_primary_summary"].fillna(False)].copy().sort_values("rate_value")
    central = df[df["rate_label"].eq("frog_cmyc_midpoint")].iloc[0]
    return {
        "ltr_history_rate_window_low": float(recommended["rate_value"].iloc[0]),
        "ltr_history_rate_window_high": float(recommended["rate_value"].iloc[-1]),
        "ltr_history_rate_central": float(central["rate_value"]),
    }


def main() -> None:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)

    species_summary = load_species_summary()
    calibrated = load_calibrated_summary()
    metadata = load_candidate_metadata()

    out = species_summary.merge(calibrated, on="species", how="left", validate="one_to_one")
    out["has_ltr_history"] = out["ltr_history_median_k2p_distance"].notna()
    out["ltr_history_source_id"] = "repo_ltr_age_species_summary"
    out["ltr_history_source_path"] = str(SPECIES_SUMMARY_FILE.relative_to(PROJECT_ROOT))
    out["ltr_history_source_sha256"] = sha256_for_file(SPECIES_SUMMARY_FILE)
    out["ltr_history_calibration_source_id"] = "repo_ltr_age_species_summary_calibrated"
    out["ltr_history_calibration_source_path"] = str(CALIBRATED_FILE.relative_to(PROJECT_ROOT))
    out["ltr_history_calibration_source_sha256"] = sha256_for_file(CALIBRATED_FILE)
    out["ltr_history_rate_source_id"] = "repo_ltr_substitution_rate_candidates"
    out["ltr_history_rate_source_path"] = str(CANDIDATES_FILE.relative_to(PROJECT_ROOT))
    out["ltr_history_rate_source_sha256"] = sha256_for_file(CANDIDATES_FILE)
    out["ltr_history_source_ids"] = (
        "repo_ltr_age_species_summary;"
        "repo_ltr_age_species_summary_calibrated;"
        "repo_ltr_substitution_rate_candidates"
    )
    for key, value in metadata.items():
        out[key] = value

    out = out.sort_values("species").reset_index(drop=True)
    out.to_csv(OUTPUT_FILE, index=False)
    print(f"Wrote {OUTPUT_FILE}")
    print(f"Rows: {len(out)}")
    print(f"Species with LTR history: {int(out['has_ltr_history'].sum())}")


if __name__ == "__main__":
    main()

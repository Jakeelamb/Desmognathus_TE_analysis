#!/usr/bin/env python3
"""Assemble corrected final-18 inputs for path-model sensitivity analysis.

No historical path-analysis table is overwritten.  The image-derived node is
named ``relative_nuclear_iod_proxy`` and is explicitly not an absolute genome
size.  Morphology and IOD variants are crossed so every downstream model fit
can name its exact measurement specification.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DIVERSITY = PROJECT_ROOT / "results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv"
COMPOSITION = PROJECT_ROOT / "results/data/corrected/diversity_pca/te_composition_matrices_analysis18_v1.csv"
ECTOPIC = PROJECT_ROOT / "results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv"
MORPHOLOGY = PROJECT_ROOT / "results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv"
IOD = PROJECT_ROOT / "results/data/corrected/microscopy/microscopy_relative_iod_sensitivity_analysis18_v1.csv"
SUPPORT = PROJECT_ROOT / "results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv"

OUTPUT_DIR = PROJECT_ROOT / "results/data/corrected/path_analysis"
INPUT_OUTPUT = OUTPUT_DIR / "corrected_path_input_sensitivity_analysis18_v1.csv"
SPECIFICATION_OUTPUT = OUTPUT_DIR / "corrected_path_measurement_specifications_analysis18_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "corrected_path_inputs_analysis18_v1.manifest.json"

FINAL_SPECIES = [
    "amphileucus", "anicetus", "apalachicolae", "auriculatus", "bairdi", "campi",
    "fuscus", "gvnigeusgwotli", "intermedius", "kanawha", "marmoratus",
    "mavrokoilius", "monticola", "ocoee", "perlapsus", "tilleyi", "valtos", "welteri",
]

IOD_COLUMNS = {
    "all_selected": "relative_iod_to_fuscus_all",
    "image_qc_pass": "relative_iod_to_fuscus_image_qc_pass",
    "high_qc": "relative_iod_to_fuscus_high_qc",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_exists(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing corrected path inputs: " + ", ".join(missing))


def order_features(composition: pd.DataFrame, diversity: pd.DataFrame) -> pd.DataFrame:
    order = composition[
        composition["te_level"].eq("order")
        & composition["composition_mode"].eq("classified_conditional")
    ].copy()
    wide = order.pivot(index="species", columns="feature", values="proportion")
    if wide.isna().any().any():
        raise ValueError("Order composition matrix contains missing values")
    if not {"LTR", "LINE"}.issubset(wide.columns):
        raise ValueError("Order composition lacks LTR or LINE")
    if (wide[["LTR", "LINE"]] <= 0).any().any():
        raise ValueError("LTR:LINE log-ratio requires positive final-panel proportions")
    wide["ltr_line_logratio"] = np.log(wide["LTR"] / wide["LINE"])
    wide = wide.reset_index()

    div = diversity[
        diversity["te_level"].eq("order")
        & diversity["composition_mode"].eq("classified_conditional")
    ][
        [
            "species", "observed_richness", "shannon_entropy", "gini_simpson",
            "hill_q1", "hill_q2", "pielou_evenness",
        ]
    ].copy()
    return wide.merge(div, on="species", how="inner", validate="one_to_one")


def morphology_wide(frame: pd.DataFrame) -> pd.DataFrame:
    wide = frame.pivot_table(
        index=["species", "estimator", "estimator_label"],
        columns="metric",
        values="estimate",
        aggfunc="first",
    ).reset_index()
    required = {"cell_area_um2", "nuc_area_um2", "nc_area_ratio"}
    if not required.issubset(wide.columns):
        raise ValueError(f"Morphology estimators lack fields: {sorted(required.difference(wide.columns))}")
    return wide.rename(columns={"estimator": "morphology_estimator"})


def iod_long(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in frame.itertuples(index=False):
        for subset, column in IOD_COLUMNS.items():
            rows.append(
                {
                    "species": row.species,
                    "iod_subset": subset,
                    "relative_nuclear_iod_proxy": float(getattr(row, column)),
                    "iod_n_images": int(row.n_images),
                    "iod_n_specimens": int(row.n_specimens),
                    "iod_n_pass_images": int(row.n_pass_images),
                    "iod_n_high_qc_nuclei": int(row.n_high_qc_nuclei),
                    "iod_proxy_is_absolute_genome_size": False,
                }
            )
    return pd.DataFrame(rows)


def cross_measurement_specifications(morphology: pd.DataFrame, iod: pd.DataFrame) -> pd.DataFrame:
    return morphology.merge(iod, on="species", how="inner", validate="many_to_many")


def main() -> None:
    source_paths = [DIVERSITY, COMPOSITION, ECTOPIC, MORPHOLOGY, IOD, SUPPORT]
    require_exists(source_paths)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    te = order_features(pd.read_csv(COMPOSITION), pd.read_csv(DIVERSITY))
    morph = morphology_wide(pd.read_csv(MORPHOLOGY))
    iod = iod_long(pd.read_csv(IOD))
    support = pd.read_csv(SUPPORT)[
        [
            "species", "n_images", "n_specimens", "n_manual_keep",
            "n_model_ranked_unlabeled", "manual_keep_fraction",
        ]
    ].rename(columns={"n_images": "morphology_n_images", "n_specimens": "morphology_n_specimens"})
    ectopic_source = pd.read_csv(ECTOPIC)
    ectopic = ectopic_source[ectopic_source["analysis_branch"].eq("all_5plus_domain_elements")][
        [
            "species", "n_elements", "ratio_median", "ratio_geometric_mean",
            "median_bootstrap_ci_low", "median_bootstrap_ci_high",
        ]
    ].rename(
        columns={
            "n_elements": "ectopic_n_elements",
            "ratio_median": "terminal_internal_ratio_median",
            "ratio_geometric_mean": "terminal_internal_ratio_geometric_mean",
            "median_bootstrap_ci_low": "terminal_internal_ratio_median_ci_low",
            "median_bootstrap_ci_high": "terminal_internal_ratio_median_ci_high",
        }
    )

    combined = cross_measurement_specifications(morph, iod)
    combined = combined.merge(te, on="species", how="left", validate="many_to_one")
    combined = combined.merge(ectopic, on="species", how="left", validate="many_to_one")
    combined = combined.merge(support, on="species", how="left", validate="many_to_one")
    combined["terminal_internal_proxy_is_ectopic_rate"] = False
    combined["analysis_interpretation"] = "exploratory_measurement_sensitivity"
    combined = combined.sort_values(["morphology_estimator", "iod_subset", "species"]).reset_index(drop=True)

    if set(combined["species"]) != set(FINAL_SPECIES):
        raise ValueError("Corrected path input does not exactly cover the final 18 species")
    expected = len(FINAL_SPECIES) * morph["morphology_estimator"].nunique() * len(IOD_COLUMNS)
    if len(combined) != expected:
        raise ValueError(f"Expected {expected} crossed path rows, observed {len(combined)}")
    if combined["relative_nuclear_iod_proxy"].isna().any():
        raise ValueError("Relative-IOD proxy is missing in a measurement specification")

    combined.to_csv(INPUT_OUTPUT, index=False, float_format="%.10g")
    specifications = (
        combined.groupby(["morphology_estimator", "estimator_label", "iod_subset"], as_index=False)
        .agg(
            n_species=("species", "nunique"),
            n_species_with_terminal_internal_proxy=("terminal_internal_ratio_median", "count"),
            minimum_morphology_specimens=("morphology_n_specimens", "min"),
            total_manual_keeps=("n_manual_keep", "sum"),
            total_model_ranked_unlabeled=("n_model_ranked_unlabeled", "sum"),
        )
    )
    specifications["absolute_genome_size_used"] = False
    specifications["confirmatory_claim_allowed"] = False
    specifications.to_csv(SPECIFICATION_OUTPUT, index=False)

    manifest = {
        "analysis_id": "corrected_path_inputs_analysis18_v1",
        "n_species": len(FINAL_SPECIES),
        "n_morphology_estimators": int(morph["morphology_estimator"].nunique()),
        "n_iod_subsets": len(IOD_COLUMNS),
        "n_measurement_specifications": int(len(specifications)),
        "n_crossed_rows": int(len(combined)),
        "absolute_genome_size_used": False,
        "historical_inputs_overwritten": False,
        "sources": [
            {"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256(path)}
            for path in source_paths
        ],
        "outputs": [
            str(INPUT_OUTPUT.relative_to(PROJECT_ROOT)),
            str(SPECIFICATION_OUTPUT.relative_to(PROJECT_ROOT)),
        ],
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {INPUT_OUTPUT} ({len(combined)} rows)")
    print(f"Wrote {SPECIFICATION_OUTPUT} ({len(specifications)} specifications)")


if __name__ == "__main__":
    main()

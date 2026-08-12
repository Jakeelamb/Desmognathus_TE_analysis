"""Replay the compact morphology release from its frozen object tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path("analyses/02_morphology/data")
IOD_PROVENANCE_PATH = Path("analyses/02_morphology/provenance/iod_source_manifest.json")
AVAILABILITY_PATH = Path("data/identity/analysis_availability.csv")
PIXEL_AREA_UM2 = 0.0144
IOD_BOOTSTRAP_REPLICATES = 2_000
IOD_BOOTSTRAP_SEED = 20260710
MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE = 0.10
MAX_KS_DISTANCE = 0.25
EXPECTED_BALANCE_FAILURES = frozenset({"D. aeneus", "D. ochrophaeus", "D. orestes", "D. wrighti"})
EXPECTED_MORPHOLOGY_SINGLE_ANIMAL_SPECIES = frozenset(
    {
        "D. anicetus",
        "D. apalachicolae",
        "D. auriculatus",
        "D. bairdi",
        "D. gvnigeusgwotli",
        "D. intermedius",
        "D. perlapsus",
    }
)


def _read(root: Path, filename: str) -> pd.DataFrame:
    return pd.read_csv(root / DATA_DIR / filename, low_memory=False)


def _numeric_mismatches(
    observed: pd.Series,
    expected: pd.Series,
    *,
    atol: float = 1e-8,
) -> list[str]:
    shared = observed.index.intersection(expected.index)
    return [
        str(species)
        for species in shared
        if not np.isclose(
            float(observed.loc[species]),
            float(expected.loc[species]),
            rtol=1e-10,
            atol=atol,
        )
    ]


def _load_conditional_anchor(root: Path) -> tuple[dict[str, object], float]:
    manifest = json.loads((root / IOD_PROVENANCE_PATH).read_text())
    anchor = manifest.get("conditional_genome_size_anchor")
    if not isinstance(anchor, dict):
        raise TypeError("IOD provenance lacks conditional_genome_size_anchor")

    required = {
        "reference_species": "D. fuscus",
        "reference_standard_processes": [
            "Process_413_green.ome.tiff",
            "Process_414_green.ome.tiff",
        ],
        "reference_standard_specimen_ids": ["32469", "32470"],
        "basis": "published_rounded_pooled_assembly_span",
        "reported_precision_gbp": 0.1,
        "ploidy_semantics": "1C_equivalent",
        "assembly_accession": "GCA_050004315.1",
        "assembly_source_doi": "10.1093/g3journal/jkaf157",
        "conversion_source_doi": "10.1002/cyto.a.10013",
        "decision_source_id": ("author_curator_decision_2026_08_12_genome_size_anchor"),
        "decision_date": "2026-08-12",
        "status": "conditional_assembly_derived_not_measured_c_value",
        "author_confirmed_same_staining_runs": True,
        "anchor_uncertainty_included": False,
    }
    mismatches = {
        key: (anchor.get(key), value) for key, value in required.items() if anchor.get(key) != value
    }
    if mismatches:
        raise ValueError(f"IOD conditional anchor metadata changed: {mismatches}")

    reported_gbp = float(anchor["reported_assembly_span_gbp"])
    gbp_per_pg = float(anchor["gbp_per_pg"])
    if reported_gbp != 16.1 or gbp_per_pg != 0.978:
        raise ValueError("IOD anchor must derive from 16.1 Gbp / 0.978 Gbp per pg")
    derived_pg = reported_gbp / gbp_per_pg
    if not np.isclose(
        float(anchor["derived_anchor_pg_1c"]),
        derived_pg,
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("Stored IOD anchor does not equal reported_gbp / gbp_per_pg")
    return anchor, derived_pg


def _hierarchical_bootstrap_draws(
    frame: pd.DataFrame,
    *,
    seed: int,
    n_bootstrap: int = IOD_BOOTSTRAP_REPLICATES,
) -> np.ndarray:
    image_values = {
        str(filename): group["nuc_iod"].dropna().to_numpy(dtype=float)
        for filename, group in frame.groupby("filename", sort=True)
    }
    image_values = {name: values for name, values in image_values.items() if len(values)}
    if not image_values:
        raise ValueError("IOD bootstrap requires at least one populated image")
    image_names = np.array(sorted(image_values), dtype=object)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        sampled_images = rng.choice(
            image_names,
            size=len(image_names),
            replace=True,
        )
        image_medians = []
        for image_name in sampled_images:
            values = image_values[str(image_name)]
            sampled_nuclei = rng.choice(values, size=len(values), replace=True)
            image_medians.append(float(np.median(sampled_nuclei)))
        draws[replicate] = float(np.mean(image_medians))
    return draws


def reconstruct_s22_calibration(root: Path | str = ROOT) -> pd.DataFrame:
    """Recompute only the explicitly conditional fuscus-anchored S22 fields."""

    resolved_root = Path(root).resolve()
    objects = _read(resolved_root, "nuclear_iod_objects.csv")
    released = _read(resolved_root, "relative_nuclear_iod_species.csv")
    anchor_metadata, anchor_pg = _load_conditional_anchor(resolved_root)
    reference_species = str(anchor_metadata["reference_species"])

    estimates = released.set_index("species")["iod_equal_image_estimate"]
    if reference_species not in estimates.index:
        raise ValueError(f"S22 lacks conditional reference species {reference_species}")
    reference_iod = float(estimates.loc[reference_species])
    if reference_iod <= 0:
        raise ValueError("Conditional reference-species IOD must be positive")

    bootstrap_draws: dict[str, np.ndarray] = {}
    for species_index, (species, group) in enumerate(objects.groupby("species", sort=True)):
        bootstrap_draws[str(species)] = _hierarchical_bootstrap_draws(
            group,
            seed=IOD_BOOTSTRAP_SEED + species_index * 1009,
        )
    reference_draws = bootstrap_draws[reference_species]
    if not np.all(reference_draws > 0):
        raise ValueError("Conditional reference bootstrap draws must be positive")

    result = released.copy()
    result["iod_ratio_to_fuscus"] = result["iod_equal_image_estimate"] / reference_iod
    calibration_factor = anchor_pg / reference_iod
    result["genome_size_pg_fuscus_anchored"] = (
        result["iod_equal_image_estimate"] * calibration_factor
    )
    for row_index, species in result["species"].items():
        calibrated_draws = bootstrap_draws[str(species)] / reference_draws * anchor_pg
        low, high = np.quantile(calibrated_draws, [0.025, 0.975])
        result.loc[row_index, "genome_size_pg_ci_low"] = float(low)
        result.loc[row_index, "genome_size_pg_ci_high"] = float(high)
    result["genome_size_reference_species"] = reference_species
    result["genome_size_reference_pg"] = anchor_pg
    result["genome_size_calibration_pg_per_iod"] = calibration_factor
    result["genome_size_rank"] = (
        result["genome_size_pg_fuscus_anchored"].rank(method="min", ascending=False).astype(int)
    )
    return result


def project_s22_availability(
    root: Path | str,
    s22: pd.DataFrame,
) -> pd.DataFrame:
    """Project the named conditional sensitivity into the shared audit table."""

    resolved_root = Path(root).resolve()
    availability = pd.read_csv(resolved_root / AVAILABILITY_PATH)
    crosswalk = pd.read_csv(resolved_root / "data/identity/species_taxonomy_crosswalk.csv")
    display_to_current = {
        f"D. {species}": str(species) for species in crosswalk["current_species"].astype(str)
    }
    values = s22[["species", "genome_size_pg_fuscus_anchored"]].copy()
    values["species"] = values["species"].map(display_to_current)
    if values["species"].isna().any() or values["species"].duplicated().any():
        raise ValueError("S22 species cannot be projected uniquely to availability")
    mapped = values.set_index("species")["genome_size_pg_fuscus_anchored"]
    result = availability.copy()
    target = "step06_genome_size_pg_fuscus_anchored"
    result[target] = result["species"].map(mapped)
    return result


def _validate_s19(root: Path, errors: list[str]) -> None:
    objects = _read(root, "cell_nucleus_objects.csv")
    released = _read(root, "cell_nucleus_species_estimates.csv")
    required = {
        "species",
        "n_size_cells",
        "n_size_images",
        "n_size_specimens",
        "largest_specimen_n",
        "largest_specimen_fraction",
        "cell_area_um2",
        "nucleus_area_um2",
    }
    missing = sorted(required - set(released.columns))
    if missing:
        errors.append(f"[S19_FIELDS] missing fields: {', '.join(missing)}")
        return
    if released["species"].duplicated().any():
        errors.append("[S19_SPECIES] species rows are not unique")
        return

    support = objects.groupby("species", sort=False).agg(
        n_size_cells=("nucleus_object_id", "size"),
        n_size_images=("filename", "nunique"),
        n_size_specimens=("specimen_id", "nunique"),
        cell_area_um2=("cell_area_um2", "median"),
        nucleus_area_um2=("nuc_area_um2", "median"),
    )
    specimen_counts = objects.groupby(["species", "specimen_id"], sort=False).size()
    support["largest_specimen_n"] = specimen_counts.groupby(level="species").max()
    support["largest_specimen_fraction"] = support["largest_specimen_n"] / support["n_size_cells"]
    observed = released.set_index("species")
    if set(observed.index) != set(support.index):
        errors.append("[S19_SPECIES] species set differs from S18")
        return

    for column in [
        "n_size_cells",
        "n_size_images",
        "n_size_specimens",
        "largest_specimen_n",
        "largest_specimen_fraction",
        "cell_area_um2",
        "nucleus_area_um2",
    ]:
        mismatches = _numeric_mismatches(observed[column], support[column])
        if mismatches:
            errors.append(f"[S19_{column.upper()}] mismatch: {', '.join(mismatches)}")


def _validate_s23(root: Path, errors: list[str]) -> None:
    objects = _read(root, "nuclear_iod_objects.csv")
    released = _read(root, "nuclear_iod_by_image.csv")
    summary_columns = [
        "specimen_id",
        "n_nuclei",
        "image_median_iod",
        "image_mean_iod",
        "image_median_nucleus_area_um2",
        "image_median_mean_od",
        "image_median_match_log_edge_sharpness",
        "image_median_match_log_relative_ring_noise",
    ]
    required = {"species", "filename", *summary_columns}
    missing = sorted(required - set(released.columns))
    if missing:
        errors.append(f"[S23_FIELDS] missing fields: {', '.join(missing)}")
        return
    if released.duplicated(["species", "filename"]).any():
        errors.append("[S23_IMAGES] species/image rows are not unique")
        return

    specimen_counts = objects.groupby(["species", "filename"])["specimen_id"].nunique()
    if not specimen_counts.eq(1).all():
        errors.append("[S23_SPECIMEN_GROUP] an image maps to multiple specimens in S21")
        return
    expected = objects.groupby(["species", "filename"], sort=False).agg(
        specimen_id=("specimen_id", "first"),
        n_nuclei=("nucleus_object_id", "size"),
        image_median_iod=("nuc_iod", "median"),
        image_mean_iod=("nuc_iod", "mean"),
        image_median_nucleus_area_um2=("nuc_area_um2", "median"),
        image_median_mean_od=("nuc_mean_od", "median"),
        image_median_match_log_edge_sharpness=("match_log_edge_sharpness", "median"),
        image_median_match_log_relative_ring_noise=(
            "match_log_relative_ring_noise",
            "median",
        ),
    )
    observed = released.set_index(["species", "filename"])
    if set(observed.index) != set(expected.index):
        errors.append("[S23_IMAGES] species/image set differs from S21")
        return
    for column in summary_columns:
        mismatches = _numeric_mismatches(observed[column], expected[column])
        if mismatches:
            errors.append(f"[S23_{column.upper()}] mismatch: {', '.join(mismatches)}")


def _validate_s22(root: Path, errors: list[str]) -> None:
    objects = _read(root, "nuclear_iod_objects.csv")
    images = _read(root, "nuclear_iod_by_image.csv")
    released = _read(root, "relative_nuclear_iod_species.csv")
    required = {
        "species",
        "n_nuclei",
        "n_images",
        "n_specimens",
        "iod_equal_image_estimate",
        "iod_equal_image_ci_low",
        "iod_equal_image_ci_high",
        "relative_iod_anchor",
        "relative_iod_index",
        "relative_iod_ci_low",
        "relative_iod_ci_high",
    }
    missing = sorted(required - set(released.columns))
    if missing:
        errors.append(f"[S22_FIELDS] missing fields: {', '.join(missing)}")
        return
    if released["species"].duplicated().any():
        errors.append("[S22_SPECIES] species rows are not unique")
        return

    expected = objects.groupby("species", sort=False).agg(
        n_nuclei=("nucleus_object_id", "size"),
        n_images=("filename", "nunique"),
        n_specimens=("specimen_id", "nunique"),
    )
    expected["iod_equal_image_estimate"] = images.groupby("species")["image_median_iod"].mean()
    anchor = float(expected["iod_equal_image_estimate"].median())
    expected["relative_iod_anchor"] = anchor
    expected["relative_iod_index"] = expected["iod_equal_image_estimate"] / anchor
    observed = released.set_index("species")
    if len(observed) != 24 or set(observed.index) != set(expected.index):
        errors.append("[S22_SPECIES] release must contain the exact 24 S21 species")
        return
    path24 = pd.read_csv(root / "data/identity/path24_panel.csv")
    crosswalk = pd.read_csv(root / "data/identity/species_taxonomy_crosswalk.csv")
    display_to_current = {
        f"D. {species}": str(species) for species in crosswalk["current_species"].astype(str)
    }
    unknown_labels = sorted(set(observed.index) - set(display_to_current))
    observed_path24 = {
        display_to_current[label] for label in observed.index if label in display_to_current
    }
    expected_path24 = set(path24["species"].astype(str))
    if unknown_labels or observed_path24 != expected_path24:
        errors.append("[S22_PATH24] S22 species do not map to the exact Path24 contract")

    for column in [
        "n_nuclei",
        "n_images",
        "n_specimens",
        "iod_equal_image_estimate",
        "relative_iod_anchor",
        "relative_iod_index",
    ]:
        mismatches = _numeric_mismatches(observed[column], expected[column])
        if mismatches:
            errors.append(f"[S22_{column.upper()}] mismatch: {', '.join(mismatches)}")

    expected_relative_low = observed["iod_equal_image_ci_low"] / anchor
    expected_relative_high = observed["iod_equal_image_ci_high"] / anchor
    for column, expected_boundary in [
        ("relative_iod_ci_low", expected_relative_low),
        ("relative_iod_ci_high", expected_relative_high),
    ]:
        mismatches = _numeric_mismatches(observed[column], expected_boundary)
        if mismatches:
            errors.append(f"[S22_{column.upper()}] mismatch: {', '.join(mismatches)}")

    iod_bounds_ok = observed["iod_equal_image_ci_low"].le(
        observed["iod_equal_image_estimate"]
    ) & observed["iod_equal_image_estimate"].le(observed["iod_equal_image_ci_high"])
    relative_bounds_ok = observed["relative_iod_ci_low"].le(
        observed["relative_iod_index"]
    ) & observed["relative_iod_index"].le(observed["relative_iod_ci_high"])
    invalid_bounds = observed.index[~(iod_bounds_ok & relative_bounds_ok)].tolist()
    if invalid_bounds:
        errors.append(
            "[S22_BOOTSTRAP_BOUNDS] point estimate outside interval: "
            + ", ".join(map(str, invalid_bounds))
        )


def _validate_s22_conditional_calibration(root: Path, errors: list[str]) -> None:
    released = _read(root, "relative_nuclear_iod_species.csv")
    calibration_columns = [
        "iod_ratio_to_fuscus",
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
        "genome_size_reference_species",
        "genome_size_reference_pg",
        "genome_size_calibration_pg_per_iod",
        "genome_size_rank",
    ]
    missing = sorted(set(calibration_columns) - set(released.columns))
    if missing:
        errors.append(f"[S22_CALIBRATION_FIELDS] missing fields: {', '.join(missing)}")
        return
    try:
        expected = reconstruct_s22_calibration(root)
    except (FileNotFoundError, OSError, TypeError, ValueError, KeyError) as error:
        errors.append(f"[S22_CALIBRATION_PROVENANCE] {type(error).__name__}: {error}")
        return
    observed = released.set_index("species")
    expected = expected.set_index("species")
    if set(observed.index) != set(expected.index):
        errors.append("[S22_CALIBRATION_SPECIES] calibrated species set changed")
        return
    for column in calibration_columns:
        if column == "genome_size_reference_species":
            mismatches = observed.index[
                observed[column].astype(str).ne(expected[column].astype(str))
            ].tolist()
        else:
            mismatches = _numeric_mismatches(
                observed[column],
                expected[column],
                atol=5.1e-9,
            )
        if mismatches:
            errors.append(f"[S22_{column.upper()}] mismatch: " + ", ".join(mismatches))

    fuscus = observed.loc["D. fuscus"]
    if not (
        np.isclose(float(fuscus["iod_ratio_to_fuscus"]), 1.0)
        and np.isclose(
            float(fuscus["genome_size_pg_fuscus_anchored"]),
            float(fuscus["genome_size_reference_pg"]),
            atol=5.1e-9,
        )
        and np.isclose(
            float(fuscus["genome_size_pg_ci_low"]),
            float(fuscus["genome_size_reference_pg"]),
            atol=5.1e-9,
        )
        and np.isclose(
            float(fuscus["genome_size_pg_ci_high"]),
            float(fuscus["genome_size_reference_pg"]),
            atol=5.1e-9,
        )
    ):
        errors.append("[S22_FUSCUS_ANCHOR] fuscus point and ratio interval must equal anchor")

    projected = project_s22_availability(root, released)
    stored_availability = pd.read_csv(root / AVAILABILITY_PATH)
    target = "step06_genome_size_pg_fuscus_anchored"
    observed_availability = stored_availability.set_index("species")[target]
    expected_availability = projected.set_index("species")[target]
    availability_mismatches = _numeric_mismatches(
        observed_availability.dropna(),
        expected_availability.dropna(),
        atol=5.1e-9,
    )
    missingness_changed = not observed_availability.isna().equals(expected_availability.isna())
    if availability_mismatches or missingness_changed:
        errors.append(
            "[S01_S22_CONDITIONAL_SCALE] availability projection differs from S22: "
            + ", ".join(availability_mismatches)
        )


def _validate_pixel_scale(root: Path, errors: list[str]) -> None:
    for release_id, filename in [
        ("S18", "cell_nucleus_objects.csv"),
        ("S21", "nuclear_iod_objects.csv"),
    ]:
        objects = _read(root, filename)
        for object_type, pixels, physical_area in [
            ("CELL", "cell_area_px", "cell_area_um2"),
            ("NUCLEUS", "nuc_area_px", "nuc_area_um2"),
        ]:
            agrees = np.isclose(
                objects[physical_area].astype(float),
                objects[pixels].astype(float) * PIXEL_AREA_UM2,
                rtol=0,
                atol=1e-8,
            )
            mismatch_count = int((~agrees).sum())
            if mismatch_count:
                errors.append(
                    f"[PIXEL_SCALE_{release_id}_{object_type}] "
                    f"{mismatch_count} rows differ from {PIXEL_AREA_UM2} um2/pixel"
                )


def _validate_s24(root: Path, errors: list[str]) -> None:
    balance = _read(root, "nuclear_iod_quality_balance.csv")
    required = {
        "species",
        "standardized_mean_difference_match_log_edge_sharpness",
        "standardized_mean_difference_match_log_relative_ring_noise",
        "ks_distance_match_log_edge_sharpness",
        "ks_distance_match_log_relative_ring_noise",
        "max_abs_standardized_mean_difference",
        "max_ks_distance",
        "passes_prespecified_balance_thresholds",
    }
    missing = sorted(required - set(balance.columns))
    if missing:
        errors.append(f"[S24_FIELDS] missing fields: {', '.join(missing)}")
        return
    if len(balance) != 24 or balance["species"].duplicated().any():
        errors.append("[S24_SPECIES] balance table must contain 24 unique species")
        return

    max_abs_smd = (
        balance[
            [
                "standardized_mean_difference_match_log_edge_sharpness",
                "standardized_mean_difference_match_log_relative_ring_noise",
            ]
        ]
        .abs()
        .max(axis=1)
    )
    max_ks = balance[
        [
            "ks_distance_match_log_edge_sharpness",
            "ks_distance_match_log_relative_ring_noise",
        ]
    ].max(axis=1)
    for tag, observed, expected in [
        (
            "S24_MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE",
            balance["max_abs_standardized_mean_difference"],
            max_abs_smd,
        ),
        ("S24_MAX_KS_DISTANCE", balance["max_ks_distance"], max_ks),
    ]:
        mismatch_count = int((~np.isclose(observed, expected, rtol=1e-10, atol=1e-8)).sum())
        if mismatch_count:
            errors.append(f"[{tag}] {mismatch_count} stored maxima disagree")

    expected_pass = max_abs_smd.le(MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE) & max_ks.le(
        MAX_KS_DISTANCE
    )
    observed_pass = balance["passes_prespecified_balance_thresholds"].astype(bool)
    mismatch_species = balance.loc[observed_pass.ne(expected_pass), "species"].tolist()
    if mismatch_species:
        errors.append("[S24_BALANCE_BOOLEAN] mismatch: " + ", ".join(map(str, mismatch_species)))
    stored_failures = set(balance.loc[~observed_pass, "species"].astype(str))
    if stored_failures != EXPECTED_BALANCE_FAILURES:
        errors.append(
            "[S24_FAILURE_SET] expected exactly: " + ", ".join(sorted(EXPECTED_BALANCE_FAILURES))
        )


def _validate_animal_support(root: Path, errors: list[str]) -> None:
    morphology = _read(root, "cell_nucleus_objects.csv")
    iod = _read(root, "nuclear_iod_objects.csv")
    morphology_animals = set(morphology["specimen_id"].dropna().astype(str))
    iod_animals = set(iod["specimen_id"].dropna().astype(str))
    if len(morphology_animals) != 42:
        errors.append(f"[MORPHOLOGY_ANIMALS] expected 42, observed {len(morphology_animals)}")
    if len(iod_animals) != 50:
        errors.append(f"[IOD_ANIMALS] expected 50, observed {len(iod_animals)}")
    union_size = len(morphology_animals | iod_animals)
    if union_size != 51:
        errors.append(f"[ANIMAL_UNION] expected 51, observed {union_size}")

    specimens_per_species = morphology.groupby("species")["specimen_id"].nunique()
    single_animal_species = frozenset(
        specimens_per_species.index[specimens_per_species.eq(1)].astype(str)
    )
    if single_animal_species != EXPECTED_MORPHOLOGY_SINGLE_ANIMAL_SPECIES:
        errors.append(
            "[MORPHOLOGY_SINGLE_ANIMAL_SPECIES] expected exactly: "
            + ", ".join(sorted(EXPECTED_MORPHOLOGY_SINGLE_ANIMAL_SPECIES))
        )


def validate_release(root: Path | str = ROOT) -> tuple[str, ...]:
    """Return every semantic mismatch in the compact morphology release."""

    resolved_root = Path(root).resolve()
    errors: list[str] = []
    _validate_s19(resolved_root, errors)
    _validate_s23(resolved_root, errors)
    _validate_s22(resolved_root, errors)
    _validate_s22_conditional_calibration(resolved_root, errors)
    _validate_pixel_scale(resolved_root, errors)
    _validate_s24(resolved_root, errors)
    _validate_animal_support(resolved_root, errors)
    return tuple(errors)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument(
        "--write-calibration",
        action="store_true",
        help=(
            "Promote the provenance-derived conditional S22 calibration and its "
            "S01 audit projection before validating the complete release."
        ),
    )
    args = parser.parse_args()
    if args.write_calibration:
        reconstructed = reconstruct_s22_calibration(args.root)
        projected = project_s22_availability(args.root, reconstructed)
        reconstructed.to_csv(
            args.root / DATA_DIR / "relative_nuclear_iod_species.csv",
            index=False,
            float_format="%.8f",
        )
        projected.to_csv(args.root / AVAILABILITY_PATH, index=False)
    errors = validate_release(args.root)
    if errors:
        print("Morphology release validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    action = "promoted and validated" if args.write_calibration else "validated"
    print(
        "PASS: morphology release replay agrees with frozen canonical tables; "
        f"conditional 16.1-Gbp assembly anchor {action}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Replay the compact morphology release from its frozen object tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path("analyses/02_morphology/data")
PIXEL_AREA_UM2 = 0.0144
MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE = 0.10
MAX_KS_DISTANCE = 0.25
EXPECTED_BALANCE_FAILURES = frozenset(
    {"D. aeneus", "D. ochrophaeus", "D. orestes", "D. wrighti"}
)
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

    support = (
        objects.groupby("species", sort=False)
        .agg(
            n_size_cells=("nucleus_object_id", "size"),
            n_size_images=("filename", "nunique"),
            n_size_specimens=("specimen_id", "nunique"),
            cell_area_um2=("cell_area_um2", "median"),
            nucleus_area_um2=("nuc_area_um2", "median"),
        )
    )
    specimen_counts = objects.groupby(["species", "specimen_id"], sort=False).size()
    support["largest_specimen_n"] = specimen_counts.groupby(level="species").max()
    support["largest_specimen_fraction"] = (
        support["largest_specimen_n"] / support["n_size_cells"]
    )
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
    expected["iod_equal_image_estimate"] = images.groupby("species")[
        "image_median_iod"
    ].mean()
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

    iod_bounds_ok = (
        observed["iod_equal_image_ci_low"].le(observed["iod_equal_image_estimate"])
        & observed["iod_equal_image_estimate"].le(observed["iod_equal_image_ci_high"])
    )
    relative_bounds_ok = (
        observed["relative_iod_ci_low"].le(observed["relative_iod_index"])
        & observed["relative_iod_index"].le(observed["relative_iod_ci_high"])
    )
    invalid_bounds = observed.index[~(iod_bounds_ok & relative_bounds_ok)].tolist()
    if invalid_bounds:
        errors.append(
            "[S22_BOOTSTRAP_BOUNDS] point estimate outside interval: "
            + ", ".join(map(str, invalid_bounds))
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

    max_abs_smd = balance[
        [
            "standardized_mean_difference_match_log_edge_sharpness",
            "standardized_mean_difference_match_log_relative_ring_noise",
        ]
    ].abs().max(axis=1)
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
        mismatch_count = int(
            (~np.isclose(observed, expected, rtol=1e-10, atol=1e-8)).sum()
        )
        if mismatch_count:
            errors.append(f"[{tag}] {mismatch_count} stored maxima disagree")

    expected_pass = max_abs_smd.le(MAX_ABS_STANDARDIZED_MEAN_DIFFERENCE) & max_ks.le(
        MAX_KS_DISTANCE
    )
    observed_pass = balance["passes_prespecified_balance_thresholds"].astype(bool)
    mismatch_species = balance.loc[observed_pass.ne(expected_pass), "species"].tolist()
    if mismatch_species:
        errors.append(
            "[S24_BALANCE_BOOLEAN] mismatch: " + ", ".join(map(str, mismatch_species))
        )
    stored_failures = set(balance.loc[~observed_pass, "species"].astype(str))
    if stored_failures != EXPECTED_BALANCE_FAILURES:
        errors.append(
            "[S24_FAILURE_SET] expected exactly: "
            + ", ".join(sorted(EXPECTED_BALANCE_FAILURES))
        )


def _validate_animal_support(root: Path, errors: list[str]) -> None:
    morphology = _read(root, "cell_nucleus_objects.csv")
    iod = _read(root, "nuclear_iod_objects.csv")
    morphology_animals = set(morphology["specimen_id"].dropna().astype(str))
    iod_animals = set(iod["specimen_id"].dropna().astype(str))
    if len(morphology_animals) != 42:
        errors.append(
            f"[MORPHOLOGY_ANIMALS] expected 42, observed {len(morphology_animals)}"
        )
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
    _validate_pixel_scale(resolved_root, errors)
    _validate_s24(resolved_root, errors)
    _validate_animal_support(resolved_root, errors)
    return tuple(errors)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    errors = validate_release(args.root)
    if errors:
        print("Morphology release validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print("PASS: morphology release replay agrees with frozen canonical tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

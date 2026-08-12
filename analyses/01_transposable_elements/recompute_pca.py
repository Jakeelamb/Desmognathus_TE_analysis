"""Reconstruct and verify the canonical TE34 diversity and CLR-PCA release."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

ANALYSIS = Path(__file__).resolve().parent
ROOT = ANALYSIS.parents[1]
DATA = ANALYSIS / "data"
PANEL = ROOT / "data" / "identity" / "te34_panel.csv"
COMPOSITION_MODE = "classified_conditional"
EXCLUSION_REASON = "not_positive_in_all_retry_averaged_te34_species"
EXPECTED_EXCLUDED_SUPERFAMILIES = {"CR1", "Chapaev", "Dada", "Ginger", "Merlin"}
DIVERSITY_COLUMNS = [
    "species",
    "te_level",
    "observed_richness",
    "shannon_entropy",
    "gini_simpson",
]


def fail(code: str, detail: str) -> None:
    raise AssertionError(f"{code}: {detail}")


def read_panel() -> list[str]:
    panel = pd.read_csv(PANEL)
    species = panel.loc[panel["te_resource_panel34_v1"], "species"].tolist()
    if len(species) != 34 or len(species) != len(set(species)):
        fail("PCA_PANEL_MISMATCH", "TE34 must contain exactly 34 unique species")
    return species


def read_compositions(panel_species: list[str]) -> dict[str, pd.DataFrame]:
    compositions: dict[str, pd.DataFrame] = {}
    for level, filename in (
        ("order", "te_order_composition.csv"),
        ("superfamily", "te_superfamily_composition.csv"),
    ):
        frame = pd.read_csv(DATA / filename)
        species = frame.pop("species").tolist()
        if species != panel_species:
            fail(
                "PCA_PANEL_MISMATCH",
                f"{filename} species or ordering differs from canonical TE34",
            )
        table = frame.astype(float)
        table.index = pd.Index(species, name="species")
        values = table.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            fail("PCA_COMPOSITION_INVALID", f"{filename} is not finite and nonnegative")
        if not np.allclose(values.sum(axis=1), 1, rtol=0, atol=5e-12):
            fail("PCA_COMPOSITION_INVALID", f"{filename} rows are not closed to one")
        compositions[level] = table
    return compositions


def read_mass_accounting(panel_species: list[str]) -> pd.DataFrame:
    mass = pd.read_csv(DATA / "dnapipete_mass_accounting.csv")
    if mass["species"].tolist() != panel_species:
        fail(
            "DIVERSITY_PANEL_MISMATCH",
            "dnapipete_mass_accounting.csv species or ordering differs from canonical TE34",
        )

    for level in ("order", "superfamily"):
        retained_fraction = mass[f"{level}_retained_fraction"].to_numpy(dtype=float)
        unresolved_fraction = mass[f"{level}_unresolved_fraction"].to_numpy(dtype=float)
        retained_bases = mass[f"{level}_retained_aligned_bases"].to_numpy(dtype=float)
        unresolved_bases = mass[f"{level}_unresolved_aligned_bases"].to_numpy(dtype=float)
        total_bases = mass["total_aligned_bases"].to_numpy(dtype=float)
        values = np.column_stack(
            (
                retained_fraction,
                unresolved_fraction,
                retained_bases,
                unresolved_bases,
                total_bases,
            )
        )
        if not np.isfinite(values).all() or (values < 0).any():
            fail(
                "DIVERSITY_MASS_INVALID",
                f"{level} mass accounting is not finite and nonnegative",
            )
        if (retained_fraction > 1).any() or (unresolved_fraction > 1).any():
            fail("DIVERSITY_MASS_INVALID", f"{level} mass fractions exceed one")
        if not np.allclose(
            retained_fraction + unresolved_fraction,
            1,
            rtol=0,
            atol=5e-12,
        ):
            fail("DIVERSITY_MASS_INVALID", f"{level} mass fractions do not close to one")
        if not np.allclose(
            retained_bases + unresolved_bases,
            total_bases,
            rtol=1e-12,
            atol=1e-3,
        ):
            fail("DIVERSITY_MASS_INVALID", f"{level} aligned-base masses do not balance")
    return mass


def diversity_metrics(probabilities: np.ndarray | pd.Series) -> dict[str, float | int]:
    values = np.asarray(probabilities, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all() or (values < 0).any():
        fail(
            "DIVERSITY_COMPOSITION_INVALID",
            "diversity probabilities must be a finite, nonnegative vector",
        )
    if not np.isclose(values.sum(), 1, rtol=0, atol=5e-12):
        fail("DIVERSITY_COMPOSITION_INVALID", "diversity probabilities do not close to one")

    positive = values[values > 0]
    observed_richness = len(positive)
    shannon_entropy = float(-np.sum(positive * np.log(positive)))
    sum_squared = float(np.sum(positive**2))
    return {
        "observed_richness": observed_richness,
        "shannon_entropy": shannon_entropy,
        "gini_simpson": 1 - sum_squared,
    }


def build_diversity(compositions: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for level, composition in compositions.items():
        for species, probabilities in composition.iterrows():
            rows.append(
                {
                    "species": species,
                    "te_level": level,
                    **diversity_metrics(probabilities),
                }
            )

    diversity = pd.DataFrame(rows, columns=DIVERSITY_COLUMNS)
    key_columns = ["species", "te_level"]
    if diversity.duplicated(key_columns).any():
        fail("DIVERSITY_RELEASE_INVALID", "diversity keys are not unique")
    strata = diversity.groupby("te_level", sort=False).size().to_dict()
    expected_strata = {level: 34 for level in compositions}
    if strata != expected_strata or len(diversity) != 68:
        fail("DIVERSITY_RELEASE_INVALID", f"unexpected diversity strata: {strata}")
    return diversity


def build_prevalence(compositions: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for level, table in compositions.items():
        for position, feature in enumerate(table.columns):
            values = table[feature]
            positive = values > 0
            n_positive = int(positive.sum())
            included = n_positive == len(table)
            rows.append(
                {
                    "te_level": level,
                    "feature": feature,
                    "n_species_total": len(table),
                    "n_species_positive": n_positive,
                    "prevalence_fraction": n_positive / len(table),
                    "absent_species": json.dumps(table.index[~positive].tolist()),
                    "minimum_proportion": values.min(),
                    "median_proportion": values.median(),
                    "mean_proportion": values.mean(),
                    "maximum_proportion": values.max(),
                    "pca_included": included,
                    "pca_exclusion_reason": None if included else EXCLUSION_REASON,
                    "_position": position,
                }
            )
    prevalence = pd.DataFrame(rows).sort_values(
        ["n_species_positive", "te_level", "_position"], kind="stable"
    )
    return prevalence.drop(columns="_position").reset_index(drop=True)


def retained_features(prevalence: pd.DataFrame, level: str) -> list[str]:
    return prevalence.loc[
        prevalence["te_level"].eq(level) & prevalence["pca_included"], "feature"
    ].tolist()


def build_clr(
    compositions: dict[str, pd.DataFrame], prevalence: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    blocks: list[pd.DataFrame] = []
    clr_tables: dict[str, pd.DataFrame] = {}
    for level, table in compositions.items():
        features = retained_features(prevalence, level)
        closed = table.loc[:, features].div(table.loc[:, features].sum(axis=1), axis=0)
        if closed.le(0).any().any():
            fail("PCA_FEATURE_CONTRACT", f"retained {level} features contain zero values")
        logged = np.log(closed)
        clr = logged.sub(logged.mean(axis=1), axis=0)
        if not np.allclose(clr.sum(axis=1), 0, rtol=0, atol=2e-12):
            fail("PCA_GEOMETRY_INVALID", f"{level} CLR rows do not sum to zero")
        clr_tables[level] = clr
        blocks.append(
            pd.DataFrame(
                {
                    "species": np.repeat(closed.index.to_numpy(), len(features)),
                    "te_level": level,
                    "composition_mode": COMPOSITION_MODE,
                    "feature": np.tile(np.asarray(features), len(closed)),
                    "closed_proportion": closed.to_numpy().ravel(),
                    "clr_value": clr.to_numpy().ravel(),
                }
            )
        )
    return pd.concat(blocks, ignore_index=True), clr_tables


def fit_pca(clr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    centered = clr.to_numpy() - clr.to_numpy().mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    n_components = min(len(clr) - 1, clr.shape[1] - 1)
    columns = [f"PC{number}" for number in range(1, n_components + 1)]
    scores = pd.DataFrame(
        u[:, :n_components] * singular[:n_components],
        index=clr.index,
        columns=columns,
    )
    loadings = pd.DataFrame(
        vt[:n_components].T,
        index=clr.columns,
        columns=columns,
    )
    for column in columns:
        anchor = loadings[column].abs().idxmax()
        if loadings.loc[anchor, column] < 0:
            loadings[column] *= -1
            scores[column] *= -1

    eigenvalues = singular[:n_components] ** 2 / (len(clr) - 1)
    variance = pd.DataFrame(
        {
            "PC": np.arange(1, n_components + 1, dtype=int),
            "eigenvalue": eigenvalues,
            "variance_explained": eigenvalues / eigenvalues.sum(),
        }
    )
    variance["cumulative_variance_explained"] = variance["variance_explained"].cumsum()

    reconstruction = scores.to_numpy() @ loadings.to_numpy().T
    if not np.allclose(reconstruction, centered, rtol=1e-10, atol=1e-11):
        fail("PCA_GEOMETRY_INVALID", "scores and loadings do not reconstruct centered CLR")
    if not np.allclose(
        loadings.to_numpy().T @ loadings.to_numpy(),
        np.eye(n_components),
        rtol=1e-10,
        atol=1e-11,
    ):
        fail("PCA_GEOMETRY_INVALID", "released loading axes are not orthonormal")
    clr_distances = np.linalg.norm(clr.to_numpy()[:, None] - clr.to_numpy()[None, :], axis=2)
    score_distances = np.linalg.norm(
        scores.to_numpy()[:, None] - scores.to_numpy()[None, :], axis=2
    )
    if not np.allclose(clr_distances, score_distances, rtol=1e-10, atol=1e-11):
        fail("PCA_GEOMETRY_INVALID", "all PC axes do not preserve CLR distances")

    scores = scores.rename_axis("species").reset_index()
    scores.insert(1, "te_level", "superfamily")
    scores.insert(2, "composition_mode", COMPOSITION_MODE)
    loadings = loadings.rename_axis("superfamily").reset_index()
    return scores, variance, loadings


def reconstruct_release() -> dict[str, pd.DataFrame]:
    panel_species = read_panel()
    compositions = read_compositions(panel_species)
    read_mass_accounting(panel_species)
    diversity = build_diversity(compositions)
    prevalence = build_prevalence(compositions)
    included_counts = (
        prevalence.loc[prevalence["pca_included"]].groupby("te_level").size().to_dict()
    )
    excluded_superfamilies = set(
        prevalence.loc[
            prevalence["te_level"].eq("superfamily") & ~prevalence["pca_included"],
            "feature",
        ]
    )
    if included_counts != {"order": 10, "superfamily": 24}:
        fail("PCA_FEATURE_CONTRACT", f"unexpected retained feature counts: {included_counts}")
    if excluded_superfamilies != EXPECTED_EXCLUDED_SUPERFAMILIES:
        fail(
            "PCA_FEATURE_CONTRACT",
            f"unexpected excluded superfamilies: {sorted(excluded_superfamilies)}",
        )

    clr_long, clr_tables = build_clr(compositions, prevalence)
    scores, variance, loadings = fit_pca(clr_tables["superfamily"])
    return {
        "diversity": diversity,
        "feature prevalence": prevalence,
        "complete CLR": clr_long,
        "scores": scores,
        "variance": variance,
        "loadings": loadings,
    }


def check_frame(
    actual: pd.DataFrame,
    expected_path: Path,
    key: str,
    error_code: str = "PCA_RELEASE_DRIFT",
) -> None:
    expected = pd.read_csv(expected_path)
    if list(actual.columns) != list(expected.columns) or len(actual) != len(expected):
        fail(
            error_code,
            f"{key} schema/row count differs from {expected_path.name}",
        )
    for column in actual.columns:
        if is_numeric_dtype(expected[column]) and not is_bool_dtype(expected[column]):
            try:
                np.testing.assert_allclose(
                    actual[column].to_numpy(dtype=float),
                    expected[column].to_numpy(dtype=float),
                    rtol=1e-9,
                    atol=5e-12,
                    equal_nan=True,
                )
            except AssertionError as error:
                fail(error_code, f"{key}.{column}: {error}")
        else:
            actual_values = actual[column].fillna("<NA>").astype(str).tolist()
            expected_values = expected[column].fillna("<NA>").astype(str).tolist()
            if actual_values != expected_values:
                first = next(
                    index
                    for index, values in enumerate(zip(actual_values, expected_values, strict=True))
                    if values[0] != values[1]
                )
                fail(
                    error_code,
                    f"{key}.{column} first differs at row {first}: "
                    f"{actual_values[first]!r} != {expected_values[first]!r}",
                )


def check_release(products: dict[str, pd.DataFrame]) -> None:
    check_frame(
        products["diversity"],
        DATA / "te_diversity.csv",
        "diversity",
        error_code="DIVERSITY_RELEASE_DRIFT",
    )
    canonical = {
        "feature prevalence": DATA / "te_feature_prevalence.csv",
        "complete CLR": DATA / "te_pca_clr_matrix.csv",
        "scores": DATA / "te_pca_scores.csv",
        "variance": DATA / "te_pca_variance.csv",
        "loadings": DATA / "te_pca_loadings.csv",
    }
    for key, expected_path in canonical.items():
        check_frame(products[key], expected_path, key)


def main() -> None:
    products = reconstruct_release()
    check_release(products)
    print(
        "PASS: compact TE34 inputs reproduce feature prevalence, complete CLR, "
        "and final superfamily PCA; all 68 canonical diversity rows also reproduce"
    )


if __name__ == "__main__":
    main()

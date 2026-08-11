#!/usr/bin/env python3
"""Rebuild and verify the compact LTR terminal:internal release.

The expensive mapping and depth-generation workflow is frozen upstream. This
script starts from the non-destructive historical element audit table, corrects
the legacy TEsorter domain-token counting error, enforces the LTR/Gypsy >=5
domain gate, reapplies the within-species two-sided Tukey IQR filter to that
eligible cohort, and reconstructs the species and resource summaries used by
the paper. By default it is read-only; ``--write`` is the deliberate
release-promotion operation.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "analyses/01_transposable_elements/data"
ELEMENT_PATH = DATA / "ltr_element_metrics.csv"
SPECIES_PATH = DATA / "ltr_species_robustness.csv"
RESOURCE_PATH = DATA / "ltr_resource_coverage.csv"
EXCLUSION_PATH = DATA / "ltr_excluded_elements.csv"
AVAILABILITY_PATH = ROOT / "data/identity/analysis_availability.csv"

RATIO_COLUMN = "ratio_terminal_internal_all_positions"
DOMAIN_GATE_RULE_ID = "ltr_gypsy_at_least_5_tesorter_domains_v1"
IQR_RULE_ID = "within_species_two_sided_tukey_1.5_iqr_v1"
PRIMARY_BRANCH = "primary_iqr_filtered_coverage_ge_80pct"
BOOTSTRAP_REPLICATES = 2_000
MIN_DOMAIN_ANNOTATIONS = 5

DOMAIN_DERIVED_COLUMNS = [
    "legacy_pipe_split_domain_count",
    "domain_gate_keep",
    "domain_gate_rule_id",
    "domain_gate_reason",
]

FILTER_COLUMNS = [
    "species_ratio_q1",
    "species_ratio_q3",
    "species_ratio_iqr",
    "species_ratio_iqr_lower_bound",
    "species_ratio_iqr_upper_bound",
    "iqr_filter_keep",
    "iqr_filter_direction",
    "iqr_filter_rule_id",
    "iqr_filter_reason",
    "primary_ltr_element",
]
RESOURCE_REBUILT_COLUMNS = [
    "resource_status",
    "n_selected_ltr_gypsy_ge5_domain_elements",
    "n_usable_elements",
    "n_selected_source_corrupt_exclusions",
    "n_coverage_ge_80pct_elements",
    "n_iqr_retained_elements",
    "n_primary_iqr_filtered_elements",
]

EXPECTED_RELEASE_COUNTS = {
    "historical_usable": 1_086,
    "eligible": 408,
    "retained": 381,
    "removed": 27,
    "lower": 1,
    "upper": 26,
    "coverage": 407,
    "primary": 380,
}


def _fail(message: str) -> None:
    raise ValueError(message)


def _read_csv(path: Path) -> pd.DataFrame:
    """Read canonical CSV floats without decimal-to-binary representation drift."""

    return pd.read_csv(path, low_memory=False, float_precision="round_trip")


def _tesorter_domain_tokens(value: object) -> tuple[str, ...]:
    """Return TEsorter DOMAIN|MODEL annotations from the cls.tsv Domains field."""

    if pd.isna(value) or not str(value).strip():
        _fail("Missing TEsorter Domains value")
    if str(value).strip().lower() == "none":
        return ()
    tokens = tuple(str(value).split())
    malformed = [
        token
        for token in tokens
        if token.count("|") != 1 or not all(token.split("|", maxsplit=1))
    ]
    if malformed:
        _fail(f"Malformed TEsorter DOMAIN|MODEL annotation tokens: {malformed}")
    return tokens


def _tesorter_domain_count(value: object) -> int:
    """Count whitespace-delimited TEsorter DOMAIN|MODEL annotations."""

    return len(_tesorter_domain_tokens(value))


def add_domain_gate(elements: pd.DataFrame) -> pd.DataFrame:
    """Correct the historical delimiter error and mark the declared LTR cohort."""

    base = elements.drop(columns=DOMAIN_DERIVED_COLUMNS, errors="ignore").copy()
    required = {
        "Domains",
        "Order",
        "Superfamily",
        "domain_count",
        "element length",
        "tesorter_annotation_found",
    }
    missing = sorted(required.difference(base.columns))
    if missing:
        _fail(f"LTR element table lacks domain-gate columns: {missing}")

    actual_count = base["Domains"].map(_tesorter_domain_count).astype(int)
    legacy_count = base["Domains"].map(
        lambda value: 0
        if pd.isna(value)
        or str(value).strip().lower() == "none"
        else len(str(value).split("|"))
    )
    eligible = (
        base["element length"].ge(3_000)
        & base["tesorter_annotation_found"].astype(bool)
        & base["Order"].eq("LTR")
        & base["Superfamily"].eq("Gypsy")
        & actual_count.ge(MIN_DOMAIN_ANNOTATIONS)
    )

    result = base.copy()
    result["legacy_pipe_split_domain_count"] = legacy_count.astype(int)
    result["domain_count"] = actual_count
    result["domain_gate_keep"] = eligible
    result["domain_gate_rule_id"] = DOMAIN_GATE_RULE_ID
    result["domain_gate_reason"] = np.select(
        [
            ~result["element length"].ge(3_000),
            ~result["tesorter_annotation_found"].astype(bool),
            ~result["Order"].eq("LTR"),
            ~result["Superfamily"].eq("Gypsy"),
            actual_count.lt(MIN_DOMAIN_ANNOTATIONS),
        ],
        [
            "element_shorter_than_3000bp",
            "missing_exact_tesorter_annotation",
            "not_ltr_order",
            "not_gypsy_superfamily",
            "fewer_than_5_tesorter_domain_annotations",
        ],
        default="eligible",
    )

    if len(result) != EXPECTED_RELEASE_COUNTS["historical_usable"]:
        _fail(f"Expected 1,086 historical usable LTR rows, found {len(result)}")
    if int(eligible.sum()) != EXPECTED_RELEASE_COUNTS["eligible"]:
        _fail(f"Corrected LTR domain gate retained {int(eligible.sum())}, expected 408")
    if result.loc[eligible, "species"].nunique() != 30:
        _fail("Corrected LTR domain gate must retain all 30 represented species")
    if actual_count.max() != 5:
        _fail(f"Frozen TEsorter domain-count maximum changed: {actual_count.max()}")
    return result


def add_iqr_filter(elements: pd.DataFrame) -> pd.DataFrame:
    """Apply the recovered two-sided 1.5-IQR rule within the eligible cohort."""

    base = elements.drop(columns=FILTER_COLUMNS, errors="ignore").copy()
    required = {
        "species",
        "element_id",
        RATIO_COLUMN,
        "domain_gate_keep",
        "coverage_ge_80pct",
        "left_ltr_positive_coverage_fraction",
        "right_ltr_positive_coverage_fraction",
        "internal_positive_coverage_fraction",
    }
    missing = sorted(required.difference(base.columns))
    if missing:
        _fail(f"LTR element table lacks required columns: {missing}")
    if base["element_id"].duplicated().any():
        _fail("LTR element table contains duplicate element_id values")
    if len(base) != EXPECTED_RELEASE_COUNTS["historical_usable"]:
        _fail(f"Expected 1,086 usable LTR elements, found {len(base)}")

    ratio = pd.to_numeric(base[RATIO_COLUMN], errors="raise")
    if not np.isfinite(ratio).all() or not ratio.gt(0).all():
        _fail("LTR ratios must be finite and strictly positive")

    eligible = base["domain_gate_keep"].astype(bool)
    q1 = pd.Series(np.nan, index=base.index, dtype=float)
    q3 = pd.Series(np.nan, index=base.index, dtype=float)
    eligible_ratios = base.loc[eligible].assign(_ratio=ratio.loc[eligible])
    grouped = eligible_ratios.groupby("species")["_ratio"]
    q1.loc[eligible] = grouped.transform(lambda values: values.quantile(0.25))
    q3.loc[eligible] = grouped.transform(lambda values: values.quantile(0.75))
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    below = eligible & ratio.lt(lower)
    above = eligible & ratio.gt(upper)
    keep = eligible & ~(below | above)

    expected_coverage = (
        base["left_ltr_positive_coverage_fraction"].ge(0.8)
        & base["right_ltr_positive_coverage_fraction"].ge(0.8)
        & base["internal_positive_coverage_fraction"].ge(0.8)
    )
    if not expected_coverage.equals(base["coverage_ge_80pct"].astype(bool)):
        _fail("Stored 80% regional positive-depth gate does not reproduce")

    result = base.copy()
    result["species_ratio_q1"] = q1
    result["species_ratio_q3"] = q3
    result["species_ratio_iqr"] = iqr
    result["species_ratio_iqr_lower_bound"] = lower
    result["species_ratio_iqr_upper_bound"] = upper
    result["iqr_filter_keep"] = keep
    result["iqr_filter_direction"] = np.select(
        [~eligible, below, above],
        ["not_eligible", "lower_outlier", "upper_outlier"],
        default="retained",
    )
    result["iqr_filter_rule_id"] = IQR_RULE_ID
    result["iqr_filter_reason"] = np.select(
        [~eligible, below, above],
        [
            "not_in_ltr_gypsy_ge5_domain_cohort",
            "below_q1_minus_1.5_iqr",
            "above_q3_plus_1.5_iqr",
        ],
        default="retained",
    )
    result["primary_ltr_element"] = keep & expected_coverage

    observed = {
        "historical_usable": len(result),
        "eligible": int(eligible.sum()),
        "retained": int(keep.sum()),
        "removed": int((eligible & ~keep).sum()),
        "lower": int(below.sum()),
        "upper": int(above.sum()),
        "coverage": int((eligible & expected_coverage).sum()),
        "primary": int(result["primary_ltr_element"].sum()),
    }
    if observed != EXPECTED_RELEASE_COUNTS:
        _fail(f"LTR IQR release counts changed: {observed}")
    if result.loc[result["primary_ltr_element"], "species"].nunique() != 30:
        _fail("Primary LTR branch must retain all 30 represented species")
    return result


def _trimmed_mean(values: np.ndarray, fraction: float = 0.1) -> float:
    ordered = np.sort(values)
    trim = int(np.floor(len(ordered) * fraction))
    kept = ordered[trim : len(ordered) - trim] if trim else ordered
    return float(kept.mean())


def _bootstrap_interval(
    values: np.ndarray, statistic: str, seed_text: str
) -> tuple[float, float]:
    seed = int(hashlib.sha256(seed_text.encode("utf-8")).hexdigest()[:16], 16)
    rng = np.random.default_rng(seed)
    sampled = rng.choice(
        values, size=(BOOTSTRAP_REPLICATES, len(values)), replace=True
    )
    if statistic == "median":
        estimates = np.median(sampled, axis=1)
    elif statistic == "geometric_mean":
        estimates = np.exp(np.mean(np.log(sampled), axis=1))
    else:
        _fail(f"Unknown bootstrap statistic: {statistic}")
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def summarize_species(elements: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct every declared branch from the row-auditable element table."""

    eligible = elements["domain_gate_keep"].astype(bool)
    coverage = eligible & elements["coverage_ge_80pct"].astype(bool)
    complete = eligible & elements["Complete"].astype(str).str.lower().eq("yes")
    branches = {
        "all_ltr_gypsy_ge5_domain_elements": eligible,
        "coverage_ge_80pct": coverage,
        "coverage_ge_80pct_and_complete_yes": coverage & complete,
        PRIMARY_BRANCH: elements["primary_ltr_element"].astype(bool),
        "tesorter_complete_yes": complete,
    }

    rows: list[dict[str, object]] = []
    for branch, mask in branches.items():
        selected = elements.loc[mask]
        for species, group in selected.groupby("species", sort=True):
            values = pd.to_numeric(group[RATIO_COLUMN], errors="raise").to_numpy()
            values = values[np.isfinite(values) & (values > 0)]
            if not len(values):
                continue
            median_ci = _bootstrap_interval(
                values, "median", f"{species}:{branch}:median"
            )
            geometric_ci = _bootstrap_interval(
                values, "geometric_mean", f"{species}:{branch}:geometric"
            )
            mean_value = float(values.mean())
            median_value = float(np.median(values))
            geometric_mean = float(np.exp(np.mean(np.log(values))))
            winsorized = np.clip(values, *np.quantile(values, [0.05, 0.95]))
            loo_means = [
                float(np.delete(values, index).mean()) for index in range(len(values))
            ]
            loo_medians = [
                float(np.median(np.delete(values, index)))
                for index in range(len(values))
            ]
            rows.append(
                {
                    "species": species,
                    "te_sra_accession": group["te_sra_accession"].iloc[0],
                    "te_assembly_accession": group["te_assembly_accession"].iloc[0],
                    "analysis_branch": branch,
                    "n_elements": len(values),
                    "ratio_arithmetic_mean": mean_value,
                    "ratio_median": median_value,
                    "ratio_geometric_mean": geometric_mean,
                    "ratio_trimmed_mean_10pct": _trimmed_mean(values),
                    "ratio_winsorized_mean_5pct": float(winsorized.mean()),
                    "ratio_q05": float(np.quantile(values, 0.05)),
                    "ratio_q95": float(np.quantile(values, 0.95)),
                    "ratio_max": float(values.max()),
                    "fraction_ratio_gt_1": float(np.mean(values > 1)),
                    "median_bootstrap_ci_low": median_ci[0],
                    "median_bootstrap_ci_high": median_ci[1],
                    "geometric_mean_bootstrap_ci_low": geometric_ci[0],
                    "geometric_mean_bootstrap_ci_high": geometric_ci[1],
                    "max_abs_leave_one_out_mean_shift": float(
                        np.max(np.abs(np.asarray(loo_means) - mean_value))
                    ),
                    "max_abs_leave_one_out_median_shift": float(
                        np.max(np.abs(np.asarray(loo_medians) - median_value))
                    ),
                    "largest_element_fraction_of_ratio_sum": float(
                        values.max() / values.sum()
                    ),
                    "median_terminal_positive_coverage_fraction": float(
                        group["terminal_positive_coverage_fraction"].median()
                    ),
                    "median_internal_positive_coverage_fraction": float(
                        group["internal_positive_coverage_fraction"].median()
                    ),
                }
            )
    result = pd.DataFrame(rows).sort_values(["analysis_branch", "species"])
    return result.reset_index(drop=True)


def summarize_resources(
    resources: pd.DataFrame, elements: pd.DataFrame
) -> pd.DataFrame:
    """Rebuild resource counts from the corrected row-level selection decisions."""

    source = resources.copy()
    historical_columns = {
        "n_historical_pipe_split_5_6_domain_elements",
        "n_historical_source_corrupt_elements",
    }
    missing = sorted(historical_columns.difference(source.columns))
    if missing:
        _fail(f"LTR resource table lacks historical audit columns: {missing}")

    result = source.drop(
        columns=RESOURCE_REBUILT_COLUMNS,
        errors="ignore",
    ).copy()
    if result["species"].duplicated().any() or len(result) != 34:
        _fail("LTR resource table must contain 34 unique species")
    selected_counts = (
        elements.loc[elements["domain_gate_keep"]].groupby("species").size()
    )
    coverage_counts = (
        elements.loc[
            elements["domain_gate_keep"] & elements["coverage_ge_80pct"]
        ]
        .groupby("species")
        .size()
    )
    iqr_counts = elements.loc[elements["iqr_filter_keep"]].groupby("species").size()
    primary_counts = (
        elements.loc[elements["primary_ltr_element"]].groupby("species").size()
    )
    result["n_selected_ltr_gypsy_ge5_domain_elements"] = (
        result["species"].map(selected_counts).fillna(0).astype(int)
    )
    result["n_usable_elements"] = result[
        "n_selected_ltr_gypsy_ge5_domain_elements"
    ]
    result["n_selected_source_corrupt_exclusions"] = 0
    result["n_coverage_ge_80pct_elements"] = (
        result["species"].map(coverage_counts).fillna(0).astype(int)
    )
    result["n_iqr_retained_elements"] = (
        result["species"].map(iqr_counts).fillna(0).astype(int)
    )
    result["n_primary_iqr_filtered_elements"] = (
        result["species"].map(primary_counts).fillna(0).astype(int)
    )
    result["resource_status"] = np.select(
        [
            ~result["tabout_exists"].astype(bool),
            result["n_selected_ltr_gypsy_ge5_domain_elements"].eq(0),
        ],
        ["missing_tabout", "no_ltr_gypsy_ge5_domain_elements"],
        default="usable",
    )
    if result["n_historical_pipe_split_5_6_domain_elements"].sum() != 1_088:
        _fail("Historical LTR resource counts do not sum to 1,088")
    if result["n_historical_source_corrupt_elements"].sum() != 2:
        _fail("Historical LTR corrupt-element counts do not sum to two")
    if result["n_selected_ltr_gypsy_ge5_domain_elements"].sum() != 408:
        _fail("Corrected LTR resource selections do not sum to 408")
    if result["n_usable_elements"].sum() != 408:
        _fail("Corrected LTR usable counts do not sum to 408")
    if result["n_selected_source_corrupt_exclusions"].sum() != 0:
        _fail("No corrected-domain eligible element may be source-corrupt")
    if result["n_coverage_ge_80pct_elements"].sum() != 407:
        _fail("LTR resource coverage-qualified counts do not sum to 407")
    if result["n_iqr_retained_elements"].sum() != 381:
        _fail("LTR resource IQR-retained counts do not sum to 381")
    if result["n_primary_iqr_filtered_elements"].sum() != 380:
        _fail("LTR resource primary counts do not sum to 380")
    return result


def correct_and_validate_exclusions(excluded: pd.DataFrame) -> pd.DataFrame:
    """Retain the two corrupt historical rows while marking them domain-ineligible."""

    expected = {
        ("aureatus", "JAUEJG010675316.1_De_4598_10648"),
        ("santeetlah", "JASANM010280042.1_De_179_6669"),
    }
    result = excluded.drop(columns=DOMAIN_DERIVED_COLUMNS, errors="ignore").copy()
    observed = set(zip(result["species"], result["element_id"], strict=True))
    if observed != expected:
        _fail(f"Approved source-corrupt LTR exclusions changed: {sorted(observed)}")
    if not result["exclusion_reason"].eq(
        "truncated_depth_file_partial_final_line"
    ).all():
        _fail("LTR exclusions no longer have the approved truncation signature")
    if "legacy_pipe_split_domain_count" not in excluded:
        _fail("Historical LTR exclusions lack the preserved legacy domain count")
    legacy_count = excluded["legacy_pipe_split_domain_count"].astype(int)
    actual_count = legacy_count - 1
    if not actual_count.eq(4).all():
        _fail("Historical source-corrupt rows must each contain four TEsorter domains")
    result["legacy_pipe_split_domain_count"] = legacy_count
    result["domain_count"] = actual_count
    result["domain_gate_keep"] = False
    result["domain_gate_rule_id"] = DOMAIN_GATE_RULE_ID
    result["domain_gate_reason"] = "fewer_than_5_tesorter_domain_annotations"
    return result


def project_ltr_availability(
    availability: pd.DataFrame, resources: pd.DataFrame
) -> pd.DataFrame:
    """Project corrected LTR resource counts into the shared availability table."""

    result = availability.copy()
    if result["species"].duplicated().any():
        _fail("Analysis availability must contain one row per species")
    resource_by_species = resources.set_index("species")
    missing_species = sorted(set(resource_by_species.index) - set(result["species"]))
    if missing_species:
        _fail(f"LTR resources are absent from analysis availability: {missing_species}")

    mapped = result["species"].isin(resource_by_species.index)
    species = result.loc[mapped, "species"]
    projections = {
        "step04_ltr_resource_status": "resource_status",
        "step04_ltr_tabout_elements": "n_tabout_elements",
        "step04_ltr_selected_elements": (
            "n_selected_ltr_gypsy_ge5_domain_elements"
        ),
        "step04_ltr_usable_elements": "n_usable_elements",
        "step04_ltr_source_corrupt_exclusions": (
            "n_selected_source_corrupt_exclusions"
        ),
    }
    for target, source in projections.items():
        result.loc[mapped, target] = species.map(resource_by_species[source]).to_numpy()
    result.loc[mapped, "step04_ltr_deletion_footprint"] = species.map(
        resource_by_species["n_selected_ltr_gypsy_ge5_domain_elements"].gt(0)
    ).to_numpy()
    return result


def _assert_frames_equal(
    expected: pd.DataFrame, observed: pd.DataFrame, label: str
) -> None:
    try:
        pd.testing.assert_frame_equal(
            expected.reset_index(drop=True),
            observed.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
    except AssertionError as error:
        raise ValueError(f"Canonical {label} does not reproduce:\n{error}") from error


def build_release() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    elements = add_domain_gate(_read_csv(ELEMENT_PATH))
    elements = add_iqr_filter(elements)
    species = summarize_species(elements)
    resources = summarize_resources(_read_csv(RESOURCE_PATH), elements)
    exclusions = correct_and_validate_exclusions(_read_csv(EXCLUSION_PATH))
    availability = project_ltr_availability(
        _read_csv(AVAILABILITY_PATH), resources
    )
    return elements, species, resources, exclusions, availability


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Promote the reconstructed LTR tables and availability projection.",
    )
    args = parser.parse_args()

    elements, species, resources, exclusions, availability = build_release()
    if args.write:
        elements.to_csv(ELEMENT_PATH, index=False)
        species.to_csv(SPECIES_PATH, index=False)
        resources.to_csv(RESOURCE_PATH, index=False)
        exclusions.to_csv(EXCLUSION_PATH, index=False)
        availability.to_csv(AVAILABILITY_PATH, index=False)
        mode = "wrote"
    else:
        _assert_frames_equal(
            elements, _read_csv(ELEMENT_PATH), "LTR element table"
        )
        _assert_frames_equal(
            species, _read_csv(SPECIES_PATH), "LTR species robustness table"
        )
        _assert_frames_equal(
            resources, _read_csv(RESOURCE_PATH), "LTR resource table"
        )
        _assert_frames_equal(
            exclusions, _read_csv(EXCLUSION_PATH), "LTR historical exclusion table"
        )
        _assert_frames_equal(
            availability,
            _read_csv(AVAILABILITY_PATH),
            "analysis availability LTR projection",
        )
        mode = "verified"

    print(
        "LTR release PASS: "
        f"{mode}; 1,086 historical usable -> 408 LTR/Gypsy >=5-domain eligible; "
        "381 IQR-retained; 407 coverage-qualified -> 380 primary elements "
        "across 30 species"
    )


if __name__ == "__main__":
    main()

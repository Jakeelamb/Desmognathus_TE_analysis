#!/usr/bin/env python3
"""
Build a literature-backed substitution-rate sensitivity layer for paired-LTR ages.

This script does not alter the canonical sequence-divergence outputs. Instead, it
reads the validated `ltr_age_pairwise_divergence.csv` table and materializes a
separate calibration layer that ties every absolute-age estimate to an explicit
primary-literature nuclear substitution rate.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT, paths  # noqa: E402


PAIRWISE_PATH = paths.results.data / "ltr_age" / "ltr_age_pairwise_divergence.csv"
OUTPUT_DIR = paths.results.data / "ltr_age"
NOTE_PATH = PROJECT_ROOT / "LTR_SUBSTITUTION_RATE_CALIBRATION.md"


RATE_CANDIDATES: List[Dict[str, object]] = [
    {
        "rate_label": "salamander_rag1_floor",
        "rate_value": 1.00e-9,
        "analysis_tier": "recommended",
        "position": "lower",
        "source_short": "Herrick & Sclavi 2014",
        "source_url": "https://doi.org/10.1016/j.crpv.2014.06.002",
        "scope": "salamander nuclear synonymous rate floor",
        "justification": (
            "Plethodontid salamander dS/Mya rates bottom out at 0.001, which maps to "
            "1.00e-9 substitutions/site/year and provides the best lineage-matched lower anchor."
        ),
    },
    {
        "rate_label": "frog_cmyc_low",
        "rate_value": 0.924e-9,
        "analysis_tier": "supporting",
        "position": "low_support",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous c-myc lower bound",
        "justification": (
            "Direct frog nuclear synonymous-site estimate from c-myc. Very close to the "
            "salamander floor, so useful as supporting corroboration rather than a separate default."
        ),
    },
    {
        "rate_label": "frog_slug",
        "rate_value": 1.03e-9,
        "analysis_tier": "supporting",
        "position": "low_support",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous slug estimate",
        "justification": (
            "Independent frog nuclear synonymous estimate from slug in Xenopus. "
            "Supports the same approximately 1e-9 lower-amphibian range."
        ),
    },
    {
        "rate_label": "frog_cmyc_midpoint",
        "rate_value": (0.924e-9 + 1.53e-9) / 2,
        "analysis_tier": "recommended",
        "position": "central",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous c-myc midpoint",
        "justification": (
            "Arithmetic midpoint of Crawford's c-myc calibration window. This is the "
            "recommended central sensitivity value because it stays within the primary "
            "frog nuclear range while avoiding the faster tyrosinase outliers."
        ),
    },
    {
        "rate_label": "frog_cmyc_high",
        "rate_value": 1.53e-9,
        "analysis_tier": "recommended",
        "position": "upper",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous c-myc upper bound",
        "justification": (
            "Upper edge of the primary frog c-myc calibration window. Used as the "
            "recommended fast-end bound for the main sensitivity set."
        ),
    },
    {
        "rate_label": "frog_tyrosinase_low",
        "rate_value": 1.69e-9,
        "analysis_tier": "exploratory",
        "position": "fast_exploratory",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous tyrosinase lower estimate",
        "justification": (
            "Faster frog nuclear estimate from tyrosinase. Retained only as exploratory "
            "because Crawford notes the ranid calibration can bias fast-end estimates upward."
        ),
    },
    {
        "rate_label": "frog_tyrosinase_high",
        "rate_value": 3.35e-9,
        "analysis_tier": "excluded",
        "position": "excluded_fast_outlier",
        "source_short": "Crawford 2003",
        "source_url": "https://doi.org/10.1007/s00239-003-2513-7",
        "scope": "frog nuclear synonymous tyrosinase upper estimate",
        "justification": (
            "Not recommended for salamander LTR dating. Crawford explicitly warns that the "
            "higher ranid estimates may be inflated by underestimated divergence ages."
        ),
    },
]


def load_pairwise() -> pd.DataFrame:
    if not PAIRWISE_PATH.exists():
        raise FileNotFoundError(f"Missing paired-LTR divergence table: {PAIRWISE_PATH}")
    df = pd.read_csv(PAIRWISE_PATH)
    required = {"species", "k2p_distance", "recommended_high_confidence", "estimation_status"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Pairwise divergence table is missing required columns: {missing}")
    df = df[df["estimation_status"].eq("estimated")].copy()
    if df.empty:
        raise ValueError("No estimated paired-LTR rows found for calibration.")
    return df


def build_candidates() -> pd.DataFrame:
    df = pd.DataFrame(RATE_CANDIDATES)
    df["rate_value"] = df["rate_value"].astype(float)
    df["included_in_pairwise_output"] = df["analysis_tier"].ne("excluded")
    df["included_in_primary_summary"] = df["analysis_tier"].eq("recommended")
    tier_order = {
        "recommended": 0,
        "supporting": 1,
        "exploratory": 2,
        "excluded": 3,
    }
    df["tier_order"] = df["analysis_tier"].map(tier_order).fillna(99)
    df = df.sort_values(["tier_order", "rate_value"]).drop(columns=["tier_order"]).reset_index(drop=True)
    return df


def calibrate_pairwise(pairwise: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for candidate in candidates.itertuples(index=False):
        if not candidate.included_in_pairwise_output:
            continue
        calibrated = pairwise.copy()
        calibrated["rate_label"] = candidate.rate_label
        calibrated["rate_value"] = candidate.rate_value
        calibrated["analysis_tier"] = candidate.analysis_tier
        calibrated["age_years_calibrated"] = calibrated["k2p_distance"] / (2 * candidate.rate_value)
        calibrated["age_mya_calibrated"] = calibrated["age_years_calibrated"] / 1_000_000
        rows.append(
            calibrated[
                [
                    "species",
                    "k2p_distance",
                    "p_distance",
                    "recommended_high_confidence",
                    "rate_label",
                    "rate_value",
                    "analysis_tier",
                    "age_years_calibrated",
                    "age_mya_calibrated",
                ]
            ]
        )
    return pd.concat(rows, ignore_index=True)


def build_summary(calibrated: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for candidate in candidates.itertuples(index=False):
        subset = calibrated[calibrated["rate_label"].eq(candidate.rate_label)].copy()
        if subset.empty:
            continue
        high_conf = subset[subset["recommended_high_confidence"].fillna(False)].copy()
        species_medians = subset.groupby("species", dropna=False)["age_mya_calibrated"].median()
        rows.append(
            {
                "rate_label": candidate.rate_label,
                "rate_value": candidate.rate_value,
                "analysis_tier": candidate.analysis_tier,
                "position": candidate.position,
                "source_short": candidate.source_short,
                "source_url": candidate.source_url,
                "n_pairs": int(len(subset)),
                "n_species": int(subset["species"].nunique()),
                "median_age_mya": float(subset["age_mya_calibrated"].median()),
                "median_age_high_confidence_mya": float(high_conf["age_mya_calibrated"].median())
                if not high_conf.empty
                else float("nan"),
                "p10_age_mya": float(subset["age_mya_calibrated"].quantile(0.10)),
                "p90_age_mya": float(subset["age_mya_calibrated"].quantile(0.90)),
                "youngest_species_median_mya": float(species_medians.min()),
                "oldest_species_median_mya": float(species_medians.max()),
            }
        )
    return pd.DataFrame(rows).sort_values("rate_value").reset_index(drop=True)


def build_species_summary(calibrated: pd.DataFrame) -> pd.DataFrame:
    species = (
        calibrated.groupby(["rate_label", "rate_value", "analysis_tier", "species"], dropna=False)
        .agg(
            n_pairs_estimated=("species", "size"),
            n_high_confidence_pairs=("recommended_high_confidence", "sum"),
            median_age_mya=("age_mya_calibrated", "median"),
            median_age_high_confidence_mya=(
                "age_mya_calibrated",
                lambda s: float("nan"),
            ),
        )
        .reset_index()
    )

    hc = calibrated[calibrated["recommended_high_confidence"].fillna(False)].copy()
    if not hc.empty:
        hc_species = (
            hc.groupby(["rate_label", "species"], dropna=False)["age_mya_calibrated"]
            .median()
            .rename("median_age_high_confidence_mya")
            .reset_index()
        )
        species = species.drop(columns=["median_age_high_confidence_mya"]).merge(
            hc_species,
            on=["rate_label", "species"],
            how="left",
        )

    return species.sort_values(["rate_value", "median_age_mya"], ascending=[True, False]).reset_index(
        drop=True
    )


def write_note(candidates: pd.DataFrame, summary: pd.DataFrame, species: pd.DataFrame) -> None:
    recommended = summary[summary["analysis_tier"].eq("recommended")].copy()
    excluded = candidates[candidates["analysis_tier"].eq("excluded")].copy()
    exploratory = candidates[candidates["analysis_tier"].eq("exploratory")].copy()
    central = recommended[recommended["position"].eq("central")].iloc[0]
    low = recommended[recommended["position"].eq("lower")].iloc[0]
    high = recommended[recommended["position"].eq("upper")].iloc[0]

    central_species = species[species["rate_label"].eq(central["rate_label"])].copy()
    youngest = central_species.nsmallest(3, "median_age_mya")[["species", "median_age_mya"]]
    oldest = central_species.nlargest(3, "median_age_mya")[["species", "median_age_mya"]]

    lines = [
        "# LTR Substitution-Rate Calibration",
        "",
        "This note records the primary-literature substitution-rate calibration used to convert",
        "paired-LTR K2P divergence into absolute age sensitivities.",
        "",
        "## Bottom Line",
        "",
        (
            f"- The recommended salamander-compatible calibration window is "
            f"`{low['rate_value']:.2e}` to `{high['rate_value']:.2e}` substitutions/site/year."
        ),
        (
            f"- The recommended central sensitivity value is `{central['rate_value']:.3e}` "
            f"substitutions/site/year (`{central['rate_label']}`)."
        ),
        (
            f"- Across all successful paired-LTR estimates, the median inferred age is "
            f"`{low['median_age_mya']:.3f}` Myr at the lower bound, "
            f"`{central['median_age_mya']:.3f}` Myr at the central value, and "
            f"`{high['median_age_mya']:.3f}` Myr at the upper bound."
        ),
        (
            f"- Across the recommended high-confidence subset, the corresponding medians are "
            f"`{low['median_age_high_confidence_mya']:.3f}`, "
            f"`{central['median_age_high_confidence_mya']:.3f}`, and "
            f"`{high['median_age_high_confidence_mya']:.3f}` Myr."
        ),
        (
            f"- The faster frog tyrosinase estimate at `{exploratory['rate_value'].min():.2e}` "
            "is retained only as exploratory sensitivity, while the very fast "
            f"`{excluded['rate_value'].min():.2e}` tyrosinase estimate is excluded from the "
            "main Desmognathus sensitivity set."
        ),
        "",
        "## Primary Sources",
        "",
        (
            "- Herrick and Sclavi (2014) used salamander `rag1` synonymous substitution rates "
            "and reported that Plethodontidae rates bottom out at `0.001 dS/Mya`, which maps "
            "to `1.00e-9` substitutions/site/year. This is the best lineage-matched lower anchor."
        ),
        (
            "- Crawford (2003) reported frog nuclear synonymous-site rates of `0.924e-9` to "
            "`1.53e-9` for `c-myc`, `1.03e-9` for `slug`, and `1.69e-9` to `3.35e-9` for "
            "`tyrosinase`."
        ),
        (
            "- Crawford also explicitly noted that the higher ranid tyrosinase estimates could "
            "be biased upward if the calibration underestimated divergence ages, so they are "
            "not used as the main Desmognathus defaults."
        ),
        "",
        "## Recommended Use",
        "",
        (
            f"- Main analysis text: report LTR ages as a sensitivity range using "
            f"`{low['rate_value']:.2e}` to `{high['rate_value']:.2e}` substitutions/site/year."
        ),
        (
            f"- Supplementary tables/figures: use `{central['rate_label']}` "
            f"(`{central['rate_value']:.3e}`) as the representative central estimate."
        ),
        "- Do not present a single absolute age as if the substitution rate were directly measured in Desmognathus.",
        "",
        "## Central-Rate Species Extremes",
        "",
    ]

    for row in oldest.itertuples(index=False):
        lines.append(f"- Oldest median species at the central rate: `{row.species}` = `{row.median_age_mya:.3f}` Myr")
    for row in youngest.itertuples(index=False):
        lines.append(f"- Youngest median species at the central rate: `{row.species}` = `{row.median_age_mya:.3f}` Myr")

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `results/data/ltr_age/ltr_substitution_rate_candidates.csv`",
            "- `results/data/ltr_age/ltr_age_calibration_summary.csv`",
            "- `results/data/ltr_age/ltr_age_species_summary_calibrated.csv`",
            "- `results/data/ltr_age/ltr_age_pairwise_divergence_calibrated.csv`",
            "",
            "## Source Links",
            "",
            "- Herrick and Sclavi 2014: https://doi.org/10.1016/j.crpv.2014.06.002",
            "- Crawford 2003: https://doi.org/10.1007/s00239-003-2513-7",
        ]
    )

    NOTE_PATH.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pairwise = load_pairwise()
    candidates = build_candidates()
    calibrated = calibrate_pairwise(pairwise, candidates)
    summary = build_summary(calibrated, candidates)
    species = build_species_summary(calibrated)

    candidates.to_csv(OUTPUT_DIR / "ltr_substitution_rate_candidates.csv", index=False)
    summary.to_csv(OUTPUT_DIR / "ltr_age_calibration_summary.csv", index=False)
    species.to_csv(OUTPUT_DIR / "ltr_age_species_summary_calibrated.csv", index=False)
    calibrated.to_csv(OUTPUT_DIR / "ltr_age_pairwise_divergence_calibrated.csv", index=False)
    write_note(candidates, summary, species)

    print(f"Wrote rate candidates to {OUTPUT_DIR / 'ltr_substitution_rate_candidates.csv'}")
    print(f"Wrote calibration summary to {OUTPUT_DIR / 'ltr_age_calibration_summary.csv'}")
    print(f"Wrote species summary to {OUTPUT_DIR / 'ltr_age_species_summary_calibrated.csv'}")
    print(f"Wrote pairwise calibrated ages to {OUTPUT_DIR / 'ltr_age_pairwise_divergence_calibrated.csv'}")
    print(f"Wrote note to {NOTE_PATH}")


if __name__ == "__main__":
    main()

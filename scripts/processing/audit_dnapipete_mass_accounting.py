#!/usr/bin/env python3
"""Preserve dnaPipeTE unresolved mass for the declared 18-species panel."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, Set

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MERGED_INPUT = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_merged_classifications.csv"
CLASS_BREAKDOWN = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_class_breakdown.csv"
ORDER_BREAKDOWN = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_order_breakdown.csv"
SUPERFAMILY_BREAKDOWN = (
    PROJECT_ROOT / "results" / "data" / "dnaPipeTE_superfamily_breakdown.csv"
)
PANEL = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "derived"
    / "panels"
    / "te_genome_primary_mediumplus.csv"
)
LOOKUP_TABLE = PROJECT_ROOT / "input_data" / "lookup_table.txt"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "corrected" / "dnapipete"
AUDIT_DIR = PROJECT_ROOT / "plans" / "publication-readiness-deep-audit"
MASS_OUTPUT = OUTPUT_DIR / "dnapipete_mass_accounting_analysis18_v2.csv"
ORDER_OUTPUT = OUTPUT_DIR / "dnaPipeTE_order_breakdown_mass_accounted_analysis18_v2.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "dnapipete_mass_accounting_analysis18_v2.manifest.json"
REPORT_OUTPUT = AUDIT_DIR / "dnapipete_mass_accounting_analysis18_v2.md"


def _canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def summarize_mass(
    merged: pd.DataFrame,
    analysis_species: Set[str],
    retained_categories: Dict[str, Set[str]],
) -> pd.DataFrame:
    """Summarize retained and unresolved aligned-base mass by species."""

    required = {"Species", "Source", "aligned_bases", *retained_categories.keys()}
    missing = sorted(required.difference(merged.columns))
    if missing:
        raise ValueError(f"dnaPipeTE merged table is missing columns: {missing}")

    frame = merged.copy()
    frame["species"] = frame["Species"].map(_canonical_species)
    frame = frame.loc[frame["species"].isin(analysis_species)].copy()
    frame["aligned_bases"] = pd.to_numeric(frame["aligned_bases"], errors="coerce")
    if frame["aligned_bases"].isna().any() or frame["aligned_bases"].lt(0).any():
        raise ValueError("dnaPipeTE merged table has invalid aligned_bases")
    observed_species = set(frame["species"].unique())
    if observed_species != set(analysis_species):
        raise ValueError(
            "dnaPipeTE species differ from the declared analysis panel: "
            f"missing={sorted(set(analysis_species) - observed_species)}"
        )

    source_counts = frame.groupby("species")["Source"].nunique(dropna=True)
    if source_counts.ne(1).any():
        bad = source_counts.loc[source_counts.ne(1)].to_dict()
        raise ValueError(f"Final-panel species do not map to exactly one SRX: {bad}")

    base = frame.groupby("species", as_index=True).agg(
        sra_accession=("Source", "first"),
        n_component_rows=("aligned_bases", "size"),
        total_aligned_bases=("aligned_bases", "sum"),
    )
    for category, retained_labels in retained_categories.items():
        prefix = category.lower()
        retained = frame[category].astype("string").isin(retained_labels)
        retained_bp = (
            frame.loc[retained]
            .groupby("species")["aligned_bases"]
            .sum()
            .reindex(base.index, fill_value=0)
        )
        retained_rows = (
            frame.loc[retained].groupby("species").size().reindex(base.index, fill_value=0)
        )
        base[f"{prefix}_retained_component_rows"] = retained_rows.astype(int)
        base[f"{prefix}_retained_aligned_bases"] = retained_bp.astype(int)
        base[f"{prefix}_unresolved_aligned_bases"] = (
            base["total_aligned_bases"] - retained_bp
        ).astype(int)
        base[f"{prefix}_retained_fraction"] = retained_bp / base["total_aligned_bases"]
        base[f"{prefix}_unresolved_fraction"] = 1 - base[f"{prefix}_retained_fraction"]

    return base.reset_index().sort_values("species").reset_index(drop=True)


def validate_analysis_resources(
    mass: pd.DataFrame,
    lookup: pd.DataFrame,
    analysis_species: Set[str],
) -> pd.DataFrame:
    """Require each observed dnaPipeTE SRX to match the active TE lookup."""

    required_mass = {"species", "sra_accession"}
    missing_mass = sorted(required_mass.difference(mass.columns))
    if missing_mass:
        raise ValueError(f"Mass ledger is missing resource columns: {missing_mass}")
    required_lookup = {"Species", "SRA_Accension", "Genome_Accension"}
    missing_lookup = sorted(required_lookup.difference(lookup.columns))
    if missing_lookup:
        raise ValueError(f"Active TE lookup is missing columns: {missing_lookup}")

    resources = lookup.loc[:, sorted(required_lookup)].copy()
    resources["species"] = resources["Species"].map(_canonical_species)
    resources = resources.loc[resources["species"].isin(analysis_species)].copy()
    if resources["species"].duplicated().any():
        duplicates = sorted(resources.loc[resources["species"].duplicated(False), "species"])
        raise ValueError(f"Active TE lookup has duplicate analysis species: {duplicates}")
    observed_lookup_species = set(resources["species"])
    if observed_lookup_species != set(analysis_species):
        raise ValueError(
            "Active TE lookup differs from the declared analysis panel: "
            f"missing={sorted(set(analysis_species) - observed_lookup_species)}"
        )

    resources = resources.rename(
        columns={
            "SRA_Accension": "te_sra_accession",
            "Genome_Accension": "te_assembly_accession",
        }
    )[["species", "te_sra_accession", "te_assembly_accession"]]
    validated = mass.merge(resources, on="species", how="left", validate="one_to_one")
    mismatched = validated["sra_accession"].ne(validated["te_sra_accession"])
    if mismatched.any():
        detail = validated.loc[
            mismatched, ["species", "sra_accession", "te_sra_accession"]
        ].to_dict("records")
        raise ValueError(
            f"Observed dnaPipeTE SRX accessions do not match the active TE lookup: {detail}"
        )
    validated = validated.drop(columns="sra_accession")
    front = ["species", "te_sra_accession", "te_assembly_accession"]
    return validated[front + [column for column in validated if column not in front]]


def recompute_relative_breakdown(
    merged: pd.DataFrame,
    analysis_species: Set[str],
    category: str,
    retained_labels: Set[str],
) -> pd.DataFrame:
    frame = merged.copy()
    frame["species"] = frame["Species"].map(_canonical_species)
    frame = frame.loc[
        frame["species"].isin(analysis_species)
        & frame[category].astype("string").isin(retained_labels)
    ].copy()
    grouped = frame.groupby(["species", category])["aligned_bases"].sum().unstack(fill_value=0)
    grouped = grouped.reindex(columns=sorted(retained_labels), fill_value=0)
    return grouped.div(grouped.sum(axis=1), axis=0) * 100


def render_report(mass: pd.DataFrame, max_order_difference: float) -> str:
    fuscus = mass.loc[mass["species"].eq("fuscus")].iloc[0]
    return f"""# dnaPipeTE mass accounting: final 18-species panel

The historical breakdown tables are closed relative compositions. This audit preserves their aligned-base denominators and the mass omitted when order or superfamily labels are unresolved.

- Species included: {len(mass)} (exactly `te_genome_primary_mediumplus`).
- Species outside that panel, including *D. orestes*, are not processed into this corrected analysis branch.
- Median order-level unresolved fraction: {mass.order_unresolved_fraction.median():.3%}.
- Order-level unresolved range: {mass.order_unresolved_fraction.min():.3%}–{mass.order_unresolved_fraction.max():.3%}.
- Median superfamily-level unresolved fraction: {mass.superfamily_unresolved_fraction.median():.3%}.
- Superfamily-level unresolved range: {mass.superfamily_unresolved_fraction.min():.3%}–{mass.superfamily_unresolved_fraction.max():.3%}.
- Recomputed retained-order composition matches the historical relative table within {max_order_difference:.3e} percentage points.

For *D. fuscus* (`{fuscus.te_sra_accession}` / `{fuscus.te_assembly_accession}`), order-level retained mass is {fuscus.order_retained_fraction:.3%} and unresolved mass is {fuscus.order_unresolved_fraction:.3%} of dnaPipeTE aligned repeat bases.

These denominators still do **not** estimate the absolute fraction of the genome occupied by TEs. They quantify representation within the dnaPipeTE aligned-repeat output. Until an absolute repeat-load estimator is validated across species, the order table must be described as relative composition.
"""


def main() -> None:
    occupied = [
        path
        for path in (MASS_OUTPUT, ORDER_OUTPUT, MANIFEST_OUTPUT, REPORT_OUTPUT)
        if path.exists()
    ]
    if occupied:
        raise FileExistsError(
            "Non-destructive mass audit requires unused outputs: "
            + ", ".join(str(path) for path in occupied)
        )
    panel_species = pd.read_csv(PANEL)["species"].map(_canonical_species)
    if panel_species.duplicated().any() or len(panel_species) != 18:
        raise ValueError("Final TE/genome primary medium-plus panel must contain 18 unique species")
    analysis_species = set(panel_species)
    breakdown_paths = {
        "Class": CLASS_BREAKDOWN,
        "Order": ORDER_BREAKDOWN,
        "Superfamily": SUPERFAMILY_BREAKDOWN,
    }
    retained_categories = {
        category: set(pd.read_csv(path, nrows=0, index_col=0).columns)
        for category, path in breakdown_paths.items()
    }
    chunks = []
    for chunk in pd.read_csv(
        MERGED_INPUT,
        usecols=["Species", "Source", "aligned_bases", "Class", "Order", "Superfamily"],
        chunksize=500_000,
        low_memory=False,
    ):
        species = chunk["Species"].map(_canonical_species)
        selected = chunk.loc[species.isin(analysis_species)]
        if not selected.empty:
            chunks.append(selected)
    merged = pd.concat(chunks, ignore_index=True)
    mass = summarize_mass(merged, analysis_species, retained_categories)
    active_lookup = pd.read_csv(LOOKUP_TABLE, sep="\t")
    mass = validate_analysis_resources(mass, active_lookup, analysis_species)

    recomputed_order = recompute_relative_breakdown(
        merged, analysis_species, "Order", retained_categories["Order"]
    )
    historical_order = pd.read_csv(ORDER_BREAKDOWN, index_col=0)
    historical_order.index = historical_order.index.to_series().map(_canonical_species)
    historical_order = historical_order.reindex(
        index=sorted(analysis_species), columns=recomputed_order.columns
    )
    recomputed_order = recomputed_order.reindex(
        index=historical_order.index, columns=historical_order.columns
    )
    max_order_difference = float((recomputed_order - historical_order).abs().max().max())
    if max_order_difference > 1e-9:
        raise RuntimeError(
            "Mass-accounted order composition does not reproduce the historical percentages"
        )

    augmented_order = historical_order.reset_index().rename(columns={"index": "species"})
    augmented_order = augmented_order.merge(
        mass[
            [
                "species",
                "te_sra_accession",
                "te_assembly_accession",
                "n_component_rows",
                "total_aligned_bases",
                "order_retained_aligned_bases",
                "order_unresolved_aligned_bases",
                "order_retained_fraction",
                "order_unresolved_fraction",
            ]
        ],
        on="species",
        how="left",
        validate="one_to_one",
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    mass.to_csv(MASS_OUTPUT, index=False)
    augmented_order.to_csv(ORDER_OUTPUT, index=False)
    REPORT_OUTPUT.write_text(render_report(mass, max_order_difference), encoding="utf-8")
    manifest = {
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "analysis_species": sorted(analysis_species),
        "n_species": len(mass),
        "merged_input": _portable(MERGED_INPUT),
        "merged_input_sha256": _sha256(MERGED_INPUT),
        "panel": _portable(PANEL),
        "panel_sha256": _sha256(PANEL),
        "active_te_lookup": _portable(LOOKUP_TABLE),
        "active_te_lookup_sha256": _sha256(LOOKUP_TABLE),
        "analysis_resources": mass[
            ["species", "te_sra_accession", "te_assembly_accession"]
        ].to_dict("records"),
        "mass_output": _portable(MASS_OUTPUT),
        "mass_output_sha256": _sha256(MASS_OUTPUT),
        "mass_augmented_order_output": _portable(ORDER_OUTPUT),
        "mass_augmented_order_output_sha256": _sha256(ORDER_OUTPUT),
        "report": _portable(REPORT_OUTPUT),
        "report_sha256": _sha256(REPORT_OUTPUT),
        "max_order_reproduction_difference_percentage_points": max_order_difference,
        "composition_interpretation": "relative_dnapipete_aligned_repeat_composition",
        "eligible_as_absolute_genome_te_load": False,
        "out_of_panel_species_processed": False,
        "supersedes_for_analysis": "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v1.csv",
        "superseded_output_preserved": True,
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {MASS_OUTPUT}")
    print(f"Wrote {ORDER_OUTPUT}")
    print(f"Wrote {MANIFEST_OUTPUT}")
    print(f"Wrote {REPORT_OUTPUT}")


if __name__ == "__main__":
    main()

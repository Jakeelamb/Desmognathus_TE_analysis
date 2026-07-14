#!/usr/bin/env python3
"""Create a TE34 sensitivity branch that averages same-species dnaPipeTE runs.

The exact-SRX TE34 cache is preserved.  For a species with rerun labels, each
run first becomes a closed within-run composition; species-level composition,
mass fractions, diversity, and PCA are then the equal-weight mean across runs.
This prevents a retry from becoming an extra species or from receiving extra
weight merely because it has more aligned bases.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/processing"))
import build_te34_descriptive_bundle as base  # noqa: E402


BASE = ROOT / "results/data/corrected/te34"
OUT = ROOT / "results/data/corrected/te34_replicate_averaged"
FIG = ROOT / "results/figures/corrected/te34_replicate_averaged"
MERGED = ROOT / "results/data/dnaPipeTE_merged_classifications.csv"
ORDER = ROOT / "results/data/dnaPipeTE_order_breakdown.csv"
SUPERFAMILY = ROOT / "results/data/dnaPipeTE_superfamily_breakdown.csv"
LOOKUP = ROOT / "input_data/lookup_table.txt"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(value: object) -> str:
    return base.canonical(value)


def active_runs(species: str, accession: str) -> tuple[pd.DataFrame, list[str]]:
    """Read exact accession plus explicitly retry-labelled runs for one species."""

    frames = []
    pattern = re.compile(rf"^{re.escape(accession)}(?:R?\d+)?$")
    for chunk in pd.read_csv(
        MERGED,
        usecols=["Species", "Source", "aligned_bases", "Order", "Superfamily"],
        chunksize=400_000,
        low_memory=False,
    ):
        chunk["species"] = chunk["Species"].map(canonical)
        selected = chunk.loc[
            chunk["species"].eq(species)
            & chunk["Source"].astype(str).map(lambda item: bool(pattern.fullmatch(item)))
        ].copy()
        if not selected.empty:
            frames.append(selected)
    result = pd.concat(frames, ignore_index=True)
    labels = sorted(result["Source"].astype(str).unique())
    if accession not in labels or not labels:
        raise ValueError(f"Missing active dnaPipeTE source {accession} for {species}")
    return result, labels


def per_run_category(frame: pd.DataFrame, category: str, valid_labels: list[str], sources: list[str]) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    total = frame.groupby("Source")["aligned_bases"].sum().reindex(sources).astype(float)
    grouped = (
        frame.loc[frame[category].astype("string").isin(valid_labels)]
        .groupby(["Source", category])["aligned_bases"].sum()
        .unstack(fill_value=0)
        .reindex(index=sources, columns=valid_labels, fill_value=0)
        .astype(float)
    )
    retained = grouped.sum(axis=1) / total
    conditional = grouped.div(grouped.sum(axis=1), axis=0)
    return conditional, retained, total


def historic_composition(path: Path, species: list[str]) -> pd.DataFrame:
    table = pd.read_csv(path, index_col=0)
    table.index = table.index.map(canonical)
    table = table.reindex(species).apply(pd.to_numeric, errors="raise")
    return table.div(table.sum(axis=1), axis=0)


def write_figure(scores: pd.DataFrame) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    data = scores.loc[scores["te_level"].eq("order")]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(data["PC1"], data["PC2"], color="#2563eb")
    for row in data.itertuples():
        ax.annotate(row.species, (row.PC1, row.PC2), fontsize=7, xytext=(2, 2), textcoords="offset points")
    ax.set(title="TE34 order CLR PCA: retry runs averaged within species", xlabel="PC1", ylabel="PC2")
    path = FIG / "te34_replicate_averaged_order_clr_pca_v1.png"
    fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig)
    return path


def main() -> None:
    files = {
        "mass": OUT / "dnapipete_mass_accounting_te34_replicate_averaged_v1.csv",
        "diversity": OUT / "te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv",
        "scores": OUT / "te_pca_scores_te34_replicate_averaged_v1.csv",
        "inventory": OUT / "repeatmasker_hit_inventory_te34_replicate_averaged_v1.csv",
        "landscape": OUT / "repeatmasker_divergence_landscape_te34_replicate_averaged_v1.csv",
        "manifest": OUT / "te34_replicate_averaged_v1.manifest.json",
    }
    if any(path.exists() for path in files.values()):
        raise FileExistsError("Replicate-averaged TE34 outputs already exist; refusing overwrite")
    mass = pd.read_csv(BASE / "dnapipete_mass_accounting_te34_v1.csv")
    species = mass["species"].tolist()
    lookup = pd.read_csv(LOOKUP, sep="\t")
    lookup["species"] = lookup["Species"].map(canonical)
    accession = lookup.set_index("species").loc["orestes", "SRA_Accension"]
    raw, sources = active_runs("orestes", accession)
    if len(sources) < 2:
        raise ValueError(f"Unexpected orestes source labels: {sources}")

    order_labels = pd.read_csv(ORDER, nrows=0, index_col=0).columns.tolist()
    super_labels = pd.read_csv(SUPERFAMILY, nrows=0, index_col=0).columns.tolist()
    order_runs, order_retained, run_total = per_run_category(raw, "Order", order_labels, sources)
    super_runs, super_retained, _ = per_run_category(raw, "Superfamily", super_labels, sources)
    order = historic_composition(ORDER, species)
    superfamily = historic_composition(SUPERFAMILY, species)
    # The historical orestes composition row is not assumed to be a particular
    # retry. Record its discrepancy and replace it with the raw-run mean.
    historical_order_vs_active_max_abs = float(
        (order.loc["orestes"] - order_runs.loc[accession]).abs().max()
    )
    historical_superfamily_vs_active_max_abs = float(
        (superfamily.loc["orestes"] - super_runs.loc[accession]).abs().max()
    )
    order.loc["orestes"] = order_runs.mean(axis=0)
    superfamily.loc["orestes"] = super_runs.mean(axis=0)

    mass["n_dnapipete_runs"] = 1
    mass["dnapipete_source_labels"] = mass["te_sra_accession"].map(lambda item: json.dumps([item]))
    mass["run_aggregation"] = "single_active_run"
    row = mass["species"].eq("orestes")
    mass.loc[row, "n_dnapipete_runs"] = len(sources)
    mass.loc[row, "dnapipete_source_labels"] = json.dumps(sources)
    mass.loc[row, "run_aggregation"] = "equal_weight_mean_of_run_level_estimates"
    mass.loc[row, "total_aligned_bases"] = float(run_total.mean())
    mass.loc[row, "n_component_rows"] = float(raw.groupby("Source").size().mean())
    for level, retained in [("order", order_retained), ("superfamily", super_retained)]:
        mass.loc[row, f"{level}_retained_fraction"] = float(retained.mean())
        mass.loc[row, f"{level}_unresolved_fraction"] = float(1 - retained.mean())
        mass.loc[row, f"{level}_retained_aligned_bases"] = float((retained * run_total).mean())
        mass.loc[row, f"{level}_unresolved_aligned_bases"] = float(((1 - retained) * run_total).mean())

    diversity_rows = []
    score_rows = []
    for level, composition in [("order", order), ("superfamily", superfamily)]:
        mass_aware = composition.mul(mass.set_index("species")[f"{level}_retained_fraction"], axis=0)
        mass_aware["Unresolved"] = mass.set_index("species")[f"{level}_unresolved_fraction"]
        for mode, table in [("classified_conditional", composition), ("mass_aware_unresolved_bin", mass_aware)]:
            for item, values in table.iterrows():
                diversity_rows.append({"species": item, "te_level": level, "composition_mode": mode, **base.diversity(values.to_numpy())})
        scores, _ = base.run_pca(composition)
        scores = scores.reset_index().rename(columns={"index": "species"})
        scores.insert(1, "te_level", level)
        scores.insert(2, "composition_mode", "classified_conditional")
        score_rows.append(scores)
    diversity = pd.DataFrame(diversity_rows)
    scores = pd.concat(score_rows, ignore_index=True)
    inventory = pd.read_csv(BASE / "repeatmasker_hit_inventory_te34_v1.csv")
    landscape = pd.read_csv(BASE / "repeatmasker_divergence_landscape_te34_v1.csv")

    OUT.mkdir(parents=True, exist_ok=True)
    for key, frame in [("mass", mass), ("diversity", diversity), ("scores", scores), ("inventory", inventory), ("landscape", landscape)]:
        frame.to_csv(files[key], index=False)
    figure = write_figure(scores)
    manifest = {
        "analysis_id": "te34_replicate_averaged_v1",
        "analysis_scope": "active_vetted_te_resource_panel34_descriptive",
        "n_te_species": 34,
        "dnapipete_run_aggregation": "equal_weight_mean_of_run_level_estimates",
        "species_with_multiple_dnapipete_runs": ["orestes"],
        "orestes_source_labels": sources,
        "retry_runs_counted_as_additional_species": False,
        "exact_srx_te34_v1_preserved": True,
        "historical_orestes_order_vs_active_run_max_abs_proportion": historical_order_vs_active_max_abs,
        "historical_orestes_superfamily_vs_active_run_max_abs_proportion": historical_superfamily_vs_active_max_abs,
        "historical_orestes_row_replaced_with_raw_run_mean": True,
        "repeatmasker_unchanged_reason": "the corrected RepeatMasker all-resource hit table contains one active resource per species; no retry-labelled duplicate was present",
        "inputs": {str(path.relative_to(ROOT)): sha256(path) for path in [MERGED, ORDER, SUPERFAMILY, *[BASE / item for item in ["dnapipete_mass_accounting_te34_v1.csv", "repeatmasker_hit_inventory_te34_v1.csv", "repeatmasker_divergence_landscape_te34_v1.csv"]]]},
        "outputs": {str(path.relative_to(ROOT)): sha256(path) for key, path in files.items() if key != "manifest"},
        "figure": {"path": str(figure.relative_to(ROOT)), "sha256": sha256(figure)},
        "expensive_upstream_tools_executed": False,
        "eligible_for_integrated_phylogenetic_path_analysis": False,
    }
    files["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"species": "orestes", "run_labels": sources, "n_runs": len(sources)}, indent=2))


if __name__ == "__main__":
    main()

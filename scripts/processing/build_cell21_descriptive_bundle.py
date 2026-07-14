#!/usr/bin/env python3
"""Build a linked-cell21 descriptive table from frozen microscopy outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv"
MORPH = ROOT / "path_analysis/data/external/derived/cellprofiler_species_morphology_summary.csv"
IOD = ROOT / "path_analysis/data/external/derived/cellprofiler_final_species_results.csv"
OUT = ROOT / "results/data/corrected/cell21"
FIG = ROOT / "results/figures/corrected/cell21"


def canonical(value: object) -> str:
    return str(value).strip().replace("Desmognathus ", "").replace("D. ", "").replace("D.", "").lower()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    output = OUT / "cell_linked_traits_cell21_v1.csv"
    manifest_output = OUT / "cell21_descriptive_bundle_v1.manifest.json"
    if output.exists() or manifest_output.exists():
        raise FileExistsError("Cell21 descriptive outputs already exist; refusing overwrite")
    panel = pd.read_csv(PANEL)
    panel["species"] = panel["species"].map(canonical)
    species = panel["species"].sort_values().tolist()
    if len(species) != 21 or panel["species"].duplicated().any():
        raise ValueError("Linked-cell panel must contain 21 unique species")
    morph = pd.read_csv(MORPH)
    morph["species"] = morph["species"].map(canonical)
    iod = pd.read_csv(IOD)
    iod["species"] = iod["species"].map(canonical)
    if morph["species"].duplicated().any() or iod["species"].duplicated().any():
        raise ValueError("Frozen microscopy summaries must be one row per species")
    table = panel[["species", "decision_basis"]].merge(
        morph[["species", "n_specimens_strict", "n_images_strict", "n_pairs_strict", "species_median_cell_area_um2", "species_median_nuc_area_um2", "linked_effective_n", "linked_support_label", "linked_support_warnings"]],
        on="species", how="left", validate="one_to_one",
    ).merge(
        iod[["species", "primary_n_images", "primary_n_specimens", "primary_n_nuclei", "primary_measurement_kind", "primary_state", "primary_support_tier", "primary_support_warnings", "primary_genome_pg", "result_status", "flag_summary"]],
        on="species", how="left", validate="one_to_one",
    )
    table = table.rename(columns={
        "species_median_cell_area_um2": "cell_area_um2",
        "species_median_nuc_area_um2": "nucleus_area_um2",
        "primary_genome_pg": "relative_nuclear_iod_raw_value",
        "primary_n_images": "iod_n_images",
        "primary_n_specimens": "iod_n_specimens",
        "primary_n_nuclei": "iod_n_nuclei",
    }).sort_values("species")
    required = ["cell_area_um2", "nucleus_area_um2", "relative_nuclear_iod_raw_value"]
    if table[required].isna().any().any():
        raise ValueError("A linked-cell species lacks a required frozen microscopy trait")
    OUT.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    scatter = ax.scatter(table["cell_area_um2"], table["nucleus_area_um2"], c=table["relative_nuclear_iod_raw_value"], cmap="viridis", s=52)
    for row in table.itertuples():
        ax.annotate(row.species, (row.cell_area_um2, row.nucleus_area_um2), fontsize=7, xytext=(2, 2), textcoords="offset points")
    ax.set(title="Linked-cell21 morphology and relative nuclear IOD", xlabel="upper-tail cell area (um²)", ylabel="corresponding nucleus area (um²)")
    fig.colorbar(scatter, ax=ax, label="relative nuclear IOD (raw image value)")
    figure = FIG / "cell21_cell_nucleus_relative_iod_v1.png"
    fig.tight_layout(); fig.savefig(figure, dpi=300); plt.close(fig)
    manifest = {
        "analysis_id": "cell21_descriptive_bundle_v1", "analysis_scope": "linked_cell_panel21_descriptive",
        "n_linked_cell_species": 21, "species": table["species"].tolist(),
        "relative_nuclear_iod_is_absolute_genome_size": False,
        "relative_nuclear_iod_interpretation": "area times mean image optical density; retained as a relative image phenotype only",
        "analysis18_outputs_preserved": True, "eligible_for_integrated_phylogenetic_path_analysis": False,
        "path_analysis_scope": "study_integrated_path_panel18_v1 only; untrusted genomic resources are not substituted for folkertsi or ochrophaeus",
        "inputs": {str(path.relative_to(ROOT)): digest(path) for path in [PANEL, MORPH, IOD]},
        "output": {"path": str(output.relative_to(ROOT)), "sha256": digest(output)},
        "figure": {"path": str(figure.relative_to(ROOT)), "sha256": digest(figure)},
    }
    manifest_output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_linked_cell_species": 21, "output": str(output.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

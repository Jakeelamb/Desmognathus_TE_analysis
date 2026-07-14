#!/usr/bin/env python3
"""Build audited inputs for the cell–nucleus–genome path-design notebook."""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Collection

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = (
    PROJECT_ROOT / "results" / "data" / "research_review" / "cell_nucleus_genome_path"
)
NOTEBOOK_PATH = (
    PROJECT_ROOT
    / "notebooks"
    / "research_review"
    / "07_cell_nucleus_genome_path_analysis.ipynb"
)
TRAIT_SUMMARY_PATH = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "research_review"
    / "genome_size_estimation"
    / "phylogeny_genome_nucleus_cell_summary.csv"
)
PRIMARY_TREE_PATH = PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"

GENOME = "genome_size"
NUCLEUS = "nucleus_size"
CELL = "cell_size"
NODES = (GENOME, NUCLEUS, CELL)
NODE_ABBREVIATIONS = {GENOME: "GS", NUCLEUS: "NS", CELL: "CS"}
DIRECTED_EDGES = tuple(
    (source, target) for source in NODES for target in NODES if source != target
)


def _dag_id(edges: Collection[tuple[str, str]]) -> str:
    edge_set = set(edges)
    if not edge_set:
        return "no_edges"
    return ";".join(
        f"{NODE_ABBREVIATIONS[source]}->{NODE_ABBREVIATIONS[target]}"
        for source, target in DIRECTED_EDGES
        if (source, target) in edge_set
    )


USER_MECHANISMS = {
    "genome_to_nucleus_to_cell": _dag_id(
        {(GENOME, NUCLEUS), (NUCLEUS, CELL)}
    ),
    "cell_to_nucleus_to_genome": _dag_id(
        {(CELL, NUCLEUS), (NUCLEUS, GENOME)}
    ),
    "nucleus_to_genome_and_cell": _dag_id(
        {(NUCLEUS, GENOME), (NUCLEUS, CELL)}
    ),
}


def _is_acyclic(edges: Collection[tuple[str, str]]) -> bool:
    remaining = set(NODES)
    edge_set = set(edges)
    while remaining:
        roots = {
            node
            for node in remaining
            if not any(target == node and source in remaining for source, target in edge_set)
        }
        if not roots:
            return False
        remaining -= roots
    return True


def _equivalence_class(
    edges: Collection[tuple[str, str]],
) -> tuple[str, str, bool]:
    edge_set = set(edges)
    skeleton = {frozenset(edge) for edge in edge_set}
    n_edges = len(skeleton)
    if n_edges == 0:
        return (
            "independent",
            "genome_size _||_ nucleus_size; genome_size _||_ cell_size; "
            "nucleus_size _||_ cell_size",
            True,
        )
    if n_edges == 1:
        pair = next(iter(skeleton))
        missing = [node for node in NODES if node not in pair][0]
        present = sorted(pair, key=NODES.index)
        class_name = {
            frozenset((GENOME, NUCLEUS)): "genome_nucleus_only",
            frozenset((NUCLEUS, CELL)): "nucleus_cell_only",
            frozenset((GENOME, CELL)): "genome_cell_only",
        }[pair]
        basis = "; ".join(f"{missing} _||_ {node}" for node in present)
        return class_name, basis, True
    if n_edges == 2:
        degrees = {
            node: sum(node in pair for pair in skeleton)
            for node in NODES
        }
        middle = next(node for node, degree in degrees.items() if degree == 2)
        endpoints = [node for node in NODES if node != middle]
        is_collider = all((endpoint, middle) in edge_set for endpoint in endpoints)
        suffix = "collider" if is_collider else "bridge"
        if is_collider:
            basis = f"{endpoints[0]} _||_ {endpoints[1]}"
        else:
            basis = f"{endpoints[0]} _||_ {endpoints[1]} | {middle}"
        middle_label = middle[:-5] if middle.endswith("_size") else middle
        return f"{middle_label}_{suffix}", basis, True
    if n_edges == 3:
        return "saturated", "none (no d-separation claim)", False
    raise AssertionError(f"Unexpected three-node skeleton with {n_edges} edges")


def enumerate_three_node_dags() -> pd.DataFrame:
    """Return all 25 labeled DAGs and their 11 Markov-equivalence classes."""
    records: list[dict[str, object]] = []
    for edge_count in range(len(DIRECTED_EDGES) + 1):
        for edges in combinations(DIRECTED_EDGES, edge_count):
            if not _is_acyclic(edges):
                continue
            class_name, basis_claim, testable = _equivalence_class(edges)
            records.append(
                {
                    "dag_id": _dag_id(edges),
                    "edges": tuple(edges),
                    "n_edges": edge_count,
                    "equivalence_class": class_name,
                    "basis_claim": basis_claim,
                    "testable_by_dsep": testable,
                }
            )
    frame = pd.DataFrame(records)
    if len(frame) != 25:
        raise AssertionError(f"Expected 25 labeled three-node DAGs, found {len(frame)}")
    frame["class_member_count"] = frame.groupby("equivalence_class")[
        "dag_id"
    ].transform("size")
    return frame.sort_values(
        ["n_edges", "equivalence_class", "dag_id"], kind="stable"
    ).reset_index(drop=True)


def load_primary_traits(path: Path = TRAIT_SUMMARY_PATH) -> pd.DataFrame:
    """Load every finalized species with all three measured traits."""
    source = pd.read_csv(path, low_memory=False)
    required = {
        "species",
        "genome_panel_status",
        "include_in_primary_genome_analysis",
        "genome_size_pg_fuscus_anchored",
        "nucleus_area_um2",
        "cell_area_um2",
        "n_genome_nuclei",
        "n_genome_images",
        "n_size_cells",
        "n_size_images",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError(f"Frozen trait summary is missing columns: {missing}")
    traits = source.loc[
        source["include_in_primary_genome_analysis"].astype(bool),
        [
            "species",
            "genome_size_pg_fuscus_anchored",
            "nucleus_area_um2",
            "cell_area_um2",
            "n_genome_nuclei",
            "n_genome_images",
            "n_size_cells",
            "n_size_images",
        ],
    ].copy()
    traits = traits.rename(
        columns={"genome_size_pg_fuscus_anchored": "genome_size_pg"}
    )
    traits.insert(1, "tree_tip", traits["species"].str.replace("D. ", "", regex=False))
    measurement_columns = ["genome_size_pg", "nucleus_area_um2", "cell_area_um2"]
    if traits.empty or traits["species"].nunique() != len(traits):
        raise ValueError("Finalized trait overlap is empty or has duplicate species.")
    if not traits[measurement_columns].notna().all().all():
        raise ValueError("Finalized overlap contains a missing trait measurement.")
    if not (traits[measurement_columns] > 0).all().all():
        raise ValueError("Finalized traits must be strictly positive before log transformation.")
    count_columns = [
        "n_genome_nuclei",
        "n_genome_images",
        "n_size_cells",
        "n_size_images",
    ]
    traits[count_columns] = traits[count_columns].astype(int)
    if (traits[count_columns] < 1).any().any():
        raise ValueError("Each finalized species must retain measured support for every trait.")
    return traits.reset_index(drop=True)


def build_measurement_bootstrap_traits(
    *, n_bootstrap: int = 250, seed: int = 20260710
) -> pd.DataFrame:
    """Build paired trait draws while preserving the shared fuscus calibration."""
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be positive")
    import build_frozen_genome_iod_phylogeny_figure as frozen_traits

    primary = load_primary_traits()
    _, long_draws, _ = frozen_traits.build_figure_data(
        n_bootstrap=n_bootstrap, seed=seed
    )
    primary_species = set(primary["species"])
    selected = long_draws.loc[long_draws["species"].isin(primary_species)].copy()
    wide = (
        selected.pivot(
            index=["bootstrap_replicate", "species"],
            columns="metric",
            values="value",
        )
        .reset_index()
        .rename_axis(columns=None)
        .rename(columns={"genome_size_pg_fuscus_anchored": "genome_size_pg"})
    )
    wide.insert(2, "tree_tip", wide["species"].str.replace("D. ", "", regex=False))
    columns = [
        "bootstrap_replicate",
        "species",
        "tree_tip",
        "genome_size_pg",
        "nucleus_area_um2",
        "cell_area_um2",
    ]
    wide = wide[columns].sort_values(
        ["bootstrap_replicate", "species"], kind="stable"
    ).reset_index(drop=True)
    counts = wide.groupby("bootstrap_replicate")["species"].nunique()
    expected_species = len(primary)
    if len(counts) != n_bootstrap or not counts.eq(expected_species).all():
        raise ValueError(
            "Every measurement-bootstrap replicate must contain the full finalized "
            f"primary panel ({expected_species} species)."
        )
    if wide.isna().any().any():
        raise ValueError("Measurement-bootstrap panel contains missing values.")
    return wide


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_outputs(*, n_bootstrap: int = 250, seed: int = 20260710) -> dict[str, object]:
    """Write the exact design tables and measurement-bootstrap input audit."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    dag_path = OUTPUT_DIR / "all_three_trait_dags.csv"
    trait_path = OUTPUT_DIR / "primary_traits_all_finalized.csv"
    bootstrap_path = OUTPUT_DIR / "measurement_bootstrap_traits_250.csv"
    manifest_path = OUTPUT_DIR / "cell_nucleus_genome_path_design_manifest.json"

    dags = enumerate_three_node_dags()
    traits = load_primary_traits()
    bootstrap = build_measurement_bootstrap_traits(
        n_bootstrap=n_bootstrap, seed=seed
    )
    dags.drop(columns="edges").to_csv(dag_path, index=False)
    traits.to_csv(trait_path, index=False)
    bootstrap.to_csv(bootstrap_path, index=False, float_format="%.8f")

    user_classes = (
        dags.set_index("dag_id")
        .loc[list(USER_MECHANISMS.values()), "equivalence_class"]
        .unique()
        .tolist()
    )
    manifest = {
        "analysis": "Cell-nucleus-genome phylogenetic path design and input audit",
        "status": "design_complete_model_fit_pending",
        "n_primary_species": int(len(traits)),
        "n_labeled_three_node_dags": int(len(dags)),
        "n_markov_equivalence_classes": int(dags["equivalence_class"].nunique()),
        "n_testable_dags": int(dags["testable_by_dsep"].sum()),
        "n_saturated_untestable_dags": int((~dags["testable_by_dsep"]).sum()),
        "user_mechanism_equivalence_classes": user_classes,
        "shared_user_mechanism_basis_claim": (
            "genome_size _||_ cell_size | nucleus_size"
        ),
        "measurement_bootstrap_replicates": int(n_bootstrap),
        "random_seed": int(seed),
        "actual_phylogenetic_path_models_fitted": False,
        "interpretation": (
            "The three proposed arrow directions are Markov-equivalent in "
            "cross-sectional three-trait data. This artifact audits candidate "
            "DAGs and inputs; it does not report a fitted causal winner."
        ),
        "outputs": {},
    }
    for label, path in (
        ("dag_table", dag_path),
        ("primary_trait_table", trait_path),
        ("measurement_bootstrap_table", bootstrap_path),
    ):
        manifest["outputs"][label] = str(path.relative_to(PROJECT_ROOT))
        manifest["outputs"][f"{label}_sha256"] = sha256_file(path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(build_outputs(), indent=2))

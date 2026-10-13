"""Build and validate the canonical genome-size-facing Path24 release tables.

The morphology release remains the authority for nuclear IOD measurements and
their conditional calibration. This builder reconstructs the Path24 traits,
measurement-bootstrap traits, candidate DAGs, and TE-overlap table directly
from canonical sources. The R runner remains the authority for fitted path
results; this module validates those frozen outputs against the canonical
genome-size schemas and identifiers.

The default mode is read-only and compares reconstructed products with the
checked-in release.  Pass ``--write`` to promote the reconstructed products.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Collection
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "analyses/03_phylogenetic_path/data"
TREE_PATH = ROOT / "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk"
IOD_SOURCE = ROOT / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
SIZE_SOURCE = ROOT / "analyses/02_morphology/data/cell_nucleus_species_estimates.csv"
IOD_OBJECT_SOURCE = ROOT / "analyses/02_morphology/data/nuclear_iod_objects.csv"
SIZE_OBJECT_SOURCE = ROOT / "analyses/02_morphology/data/cell_nucleus_objects.csv"
IOD_MANIFEST = ROOT / "analyses/02_morphology/provenance/iod_source_manifest.json"
TE_DIVERSITY_SOURCE = ROOT / "analyses/01_transposable_elements/data/te_diversity.csv"
TE_PANEL_SOURCE = ROOT / "data/identity/te34_panel.csv"
PATH_PANEL_SOURCE = ROOT / "data/identity/path24_panel.csv"

PATH_TRAIT_DEFINITION = (
    "conditional_fuscus_assembly_anchored_genome_size_pg_1c; "
    "iod_derived_not_measured_c_value"
)
GENOME_SIZE_DEFINITION = (
    "conditional D. fuscus assembly-anchored genome-size estimate (pg/1C) "
    "derived from nuclear IOD; anchor 16.1 Gbp / 0.978 Gbp per pg = "
    "16.462167689 pg/1C; not a measured species C-value; anchor uncertainty excluded"
)
STANDARDIZED_DEFINITION = (
    "z score of log10 conditional D. fuscus assembly-anchored genome size (pg/1C); "
    "IOD-derived, not a measured species C-value"
)
TE_DIVERSITY_DEFINITION = (
    "each species' ten most abundant positive classified TE superfamilies; "
    "selected species-specifically and reclosed to one; Shannon uses natural logs; "
    "Gini-Simpson is 1 - sum(p_i^2)"
)

PATH24_SCHEMA = [
    "species",
    "tree_tip",
    "genome_size_pg_fuscus_anchored",
    "genome_size_pg_ci_low",
    "genome_size_pg_ci_high",
    "nucleus_area_um2",
    "cell_area_um2",
    "n_genome_size_nuclei",
    "n_genome_size_images",
    "n_genome_size_specimens",
    "n_size_cells",
    "n_size_images",
    "genome_size_definition",
]
TE_OVERLAP_SCHEMA = [
    "species",
    "scientific_name",
    "manual_review_inclusion",
    "te_sra_accession",
    "te_assembly_accession",
    "n_superfamilies_included",
    "shannon_entropy",
    "gini_simpson",
    "genome_size_pg_fuscus_anchored",
    "genome_size_pg_ci_low",
    "genome_size_pg_ci_high",
    "n_nuclei",
    "n_images",
    "n_specimens",
    "estimate_support",
    "sample_size_support",
    "quality_match_status",
    "te_diversity_definition",
    "genome_size_definition",
]

FROZEN_RESULT_FILES = (
    "evolutionary_model_sensitivity.csv",
    "full_simulation_rankings.csv.gz",
    "pairwise_pgls.csv",
    "pairwise_phylogenetic_regressions.csv",
    "path_dsep_tests.csv",
    "path_edge_estimates.csv",
    "path_leave_one_species_out.csv",
    "path_measurement_sensitivity.csv",
    "path_model_comparison.csv",
    "path_model_stability.csv",
    "path_simulation_calibration.csv",
    "path_standardized_traits.csv",
    "path_tree_sensitivity.csv",
)

GENOME = "genome_size"
NUCLEUS = "nucleus_size"
CELL = "cell_size"
NODES = (GENOME, NUCLEUS, CELL)
NODE_ABBREVIATIONS = {GENOME: "GS", NUCLEUS: "NS", CELL: "CS"}
DIRECTED_EDGES = tuple(
    (source, target) for source in NODES for target in NODES if source != target
)
MODEL_IDS = {
    "independent",
    "genome_nucleus_only",
    "nucleus_cell_only",
    "genome_cell_only",
    "nucleus_bridge",
    "nucleus_collider",
    "cell_bridge",
    "cell_collider",
    "genome_bridge",
    "genome_collider",
}
RELATION_IDS = {
    "genome_to_nucleus",
    "nucleus_to_genome",
    "genome_to_cell",
    "cell_to_genome",
    "nucleus_to_cell",
    "cell_to_nucleus",
    "genome_and_nucleus_to_cell",
    "genome_and_cell_to_nucleus",
    "nucleus_and_cell_to_genome",
}
MECHANISM_IDS = {
    "genome_to_nucleus_to_cell",
    "cell_to_nucleus_to_genome",
    "nucleus_to_genome_and_cell",
}


def read_csv(path: Path) -> pd.DataFrame:
    """Read deposited decimals without the default parser's one-ULP drift."""
    return pd.read_csv(path, float_precision="round_trip")


def path_tip_order() -> list[str]:
    tree = Phylo.read(TREE_PATH, "newick")
    tips = [str(tip.name) for tip in tree.get_terminals()]
    if len(tips) != 24 or len(set(tips)) != 24:
        raise ValueError("Path24 tree must contain exactly 24 unique terminal labels")
    return tips


def build_path24_traits() -> pd.DataFrame:
    iod = read_csv(IOD_SOURCE)
    size = read_csv(SIZE_SOURCE)
    required_iod = {
        "species",
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
        "n_nuclei",
        "n_images",
        "n_specimens",
    }
    if missing := required_iod - set(iod.columns):
        raise ValueError(f"IOD/calibration source is missing columns: {sorted(missing)}")
    required_size = {
        "species",
        "nucleus_area_um2",
        "cell_area_um2",
        "n_size_cells",
        "n_size_images",
    }
    if missing := required_size - set(size.columns):
        raise ValueError(f"Morphology source is missing columns: {sorted(missing)}")

    traits = iod[
        [
            "species",
            "genome_size_pg_fuscus_anchored",
            "genome_size_pg_ci_low",
            "genome_size_pg_ci_high",
            "n_nuclei",
            "n_images",
            "n_specimens",
        ]
    ].merge(
        size[
            [
                "species",
                "nucleus_area_um2",
                "cell_area_um2",
                "n_size_cells",
                "n_size_images",
            ]
        ],
        on="species",
        how="inner",
        validate="one_to_one",
    )
    traits["tree_tip"] = traits["species"].str.removeprefix("D. ")
    order = path_tip_order()
    if set(traits["tree_tip"]) != set(order):
        raise ValueError("Morphology/calibration join does not equal the frozen Path24 tree panel")
    traits = traits.set_index("tree_tip").loc[order].reset_index()
    traits = traits.rename(
        columns={
            "n_nuclei": "n_genome_size_nuclei",
            "n_images": "n_genome_size_images",
            "n_specimens": "n_genome_size_specimens",
        }
    )
    traits["genome_size_definition"] = GENOME_SIZE_DEFINITION
    traits = traits[PATH24_SCHEMA]
    numeric = traits.select_dtypes(include="number")
    if numeric.isna().any().any() or (numeric <= 0).any().any():
        raise ValueError("Path24 modeled traits and support counts must be finite and positive")
    return traits


def build_te_genome_size_overlap() -> pd.DataFrame:
    te_panel = read_csv(TE_PANEL_SOURCE)
    path_panel = read_csv(PATH_PANEL_SOURCE)
    diversity = read_csv(TE_DIVERSITY_SOURCE)
    iod = read_csv(IOD_SOURCE).assign(
        species=lambda frame: frame["species"].str.removeprefix("D. ")
    )

    te = te_panel.loc[
        te_panel["te_resource_panel34_v1"].eq(True),
        [
            "species",
            "manual_review_inclusion",
            "te_sra_accession",
            "te_assembly_accession",
        ],
    ]
    path_species = set(
        path_panel.loc[path_panel["path_analysis_panel24_v1"].eq(True), "species"]
    )
    te = te.loc[te["species"].isin(path_species)]
    diversity_columns = [
        "species",
        "n_superfamilies_included",
        "shannon_entropy",
        "gini_simpson",
    ]
    if list(diversity.columns) != diversity_columns:
        raise ValueError(
            "TE diversity must be the canonical 34-row Wang-style top-10-superfamily release"
        )
    if len(diversity) != 34 or not diversity["species"].is_unique:
        raise ValueError("TE diversity must contain exactly 34 unique TE34 species")
    if not diversity["n_superfamilies_included"].eq(10).all():
        raise ValueError("every canonical diversity row must contain exactly 10 superfamilies")
    overlap = te.merge(diversity, on="species", how="inner", validate="one_to_one")
    overlap = overlap.merge(
        iod[
            [
                "species",
                "genome_size_pg_fuscus_anchored",
                "genome_size_pg_ci_low",
                "genome_size_pg_ci_high",
                "n_nuclei",
                "n_images",
                "n_specimens",
                "estimate_support",
                "sample_size_support",
                "quality_match_status",
            ]
        ],
        on="species",
        how="inner",
        validate="one_to_one",
    )
    overlap.insert(1, "scientific_name", "D. " + overlap["species"])
    overlap["te_diversity_definition"] = TE_DIVERSITY_DEFINITION
    overlap["genome_size_definition"] = GENOME_SIZE_DEFINITION
    overlap = overlap[TE_OVERLAP_SCHEMA].sort_values("species").reset_index(drop=True)
    if len(overlap) != 21 or not overlap["species"].is_unique:
        raise ValueError("TE-genome-size overlap must contain 21 unique species")
    return overlap


def dag_id(edges: Collection[tuple[str, str]]) -> str:
    edge_set = set(edges)
    if not edge_set:
        return "no_edges"
    return ";".join(
        f"{NODE_ABBREVIATIONS[source]}->{NODE_ABBREVIATIONS[target]}"
        for source, target in DIRECTED_EDGES
        if (source, target) in edge_set
    )


def is_acyclic(edges: Collection[tuple[str, str]]) -> bool:
    remaining = set(NODES)
    edge_set = set(edges)
    while remaining:
        roots = {
            node
            for node in remaining
            if not any(
                target == node and source in remaining for source, target in edge_set
            )
        }
        if not roots:
            return False
        remaining -= roots
    return True


def equivalence_class(
    edges: Collection[tuple[str, str]],
) -> tuple[str, str, bool]:
    edge_set = set(edges)
    skeleton = {frozenset(edge) for edge in edge_set}
    n_edges = len(skeleton)
    if n_edges == 0:
        return (
            "independent",
            (
                "genome_size _||_ nucleus_size; genome_size _||_ cell_size; "
                "nucleus_size _||_ cell_size"
            ),
            True,
        )
    if n_edges == 1:
        pair = next(iter(skeleton))
        missing = next(node for node in NODES if node not in pair)
        present = sorted(pair, key=NODES.index)
        class_name = {
            frozenset((GENOME, NUCLEUS)): "genome_nucleus_only",
            frozenset((NUCLEUS, CELL)): "nucleus_cell_only",
            frozenset((GENOME, CELL)): "genome_cell_only",
        }[pair]
        basis = "; ".join(f"{missing} _||_ {node}" for node in present)
        return class_name, basis, True
    if n_edges == 2:
        degrees = {node: sum(node in pair for pair in skeleton) for node in NODES}
        middle = next(node for node, degree in degrees.items() if degree == 2)
        endpoints = [node for node in NODES if node != middle]
        is_collider = all((endpoint, middle) in edge_set for endpoint in endpoints)
        suffix = "collider" if is_collider else "bridge"
        if is_collider:
            basis = f"{endpoints[0]} _||_ {endpoints[1]}"
        else:
            basis = f"{endpoints[0]} _||_ {endpoints[1]} | {middle}"
        middle_label = middle.removesuffix("_size")
        return f"{middle_label}_{suffix}", basis, True
    if n_edges == 3:
        return "saturated", "none (no d-separation claim)", False
    raise AssertionError(f"Unexpected three-node skeleton with {n_edges} edges")


def build_candidate_dags() -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for edge_count in range(len(DIRECTED_EDGES) + 1):
        for edges in combinations(DIRECTED_EDGES, edge_count):
            if not is_acyclic(edges):
                continue
            class_name, basis_claim, testable = equivalence_class(edges)
            records.append(
                {
                    "dag_id": dag_id(edges),
                    "n_edges": edge_count,
                    "equivalence_class": class_name,
                    "basis_claim": basis_claim,
                    "testable_by_dsep": testable,
                }
            )
    frame = pd.DataFrame(records)
    if len(frame) != 25:
        raise AssertionError(f"Expected 25 labeled DAGs, found {len(frame)}")
    frame["class_member_count"] = frame.groupby("equivalence_class")[
        "dag_id"
    ].transform("size")
    return frame.sort_values(
        ["n_edges", "equivalence_class", "dag_id"], kind="stable"
    ).reset_index(drop=True)


def bootstrap_equal_image_intensity(
    frame: pd.DataFrame, *, n_bootstrap: int, seed: int
) -> np.ndarray:
    images = {
        str(filename): group["nuc_iod"].dropna().to_numpy(dtype=float)
        for filename, group in frame.groupby("filename", sort=True)
    }
    images = {name: values for name, values in images.items() if len(values)}
    if not images:
        raise ValueError("No usable image-specific intensity values")
    names = np.array(sorted(images), dtype=object)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        sampled_names = rng.choice(names, size=len(names), replace=True)
        image_medians = []
        for name in sampled_names:
            values = images[str(name)]
            resampled = rng.choice(values, size=len(values), replace=True)
            image_medians.append(float(np.median(resampled)))
        draws[replicate] = float(np.mean(image_medians))
    return draws


def bootstrap_paired_size_medians(
    frame: pd.DataFrame, *, n_bootstrap: int, seed: int
) -> np.ndarray:
    values = frame[["cell_area_um2", "nuc_area_um2"]].to_numpy(dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Size panel contains no usable paired measurements")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(n_bootstrap, len(values)))
    return np.median(values[indices], axis=1)


def conditional_anchor_pg() -> float:
    manifest = json.loads(IOD_MANIFEST.read_text(encoding="utf-8"))
    anchor = manifest["conditional_genome_size_anchor"]
    if float(anchor["reported_assembly_span_gbp"]) != 16.1:
        raise ValueError("Expected the current 16.1-Gbp D. fuscus assembly anchor")
    if float(anchor["gbp_per_pg"]) != 0.978:
        raise ValueError("Expected the audited 0.978-Gbp-per-pg conversion")
    return float(anchor["derived_anchor_pg_1c"])


def build_measurement_bootstrap(
    *, n_bootstrap: int = 250, seed: int = 20260710
) -> pd.DataFrame:
    intensity = read_csv(IOD_OBJECT_SOURCE)
    sizes = read_csv(SIZE_OBJECT_SOURCE)
    for label, source in (("intensity", intensity), ("size", sizes)):
        if not source["review_status"].eq("reviewed_keep").all():
            raise ValueError(f"{label} bootstrap source contains a non-finalized row")
    path_species = {f"D. {tip}" for tip in path_tip_order()}
    if set(intensity["species"]) != path_species or set(sizes["species"]) != path_species:
        raise ValueError("Measurement bootstrap sources do not equal Path24")

    species_order = sorted(path_species)
    intensity_seed_index = {
        species: index for index, species in enumerate(species_order)
    }
    reference_species = "D. fuscus"
    reference_draws = bootstrap_equal_image_intensity(
        intensity.loc[intensity["species"].eq(reference_species)],
        n_bootstrap=n_bootstrap,
        seed=seed + intensity_seed_index[reference_species] * 1009,
    )
    anchor_pg = conditional_anchor_pg()
    frames: list[pd.DataFrame] = []
    for species_index, species in enumerate(species_order):
        size_draws = bootstrap_paired_size_medians(
            sizes.loc[sizes["species"].eq(species)],
            n_bootstrap=n_bootstrap,
            seed=seed + species_index * 2017 + 1,
        )
        species_draws = bootstrap_equal_image_intensity(
            intensity.loc[intensity["species"].eq(species)],
            n_bootstrap=n_bootstrap,
            seed=seed + intensity_seed_index[species] * 1009,
        )
        frames.append(
            pd.DataFrame(
                {
                    "bootstrap_replicate": np.arange(1, n_bootstrap + 1),
                    "species": species,
                    "tree_tip": species.removeprefix("D. "),
                    "nucleus_area_um2": size_draws[:, 1],
                    "cell_area_um2": size_draws[:, 0],
                    "genome_size_pg_fuscus_anchored": (
                        species_draws / reference_draws * anchor_pg
                    ),
                }
            )
        )
    result = pd.concat(frames, ignore_index=True).sort_values(
        ["bootstrap_replicate", "species"], kind="stable"
    ).reset_index(drop=True)
    counts = result.groupby("bootstrap_replicate")["species"].nunique()
    if len(counts) != n_bootstrap or not counts.eq(24).all():
        raise ValueError("Each measurement bootstrap replicate must contain Path24")
    fuscus = result.loc[
        result["species"].eq(reference_species),
        "genome_size_pg_fuscus_anchored",
    ]
    if not np.array_equal(fuscus.to_numpy(), np.full(n_bootstrap, anchor_pg)):
        raise ValueError("D. fuscus bootstrap draws must equal the conditional anchor")
    return result


def validate_frozen_results(path_traits: pd.DataFrame) -> None:
    tables = {name: read_csv(DATA_DIR / name) for name in FROZEN_RESULT_FILES}
    ranking_names = (
        "evolutionary_model_sensitivity.csv",
        "full_simulation_rankings.csv.gz",
        "path_leave_one_species_out.csv",
        "path_measurement_sensitivity.csv",
        "path_model_comparison.csv",
        "path_tree_sensitivity.csv",
    )
    for name in ranking_names:
        if set(tables[name]["model"]) != MODEL_IDS:
            raise ValueError(f"{name} does not contain the ten canonical model IDs")
    if set(tables["path_model_stability.csv"]["model"]) != MODEL_IDS:
        raise ValueError("Path stability table has noncanonical model IDs")
    if set(tables["path_dsep_tests.csv"]["candidate_model"]) != MODEL_IDS:
        raise ValueError("D-separation table has noncanonical candidate IDs")
    regressions = tables["pairwise_phylogenetic_regressions.csv"]
    if set(regressions["relation_id"]) != RELATION_IDS:
        raise ValueError("Pairwise regressions have noncanonical relation IDs")
    if set(tables["path_edge_estimates.csv"]["mechanism"]) != MECHANISM_IDS:
        raise ValueError("Path edges have noncanonical mechanism IDs")

    for name, table in tables.items():
        if "path_trait_definition" in table:
            definitions = set(table["path_trait_definition"].dropna())
            if definitions != {PATH_TRAIT_DEFINITION}:
                raise ValueError(f"{name} has a noncanonical trait definition")

    standardized = tables["path_standardized_traits.csv"]
    if list(standardized.columns) != [
        "species",
        "tree_tip",
        "genome_size_log10_z",
        "nucleus_area_log10_z",
        "cell_area_log10_z",
        "genome_size_definition",
    ]:
        raise ValueError("Standardized traits have a noncanonical schema")
    if set(standardized["genome_size_definition"]) != {STANDARDIZED_DEFINITION}:
        raise ValueError("Standardized traits have a noncanonical definition")

    pairwise = tables["pairwise_pgls.csv"]
    comparisons = {
        "genome_size_vs_nucleus_area",
        "genome_size_vs_cell_area",
        "nucleus_area_vs_cell_area",
    }
    if set(pairwise["comparison"]) != comparisons:
        raise ValueError("Pairwise PGLS table has noncanonical comparison IDs")
    genome_metric = "genome_size_pg_fuscus_anchored"
    genome_rows = pairwise["left_metric"].eq(genome_metric)
    if genome_rows.sum() != 2:
        raise ValueError("Pairwise PGLS table must have two genome-size predictors")
    expected_mean = round(np.log10(path_traits[genome_metric]).mean(), 8)
    if not np.allclose(
        pairwise.loc[genome_rows, "predictor_log10_mean"],
        expected_mean,
        rtol=0,
        atol=5e-9,
    ):
        raise ValueError("Pairwise PGLS genome-size scale metadata is stale")


def reconstruct_release() -> dict[str, pd.DataFrame]:
    return {
        "candidate_dags.csv": build_candidate_dags(),
        "path24_traits.csv": build_path24_traits(),
        "measurement_bootstrap_traits.csv": build_measurement_bootstrap(),
        "te_genome_size_overlap.csv": build_te_genome_size_overlap(),
    }


def write_release(products: dict[str, pd.DataFrame]) -> None:
    for name, table in products.items():
        path = DATA_DIR / name
        if path.suffix == ".gz":
            table.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
        else:
            table.to_csv(path, index=False)


def check_release(products: dict[str, pd.DataFrame]) -> None:
    failures: list[str] = []
    for name, expected in products.items():
        path = DATA_DIR / name
        if not path.exists():
            failures.append(f"missing {path.relative_to(ROOT)}")
            continue
        observed = read_csv(path)
        try:
            pd.testing.assert_frame_equal(observed, expected, check_exact=True, check_dtype=True)
        except AssertionError as error:
            failures.append(f"{name}: {error}")
    if failures:
        raise AssertionError("Path24 release differs from reconstruction:\n" + "\n".join(failures))
    validate_frozen_results(products["path24_traits.csv"])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="promote reconstructed tables into analyses/03_phylogenetic_path/data",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    products = reconstruct_release()
    if args.write:
        write_release(products)
    check_release(products)
    mode = "wrote and validated" if args.write else "validated"
    print(
        f"{mode} {len(products)} reconstructed Path24 tables and validated "
        f"{len(FROZEN_RESULT_FILES)} frozen R result tables"
    )


if __name__ == "__main__":
    main()

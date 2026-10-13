"""Prepare and audit an exploratory Bonsai comparison against the canonical TE34 PCA."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
ROOT = ANALYSIS.parents[1]
TE_DATA = ANALYSIS / "data"
TRACKED_DATA = HERE / "data"

PANEL = ROOT / "data" / "identity" / "te34_panel.csv"
CLR = TE_DATA / "te_pca_clr_matrix.csv"
PCA_SCORES = TE_DATA / "te_pca_scores.csv"
TRAITS = ROOT / "analyses" / "03_phylogenetic_path" / "data" / "organismal_traits.csv"
TRACKED_TREE = TRACKED_DATA / "te34_bonsai_tree.nwk"
TRACKED_METADATA = TRACKED_DATA / "te34_bonsai_upstream_metadata.csv"

UPSTREAM_URL = "https://github.com/dhdegroot/Bonsai-data-representation"
UPSTREAM_COMMIT = "b16ef6f0a164610327022c7d0b92de82b3596aad"
PAPER_DOI = "10.1038/s41587-026-03220-2"
POST_URL = "https://x.com/bravo_abad/status/2091911661183451532"
N_SPECIES = 34
N_FEATURES = 24
N_PCS = 23
K_NEIGHBORS = 5

SUMMARY_FILENAME = "te34_bonsai_distance_fidelity.csv"
BY_SPECIES_FILENAME = "te34_bonsai_distance_fidelity_by_species.csv"
MANIFEST_FILENAME = "te34_bonsai_comparison_manifest.json"
TREE_FILENAME = "te34_bonsai_tree.nwk"
METADATA_FILENAME = "te34_bonsai_upstream_metadata.csv"

REPRESENTATIONS = (
    ("pca_2d", "PCA PC1-PC2", 2),
    ("pca_6d", "PCA PC1-PC6", 6),
    ("pca_full_23d", "PCA all 23 PCs", 23),
)


def fail(code: str, detail: str) -> None:
    raise AssertionError(f"{code}: {detail}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_species() -> list[str]:
    panel = pd.read_csv(PANEL)
    species = panel.loc[panel["te_resource_panel34_v1"], "species"].tolist()
    if len(species) != N_SPECIES or len(species) != len(set(species)):
        fail("BONSAI_PANEL_MISMATCH", "TE34 must contain 34 unique species")
    return species


def load_clr_matrix() -> pd.DataFrame:
    long = pd.read_csv(CLR, float_precision="round_trip")
    long = long.loc[
        long["te_level"].eq("superfamily")
        & long["composition_mode"].eq("classified_conditional")
    ].copy()
    keys = ["species", "feature"]
    if long.duplicated(keys).any():
        fail("BONSAI_CLR_INVALID", "superfamily CLR keys are not unique")

    species = canonical_species()
    features = long["feature"].drop_duplicates().tolist()
    matrix = long.pivot(index="species", columns="feature", values="clr_value")
    matrix = matrix.reindex(index=species, columns=features)
    values = matrix.to_numpy(dtype=float)
    if matrix.shape != (N_SPECIES, N_FEATURES):
        fail("BONSAI_CLR_INVALID", f"expected 34 x 24 CLR matrix, found {matrix.shape}")
    if not np.isfinite(values).all():
        fail("BONSAI_CLR_INVALID", "CLR matrix contains missing or non-finite values")
    if not np.allclose(values.sum(axis=1), 0, rtol=0, atol=2e-12):
        fail("BONSAI_CLR_INVALID", "CLR rows do not sum to zero")
    return matrix


def load_pca_scores(species: list[str]) -> pd.DataFrame:
    scores = pd.read_csv(PCA_SCORES, float_precision="round_trip")
    scores = scores.loc[
        scores["te_level"].eq("superfamily")
        & scores["composition_mode"].eq("classified_conditional")
    ].set_index("species")
    pc_columns = [f"PC{index}" for index in range(1, N_PCS + 1)]
    scores = scores.reindex(index=species).loc[:, pc_columns]
    if scores.shape != (N_SPECIES, N_PCS) or not np.isfinite(
        scores.to_numpy(dtype=float)
    ).all():
        fail("BONSAI_PCA_INVALID", "canonical PCA scores are incomplete")
    return scores


def load_tree(tree_path: Path, species: list[str]) -> Phylo.BaseTree.Tree:
    tree = Phylo.read(tree_path, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != N_SPECIES or set(tips) != set(species) or len(tips) != len(set(tips)):
        fail("BONSAI_TREE_INVALID", "Bonsai tree tips must match TE34 exactly")
    for clade in tree.find_clades():
        length = clade.branch_length
        if length is not None and (not np.isfinite(length) or length < 0):
            fail("BONSAI_TREE_INVALID", "branch lengths must be finite and nonnegative")
    return tree


def squared_euclidean(values: np.ndarray) -> np.ndarray:
    differences = values[:, None, :] - values[None, :, :]
    return np.einsum("ijk,ijk->ij", differences, differences)


def tree_path_distances(tree: Phylo.BaseTree.Tree, species: list[str]) -> np.ndarray:
    result = np.zeros((len(species), len(species)), dtype=float)
    for row, first in enumerate(species):
        for column in range(row + 1, len(species)):
            result[row, column] = result[column, row] = tree.distance(
                first, species[column]
            )
    return result


def pearson_r(first: np.ndarray, second: np.ndarray) -> float:
    if first.ndim != 1 or second.ndim != 1 or first.shape != second.shape:
        fail("BONSAI_METRIC_INVALID", "Pearson inputs must be equal-length vectors")
    if np.ptp(first) == 0 or np.ptp(second) == 0:
        fail("BONSAI_METRIC_INVALID", "Pearson inputs must vary")
    return float(np.corrcoef(first, second)[0, 1])


def nearest_neighbor_recall(
    source: np.ndarray, representation: np.ndarray, k: int = K_NEIGHBORS
) -> np.ndarray:
    n_species = source.shape[0]
    recalls = np.empty(n_species, dtype=float)
    tie_break = np.arange(n_species)
    for index in range(n_species):
        source_order = [
            item
            for item in np.lexsort((tie_break, source[index]))
            if item != index
        ][:k]
        representation_order = [
            item
            for item in np.lexsort((tie_break, representation[index]))
            if item != index
        ][:k]
        recalls[index] = len(set(source_order) & set(representation_order)) / k
    return recalls


def normalized_stress(source: np.ndarray, representation: np.ndarray) -> tuple[float, float]:
    upper = np.triu_indices_from(source, k=1)
    source_vector = source[upper]
    representation_vector = representation[upper]
    denominator = float(np.dot(representation_vector, representation_vector))
    if denominator <= 0:
        fail("BONSAI_METRIC_INVALID", "representation distances have zero norm")
    scale = float(np.dot(source_vector, representation_vector) / denominator)
    residual = source_vector - scale * representation_vector
    stress = float(np.sqrt(np.dot(residual, residual) / np.dot(source_vector, source_vector)))
    return scale, stress


def evaluate_representation(
    source: np.ndarray,
    representation: np.ndarray,
    species: list[str],
    representation_id: str,
    label: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    upper = np.triu_indices_from(source, k=1)
    per_species_r = np.asarray(
        [
            pearson_r(np.delete(source[index], index), np.delete(representation[index], index))
            for index in range(len(species))
        ]
    )
    recalls = nearest_neighbor_recall(source, representation)
    scale, stress = normalized_stress(source, representation)
    summary = {
        "representation_id": representation_id,
        "representation_label": label,
        "n_species": len(species),
        "n_pairwise_distances": len(source[upper]),
        "global_distance_pearson_r": pearson_r(source[upper], representation[upper]),
        "median_species_distance_pearson_r": float(np.median(per_species_r)),
        "minimum_species_distance_pearson_r": float(per_species_r.min()),
        "maximum_species_distance_pearson_r": float(per_species_r.max()),
        "mean_five_neighbor_recall": float(recalls.mean()),
        "least_squares_distance_scale": scale,
        "normalized_stress_after_scale": stress,
    }
    by_species = [
        {
            "species": name,
            "representation_id": representation_id,
            "representation_label": label,
            "distance_pearson_r": per_species_r[index],
            "five_neighbor_recall": recalls[index],
        }
        for index, name in enumerate(species)
    ]
    return summary, by_species


def read_upstream_metadata(metadata_path: Path) -> dict[str, object]:
    metadata = pd.read_csv(metadata_path)
    if len(metadata) != 1:
        fail("BONSAI_METADATA_INVALID", "expected one upstream metadata row")
    row = metadata.iloc[0]
    required = ["nGenes", "nCells", "reconsTime", "zscoreCutoff", "coresUsed", "finalLoglik"]
    if any(column not in metadata.columns for column in required):
        fail("BONSAI_METADATA_INVALID", "upstream metadata columns are incomplete")
    if int(row["nGenes"]) != N_FEATURES or int(row["nCells"]) != N_SPECIES:
        fail("BONSAI_METADATA_INVALID", "upstream metadata dimensions do not match TE34")
    return {
        "n_features": int(row["nGenes"]),
        "n_species": int(row["nCells"]),
        "reconstruction_seconds_reported_by_bonsai": float(row["reconsTime"]),
        "zscore_cutoff": float(row["zscoreCutoff"]),
        "cores_reported_by_bonsai": int(row["coresUsed"]),
        "final_log_likelihood": float(row["finalLoglik"]),
    }


def build_analysis(tree_path: Path, metadata_path: Path) -> dict[str, object]:
    clr = load_clr_matrix()
    species = clr.index.tolist()
    pca = load_pca_scores(species)
    tree = load_tree(tree_path, species)

    source = squared_euclidean(clr.to_numpy(dtype=float))
    summaries: list[dict[str, object]] = []
    by_species: list[dict[str, object]] = []
    for representation_id, label, n_components in REPRESENTATIONS:
        distances = squared_euclidean(pca.iloc[:, :n_components].to_numpy(dtype=float))
        summary, species_rows = evaluate_representation(
            source, distances, species, representation_id, label
        )
        summaries.append(summary)
        by_species.extend(species_rows)

    bonsai_distances = tree_path_distances(tree, species)
    summary, species_rows = evaluate_representation(
        source, bonsai_distances, species, "bonsai_tree", "Bonsai tree"
    )
    summaries.append(summary)
    by_species.extend(species_rows)

    full_pca = squared_euclidean(pca.to_numpy(dtype=float))
    maximum_full_pca_difference = float(np.max(np.abs(source - full_pca)))
    if maximum_full_pca_difference >= 1e-12:
        fail("BONSAI_GEOMETRY_MISMATCH", "all released PCs do not reproduce CLR distances")

    summary_frame = pd.DataFrame(summaries)
    species_frame = pd.DataFrame(by_species)
    metadata = read_upstream_metadata(metadata_path)
    zero_length_branches = sum(
        clade.branch_length == 0 for clade in tree.find_clades() if clade.branch_length is not None
    )

    return {
        "summary": summary_frame,
        "by_species": species_frame,
        "metadata": metadata,
        "maximum_full_pca_squared_distance_difference": maximum_full_pca_difference,
        "zero_length_branches": zero_length_branches,
    }


def manifest(
    analysis: dict[str, object], tree_path: Path, metadata_path: Path
) -> dict[str, object]:
    summary = analysis["summary"]
    assert isinstance(summary, pd.DataFrame)
    result_summary = {
        row["representation_id"]: {
            "global_distance_pearson_r": row["global_distance_pearson_r"],
            "median_species_distance_pearson_r": row[
                "median_species_distance_pearson_r"
            ],
            "mean_five_neighbor_recall": row["mean_five_neighbor_recall"],
            "normalized_stress_after_scale": row["normalized_stress_after_scale"],
        }
        for row in summary.to_dict(orient="records")
    }
    input_paths = {
        "te34_panel": PANEL,
        "superfamily_clr": CLR,
        "pca_scores": PCA_SCORES,
        "organismal_traits": TRAITS,
        "bonsai_tree": tree_path,
        "bonsai_metadata": metadata_path,
        "comparison_script": Path(__file__).resolve(),
        "plot_script": HERE / "plot_bonsai_comparison.R",
        "upstream_runner": HERE / "run_upstream_bonsai.sh",
        "python_project": ROOT / "pyproject.toml",
        "python_lock": ROOT / "uv.lock",
    }
    inputs = {
        label: {"path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path), "sha256": sha256(path)}
        for label, path in input_paths.items()
    }
    return {
        "analysis_id": "te34_bonsai_vs_pca_exploratory_v1",
        "analysis_date": "2026-08-25",
        "analysis_status": "exploratory_method_comparison_with_validated_publication_representation",
        "source": {
            "request_post": POST_URL,
            "paper_doi": PAPER_DOI,
            "paper_title": "Bonsai reconstructs tree representations for distortion-free visualization and exploration of high-dimensional data",
            "upstream_repository": UPSTREAM_URL,
            "upstream_commit": UPSTREAM_COMMIT,
            "upstream_license": "CC-BY-NC-4.0",
        },
        "inputs": inputs,
        "method_adaptation": {
            "objects": "34 TE34 species",
            "features": "24 superfamilies positive in every TE34 species",
            "feature_values": "canonical classified-conditional CLR values",
            "measurement_standard_deviations": None,
            "measurement_uncertainty_treatment": "omitted; upstream assumes a very small error bar",
            "input_is_sanity_output": False,
            "zscore_cutoff": -1,
            "rescale_by_var": False,
            "rescale_reason": "preserve the same unscaled CLR/Aitchison geometry used by the canonical PCA",
            "use_knn": -1,
            "use_knn_reason": "exact all-pairs candidate search is tractable for 34 species",
            "nnn_n_randommoves": 1000,
            "nnn_n_randomtrees": 10,
            "UB_ellipsoid_size": 1.0,
        },
        "distance_fidelity": {
            "source_distance": "squared Euclidean distance in the 24-feature CLR matrix",
            "pca_distance": "squared Euclidean distance in the stated PC score subspace",
            "bonsai_distance": "sum of fitted diffusion-time branch lengths along the tree path",
            "primary_summary": "Pearson correlation across distances from each species to the other 33 species",
            "secondary_summaries": [
                "global pairwise Pearson correlation",
                "five-nearest-neighbor recall",
                "normalized least-squares stress after one global scale factor",
            ],
            "maximum_full_pca_squared_distance_difference": analysis[
                "maximum_full_pca_squared_distance_difference"
            ],
        },
        "upstream_run": analysis["metadata"],
        "tree_diagnostics": {
            "tips": N_SPECIES,
            "zero_length_branches": analysis["zero_length_branches"],
            "bootstrap_or_measurement_error_sensitivity": (
                "not part of this fit; see the companion bounded feature-bootstrap "
                "strategy validation"
            ),
        },
        "results": result_summary,
        "interpretation_boundaries": [
            "The Bonsai tree is a representation of TE-composition similarity, not an inferred species phylogeny.",
            "This 34-species use is far smaller than the single-cell applications emphasized by the paper.",
            "No defensible per-species feature-error matrix exists, so the run does not exercise Bonsai's uncertainty model.",
            "The direct visual comparison is Bonsai versus PC1-PC2; all 23 PCs remain the exact distance-preserving PCA control.",
            "A companion 200-replicate feature-resampling analysis measures topology sensitivity but is not a biological-specimen bootstrap.",
            "The method comparison remains exploratory and does not alter canonical PCA files; compact validated Bonsai outputs are promoted separately to the controlled publication release.",
        ],
    }


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def write_analysis(
    analysis: dict[str, object], tree_path: Path, metadata_path: Path, output_dir: Path
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = analysis["summary"]
    by_species = analysis["by_species"]
    assert isinstance(summary, pd.DataFrame)
    assert isinstance(by_species, pd.DataFrame)
    write_frame(summary, output_dir / SUMMARY_FILENAME)
    write_frame(by_species, output_dir / BY_SPECIES_FILENAME)

    tree_destination = output_dir / TREE_FILENAME
    metadata_destination = output_dir / METADATA_FILENAME
    if tree_path.resolve() != tree_destination.resolve():
        shutil.copyfile(tree_path, tree_destination)
    if metadata_path.resolve() != metadata_destination.resolve():
        shutil.copyfile(metadata_path, metadata_destination)

    output_manifest = manifest(analysis, tree_destination, metadata_destination)
    (output_dir / MANIFEST_FILENAME).write_text(
        json.dumps(output_manifest, indent=2, sort_keys=False) + "\n"
    )


def check_frame(observed: pd.DataFrame, expected_path: Path) -> None:
    expected = pd.read_csv(expected_path, float_precision="round_trip")
    pd.testing.assert_frame_equal(
        observed,
        expected,
        check_exact=False,
        rtol=2e-13,
        atol=2e-15,
    )


def check_tracked_analysis(analysis: dict[str, object]) -> None:
    summary = analysis["summary"]
    by_species = analysis["by_species"]
    assert isinstance(summary, pd.DataFrame)
    assert isinstance(by_species, pd.DataFrame)
    check_frame(summary, TRACKED_DATA / SUMMARY_FILENAME)
    check_frame(by_species, TRACKED_DATA / BY_SPECIES_FILENAME)
    observed_manifest = manifest(analysis, TRACKED_TREE, TRACKED_METADATA)
    expected_manifest = json.loads((TRACKED_DATA / MANIFEST_FILENAME).read_text())
    if observed_manifest != expected_manifest:
        fail("BONSAI_MANIFEST_MISMATCH", "tracked comparison manifest is stale")


def prepare_run_directory(run_dir: Path) -> None:
    generated_paths = [
        run_dir / "input" / "features.txt",
        run_dir / "input" / "cellID.txt",
        run_dir / "input" / "geneID.txt",
        run_dir / "bonsai_config.yaml",
    ]
    existing = [path for path in generated_paths if path.exists()]
    if existing:
        fail(
            "BONSAI_RUN_DIR_NOT_EMPTY",
            "refusing to overwrite: " + ", ".join(str(path) for path in existing),
        )

    clr = load_clr_matrix()
    input_dir = run_dir / "input"
    results_dir = run_dir / "results"
    input_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    clr.T.to_csv(
        input_dir / "features.txt",
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )
    (input_dir / "cellID.txt").write_text("\n".join(clr.index) + "\n")
    (input_dir / "geneID.txt").write_text("\n".join(clr.columns) + "\n")
    config = {
        "dataset": "te34_superfamily_clr",
        "data_folder": str(input_dir.resolve()),
        "filenames_data": "features.txt",
        "results_folder": str(results_dir.resolve()),
        "verbose": False,
        "input_is_sanity_output": False,
        "zscore_cutoff": -1.0,
        "rescale_by_var": False,
        "nnn_n_randommoves": 1000,
        "nnn_n_randomtrees": 10,
        "use_knn": -1,
        "UB_ellipsoid_size": 1.0,
        "spr_strategy": "large_tree",
        "pickup_intermediate": False,
        "tmp_folder": "",
    }
    lines = [
        "# Generated for the exploratory TE34 Bonsai comparison.",
        f"# Upstream commit: {UPSTREAM_COMMIT}",
    ]
    for key, value in config.items():
        if isinstance(value, bool):
            rendered = str(value).lower()
        elif isinstance(value, str):
            rendered = json.dumps(value)
        else:
            rendered = str(value)
        lines.append(f"{key}: {rendered}")
    (run_dir / "bonsai_config.yaml").write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="write Bonsai input and pinned run config")
    prepare.add_argument("--run-dir", type=Path, required=True)

    analyze = subparsers.add_parser("analyze", help="analyze a completed upstream tree")
    analyze.add_argument("--tree", type=Path, required=True)
    analyze.add_argument("--metadata", type=Path, required=True)
    analyze.add_argument("--output-dir", type=Path, required=True)
    analyze.add_argument("--write", action="store_true")

    subparsers.add_parser("check", help="replay and verify the tracked exploratory result")
    return parser.parse_args()


def print_summary(summary: pd.DataFrame) -> None:
    columns = [
        "representation_label",
        "global_distance_pearson_r",
        "median_species_distance_pearson_r",
        "mean_five_neighbor_recall",
        "normalized_stress_after_scale",
    ]
    print(summary.loc[:, columns].to_string(index=False, float_format=lambda value: f"{value:.4f}"))


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare_run_directory(args.run_dir)
        print(f"PASS: wrote exact TE34 Bonsai inputs and config under {args.run_dir}")
        return

    if args.command == "check":
        analysis = build_analysis(TRACKED_TREE, TRACKED_METADATA)
        check_tracked_analysis(analysis)
        print_summary(analysis["summary"])
        print("PASS: tracked exploratory Bonsai comparison reproduces exactly")
        return

    analysis = build_analysis(args.tree, args.metadata)
    print_summary(analysis["summary"])
    if args.write:
        write_analysis(analysis, args.tree, args.metadata, args.output_dir)
        print(f"PASS: wrote exploratory Bonsai comparison under {args.output_dir}")
    else:
        print("PASS: analysis completed read-only; pass --write to store outputs")


if __name__ == "__main__":
    main()

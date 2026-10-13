"""Validate reproductive-strategy localization on the exploratory TE34 Bonsai tree."""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import io
import json
import shutil
from pathlib import Path
from typing import Any

import compare_bonsai_to_pca as comparison
import numpy as np
import pandas as pd
from Bio import Phylo

HERE = Path(__file__).resolve().parent
ROOT = comparison.ROOT
TRACKED_DATA = HERE / "data"
REFERENCE_TREE = TRACKED_DATA / comparison.TREE_FILENAME

STRATEGIES = (
    "aquatic_eggs",
    "direct_development",
    "terrestrial_eggs_aquatic_larvae",
)
STRATEGY_LABELS = {
    "aquatic_eggs": "Aquatic eggs",
    "direct_development": "Direct development",
    "terrestrial_eggs_aquatic_larvae": "Terrestrial eggs + aquatic larvae",
}
STRATEGY_CODES = {strategy: index + 1 for index, strategy in enumerate(STRATEGIES)}

PERMUTATIONS = 99_999
PERMUTATION_SEED = 20260825
BOOTSTRAP_REPLICATES = 200
BOOTSTRAP_SEED = 20260826
BOOTSTRAP_RANDOM_MOVES = 100
BOOTSTRAP_RANDOM_TREES = 1

PERMUTATION_FILENAME = "te34_bonsai_strategy_label_permutation_tests.csv"
PERMUTATION_HISTOGRAM_FILENAME = "te34_bonsai_strategy_label_permutation_fitch_null.csv"
BOOTSTRAP_SUMMARY_FILENAME = "te34_bonsai_feature_bootstrap_summary.csv"
BOOTSTRAP_REPLICATE_FILENAME = "te34_bonsai_feature_bootstrap_replicates.csv"
BOOTSTRAP_STRATEGY_FILENAME = "te34_bonsai_feature_bootstrap_strategy_stability.csv"
BOOTSTRAP_SPLIT_FILENAME = "te34_bonsai_feature_bootstrap_split_support.csv"
BOOTSTRAP_DRAWS_FILENAME = "te34_bonsai_feature_bootstrap_draws.csv"
BOOTSTRAP_PLAN_FILENAME = "te34_bonsai_feature_bootstrap_plan.json"
BOOTSTRAP_ARCHIVE_FILENAME = "te34_bonsai_feature_bootstrap_trees.nwk.gz"
BOOTSTRAP_INDEX_FILENAME = "te34_bonsai_feature_bootstrap_tree_index.csv"
SUPPORT_TREE_FILENAME = "te34_bonsai_tree_with_feature_bootstrap_support.nwk"
TIP_ANNOTATION_FILENAME = "te34_bonsai_tip_annotations.csv"
PUBLICATION_SUMMARY_FILENAME = "te34_bonsai_strategy_publication_summary.csv"
VALIDATION_MANIFEST_FILENAME = "te34_bonsai_strategy_validation_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_tip_names(tree: Phylo.BaseTree.Tree) -> list[str]:
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != comparison.N_SPECIES or len(tips) != len(set(tips)):
        comparison.fail("BONSAI_VALIDATION_TREE_INVALID", "tree must have 34 unique tips")
    return tips


def validate_tree(tree: Phylo.BaseTree.Tree, species: list[str]) -> None:
    tips = tree_tip_names(tree)
    if set(tips) != set(species):
        comparison.fail("BONSAI_VALIDATION_TREE_INVALID", "tree tips do not match TE34")
    for clade in tree.find_clades():
        length = clade.branch_length
        if length is not None and (not np.isfinite(length) or length < 0):
            comparison.fail(
                "BONSAI_VALIDATION_TREE_INVALID", "branch lengths must be finite and nonnegative"
            )


def load_strategy_codes(species: list[str]) -> tuple[np.ndarray, dict[str, set[str]]]:
    traits = pd.read_csv(comparison.TRAITS, dtype=str, keep_default_na=False).set_index("species")
    codes = np.zeros(len(species), dtype=np.int8)
    groups: dict[str, set[str]] = {}
    for strategy, code in STRATEGY_CODES.items():
        group = set(traits.index[traits["reproductive_strategy"].eq(strategy)]) & set(species)
        groups[strategy] = group
        for member in group:
            codes[species.index(member)] = code
    observed_counts = {strategy: len(groups[strategy]) for strategy in STRATEGIES}
    if observed_counts != {
        "aquatic_eggs": 7,
        "direct_development": 3,
        "terrestrial_eggs_aquatic_larvae": 24,
    }:
        comparison.fail(
            "BONSAI_VALIDATION_LABEL_INVALID", f"unexpected TE34 strategy counts: {observed_counts}"
        )
    if int(np.count_nonzero(codes)) != comparison.N_SPECIES:
        comparison.fail("BONSAI_VALIDATION_LABEL_INVALID", "all 34 TE34 tips must be assigned")
    return codes, groups


def tip_annotation_frame(species: list[str], codes: np.ndarray) -> pd.DataFrame:
    code_to_strategy = {code: strategy for strategy, code in STRATEGY_CODES.items()}
    rows = []
    for name, code_value in zip(species, codes, strict=True):
        code = int(code_value)
        strategy = code_to_strategy.get(code, "unassigned")
        rows.append(
            {
                "species": name,
                "display_species": f"D. {name}",
                "reproductive_strategy": strategy,
                "reproductive_strategy_label": STRATEGY_LABELS.get(
                    strategy, "Unassigned"
                ),
                "strategy_assigned": code != 0,
                "annotation_timing": "expert_post_hoc_after_initial_pca_review",
            }
        )
    frame = pd.DataFrame(rows)
    if frame.shape != (comparison.N_SPECIES, 6) or not frame["species"].is_unique:
        comparison.fail(
            "BONSAI_VALIDATION_ANNOTATION_INVALID",
            "publication tip annotations must contain 34 unique species",
        )
    return frame


def publication_summary_frame(
    label_products: dict[str, pd.DataFrame], bootstrap: dict[str, object]
) -> pd.DataFrame:
    tests = label_products["tests"].set_index("test_id")
    strategies = bootstrap["strategies"].set_index("strategy")
    global_test = tests.loc["global_fitch_parsimony_steps"]
    rows = []
    for strategy in STRATEGIES:
        strategy_test = tests.loc[f"{strategy}_smallest_containing_side_tips"]
        stability = strategies.loc[strategy]
        rows.append(
            {
                "strategy": strategy,
                "strategy_label": STRATEGY_LABELS[strategy],
                "group_n": int(stability["group_n"]),
                "full_fit_smallest_containing_side_tips": int(
                    stability["full_fit_smallest_containing_side_tips"]
                ),
                "full_fit_edge_delimited": bool(
                    stability["full_fit_edge_delimited"]
                ),
                "fixed_tree_smallest_side_p_value": float(strategy_test["p_value"]),
                "fixed_tree_smallest_side_p_value_holm": float(
                    strategy_test["p_value_holm_all_nine_tests"]
                ),
                "feature_bootstrap_edge_delimited_count": int(
                    stability["bootstrap_edge_delimited_count"]
                ),
                "feature_bootstrap_replicates": int(
                    stability["bootstrap_replicates"]
                ),
                "feature_bootstrap_edge_delimited_fraction": float(
                    stability["bootstrap_edge_delimited_fraction"]
                ),
                "feature_bootstrap_smallest_side_median": float(
                    stability["bootstrap_smallest_side_median"]
                ),
                "global_fitch_observed_steps": int(global_test["observed"]),
                "global_fitch_null_mean_steps": float(global_test["null_mean"]),
                "global_fitch_p_value": float(global_test["p_value"]),
                "global_fitch_p_value_holm": float(
                    global_test["p_value_holm_all_nine_tests"]
                ),
            }
        )
    return pd.DataFrame(rows)


def canonical_split(side: set[str], all_tips: set[str]) -> frozenset[str]:
    complement = all_tips - side
    if len(side) < len(complement):
        selected = side
    elif len(complement) < len(side):
        selected = complement
    else:
        selected = side if tuple(sorted(side)) <= tuple(sorted(complement)) else complement
    return frozenset(selected)


def unrooted_splits(tree: Phylo.BaseTree.Tree, all_tips: set[str]) -> set[frozenset[str]]:
    splits: set[frozenset[str]] = set()
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        side = {tip.name for tip in clade.get_terminals()}
        if 1 < len(side) < len(all_tips) - 1:
            splits.add(canonical_split(side, all_tips))
    return splits


def split_side_masks(splits: set[frozenset[str]], species: list[str]) -> list[int]:
    index = {name: position for position, name in enumerate(species)}
    full_mask = (1 << len(species)) - 1
    masks = {full_mask}
    for split in splits:
        mask = sum(1 << index[name] for name in split)
        masks.add(mask)
        masks.add(full_mask ^ mask)
    return sorted(masks, key=lambda value: (value.bit_count(), value))


def smallest_containing_side(group_indices: np.ndarray, side_masks: list[int]) -> int:
    group_mask = sum(1 << int(index) for index in group_indices)
    return min(mask.bit_count() for mask in side_masks if group_mask & ~mask == 0)


def tree_distance_matrix(tree: Phylo.BaseTree.Tree, species: list[str]) -> np.ndarray:
    result = np.zeros((len(species), len(species)), dtype=float)
    for first in range(len(species)):
        for second in range(first + 1, len(species)):
            result[first, second] = result[second, first] = tree.distance(
                species[first], species[second]
            )
    return result


def compile_fitch_tree(
    clade: Phylo.BaseTree.Clade, species_index: dict[str, int]
) -> int | tuple[Any, ...]:
    if clade.is_terminal():
        return species_index[clade.name]
    return tuple(compile_fitch_tree(child, species_index) for child in clade.clades)


def fitch_score(node: int | tuple[Any, ...], codes: np.ndarray) -> tuple[int, int]:
    if isinstance(node, int):
        code = int(codes[node])
        return (0 if code == 0 else 1 << (code - 1)), 0
    child_results = [fitch_score(child, codes) for child in node]
    score = sum(result[1] for result in child_results)
    states = [result[0] for result in child_results if result[0] != 0]
    if not states:
        return 0, score
    combined = states[0]
    for state in states[1:]:
        intersection = combined & state
        if intersection:
            combined = intersection
        else:
            combined |= state
            score += 1
    return combined, score


def mean_within_distance(distance: np.ndarray, indices: np.ndarray) -> float:
    if len(indices) < 2:
        comparison.fail("BONSAI_VALIDATION_METRIC_INVALID", "strategy has fewer than two tips")
    block = distance[np.ix_(indices, indices)]
    return float(block[np.triu_indices(len(indices), k=1)].mean())


class LocalizationContext:
    def __init__(
        self,
        tree: Phylo.BaseTree.Tree,
        species: list[str],
        observed_codes: np.ndarray,
    ) -> None:
        self.species = species
        self.observed_codes = observed_codes
        self.distance = tree_distance_matrix(tree, species)
        self.fitch_tree = compile_fitch_tree(tree.root, {name: i for i, name in enumerate(species)})
        self.side_masks = split_side_masks(unrooted_splits(tree, set(species)), species)
        self.assigned_indices = np.flatnonzero(observed_codes)
        nearest: list[int] = []
        for index in self.assigned_indices:
            candidates = [item for item in self.assigned_indices if item != index]
            nearest.append(min(candidates, key=lambda item: (self.distance[index, item], species[item])))
        self.nearest_assigned = np.asarray(nearest, dtype=int)

    def metrics(self, codes: np.ndarray) -> dict[str, float]:
        _, parsimony_steps = fitch_score(self.fitch_tree, codes)
        metrics: dict[str, float] = {
            "global_fitch_parsimony_steps": float(parsimony_steps),
            "global_one_nearest_neighbor_agreement": float(
                np.mean(codes[self.assigned_indices] == codes[self.nearest_assigned])
            ),
        }
        strategy_means: list[float] = []
        for strategy, code in STRATEGY_CODES.items():
            indices = np.flatnonzero(codes == code)
            within = mean_within_distance(self.distance, indices)
            strategy_means.append(within)
            metrics[f"{strategy}_smallest_containing_side_tips"] = float(
                smallest_containing_side(indices, self.side_masks)
            )
            metrics[f"{strategy}_mean_within_patristic_distance"] = within
        metrics["global_mean_of_strategy_mean_within_patristic_distance"] = float(
            np.mean(strategy_means)
        )
        return metrics


def permutation_test_definitions() -> list[dict[str, str]]:
    definitions = [
        {
            "test_id": "global_fitch_parsimony_steps",
            "scope": "global",
            "strategy": "all_assigned",
            "statistic": "Fitch parsimony steps",
            "alternative": "lower",
        },
        {
            "test_id": "global_one_nearest_neighbor_agreement",
            "scope": "global",
            "strategy": "all_assigned",
            "statistic": "fraction whose nearest assigned tree neighbor has the same strategy",
            "alternative": "higher",
        },
        {
            "test_id": "global_mean_of_strategy_mean_within_patristic_distance",
            "scope": "global",
            "strategy": "all_assigned",
            "statistic": "unweighted mean of the three within-strategy mean patristic distances",
            "alternative": "lower",
        },
    ]
    for strategy in STRATEGIES:
        definitions.extend(
            [
                {
                    "test_id": f"{strategy}_smallest_containing_side_tips",
                    "scope": "strategy",
                    "strategy": strategy,
                    "statistic": "smallest edge-delimited tree side containing every group member",
                    "alternative": "lower",
                },
                {
                    "test_id": f"{strategy}_mean_within_patristic_distance",
                    "scope": "strategy",
                    "strategy": strategy,
                    "statistic": "mean within-group patristic distance",
                    "alternative": "lower",
                },
            ]
        )
    return definitions


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ordered = p_values[order]
    adjusted_ordered = np.maximum.accumulate((len(p_values) - np.arange(len(p_values))) * ordered)
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(adjusted_ordered, 1.0)
    return adjusted


def run_label_permutations(
    tree: Phylo.BaseTree.Tree,
    species: list[str],
    observed_codes: np.ndarray,
    permutations: int = PERMUTATIONS,
    seed: int = PERMUTATION_SEED,
) -> dict[str, pd.DataFrame]:
    context = LocalizationContext(tree, species, observed_codes)
    observed = context.metrics(observed_codes)
    definitions = permutation_test_definitions()
    test_ids = [definition["test_id"] for definition in definitions]
    null = np.empty((permutations, len(definitions)), dtype=float)
    rng = np.random.default_rng(seed)
    assigned_labels = observed_codes[context.assigned_indices].copy()
    permuted_codes = observed_codes.copy()
    for replicate in range(permutations):
        permuted_codes[context.assigned_indices] = rng.permutation(assigned_labels)
        metrics = context.metrics(permuted_codes)
        null[replicate] = [metrics[test_id] for test_id in test_ids]

    rows: list[dict[str, object]] = []
    raw_p_values: list[float] = []
    for index, definition in enumerate(definitions):
        observed_value = observed[definition["test_id"]]
        null_values = null[:, index]
        if definition["alternative"] == "lower":
            extreme = int(np.count_nonzero(null_values <= observed_value))
        else:
            extreme = int(np.count_nonzero(null_values >= observed_value))
        p_value = (extreme + 1) / (permutations + 1)
        raw_p_values.append(p_value)
        rows.append(
            {
                **definition,
                "observed": observed_value,
                "null_mean": float(null_values.mean()),
                "null_standard_deviation": float(null_values.std(ddof=1)),
                "null_q025": float(np.quantile(null_values, 0.025)),
                "null_q975": float(np.quantile(null_values, 0.975)),
                "permutations": permutations,
                "seed": seed,
                "extreme_permutations": extreme,
                "p_value": p_value,
            }
        )
    adjusted = holm_adjust(np.asarray(raw_p_values))
    for row, adjusted_value in zip(rows, adjusted, strict=True):
        row["p_value_holm_all_nine_tests"] = float(adjusted_value)
        strategy = row["strategy"]
        row["group_n"] = (
            int(np.count_nonzero(observed_codes == STRATEGY_CODES[strategy]))
            if strategy in STRATEGY_CODES
            else int(np.count_nonzero(observed_codes))
        )

    fitch_values, counts = np.unique(null[:, test_ids.index("global_fitch_parsimony_steps")], return_counts=True)
    histogram = pd.DataFrame(
        {
            "fitch_parsimony_steps": fitch_values.astype(int),
            "permutation_count": counts.astype(int),
            "permutation_fraction": counts / permutations,
        }
    )
    return {"tests": pd.DataFrame(rows), "fitch_histogram": histogram}


def bootstrap_plan(replicates: int, seed: int) -> dict[str, object]:
    return {
        "analysis_id": "te34_bonsai_feature_bootstrap_v1",
        "reference_tree": str(REFERENCE_TREE.relative_to(ROOT)),
        "reference_tree_sha256": sha256(REFERENCE_TREE),
        "source_clr": str(comparison.CLR.relative_to(ROOT)),
        "source_clr_sha256": sha256(comparison.CLR),
        "replicates": replicates,
        "seed": seed,
        "sampling_unit": "24 CLR-superfamily feature coordinates sampled with replacement",
        "calibration_replicate": "replicate_000 uses all 24 features exactly once",
        "bounded_bonsai_config": {
            "input_is_sanity_output": False,
            "measurement_standard_deviations": None,
            "zscore_cutoff": -1.0,
            "rescale_by_var": False,
            "use_knn": -1,
            "nnn_n_randommoves": BOOTSTRAP_RANDOM_MOVES,
            "nnn_n_randomtrees": BOOTSTRAP_RANDOM_TREES,
            "UB_ellipsoid_size": 1.0,
        },
        "upstream_repository": comparison.UPSTREAM_URL,
        "upstream_commit": comparison.UPSTREAM_COMMIT,
        "upstream_license": "CC-BY-NC-4.0",
    }


def render_yaml(config: dict[str, object]) -> str:
    lines = [
        "# Generated for the bounded TE34 Bonsai feature bootstrap.",
        f"# Upstream commit: {comparison.UPSTREAM_COMMIT}",
    ]
    for key, value in config.items():
        if isinstance(value, bool):
            rendered = str(value).lower()
        elif isinstance(value, str):
            rendered = json.dumps(value)
        else:
            rendered = str(value)
        lines.append(f"{key}: {rendered}")
    return "\n".join(lines) + "\n"


def prepare_feature_bootstrap(
    run_dir: Path,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
) -> None:
    plan = bootstrap_plan(replicates, seed)
    plan_path = run_dir / BOOTSTRAP_PLAN_FILENAME
    draws_path = run_dir / BOOTSTRAP_DRAWS_FILENAME
    if plan_path.exists():
        if json.loads(plan_path.read_text()) != plan:
            comparison.fail("BONSAI_BOOTSTRAP_PLAN_MISMATCH", "existing bootstrap plan differs")
        if not draws_path.exists():
            comparison.fail("BONSAI_BOOTSTRAP_PLAN_MISMATCH", "existing feature draws are missing")
        print(f"PASS: existing {replicates}-replicate bootstrap plan is resumable")
        return

    clr = comparison.load_clr_matrix()
    features = np.asarray(clr.columns, dtype=object)
    rng = np.random.default_rng(seed)
    draws = np.vstack(
        [np.arange(len(features), dtype=int), rng.integers(0, len(features), size=(replicates, len(features)))]
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    draw_rows: list[dict[str, object]] = []
    for replicate_id, draw in enumerate(draws):
        replicate_dir = run_dir / "bootstrap" / f"replicate_{replicate_id:03d}"
        input_dir = replicate_dir / "input"
        results_dir = replicate_dir / "results"
        input_dir.mkdir(parents=True, exist_ok=True)
        results_dir.mkdir(parents=True, exist_ok=True)
        selected = clr.iloc[:, draw]
        selected.T.to_csv(
            input_dir / "features.txt",
            sep="\t",
            header=False,
            index=False,
            lineterminator="\n",
            float_format="%.17g",
        )
        (input_dir / "cellID.txt").write_text("\n".join(clr.index) + "\n")
        feature_ids = [f"{features[index]}__draw_{position + 1:02d}" for position, index in enumerate(draw)]
        (input_dir / "geneID.txt").write_text("\n".join(feature_ids) + "\n")
        config = {
            "dataset": f"te34_feature_bootstrap_{replicate_id:03d}",
            "data_folder": str(input_dir.resolve()),
            "filenames_data": "features.txt",
            "results_folder": str(results_dir.resolve()),
            "verbose": False,
            "input_is_sanity_output": False,
            "zscore_cutoff": -1.0,
            "rescale_by_var": False,
            "nnn_n_randommoves": BOOTSTRAP_RANDOM_MOVES,
            "nnn_n_randomtrees": BOOTSTRAP_RANDOM_TREES,
            "use_knn": -1,
            "UB_ellipsoid_size": 1.0,
            "spr_strategy": "large_tree",
            "pickup_intermediate": False,
            "tmp_folder": "",
        }
        (replicate_dir / "bonsai_config.yaml").write_text(render_yaml(config))
        draw_rows.append(
            {
                "replicate_id": replicate_id,
                "is_full_feature_calibration": replicate_id == 0,
                "n_unique_features": len(set(draw)),
                **{
                    f"draw_{position + 1:02d}": features[index]
                    for position, index in enumerate(draw)
                },
            }
        )
    pd.DataFrame(draw_rows).to_csv(draws_path, index=False, lineterminator="\n")
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    print(f"PASS: prepared calibration plus {replicates} feature-bootstrap replicates")


def read_upstream_runtime(metadata_path: Path) -> float:
    metadata = pd.read_csv(metadata_path)
    if len(metadata) != 1 or "reconsTime" not in metadata:
        comparison.fail("BONSAI_BOOTSTRAP_METADATA_INVALID", str(metadata_path))
    return float(metadata.iloc[0]["reconsTime"])


def collect_run_trees(
    run_dir: Path, replicates: int
) -> tuple[dict[int, Phylo.BaseTree.Tree], dict[int, float]]:
    trees: dict[int, Phylo.BaseTree.Tree] = {}
    runtimes: dict[int, float] = {}
    species = comparison.canonical_species()
    missing: list[int] = []
    for replicate_id in range(replicates + 1):
        result_dir = run_dir / "bootstrap" / f"replicate_{replicate_id:03d}" / "results" / "final_bonsai"
        tree_path = result_dir / "tree.nwk"
        metadata_path = result_dir / "metadata_bonsai.txt"
        if not tree_path.exists() or not metadata_path.exists():
            missing.append(replicate_id)
            continue
        tree = Phylo.read(tree_path, "newick")
        validate_tree(tree, species)
        trees[replicate_id] = tree
        runtimes[replicate_id] = read_upstream_runtime(metadata_path)
    if missing:
        comparison.fail(
            "BONSAI_BOOTSTRAP_INCOMPLETE",
            f"missing {len(missing)} results; first IDs: {missing[:10]}",
        )
    return trees, runtimes


def normalized_rf(
    first: set[frozenset[str]], second: set[frozenset[str]]
) -> tuple[int, int, float]:
    shared = len(first & second)
    rf = len(first ^ second)
    denominator = len(first) + len(second)
    return shared, rf, rf / denominator


def strategy_topology_metrics(
    tree: Phylo.BaseTree.Tree,
    species: list[str],
    groups: dict[str, set[str]],
) -> dict[str, float | bool]:
    splits = unrooted_splits(tree, set(species))
    side_masks = split_side_masks(splits, species)
    index = {name: position for position, name in enumerate(species)}
    metrics: dict[str, float | bool] = {}
    for strategy, group in groups.items():
        indices = np.asarray([index[name] for name in group], dtype=int)
        side_size = smallest_containing_side(indices, side_masks)
        metrics[f"{strategy}_smallest_containing_side_tips"] = side_size
        metrics[f"{strategy}_edge_delimited"] = side_size == len(group)
    return metrics


def bootstrap_products(
    trees: dict[int, Phylo.BaseTree.Tree],
    runtimes: dict[int, float],
    draws: pd.DataFrame,
    reference_tree: Phylo.BaseTree.Tree,
) -> dict[str, object]:
    species = comparison.canonical_species()
    observed_codes, groups = load_strategy_codes(species)
    del observed_codes
    all_tips = set(species)
    reference_splits = unrooted_splits(reference_tree, all_tips)
    calibration_splits = unrooted_splits(trees[0], all_tips)
    reference_strategy = strategy_topology_metrics(reference_tree, species, groups)
    calibration_strategy = strategy_topology_metrics(trees[0], species, groups)

    replicate_rows: list[dict[str, object]] = []
    bootstrap_split_sets: list[set[frozenset[str]]] = []
    for replicate_id in sorted(trees):
        if replicate_id == 0:
            continue
        splits = unrooted_splits(trees[replicate_id], all_tips)
        bootstrap_split_sets.append(splits)
        shared_full, rf_full, normalized_full = normalized_rf(reference_splits, splits)
        shared_calibration, rf_calibration, normalized_calibration = normalized_rf(
            calibration_splits, splits
        )
        strategy_metrics = strategy_topology_metrics(trees[replicate_id], species, groups)
        draw = draws.loc[draws["replicate_id"].eq(replicate_id)].iloc[0]
        replicate_rows.append(
            {
                "replicate_id": replicate_id,
                "n_unique_features": int(draw["n_unique_features"]),
                "reconstruction_seconds": runtimes[replicate_id],
                "shared_splits_with_full_fit": shared_full,
                "robinson_foulds_to_full_fit": rf_full,
                "normalized_robinson_foulds_to_full_fit": normalized_full,
                "shared_splits_with_matched_calibration": shared_calibration,
                "robinson_foulds_to_matched_calibration": rf_calibration,
                "normalized_robinson_foulds_to_matched_calibration": normalized_calibration,
                **strategy_metrics,
            }
        )
    replicate_frame = pd.DataFrame(replicate_rows)

    split_rows: list[dict[str, object]] = []
    sorted_splits = sorted(reference_splits, key=lambda value: (len(value), tuple(sorted(value))))
    for split_number, split in enumerate(sorted_splits, start=1):
        complement = all_tips - set(split)
        focal = [
            strategy
            for strategy, group in groups.items()
            if group == set(split) or group == complement
        ]
        support_count = sum(split in bootstrap_splits for bootstrap_splits in bootstrap_split_sets)
        split_rows.append(
            {
                "split_id": f"split_{split_number:03d}",
                "smaller_side_n": len(split),
                "smaller_side_species": json.dumps(sorted(split)),
                "focal_strategy": focal[0] if focal else "",
                "present_in_matched_calibration": split in calibration_splits,
                "bootstrap_support_count": support_count,
                "bootstrap_replicates": len(bootstrap_split_sets),
                "bootstrap_support_fraction": support_count / len(bootstrap_split_sets),
            }
        )
    split_frame = pd.DataFrame(split_rows)

    strategy_rows: list[dict[str, object]] = []
    for strategy in STRATEGIES:
        size_column = f"{strategy}_smallest_containing_side_tips"
        edge_column = f"{strategy}_edge_delimited"
        sizes = replicate_frame[size_column].to_numpy(dtype=float)
        edge_delimited = replicate_frame[edge_column].to_numpy(dtype=bool)
        strategy_rows.append(
            {
                "strategy": strategy,
                "strategy_label": STRATEGY_LABELS[strategy],
                "group_n": len(groups[strategy]),
                "full_fit_smallest_containing_side_tips": reference_strategy[
                    size_column
                ],
                "full_fit_edge_delimited": reference_strategy[edge_column],
                "matched_calibration_smallest_containing_side_tips": calibration_strategy[
                    size_column
                ],
                "matched_calibration_edge_delimited": calibration_strategy[edge_column],
                "bootstrap_edge_delimited_count": int(edge_delimited.sum()),
                "bootstrap_replicates": len(edge_delimited),
                "bootstrap_edge_delimited_fraction": float(edge_delimited.mean()),
                "bootstrap_smallest_side_median": float(np.median(sizes)),
                "bootstrap_smallest_side_q05": float(np.quantile(sizes, 0.05)),
                "bootstrap_smallest_side_q95": float(np.quantile(sizes, 0.95)),
            }
        )
    strategy_frame = pd.DataFrame(strategy_rows)

    full_cal_shared, full_cal_rf, full_cal_normalized = normalized_rf(
        reference_splits, calibration_splits
    )
    summary_rows = [
        {
            "metric_id": "matched_calibration_vs_full_fit_normalized_rf",
            "n": 1,
            "value": full_cal_normalized,
            "q05": full_cal_normalized,
            "median": full_cal_normalized,
            "q95": full_cal_normalized,
            "detail": f"shared={full_cal_shared}; RF={full_cal_rf}; splits={len(reference_splits)}+{len(calibration_splits)}",
        },
        {
            "metric_id": "bootstrap_normalized_rf_to_full_fit",
            "n": len(replicate_frame),
            "value": float(replicate_frame["normalized_robinson_foulds_to_full_fit"].mean()),
            "q05": float(replicate_frame["normalized_robinson_foulds_to_full_fit"].quantile(0.05)),
            "median": float(replicate_frame["normalized_robinson_foulds_to_full_fit"].median()),
            "q95": float(replicate_frame["normalized_robinson_foulds_to_full_fit"].quantile(0.95)),
            "detail": "value is mean; quantiles are across feature-bootstrap trees",
        },
        {
            "metric_id": "bootstrap_normalized_rf_to_matched_calibration",
            "n": len(replicate_frame),
            "value": float(
                replicate_frame["normalized_robinson_foulds_to_matched_calibration"].mean()
            ),
            "q05": float(
                replicate_frame["normalized_robinson_foulds_to_matched_calibration"].quantile(0.05)
            ),
            "median": float(
                replicate_frame["normalized_robinson_foulds_to_matched_calibration"].median()
            ),
            "q95": float(
                replicate_frame["normalized_robinson_foulds_to_matched_calibration"].quantile(0.95)
            ),
            "detail": "value is mean; matched bounded-search calibration uses all features once",
        },
        {
            "metric_id": "bootstrap_reference_split_support",
            "n": len(split_frame),
            "value": float(split_frame["bootstrap_support_fraction"].mean()),
            "q05": float(split_frame["bootstrap_support_fraction"].quantile(0.05)),
            "median": float(split_frame["bootstrap_support_fraction"].median()),
            "q95": float(split_frame["bootstrap_support_fraction"].quantile(0.95)),
            "detail": "value is mean support across the 31 full-fit internal splits",
        },
        {
            "metric_id": "bootstrap_reconstruction_seconds",
            "n": len(replicate_frame),
            "value": float(replicate_frame["reconstruction_seconds"].mean()),
            "q05": float(replicate_frame["reconstruction_seconds"].quantile(0.05)),
            "median": float(replicate_frame["reconstruction_seconds"].median()),
            "q95": float(replicate_frame["reconstruction_seconds"].quantile(0.95)),
            "detail": "Bonsai-reported reconstruction time per replicate",
        },
    ]
    return {
        "summary": pd.DataFrame(summary_rows),
        "replicates": replicate_frame,
        "strategies": strategy_frame,
        "splits": split_frame,
    }


def support_tree_text(
    reference_tree: Phylo.BaseTree.Tree, split_support: pd.DataFrame
) -> str:
    tree = copy.deepcopy(reference_tree)
    all_tips = {tip.name for tip in tree.get_terminals()}
    support_by_split = {
        frozenset(json.loads(row["smaller_side_species"])): float(row["bootstrap_support_fraction"])
        for row in split_support.to_dict(orient="records")
    }
    seen: set[frozenset[str]] = set()
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or clade is tree.root:
            continue
        clade.confidence = None
        clade.name = None
        side = {tip.name for tip in clade.get_terminals()}
        if not 1 < len(side) < len(all_tips) - 1:
            continue
        split = canonical_split(side, all_tips)
        if split in seen:
            continue
        seen.add(split)
        clade.name = f"{100 * support_by_split[split]:.0f}"
    tree.root.name = None
    buffer = io.StringIO()
    Phylo.write(tree, buffer, "newick")
    return buffer.getvalue()


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, lineterminator="\n", float_format="%.17g")


def write_tree_archive(
    trees: dict[int, Phylo.BaseTree.Tree], runtimes: dict[int, float], output_dir: Path
) -> pd.DataFrame:
    index_rows: list[dict[str, object]] = []
    with gzip.open(output_dir / BOOTSTRAP_ARCHIVE_FILENAME, "wt") as handle:
        for line_number, replicate_id in enumerate(sorted(trees), start=1):
            Phylo.write(trees[replicate_id], handle, "newick")
            index_rows.append(
                {
                    "archive_line": line_number,
                    "replicate_id": replicate_id,
                    "is_full_feature_calibration": replicate_id == 0,
                    "reconstruction_seconds": runtimes[replicate_id],
                }
            )
    index_frame = pd.DataFrame(index_rows)
    write_frame(index_frame, output_dir / BOOTSTRAP_INDEX_FILENAME)
    return index_frame


def read_tree_archive(data_dir: Path) -> tuple[dict[int, Phylo.BaseTree.Tree], dict[int, float]]:
    index = pd.read_csv(data_dir / BOOTSTRAP_INDEX_FILENAME)
    with gzip.open(data_dir / BOOTSTRAP_ARCHIVE_FILENAME, "rt") as handle:
        parsed = list(Phylo.parse(handle, "newick"))
    if len(parsed) != len(index):
        comparison.fail("BONSAI_BOOTSTRAP_ARCHIVE_INVALID", "tree archive and index differ")
    trees = {
        int(row.replicate_id): parsed[position]
        for position, row in enumerate(index.itertuples(index=False))
    }
    runtimes = {
        int(row.replicate_id): float(row.reconstruction_seconds)
        for row in index.itertuples(index=False)
    }
    return trees, runtimes


def validation_manifest(
    label_products: dict[str, pd.DataFrame],
    bootstrap: dict[str, object],
    data_dir: Path,
) -> dict[str, object]:
    tests = label_products["tests"].set_index("test_id")
    strategies = bootstrap["strategies"].set_index("strategy")
    summary = bootstrap["summary"].set_index("metric_id")
    input_paths = {
        "reference_bonsai_tree": REFERENCE_TREE,
        "superfamily_clr": comparison.CLR,
        "organismal_traits": comparison.TRAITS,
        "bootstrap_plan": data_dir / BOOTSTRAP_PLAN_FILENAME,
        "bootstrap_draws": data_dir / BOOTSTRAP_DRAWS_FILENAME,
        "bootstrap_tree_archive": data_dir / BOOTSTRAP_ARCHIVE_FILENAME,
        "bootstrap_tree_index": data_dir / BOOTSTRAP_INDEX_FILENAME,
        "validation_script": Path(__file__).resolve(),
        "bootstrap_runner": HERE / "run_feature_bootstrap.sh",
        "validation_plot_script": HERE / "plot_bonsai_validation.R",
        "python_project": ROOT / "pyproject.toml",
        "python_lock": ROOT / "uv.lock",
    }
    inputs = {
        label: {
            "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
            "sha256": sha256(path),
        }
        for label, path in input_paths.items()
    }
    return {
        "analysis_id": "te34_bonsai_reproductive_strategy_validation_v1",
        "analysis_date": "2026-08-25",
        "analysis_status": "exploratory_post_hoc_label_localization_and_bounded_feature_bootstrap",
        "inputs": inputs,
        "label_permutation": {
            "assigned_tips": 34,
            "unassigned_tips": [],
            "group_sizes": {
                strategy: int(tests.loc[f"{strategy}_smallest_containing_side_tips", "group_n"])
                for strategy in STRATEGIES
            },
            "permutations": PERMUTATIONS,
            "seed": PERMUTATION_SEED,
            "primary_endpoint": "global_fitch_parsimony_steps",
            "primary_observed": float(tests.loc["global_fitch_parsimony_steps", "observed"]),
            "primary_p_value": float(tests.loc["global_fitch_parsimony_steps", "p_value"]),
            "holm_adjustment": "all nine reported label-localization tests",
        },
        "feature_bootstrap": {
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "sampling_unit": "24 CLR-superfamily coordinates sampled with replacement",
            "bounded_search": {
                "nnn_n_randommoves": BOOTSTRAP_RANDOM_MOVES,
                "nnn_n_randomtrees": BOOTSTRAP_RANDOM_TREES,
                "full_fit_nnn_n_randommoves": 1000,
                "full_fit_nnn_n_randomtrees": 10,
            },
            "matched_calibration_normalized_rf_to_full_fit": float(
                summary.loc["matched_calibration_vs_full_fit_normalized_rf", "value"]
            ),
            "strategy_edge_delimited_support": {
                strategy: float(strategies.loc[strategy, "bootstrap_edge_delimited_fraction"])
                for strategy in STRATEGIES
            },
        },
        "interpretation_boundaries": [
            "Reproductive-strategy labels were interpreted after review of the initial PCA, so all tests are post hoc and exploratory.",
            "Permuting labels on the fixed TE-composition tree tests localization, not causation or phylogenetically independent association.",
            "The two focal exact groups are edge-delimited TE-composition clusters, not inferred species clades.",
            "Feature coordinates are compositional and correlated; resampling the 24 CLR coordinates is a perturbation analysis, not a biological-specimen bootstrap.",
            "Bootstrap runs use one random tree and 100 random NNI moves versus ten trees and 1,000 moves in the full fit.",
            "The matched full-feature calibration differs from the full fit, so support combines feature-resampling and bounded-search sensitivity.",
            "No per-species feature-error matrix is available; canonical PCA outputs remain as a supplementary audit while the validated Bonsai representation supplies the publication figure.",
        ],
    }


def write_validation(
    label_products: dict[str, pd.DataFrame],
    bootstrap: dict[str, object],
    trees: dict[int, Phylo.BaseTree.Tree],
    runtimes: dict[int, float],
    run_dir: Path,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_frame(label_products["tests"], output_dir / PERMUTATION_FILENAME)
    write_frame(
        label_products["fitch_histogram"], output_dir / PERMUTATION_HISTOGRAM_FILENAME
    )
    write_frame(bootstrap["summary"], output_dir / BOOTSTRAP_SUMMARY_FILENAME)
    write_frame(bootstrap["replicates"], output_dir / BOOTSTRAP_REPLICATE_FILENAME)
    write_frame(bootstrap["strategies"], output_dir / BOOTSTRAP_STRATEGY_FILENAME)
    write_frame(bootstrap["splits"], output_dir / BOOTSTRAP_SPLIT_FILENAME)
    species = comparison.canonical_species()
    codes, _ = load_strategy_codes(species)
    write_frame(
        tip_annotation_frame(species, codes), output_dir / TIP_ANNOTATION_FILENAME
    )
    write_frame(
        publication_summary_frame(label_products, bootstrap),
        output_dir / PUBLICATION_SUMMARY_FILENAME,
    )
    shutil.copyfile(run_dir / BOOTSTRAP_DRAWS_FILENAME, output_dir / BOOTSTRAP_DRAWS_FILENAME)
    shutil.copyfile(run_dir / BOOTSTRAP_PLAN_FILENAME, output_dir / BOOTSTRAP_PLAN_FILENAME)
    write_tree_archive(trees, runtimes, output_dir)
    reference_tree = Phylo.read(REFERENCE_TREE, "newick")
    (output_dir / SUPPORT_TREE_FILENAME).write_text(
        support_tree_text(reference_tree, bootstrap["splits"])
    )
    output_manifest = validation_manifest(label_products, bootstrap, output_dir)
    (output_dir / VALIDATION_MANIFEST_FILENAME).write_text(
        json.dumps(output_manifest, indent=2) + "\n"
    )


def check_frame(observed: pd.DataFrame, expected_path: Path) -> None:
    expected = pd.read_csv(
        expected_path, float_precision="round_trip", keep_default_na=False
    )
    pd.testing.assert_frame_equal(
        observed,
        expected,
        check_dtype=False,
        check_exact=False,
        rtol=2e-13,
        atol=2e-15,
    )


def check_tracked() -> None:
    species = comparison.canonical_species()
    reference_tree = Phylo.read(REFERENCE_TREE, "newick")
    validate_tree(reference_tree, species)
    codes, _ = load_strategy_codes(species)
    label_products = run_label_permutations(reference_tree, species, codes)
    trees, runtimes = read_tree_archive(TRACKED_DATA)
    draws = pd.read_csv(TRACKED_DATA / BOOTSTRAP_DRAWS_FILENAME)
    bootstrap = bootstrap_products(trees, runtimes, draws, reference_tree)
    for frame, filename in (
        (label_products["tests"], PERMUTATION_FILENAME),
        (label_products["fitch_histogram"], PERMUTATION_HISTOGRAM_FILENAME),
        (bootstrap["summary"], BOOTSTRAP_SUMMARY_FILENAME),
        (bootstrap["replicates"], BOOTSTRAP_REPLICATE_FILENAME),
        (bootstrap["strategies"], BOOTSTRAP_STRATEGY_FILENAME),
        (bootstrap["splits"], BOOTSTRAP_SPLIT_FILENAME),
        (tip_annotation_frame(species, codes), TIP_ANNOTATION_FILENAME),
        (
            publication_summary_frame(label_products, bootstrap),
            PUBLICATION_SUMMARY_FILENAME,
        ),
    ):
        check_frame(frame, TRACKED_DATA / filename)
    observed_support_tree = support_tree_text(reference_tree, bootstrap["splits"])
    if observed_support_tree != (TRACKED_DATA / SUPPORT_TREE_FILENAME).read_text():
        comparison.fail("BONSAI_BOOTSTRAP_SUPPORT_TREE_MISMATCH", "support tree is stale")
    observed_manifest = validation_manifest(label_products, bootstrap, TRACKED_DATA)
    expected_manifest = json.loads((TRACKED_DATA / VALIDATION_MANIFEST_FILENAME).read_text())
    if observed_manifest != expected_manifest:
        comparison.fail("BONSAI_VALIDATION_MANIFEST_MISMATCH", "validation manifest is stale")
    print_label_summary(label_products["tests"])
    print_bootstrap_summary(bootstrap["strategies"])
    print("PASS: tracked label permutations and feature bootstrap reproduce exactly")


def print_label_summary(tests: pd.DataFrame) -> None:
    columns = ["test_id", "observed", "null_mean", "p_value", "p_value_holm_all_nine_tests"]
    print(tests.loc[:, columns].to_string(index=False, float_format=lambda value: f"{value:.5g}"))


def print_bootstrap_summary(strategies: pd.DataFrame) -> None:
    columns = [
        "strategy_label",
        "group_n",
        "full_fit_edge_delimited",
        "matched_calibration_edge_delimited",
        "bootstrap_edge_delimited_fraction",
        "bootstrap_smallest_side_median",
    ]
    print(strategies.loc[:, columns].to_string(index=False, float_format=lambda value: f"{value:.4f}"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    permutation = subparsers.add_parser("permutation", help="run fixed-tree label permutations")
    permutation.add_argument("--permutations", type=int, default=PERMUTATIONS)
    permutation.add_argument("--seed", type=int, default=PERMUTATION_SEED)

    prepare = subparsers.add_parser("prepare-bootstrap", help="prepare resumable feature draws")
    prepare.add_argument("--run-dir", type=Path, required=True)
    prepare.add_argument("--replicates", type=int, default=BOOTSTRAP_REPLICATES)
    prepare.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)

    analyze = subparsers.add_parser("analyze-bootstrap", help="summarize completed upstream runs")
    analyze.add_argument("--run-dir", type=Path, required=True)
    analyze.add_argument("--output-dir", type=Path, required=True)
    analyze.add_argument("--write", action="store_true")

    subparsers.add_parser("check", help="replay all tracked validation artifacts")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare-bootstrap":
        prepare_feature_bootstrap(args.run_dir, args.replicates, args.seed)
        return
    if args.command == "check":
        check_tracked()
        return

    species = comparison.canonical_species()
    reference_tree = Phylo.read(REFERENCE_TREE, "newick")
    validate_tree(reference_tree, species)
    codes, _ = load_strategy_codes(species)
    label_products = run_label_permutations(
        reference_tree,
        species,
        codes,
        permutations=args.permutations if args.command == "permutation" else PERMUTATIONS,
        seed=args.seed if args.command == "permutation" else PERMUTATION_SEED,
    )
    print_label_summary(label_products["tests"])
    if args.command == "permutation":
        print("PASS: fixed-tree reproductive-strategy label permutation completed")
        return

    plan = json.loads((args.run_dir / BOOTSTRAP_PLAN_FILENAME).read_text())
    if plan != bootstrap_plan(BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED):
        comparison.fail("BONSAI_BOOTSTRAP_PLAN_MISMATCH", "run plan is not the declared campaign")
    trees, runtimes = collect_run_trees(args.run_dir, BOOTSTRAP_REPLICATES)
    draws = pd.read_csv(args.run_dir / BOOTSTRAP_DRAWS_FILENAME)
    bootstrap = bootstrap_products(trees, runtimes, draws, reference_tree)
    print_bootstrap_summary(bootstrap["strategies"])
    if args.write:
        write_validation(
            label_products, bootstrap, trees, runtimes, args.run_dir, args.output_dir
        )
        print(f"PASS: wrote strategy validation under {args.output_dir}")
    else:
        print("PASS: bootstrap analysis completed read-only; pass --write to store outputs")


if __name__ == "__main__":
    main()

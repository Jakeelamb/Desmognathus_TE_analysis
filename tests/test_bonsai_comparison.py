from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "analyses"
    / "01_transposable_elements"
    / "bonsai_comparison"
    / "compare_bonsai_to_pca.py"
)
SPEC = importlib.util.spec_from_file_location("bonsai_comparison", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
BONSAI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BONSAI)


def test_bonsai_uses_the_exact_canonical_clr_and_pca_geometry() -> None:
    clr = BONSAI.load_clr_matrix()
    pca = BONSAI.load_pca_scores(clr.index.tolist())

    assert clr.shape == (34, 24)
    assert pca.shape == (34, 23)
    source = BONSAI.squared_euclidean(clr.to_numpy())
    full_pca = BONSAI.squared_euclidean(pca.to_numpy())
    assert np.max(np.abs(source - full_pca)) < 1e-12


def test_tracked_bonsai_comparison_replays_and_preserves_claim_boundaries() -> None:
    analysis = BONSAI.build_analysis(BONSAI.TRACKED_TREE, BONSAI.TRACKED_METADATA)
    BONSAI.check_tracked_analysis(analysis)

    summary = analysis["summary"].set_index("representation_id")
    assert set(summary.index) == {"pca_2d", "pca_6d", "pca_full_23d", "bonsai_tree"}
    assert summary.loc["pca_full_23d", "global_distance_pearson_r"] > 1 - 1e-14
    assert summary.loc["pca_full_23d", "normalized_stress_after_scale"] < 1e-13
    assert summary.loc["bonsai_tree", "global_distance_pearson_r"] < summary.loc[
        "pca_2d", "global_distance_pearson_r"
    ]
    assert summary.loc["bonsai_tree", "mean_five_neighbor_recall"] > summary.loc[
        "pca_2d", "mean_five_neighbor_recall"
    ]
    assert summary.loc["bonsai_tree", "normalized_stress_after_scale"] < summary.loc[
        "pca_2d", "normalized_stress_after_scale"
    ]
    assert summary.loc["pca_6d", "global_distance_pearson_r"] > summary.loc[
        "bonsai_tree", "global_distance_pearson_r"
    ]

    manifest = json.loads((BONSAI.TRACKED_DATA / BONSAI.MANIFEST_FILENAME).read_text())
    boundaries = " ".join(manifest["interpretation_boundaries"])
    assert (
        manifest["analysis_status"]
        == "exploratory_method_comparison_with_validated_publication_representation"
    )
    assert manifest["method_adaptation"]["measurement_standard_deviations"] is None
    assert manifest["method_adaptation"]["rescale_by_var"] is False
    assert "not an inferred species phylogeny" in boundaries
    assert "feature-resampling analysis" in boundaries

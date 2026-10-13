from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "Publication/datasets"
FIGURE_DIR = ROOT / "Publication/figures"
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
NUMBERED_MANIFEST = FIGURE_DIR / "FIGURE_MANIFEST.csv"
AUXILIARY_MANIFEST = FIGURE_DIR / "AUXILIARY_FIGURE_MANIFEST.csv"
BASENAME = "Path24_development_and_trait_phylogeny_reference"
AUXILIARY_ID = "Auxiliary_Path24_development_traits_reference"


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _path24_development_traits() -> pd.DataFrame:
    path24 = pd.read_csv(DATASETS / "Supplementary_Data_S04_path_analysis_traits24.csv")[
        [
            "species",
            "tree_tip",
            "genome_size_pg_fuscus_anchored",
            "genome_size_pg_ci_low",
            "genome_size_pg_ci_high",
        ]
    ]
    morphology = pd.read_csv(
        DATASETS / "Supplementary_Data_S19_cell_nucleus_species_estimates.csv"
    )[
        [
            "species",
            "cell_area_um2",
            "cell_area_um2_ci_low",
            "cell_area_um2_ci_high",
            "nucleus_area_um2",
            "nucleus_area_um2_ci_low",
            "nucleus_area_um2_ci_high",
        ]
    ]
    development = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/organismal_traits.csv",
        keep_default_na=False,
    )[["species", "development_mode"]].rename(columns={"species": "tree_tip"})

    joined = (
        path24.merge(
            morphology,
            on="species",
            how="left",
            validate="one_to_one",
            indicator="morphology_join",
        )
        .merge(
            development,
            on="tree_tip",
            how="left",
            validate="one_to_one",
            indicator="development_join",
        )
    )
    return joined


def test_path24_trait_reference_uses_the_exact_complete_join_without_imputation() -> None:
    joined = _path24_development_traits()
    tree_nodes = pd.read_csv(
        DATASETS / "Supplementary_Data_S39_phylogeny_nodes.csv",
        keep_default_na=False,
    )
    tree_tips = set(tree_nodes.loc[tree_nodes["is_tip"], "tip_label"])
    organismal_traits = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/organismal_traits.csv",
        keep_default_na=False,
    )

    assert (len(joined), joined["species"].nunique(), joined["tree_tip"].nunique()) == (
        24,
        24,
        24,
    )
    for indicator in ("morphology_join", "development_join"):
        assert joined[indicator].eq("both").all()
    assert set(joined["tree_tip"]) == tree_tips

    direct = joined.loc[joined["development_mode"].eq("direct_development")]
    assert set(direct["tree_tip"]) == {"aeneus", "wrighti"}
    assert len(direct) == 2
    assert joined["development_mode"].eq("aquatic_larva").sum() == 22
    all_direct = organismal_traits.loc[
        organismal_traits["development_mode"].eq("direct_development"), "species"
    ]
    assert set(all_direct) == {"aeneus", "organi", "wrighti"}
    assert "organi" not in set(joined["tree_tip"])
    assert "D. organi" not in set(joined["species"])


def test_path24_trait_reference_intervals_are_finite_and_bracket_estimates() -> None:
    joined = _path24_development_traits()
    intervals = {
        "genome_size_pg_fuscus_anchored": (
            "genome_size_pg_ci_low",
            "genome_size_pg_ci_high",
        ),
        "nucleus_area_um2": (
            "nucleus_area_um2_ci_low",
            "nucleus_area_um2_ci_high",
        ),
        "cell_area_um2": ("cell_area_um2_ci_low", "cell_area_um2_ci_high"),
    }

    for estimate, (low, high) in intervals.items():
        values = joined[[estimate, low, high]]
        assert values.notna().all().all(), estimate
        assert (values > 0).all().all(), estimate
        assert values[low].le(values[estimate]).all(), estimate
        assert values[estimate].le(values[high]).all(), estimate


def test_path24_trait_reference_is_auxiliary_and_not_a_numbered_figure() -> None:
    source = FIGURE_SOURCE.read_text(encoding="utf-8")
    block = source.split("## ---- path24-development-traits-reference", maxsplit=1)[1].split(
        "## ---- path-dag-reference", maxsplit=1
    )[0]

    assert "record_auxiliary_figure(" in block
    assert "record_figure(" not in block
    assert BASENAME in block
    assert "D. organi is absent" in block
    assert "not imputed" in block
    assert "genome_size_path24 <- path_traits" in block
    assert "Supplementary_Data_S22_relative_nuclear_iod_species.csv" not in block

    numbered_rows = _read_manifest(NUMBERED_MANIFEST)
    assert all(row["basename"] != BASENAME for row in numbered_rows)

    auxiliary_rows = _read_manifest(AUXILIARY_MANIFEST)
    matching = [row for row in auxiliary_rows if row["basename"] == BASENAME]
    assert len(matching) == 1
    assert matching[0]["auxiliary_id"] == AUXILIARY_ID


def test_path24_trait_reference_assets_exist_in_all_release_formats() -> None:
    matching = [
        row for row in _read_manifest(AUXILIARY_MANIFEST) if row["auxiliary_id"] == AUXILIARY_ID
    ]
    assert len(matching) == 1
    row = matching[0]

    for column in ("pdf", "png", "tiff"):
        path = FIGURE_DIR / row[column]
        assert path.is_file(), path
        assert path.stat().st_size > 10_000, path


def test_path_dag_node_is_genome_size_not_relative_iod() -> None:
    source = FIGURE_SOURCE.read_text(encoding="utf-8")
    block = source.split("## ---- path-dag-reference", maxsplit=1)[1].split(
        "## ---- figure-manifest", maxsplit=1
    )[0]

    assert 'node_label = c("Genome\\nsize", "Nucleus\\narea", "Cell\\narea")' in block
    assert '"Relative\\nIOD"' not in block
    assert "Genome size is the D. fuscus assembly-anchored pg/1C estimate" in block
    assert "not a direct species C-value" in block

from __future__ import annotations

from pathlib import Path
from runpy import run_path
from xml.etree import ElementTree
from zipfile import ZipFile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
METHODS = ROOT / "Publication/methods"


def _docx_text(path: Path) -> str:
    with ZipFile(path) as archive:
        document_xml = archive.read("word/document.xml")
    root = ElementTree.fromstring(document_xml)
    word_text = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t"
    return " ".join((node.text or "") for node in root.iter(word_text))


def test_methods_evidence_matrix_uses_only_active_sources() -> None:
    evidence = pd.read_csv(METHODS / "METHODS_EVIDENCE_MATRIX.csv")
    assert evidence["section_id"].tolist() == [f"M{number:02d}" for number in range(1, 21)]

    forbidden = ("input_data/", "path_analysis/", "results/", "plans/")
    for source_list in evidence["primary_sources"].astype(str):
        for source in (item.strip() for item in source_list.split(";")):
            assert not source.startswith(forbidden), source
            if source == "none" or source.startswith("doi:"):
                continue
            assert (ROOT / source).is_file(), source


def test_methods_scaffold_records_finalized_release_boundaries() -> None:
    scaffold = (METHODS / "METHODS_BULLET_SCAFFOLD.md").read_text(encoding="utf-8")
    required_phrases = {
        "one 68-row, five-column table with 34 order-level and 34 superfamily-level rows",
        "order-level stratum is the primary paper-facing diversity analysis",
        "superfamily-level stratum is a sensitivity/audit view",
        "classified categories were reclosed to one",
        "Natural-log Shannon entropy was `H′ = −Σpᵢ ln(pᵢ)`",
        "explicitly named Gini-Simpson index was `1 − Σpᵢ²`",
        "These are the only paper-facing diversity indices",
        "Observed richness is retained as an audit/support count",
        "Unresolved aligned-base mass remains only in S06-S07 dnaPipeTE quality-control tables",
        "Rare classified features excluded from the shared-feature PCA remain in both order- and superfamily-level diversity calculations",
        "left LTR, right LTR, and internal region",
        "inclusive, within-species, two-sided Tukey filter",
        "380 primary elements",
        "whitespace-delimited `DOMAIN|MODEL` annotations",
        "0.12 µm/pixel",
        "Seven species are represented by one animal",
        "`nuclear_iod_by_image.csv`",
        "fail this gate and remain visibly flagged",
    }
    assert all(phrase in scaffold for phrase in required_phrases)
    assert "Hill q1" not in scaffold
    assert "Hill q2" not in scaffold
    assert "Pielou" not in scaffold
    assert "mass-aware" not in scaffold
    assert "mass_aware_unresolved_bin" not in scaffold
    assert "136 rows" not in scaffold
    assert "four exact 34-species strata" not in scaffold


def test_s08_metadata_declares_one_classified_only_diversity_table() -> None:
    metadata = run_path(ROOT / "scripts/publication_metadata.py")
    dependencies = metadata["SUPPLEMENT_DEPENDENCIES"]
    column_definitions = metadata["SUPPLEMENT_COLUMN_DEFINITIONS"]
    supplement_metadata = metadata["SUPPLEMENT_METADATA"]

    assert dependencies["S08"] == (
        "analyses/01_transposable_elements/data/te_order_composition.csv",
        "analyses/01_transposable_elements/data/te_superfamily_composition.csv",
    )
    assert "composition_mode" not in column_definitions["S08"]

    description = supplement_metadata["S08"]["description"]
    manuscript_use = supplement_metadata["S08"]["manuscript_use"]
    assert "68-row, five-column classified-only diversity table" in description
    assert "34 order-level and 34 superfamily-level" in description
    assert "Unresolved aligned-base mass is excluded from diversity" in manuscript_use
    assert "S06-S07" in manuscript_use


def test_conditional_assembly_anchor_metadata_is_traceable_and_limited() -> None:
    metadata = run_path(ROOT / "scripts/publication_metadata.py")
    dependencies = metadata["SUPPLEMENT_DEPENDENCIES"]
    column_definitions = metadata["SUPPLEMENT_COLUMN_DEFINITIONS"]
    supplement_metadata = metadata["SUPPLEMENT_METADATA"]

    anchor_dependencies = {
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
        "analyses/02_morphology/validate_release.py",
    }
    assert anchor_dependencies <= set(dependencies["S01"])
    assert anchor_dependencies - {
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
    } <= set(dependencies["S22"])

    for supplement_id in ["S01", "S22"]:
        description = supplement_metadata[supplement_id]["description"]
        manuscript_use = supplement_metadata[supplement_id]["manuscript_use"]
        assert "conditional" in description.casefold()
        assert "16.1-Gbp" in description
        assert "0.978 Gbp per pg" in description
        assert "not a measured C-value" in manuscript_use

    s01_definitions = column_definitions["S01"]
    assert "Exact S22 projection" in s01_definitions["step06_genome_size_pg_fuscus_anchored"]
    assert "not a measured C-value" in s01_definitions["step06_genome_size_pg_fuscus_anchored"]

    s22_definitions = column_definitions["S22"]
    expected_calibration_columns = {
        "iod_ratio_to_fuscus",
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
        "genome_size_reference_species",
        "genome_size_reference_pg",
        "genome_size_calibration_pg_per_iod",
        "genome_size_rank",
    }
    assert expected_calibration_columns <= set(s22_definitions)
    assert "16.1 Gbp / 0.978 Gbp per pg" in s22_definitions["genome_size_reference_pg"]
    assert "anchor-value uncertainty is excluded" in s22_definitions["genome_size_pg_ci_low"]
    assert "not a measured C-value" in s22_definitions["genome_size_pg_fuscus_anchored"]


def test_word_methods_exports_do_not_restore_removed_diversity_mode() -> None:
    word_exports = (
        METHODS / "MATERIALS_AND_METHODS_NARRATIVE_DRAFT.docx",
        METHODS / "METHODS_BULLET_SCAFFOLD_WORD_READY.docx",
    )
    forbidden = (
        "mass-aware",
        "mass aware",
        "mass_aware_unresolved_bin",
        "136 rows",
        "four exact 34-species strata",
        "values that included unresolved repeat mass",
    )
    for path in word_exports:
        assert path.is_file()
        text = _docx_text(path).casefold()
        assert not any(phrase.casefold() in text for phrase in forbidden), path

    narrative = _docx_text(word_exports[0])
    assert "Unresolved aligned-base mass was excluded from all diversity values" in narrative
    assert "34 superfamily-level values" in narrative
    assert "Process_413/specimen 32469" in narrative
    assert "16.462167689 pg/1C" in narrative
    assert "relative IOD was retained as the primary phenotype" in narrative

    scaffold = _docx_text(word_exports[1])
    assert "12 August 2026" in scaffold
    assert "one 68-row, five-column table" in scaffold
    assert "does not enter S08 diversity" in scaffold
    assert "Process_414/specimen 32470" in scaffold
    assert "16.462167689 pg/1C" in scaffold

    for path in word_exports:
        assert "16.36-pg" not in _docx_text(path)


def test_unresolved_methods_are_recorded_without_recreating_old_paths() -> None:
    required = pd.read_csv(METHODS / "METHODS_REQUIRED_AUTHOR_INPUT.csv").set_index("item_id")
    assert required.loc["W12", "status"] == "partial"
    assert required.loc["W16", "status"] == "partial"
    assert "multimapper policy" in required.loc["W12", "required_detail"]
    assert "checkpoint-to-mask-release mapping" in required.loc["W16", "required_detail"]

    old_image_table = ROOT / "analyses/02_morphology/data/relative_nuclear_iod_by_image.csv"
    new_image_table = ROOT / "analyses/02_morphology/data/nuclear_iod_by_image.csv"
    assert not old_image_table.exists()
    assert new_image_table.is_file()

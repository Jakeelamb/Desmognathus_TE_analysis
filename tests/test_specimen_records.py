from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/publication/reconcile_specimen_records.py"
SPEC = importlib.util.spec_from_file_location("specimen_records", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
SOURCE = REPO / MODULE.RECORD_DIR


def source_records():
    tables = MODULE.parse_tables(SOURCE / "AUTHOR_SUPPLIED_RECORDS.txt")
    return tables, tables["MVZ"] + tables["RAP"]


def inventory(specimen="3372", name="Process_314_raw_green.ome.tiff", species="D. fuscus"):
    return [{"assay": "nuclear_iod", "specimen_id": specimen, "filename": name, "species": species}]


def test_original_dates_names_and_comments_survive():
    tables, _ = source_records()
    assert len(tables["MVZ"]) == 97
    assert len(tables["RAP"]) == 47
    mvz = next(r for r in tables["MVZ"] if r["numeric_identifier"] == "109787")
    assert mvz["Date"] == "1/12/23"
    assert mvz["Verbatim Date"] == "1 Jul 1972"
    assert mvz["Species"] == "D. ochrophaeus"
    assert mvz["Suggested Species"] == "D. orestes"
    rap = next(r for r in tables["RAP"] if r["source_identifier"] == "RAP3518")
    assert rap["Source comment"] == "Contaminated???"
    assert json.loads(tables["RAP"][-1]["raw_fields_json"])[-1] == "Whole + Tissue"
    assert len(json.loads(tables["RAP"][-1]["raw_fields_json"])) == 16


def test_numeric_only_identity_remains_candidate():
    _, records = source_records()
    row = MODULE.build_associations(inventory(), records, [])[0]
    assert row["source_identifier"] == "RAP3372"
    assert row["association_status"] == "unique_numeric_candidate"
    assert row["namespace_evidence"] == "namespace_not_recorded"
    assert row["permit_verbatim"] == "WVDNR 2023.001"


def test_explicit_namespace_never_crosses_registry():
    tables, _ = source_records()
    mvz = {**tables["MVZ"][0], "numeric_identifier": "3372", "source_identifier": "MVZ:3372"}
    rap = next(r for r in tables["RAP"] if r["source_identifier"] == "RAP3372")
    rows = MODULE.build_associations(inventory("MVZ:3372"), [mvz, rap], [])
    assert len(rows) == 1
    assert rows[0]["source_namespace"] == "MVZ"
    assert rows[0]["association_status"] == "confirmed_namespace"
    ambiguous = MODULE.build_associations(inventory(), [mvz, rap], [])
    assert len(ambiguous) == 2
    assert {r["association_status"] for r in ambiguous} == {"ambiguous_source_rows"}


def test_duplicate_source_rows_are_preserved():
    _, records = source_records()
    rap = next(r for r in records if r["source_identifier"] == "RAP3372")
    rows = MODULE.build_associations(inventory(), [rap, {**rap, "source_row": "999"}], [])
    assert len(rows) == 2
    assert {r["source_row"] for r in rows} == {rap["source_row"], "999"}
    assert {r["association_status"] for r in rows} == {"ambiguous_source_rows"}


def test_suggested_taxon_never_overrides_analysis():
    _, records = source_records()
    row = MODULE.build_associations(
        inventory("109787", "Process_409_green.ome.tiff", "D. orestes"), records, []
    )[0]
    assert row["analysis_species"] == "D. orestes"
    assert row["source_species"] == "D. ochrophaeus"
    assert row["suggested_species"] == "D. orestes"
    assert row["species_exact_agreement"] == "false"
    assert row["collection_date_verbatim"] == "1 Jul 1972"
    assert row["ambiguous_spreadsheet_date"] == "1/12/23"


def test_missing_3463_is_not_replaced_by_3462():
    _, records = source_records()
    row = MODULE.build_associations(inventory("3463", species="D. auriculatus"), records, [])[0]
    assert row["association_status"] == "unmatched"
    assert row["source_identifier"] == ""
    assert row["analysis_specimen_id"] == "3463"


def test_changed_source_is_rejected(tmp_path):
    path = tmp_path / "source.txt"
    path.write_bytes((SOURCE / "AUTHOR_SUPPLIED_RECORDS.txt").read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        MODULE.parse_tables(path)


def test_real_cli_and_repeat_are_deterministic(tmp_path):
    destinations = [tmp_path / "first", tmp_path / "second"]
    for output in destinations:
        completed = subprocess.run(
            [sys.executable, str(SCRIPT), "--output-dir", str(output)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
    for path in destinations[0].iterdir():
        assert path.read_bytes() == (destinations[1] / path.name).read_bytes()
    receipt = json.loads((destinations[0] / "receipt.json").read_text())
    assert receipt["assay_counts"] == {
        "morphometry": {"objects": 1152, "animals": 42},
        "nuclear_iod": {"objects": 805, "animals": 50},
    }
    assert receipt["association_counts"]["reference_standard"] == {"confirmed_namespace": 2}
    repeated = subprocess.run(
        [sys.executable, str(SCRIPT), "--output-dir", str(destinations[0])],
        capture_output=True,
        text=True,
        check=False,
    )
    assert repeated.returncode == 2
    assert "already exists" in repeated.stderr
    for path in destinations[0].iterdir():
        assert path.read_bytes() == (destinations[1] / path.name).read_bytes()


def test_no_write_to_retained_science_directories(tmp_path):
    for name in ("Publication", "analyses", "data", "outputs"):
        path = REPO / name / "specimen_records_refusal_control"
        with pytest.raises(ValueError, match="retained scientific"):
            MODULE.run(REPO, SOURCE, path)
        assert not path.exists()


def test_changed_metadata_is_rejected_before_output(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for filename in (
        "AUTHOR_SUPPLIED_RECORDS.txt",
        "MASTER_IMAGE_METADATA.csv",
        "METADATA_SOURCE.json",
        "SOURCE.json",
    ):
        (source / filename).write_bytes((SOURCE / filename).read_bytes())
    metadata = source / "MASTER_IMAGE_METADATA.csv"
    metadata.write_bytes(metadata.read_bytes() + b"\n")
    output = tmp_path / "audit"
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--source-dir", str(source), "--output-dir", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 2
    assert "Metadata snapshot checksum" in completed.stderr
    assert not output.exists()

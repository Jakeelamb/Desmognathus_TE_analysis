import csv
import hashlib
import json
import re
import shutil
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


root = Path(__file__).resolve().parent
project = root.parent.parent
template = project / "outputs/submission_preparation_20261013/STAINING_RUN_MAP_TO_COMPLETE.csv"
images = project / "analyses/02_morphology/data/nuclear_iod_by_image.csv"
metadata = Path("/home/jake/Projects/cellprofiler_test/master_image_metadata.csv")
glasses = metadata.with_name("slide_glass_mapping.csv")
inputs = [template, images, metadata, glasses]
before = {str(path): digest(path) for path in inputs}
attestation = json.loads((root / "AUTHOR_STAINING_CONFIRMATION.json").read_text())

with template.open(newline="") as handle:
    reader = csv.DictReader(handle)
    fields = reader.fieldnames + [
        "confirmation_basis", "confirmation_date", "run_id_is_internal_label",
        "physical_slide_glass_id", "scan_date", "original_vsi_filename",
    ]
    rows = list(reader)
with images.open(newline="") as handle:
    retained = list(csv.DictReader(handle))
with metadata.open(newline="") as handle:
    by_filename = {row["filename"]: (line, row) for line, row in enumerate(csv.DictReader(handle), 2)}
with glasses.open(newline="") as handle:
    by_process = {row["process_id"]: (line, row) for line, row in enumerate(csv.DictReader(handle), 2)}

assert len(rows) == 53 and len(retained) == 51
target_rows = [row for row in rows if row["role"] in {"target", "current species denominator"}]
standard_rows = [row for row in rows if row["role"] == "designated standard; QC unresolved"]
assert len(target_rows) == 51 and len(standard_rows) == 2
assert {row["filename"] for row in target_rows} == {row["filename"] for row in retained}
assert {row["specimen_id"] for row in standard_rows} == {"32469", "32470"}
for row in rows:
    line, record = by_filename[row["filename"]]
    assert record["specimen_id"] == row["specimen_id"]
    process = re.fullmatch(r"Process_(\d+).*", row["filename"]).group(1)
    glass_line, glass = by_process.get(process, (None, {}))
    row.update({
        "staining_run_id": attestation["staining_run_id"],
        "staining_date": "",
        "standard_process_for_run": "Process_413;Process_414",
        "record_location": "AUTHOR_STAINING_CONFIRMATION.json",
        "acquisition_record_location": f"{metadata}:{line}" + (f";{glasses}:{glass_line}" if glass_line else ""),
        "notes": "Run membership confirmed by author. Run ID is an internal label. Staining date unknown. Standard QC and calibration remain unresolved.",
        "confirmation_basis": "author attestation",
        "confirmation_date": attestation["confirmation_date"],
        "run_id_is_internal_label": "true",
        "physical_slide_glass_id": glass.get("slide_glass_id", ""),
        "scan_date": glass.get("scan_date", ""),
        "original_vsi_filename": glass.get("vsi_filename", ""),
    })

mapping = root / "STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv"
with mapping.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

exports = [
    ("/tmp/browser-use/exports/Desmognathus dataset-52841d9f-57d0-46ab-8166-49d734618bbf.xlsx", "Desmognathus_dataset_source_snapshot.xlsx", "https://docs.google.com/spreadsheets/d/1CaMYntCP1eDSa1kX297Y36b6hhb9Jcps_1SQqnE2M6E/edit"),
    ("/tmp/browser-use/exports/Desmognathus GBE manuscript-6d134466-5ba1-47a5-a600-0addf1550fc7.md", "Desmognathus_GBE_manuscript_source_snapshot.md", "https://docs.google.com/document/d/1IKIZmB_hmDO3XcjcJsQnZ_HNOGnx8SBy6rlRPuDXxVc/edit"),
    ("/tmp/browser-use/exports/Cell Stain and image methodology-6dca3a4b-8e96-4336-8754-b6dae2ce2e8a.md", "Cell_Stain_and_image_methodology_source_snapshot.md", "https://docs.google.com/document/d/1uTvZyVQs7ZlLS_MY44zNUgSO0yTtLN1wliRnJw7lfI0/edit"),
    ("/tmp/browser-use/exports/Staining Procedure-014ec278-72a6-4bf3-b041-363a727c8cf4.md", "Staining_Procedure_source_snapshot.md", "https://docs.google.com/document/d/1JYC3_nhW3JNp66Th7hJNtBI5jL3G54xA4pPLNIQMd3Y/edit"),
]
snapshots = []
for source, name, url in exports:
    source = Path(source)
    destination = root / name
    shutil.copyfile(source, destination)
    assert digest(source) == digest(destination)
    snapshots.append({"path": name, "source_url": url, "sha256": digest(destination), "format": destination.suffix[1:]})
lab_work = root / "Salamander_slides_Lab_Work_source_snapshot.xlsx"
snapshots.append({"path": lab_work.name, "source_url": "https://docs.google.com/spreadsheets/d/1qmsaXn9pffPq68fphp9Wz5Tcvup43oas2AN8DKQh8QY/edit", "sha256": digest(lab_work), "format": "xlsx"})

assert before == {str(path): digest(path) for path in inputs}
receipt = {
    "verified_browser_account": "jakelamb74@gmail.com",
    "inspection_date": "2026-10-13",
    "throughput_checkpoint": "n/a, read-only investigation",
    "target_image_rows": 51,
    "designated_standard_rows": 2,
    "retained_nuclei": sum(int(row["n_nuclei"]) for row in retained),
    "images_with_glass_records": sum(bool(row["physical_slide_glass_id"]) for row in rows),
    "mapping_sha256": digest(mapping),
    "author_confirmation_sha256": digest(root / "AUTHOR_STAINING_CONFIRMATION.json"),
    "local_inputs_unchanged": before,
    "read_only_source_snapshots": snapshots,
    "export_limit": "Browser exports do not independently establish complete native comments, notes, formatting, or revision history.",
    "cloud_writes_performed": False,
    "calibration_validated": False,
}
(root / "ONLINE_SOURCE_PINS.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps({key: receipt[key] for key in ["target_image_rows", "designated_standard_rows", "retained_nuclei", "images_with_glass_records", "cloud_writes_performed", "calibration_validated"]}))

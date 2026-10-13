"""Inventory selected artifacts and verify preserved source bytes without refits."""

import csv
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
RUN = Path(__file__).resolve().parent
ART = ROOT / "outputs/manuscript_polished_20260924_01"
TABLES = ROOT / "outputs/manuscript_story_20260924_01/tables/integrated"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


baseline = json.loads((RUN / "baseline.json").read_text())
changed = [
    row["path"] for row in baseline["protected_figure_run"]
    if not (ROOT / row["path"]).is_file()
    or digest(ROOT / row["path"]) != row["sha256"]
]
assert not changed, changed
with ZipFile(ART / "Manuscript_figures_1-8.zip") as archive:
    assert archive.testzip() is None
    package = json.loads(archive.read("PACKAGE_MANIFEST.json"))
    for row in package["files"]:
        assert digest(ART / row["path"]) == row["sha256"]
        assert hashlib.sha256(archive.read(row["path"])).hexdigest() == row["sha256"]
    assert len(package["files"]) == 26

inventory = []


def add(category, item, path, disposition):
    inventory.append({
        "category": category, "item": item,
        "path": str(path.relative_to(ROOT)), "sha256": digest(path),
        "disposition": disposition,
    })


for number in range(1, 9):
    folder = ART / "figures" / f"figure_{number}"
    with (folder / "FIGURE_MANIFEST.csv").open() as handle:
        row = next(csv.DictReader(handle))
    for extension in ("pdf", "png", "tiff"):
        add("main_figure", f"Figure {number}: {row['title']}",
            folder / row[extension], "review candidate; measurements conditional")
for number, name in enumerate((
    "Table_1_data_coverage.csv", "Table_2_focal_associations.csv",
    "Table_3_path_comparisons.csv", "Table_4_species_measurements.csv",
    "Table_5_sampling_sensitivity.csv",
), 1):
    add("supporting_table", f"Proposed supplementary table S{number}",
        TABLES / name, "numbering and manuscript placement proposed")
add("supporting_workbook", "Five integrated tables", TABLES / "Manuscript_tables_integrated.xlsx",
    "source workbook; supplement PDF not assembled")
for extension in ("md", "docx"):
    add("methods", "Reconciled researcher Methods",
        ROOT / "Publication/methods" / f"MATERIALS_AND_METHODS_RESEARCHER_DRAFT.{extension}",
        "synchronized and reviewed; laboratory records remain open")
add("captions", "Matched findings, captions, and alt text", ART / "FINDINGS_AND_CAPTIONS.md",
    "current eight-figure source")
add("figure_package", "Eight collages and captions", ART / "Manuscript_figures_1-8.zip",
    "all 26 manifested entries verified against current masters")
with (RUN / "SUBMISSION_INVENTORY.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(inventory[0]))
    writer.writeheader()
    writer.writerows(inventory)

with (ROOT / "analyses/02_morphology/data/nuclear_iod_by_image.csv").open() as handle:
    image_rows = list(csv.DictReader(handle))
assert len(image_rows) == 51
fields = ["species", "filename", "specimen_id", "role", "staining_run_id",
          "staining_date", "standard_process_for_run", "record_location",
          "acquisition_record_location", "notes"]
with (RUN / "STAINING_RUN_MAP_TO_COMPLETE.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for row in image_rows:
        writer.writerow({
            **{key: row[key] for key in ("species", "filename", "specimen_id")},
            "role": "current species denominator" if row["species"] == "D. fuscus" else "target",
        })
    for process, specimen in ((413, 32469), (414, 32470)):
        writer.writerow({"species": "D. fuscus", "filename": f"Process_{process}_green.ome.tiff",
                         "specimen_id": specimen, "role": "designated standard; QC unresolved"})

receipt = {
    "recorded_at": datetime.now(UTC).isoformat(), "status": "PASS",
    "protected_figure_files_unchanged": len(baseline["protected_figure_run"]),
    "package_manifest_entries_verified": len(package["files"]),
    "selected_inventory_entries": len(inventory), "target_images_in_mapping_template": 51,
    "designated_standards_in_mapping_template": 2,
    "scientific_models_refitted": False,
    "qualification": "artifact integrity and local inventory; submission readiness remains open",
}
(RUN / "artifact_integrity.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt, indent=2))

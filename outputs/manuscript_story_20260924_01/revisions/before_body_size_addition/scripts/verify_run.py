"""Record preservation evidence and check this isolated publication run."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
ROOT = RUN.parents[1]
QA = RUN / "qa"

def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

def original_files():
    for directory in ("Publication", "analyses", "scripts/publication", "tests"):
        for path in sorted((ROOT / directory).rglob("*")):
            if path.is_file() and "__pycache__" not in path.parts:
                yield path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("baseline", "verify"))
    args = parser.parse_args()
    QA.mkdir(exist_ok=True)
    baseline = QA / "original_checksums.json"
    if args.mode == "baseline":
        if baseline.exists():
            raise RuntimeError("Preservation baseline already exists; refusing replacement")
        records = {str(path.relative_to(ROOT)): digest(path) for path in original_files()}
        baseline.write_text(json.dumps(records, indent=2) + "\n")
        print(f"Recorded {len(records)} original-file SHA-256 hashes")
        return
    records = json.loads(baseline.read_text())
    changed = [name for name, value in records.items() if not (ROOT / name).is_file() or digest(ROOT / name) != value]
    preserved_path = QA / "pre_path_addition_assets.json"
    preserved = json.loads(preserved_path.read_text()) if preserved_path.exists() else {}
    # The user-authorized label revision retains the earlier exports verbatim
    # in this archive while the current gallery uses the revised labels.
    label_archive = RUN / "revisions/before_genome_size_labels"
    changed_existing_assets = []
    archived_existing_assets = []
    for name, value in preserved.items():
        if (RUN / name).is_file() and digest(RUN / name) == value:
            continue
        if (label_archive / name).is_file() and digest(label_archive / name) == value:
            archived_existing_assets.append(name)
        else:
            changed_existing_assets.append(name)
    manifest_paths = sorted((RUN / "figures").glob("figure_*/FIGURE_MANIFEST.csv"))
    assert len(manifest_paths) == 6, f"Expected 6 main figures, found {len(manifest_paths)}"
    assets = []
    figures = []
    panel_count = 0
    for path in manifest_paths:
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 1
        figures.append(rows[0]["figure_id"])
        with (path.parent / "panels/PANEL_MANIFEST.csv").open() as handle:
            panels = list(csv.DictReader(handle))
        panel_count += len(panels)
        for base, row in [(path.parent, rows[0])] + [(path.parent / "panels", p) for p in panels]:
            assert row["final_font_family"] == "Arial"
            for kind in ("pdf", "png", "tiff"):
                asset = base / row[kind]
                assert asset.is_file() and asset.stat().st_size > 1000, asset
                assets.append(str(asset.relative_to(RUN)))
    assert sorted(figures) == [f"Figure_{i}" for i in range(1, 7)]
    font_audit = []
    for relative in assets:
        if relative.endswith(".pdf"):
            output = subprocess.check_output(["pdffonts", str(RUN / relative)], text=True)
            font_lines = output.splitlines()[2:]
            assert font_lines, relative
            for line in font_lines:
                assert "Arial" in line, (relative, line)
                assert "yes" in line.split(), (relative, line)
            font_audit.append({"file": relative, "pdffonts": output})
    (QA / "font_audit.json").write_text(json.dumps(font_audit, indent=2) + "\n")
    result = {"original_files_checked": len(records), "changed_original_files": changed,
              "earlier_candidate_assets_checked": len(preserved),
              "earlier_candidate_assets_preserved_in_label_archive": len(archived_existing_assets),
              "changed_earlier_candidate_assets": changed_existing_assets,
              "main_figures": len(figures), "isolated_panels": panel_count,
              "figure_assets": len(assets), "all_figure_fonts_Arial": True}
    (QA / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    input_records = {str(p.relative_to(RUN)): digest(p) for p in sorted((RUN / "inputs").rglob("*")) if p.is_file() and "__pycache__" not in p.parts}
    (QA / "input_checksums.json").write_text(json.dumps(input_records, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    assert not changed, "Original files changed; inspect preservation report"
    assert not changed_existing_assets, "Earlier candidate asset is not preserved live or in the label-revision archive"

if __name__ == "__main__":
    main()

"""Preserve the supplied specimen tables and audit their released assay associations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

RAW_SHA256 = "93dae07e9ca1ac0566e642634e8aac3a6430bc7f2302bfd1b428ac88f3447963"
RECORD_DIR = Path("docs/submission_work_20261013/specimen_records")
ASSAYS = {
    "morphometry": ("cell_nucleus_objects.csv", 1152, 42),
    "nuclear_iod": ("nuclear_iod_objects.csv", 805, 50),
}
STANDARDS = {"413": "32469", "414": "32470"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_tables(raw: Path) -> dict[str, list[dict[str, str]]]:
    if sha256(raw) != RAW_SHA256:
        raise ValueError("Author-supplied source checksum differs from the preserved answer")
    text = raw.read_text(encoding="utf-8")
    first, rap = text.split(" Field Series\t", 1)
    _, mvz = first.split(" Date\t", 1)
    result = {}
    for namespace, body in {"MVZ": "Date\t" + mvz, "RAP": "Field Series\t" + rap}.items():
        table = list(csv.reader(body.splitlines(), delimiter="\t"))
        header = table[0]
        # The unnamed final RAP column contains the supplied contamination comment.
        names = [name.strip() or "Source comment" for name in header]
        records = []
        for index, fields in enumerate(table[1:], 1):
            if len(fields) > len(header) or len(fields) < len(header) - 1:
                raise ValueError(f"Unexpected column count in {namespace} row {index}")
            values = dict(zip(names, fields, strict=False))
            identifier = values["MVZ:ID" if namespace == "MVZ" else "Field Series"]
            if not re.fullmatch(r"\d+" if namespace == "MVZ" else r"RAP\d+", identifier):
                raise ValueError(f"Invalid {namespace} identity in row {index}")
            numeric = identifier if namespace == "MVZ" else identifier.removeprefix("RAP")
            records.append(
                {
                    "source_row": str(index),
                    "namespace": namespace,
                    "source_identifier": f"MVZ:{identifier}" if namespace == "MVZ" else identifier,
                    "numeric_identifier": numeric,
                    "raw_fields_json": json.dumps(fields, ensure_ascii=False),
                    **{name: values.get(name, "") for name in names},
                }
            )
        result[namespace] = records
    if {key: len(value) for key, value in result.items()} != {"MVZ": 97, "RAP": 47}:
        raise ValueError("Preserved answer must contain 97 MVZ and 47 RAP rows")
    return result


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def species_name(row: dict[str, str]) -> str:
    return row["Species"] if row["namespace"] == "MVZ" else f"D. {row['Species']}"


def namespace_evidence(specimen: str, filenames: set[str]) -> tuple[str, str]:
    explicit = re.fullmatch(r"(MVZ:|RAP)(\d+)", specimen)
    if explicit:
        return ("MVZ" if explicit[1] == "MVZ:" else "RAP", "explicit_specimen_identifier")
    labels = {
        ns
        for name in filenames
        for ns in ("MVZ", "RAP")
        if re.search(rf"(?:^|[^A-Za-z0-9]){ns}[:_ -]?{re.escape(specimen)}(?:[^0-9]|$)", name)
    }
    if len(labels) == 1:
        return next(iter(labels)), "explicit_source_filename"
    return "", "namespace_not_recorded"


def build_associations(
    inventory: list[dict[str, str]], records: list[dict[str, str]], metadata: list[dict[str, str]]
) -> list[dict[str, str]]:
    source_index = defaultdict(list)
    for row in records:
        source_index[row["numeric_identifier"]].append(row)
    metadata_index = defaultdict(list)
    for row in metadata:
        metadata_index[(row["filename"], row["specimen_id"])].append(row)
    groups = defaultdict(list)
    for row in inventory:
        groups[(row["assay"], row["specimen_id"], row["species"])].append(row)
    result = []
    for (assay, specimen, species), objects in sorted(groups.items()):
        filenames = {row["filename"] for row in objects}
        namespace, evidence = namespace_evidence(specimen, filenames)
        numeric = re.sub(r"^(MVZ:|RAP)", "", specimen)
        candidates = source_index[numeric]
        if namespace:
            candidates = [row for row in candidates if row["namespace"] == namespace]
        if assay == "reference_standard":
            namespace, evidence = "MVZ", "author_designated_standard_and_supplied_MVZ_record"
            candidates = [row for row in candidates if row["namespace"] == namespace]
        status = (
            "unmatched"
            if not candidates
            else "ambiguous_source_rows"
            if len(candidates) > 1
            else "confirmed_namespace"
            if namespace
            else "unique_numeric_candidate"
        )
        preserved = set()
        all_metadata_present = True
        for name in filenames:
            rows = metadata_index[(name, numeric)]
            if len(rows) != 1:
                all_metadata_present = False
            preserved.update(row["preservation"] for row in rows)
        base = {
            "assay": assay,
            "analysis_specimen_id": specimen,
            "analysis_species": species,
            "object_count": str(len(objects)) if assay != "reference_standard" else "",
            "image_count": str(len(filenames)),
            "image_filenames": "|".join(sorted(filenames)),
            "namespace_evidence": evidence,
            "association_status": status,
            "candidate_row_count": str(len(candidates)),
            "metadata_image_identity_complete": str(all_metadata_present).lower(),
            "historical_preservation_labels": "|".join(sorted(preserved)),
        }
        for source in candidates or [None]:
            result.append(
                {
                    **base,
                    "source_namespace": source["namespace"] if source else "",
                    "source_identifier": source["source_identifier"] if source else "",
                    "source_row": source["source_row"] if source else "",
                    "source_species": species_name(source) if source else "",
                    "suggested_species": source.get("Suggested Species", "") if source else "",
                    "species_exact_agreement": str(species_name(source) == species).lower()
                    if source
                    else "",
                    "collection_date_verbatim": source.get("Verbatim Date", "") if source else "",
                    "collection_day": source.get("Day", "") if source else "",
                    "collection_month": source.get("Month", "") if source else "",
                    "collection_year": source.get("Year", "") if source else "",
                    "ambiguous_spreadsheet_date": source.get("Date", "") if source else "",
                    "locality": source.get("Specific Locality", source.get("Detailed Locality", ""))
                    if source
                    else "",
                    "permit_verbatim": source.get("Permit", "") if source else "",
                    "source_comment": source.get("Source comment", "") if source else "",
                }
            )
    return result


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(repo: Path, source_dir: Path, output: Path) -> None:
    if output.exists():
        raise ValueError(f"Output directory already exists: {output}")
    for protected in ("Publication", "data", "analyses", "outputs"):
        if output.resolve().is_relative_to((repo / protected).resolve()):
            raise ValueError("Write the audit outside retained scientific artifact directories")
    raw = source_dir / "AUTHOR_SUPPLIED_RECORDS.txt"
    metadata_path = source_dir / "MASTER_IMAGE_METADATA.csv"
    tables = parse_tables(raw)
    metadata_source = json.loads((source_dir / "METADATA_SOURCE.json").read_text())
    if sha256(metadata_path) != metadata_source["sha256"]:
        raise ValueError("Metadata snapshot checksum differs from its source record")
    metadata = read_csv(metadata_path)
    inventory = []
    pins = {
        "author_supplied_records": sha256(raw),
        "master_image_metadata_snapshot": sha256(metadata_path),
        "metadata_source_record": sha256(source_dir / "METADATA_SOURCE.json"),
        "author_source_record": sha256(source_dir / "SOURCE.json"),
    }
    assay_counts = {}
    for assay, (filename, expected_objects, expected_animals) in ASSAYS.items():
        path = repo / "analyses/02_morphology/data" / filename
        rows = read_csv(path)
        animals = {row["specimen_id"] for row in rows}
        if len(rows) != expected_objects or len(animals) != expected_animals:
            raise ValueError(f"Unexpected retained {assay} inventory; no output written")
        pins[filename] = sha256(path)
        assay_counts[assay] = {"objects": len(rows), "animals": len(animals)}
        inventory.extend({"assay": assay, **row} for row in rows)
    for slide, specimen in STANDARDS.items():
        rows = [
            row for row in metadata if row["slide_id"] == slide and row["specimen_id"] == specimen
        ]
        if len(rows) != 1:
            raise ValueError(f"Standard {slide}/{specimen} metadata is missing or ambiguous")
        inventory.append({"assay": "reference_standard", **rows[0]})
    associations = build_associations(inventory, tables["MVZ"] + tables["RAP"], metadata)
    summary = {}
    for assay in (*ASSAYS, "reference_standard"):
        subset = [row for row in associations if row["assay"] == assay]
        grouped = {
            (r["analysis_specimen_id"], r["analysis_species"]): r["association_status"]
            for r in subset
        }
        summary[assay] = dict(Counter(grouped.values()))
    output.mkdir(parents=True)
    for namespace, rows in tables.items():
        write_csv(output / f"{namespace.lower()}_source_records.csv", rows)
    write_csv(output / "assay_specimen_associations.csv", associations)
    receipt = {
        "raw_rows": {key: len(rows) for key, rows in tables.items()},
        "assay_counts": assay_counts,
        "association_counts": summary,
        "input_sha256": pins,
        "script_sha256": sha256(Path(__file__)),
        "source_kind": "author-supplied paste; Google Sheets metadata access denied",
        "release_modified": False,
    }
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    unmatched = sorted(
        {
            row["analysis_specimen_id"]
            for row in associations
            if row["association_status"] == "unmatched"
        }
    )
    differences = [row for row in associations if row["species_exact_agreement"] == "false"]
    permits = Counter(
        row["permit_verbatim"]
        for row in associations
        if row["assay"] == "nuclear_iod" and row["source_namespace"] == "RAP"
    )
    lines = [
        "# Specimen record associations",
        "",
        "The supplied answer contains 97 MVZ records and 47 RAP records. The CSVs retain original species names, suggested names, comments and all source fields. `raw_fields_json` preserves the original field counts, including absent trailing empty fields.",
        "",
        "Released specimen IDs have no MVZ/RAP namespace. A unique numeric correspondence is a candidate association, pending an explicit namespace record or author confirmation. Species agreement is reported separately and never changes the analysis taxonomy. The two reference standards have an author-designated MVZ association.",
        "",
        "| Assay | Confirmed namespace | Unique numeric candidate | Unmatched | Ambiguous |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for assay, counts in summary.items():
        lines.append(
            f"| {assay} | {counts.get('confirmed_namespace', 0)} | {counts.get('unique_numeric_candidate', 0)} | {counts.get('unmatched', 0)} | {counts.get('ambiguous_source_rows', 0)} |"
        )
    lines.extend(
        [
            "",
            f"Unmatched analysis specimen IDs: {', '.join(unmatched) or 'none'}. Do not substitute RAP3462 for analysis specimen 3463.",
            "",
            "MVZ `Date` is an ambiguous spreadsheet/lab date. Collection information comes from `Verbatim Date`; incomplete dates remain verbatim. RAP dates retain separate day, month and year fields. Collection records and historical `Dried Blood`/`Fixed` labels do not establish lifetime fixation and preservation.",
            "",
            "Process_413/MVZ32469 and Process_414/MVZ32470 are the intended standards. Both source records report 7 Apr 1937, 6 mi E Calderwood, Tennessee. Released target D. fuscus measurements use other animals, including numeric IDs 3372/3373 associated provisionally with RAP3372/RAP3373. The standards are not members of the released 42/50-animal target panels.",
            "",
            "Candidate RAP permit counts for the 50-animal IOD inventory:",
            "",
            "| Permit verbatim | Candidate animals |",
            "| --- | ---: |",
        ]
    )
    lines.extend(f"| {permit} | {count} |" for permit, count in sorted(permits.items()))
    lines.extend(
        [
            "",
            "Exact species-name disagreements:",
            "",
            "| Assay | Analysis ID | Analysis species | Source species | Suggested species |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    lines.extend(
        f"| {r['assay']} | {r['analysis_specimen_id']} | {r['analysis_species']} | {r['source_species']} | {r['suggested_species']} |"
        for r in differences
    )
    lines.extend(
        [
            "",
            "RAP3518 retains the supplied `Contaminated???` comment. This is a source flag, not an automatic exclusion. MVZ suggested reidentifications remain in the source CSV and association table.",
            "",
            "The receipt pins the source answer, metadata snapshot, released inventories and script. The script refuses existing output directories and writes no retained analysis artifact.",
        ]
    )
    (output / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.repo_root, args.source_dir or args.repo_root / RECORD_DIR, args.output_dir)
    except (ValueError, KeyError) as error:
        parser.exit(2, f"Error: {error}\n")


if __name__ == "__main__":
    main()

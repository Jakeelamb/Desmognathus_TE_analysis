#!/usr/bin/env python3
"""
Resolve and optionally download the genome assemblies referenced in the lookup table.

This script uses the canonical `Genome_Accension` column in `input_data/lookup_table.txt`
to build a verified local genome FASTA cache under `input_data/genomes/`. It writes a
manifest to `results/data/ltr_age/` so the LTR-age workflow can distinguish resolved,
downloaded, and checksum-verified assemblies from unresolved or missing ones.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT, load_lookup_table, paths  # noqa: E402


MANIFEST_PATH = paths.results.data / "ltr_age" / "genome_assembly_manifest.csv"
DOWNLOAD_ROOT = getattr(paths.input_data, "genomes", PROJECT_ROOT / "input_data" / "genomes")
DIRECTORY_RE_TEMPLATE = r'href="({accession}[^\"]*/)"'
CONTENT_LENGTH_RE = re.compile(r"(?im)^content-length:\s*(\d+)\s*$")
MD5_RE_TEMPLATE = r"^([0-9a-f]{{32}})\s+\./{filename}$"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve and download genome FASTA assemblies referenced in the lookup table."
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download and checksum-verify the resolved assembly FASTA files.",
    )
    parser.add_argument(
        "--species",
        nargs="*",
        default=None,
        help="Optional species subset matching the `Species` column in the lookup table.",
    )
    parser.add_argument(
        "--accessions",
        nargs="*",
        default=None,
        help="Optional accession subset matching the `Genome_Accension` column.",
    )
    parser.add_argument(
        "--no-head",
        action="store_true",
        help="Skip remote HEAD requests for content-length metadata.",
    )
    return parser.parse_args()


def accession_to_parent_url(gca_accession: str) -> str:
    digits = gca_accession.split("_", 1)[1].split(".", 1)[0]
    triplets = [digits[i : i + 3] for i in range(0, len(digits), 3)]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{'/'.join(triplets)}/"


def run_curl(args: List[str], text: bool = True) -> str:
    result = subprocess.run(
        ["curl", "--retry", "3", "--retry-delay", "2", *args],
        check=True,
        capture_output=True,
        text=text,
    )
    return result.stdout


def resolve_accession(gca_accession: str, include_head: bool = True) -> Dict[str, str]:
    parent_url = accession_to_parent_url(gca_accession)
    directory_html = run_curl(["-fsSL", parent_url])
    directory_re = re.compile(DIRECTORY_RE_TEMPLATE.format(accession=re.escape(gca_accession)))
    match = directory_re.search(directory_html)
    if not match:
        raise ValueError(f"Could not resolve assembly directory for {gca_accession} at {parent_url}")

    assembly_dir = match.group(1).rstrip("/")
    fasta_name = f"{assembly_dir}_genomic.fna.gz"
    fasta_url = f"{parent_url}{assembly_dir}/{fasta_name}"
    md5_url = f"{parent_url}{assembly_dir}/md5checksums.txt"

    md5_text = run_curl(["-fsSL", md5_url])
    md5_re = re.compile(MD5_RE_TEMPLATE.format(filename=re.escape(fasta_name)), re.MULTILINE)
    md5_match = md5_re.search(md5_text)
    if not md5_match:
        raise ValueError(f"Could not find FASTA MD5 for {gca_accession} in {md5_url}")

    remote_bytes = ""
    if include_head:
        head_text = run_curl(["-fsSLI", fasta_url])
        head_match = CONTENT_LENGTH_RE.search(head_text)
        if head_match:
            remote_bytes = head_match.group(1)

    return {
        "parent_url": parent_url,
        "assembly_dir": assembly_dir,
        "fasta_name": fasta_name,
        "fasta_url": fasta_url,
        "md5_url": md5_url,
        "expected_md5": md5_match.group(1),
        "remote_bytes": remote_bytes,
    }


def compute_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "curl",
            "--retry",
            "3",
            "--retry-delay",
            "2",
            "-fL",
            "-C",
            "-",
            "-o",
            str(destination),
            url,
        ],
        check=True,
    )


def rotate_bad_file(path: Path) -> Path:
    candidate = path.with_name(path.name + ".md5_mismatch")
    counter = 2
    while candidate.exists():
        candidate = path.with_name(path.name + f".md5_mismatch_{counter}")
        counter += 1
    path.rename(candidate)
    return candidate


def filter_lookup(lookup, species_filter: Optional[Iterable[str]], accession_filter: Optional[Iterable[str]]):
    filtered = lookup.copy()
    if species_filter:
        wanted = set(species_filter)
        filtered = filtered[filtered["Species"].isin(wanted)].copy()
    if accession_filter:
        wanted = set(accession_filter)
        filtered = filtered[filtered["Genome_Accension"].isin(wanted)].copy()
    return filtered.sort_values(["Species", "Genome_Accension"]).reset_index(drop=True)


def local_download_path(assembly_dir: str, fasta_name: str) -> Path:
    return Path(DOWNLOAD_ROOT) / assembly_dir / fasta_name


def write_manifest(rows: List[Dict[str, str]]) -> None:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "species",
        "gca_accession",
        "parent_url",
        "assembly_dir",
        "fasta_name",
        "fasta_url",
        "md5_url",
        "expected_md5",
        "remote_bytes",
        "local_path",
        "local_exists",
        "local_bytes",
        "local_md5",
        "md5_match",
        "status",
        "notes",
    ]
    with MANIFEST_PATH.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    lookup = filter_lookup(load_lookup_table(), args.species, args.accessions)
    if lookup.empty:
        raise SystemExit("No lookup rows matched the requested filters.")

    Path(DOWNLOAD_ROOT).mkdir(parents=True, exist_ok=True)
    manifest_rows: List[Dict[str, str]] = []
    n_verified = 0
    n_downloaded = 0
    n_failed = 0

    for row in lookup.to_dict("records"):
        species = row["Species"]
        gca_accession = row["Genome_Accension"]
        manifest_row: Dict[str, str] = {
            "species": species,
            "gca_accession": gca_accession,
            "parent_url": "",
            "assembly_dir": "",
            "fasta_name": "",
            "fasta_url": "",
            "md5_url": "",
            "expected_md5": "",
            "remote_bytes": "",
            "local_path": "",
            "local_exists": "False",
            "local_bytes": "",
            "local_md5": "",
            "md5_match": "",
            "status": "",
            "notes": "",
        }
        try:
            resolved = resolve_accession(gca_accession, include_head=not args.no_head)
            manifest_row.update(resolved)
            destination = local_download_path(resolved["assembly_dir"], resolved["fasta_name"])
            manifest_row["local_path"] = str(destination.relative_to(PROJECT_ROOT))

            local_md5 = ""
            local_exists = destination.exists()
            if local_exists:
                local_md5 = compute_md5(destination)
                manifest_row["local_exists"] = "True"
                manifest_row["local_bytes"] = str(destination.stat().st_size)
                manifest_row["local_md5"] = local_md5
                manifest_row["md5_match"] = str(local_md5 == resolved["expected_md5"])

            if args.download:
                if local_exists and local_md5 == resolved["expected_md5"]:
                    manifest_row["status"] = "verified_existing"
                    n_verified += 1
                else:
                    if local_exists:
                        rotated = rotate_bad_file(destination)
                        manifest_row["notes"] = f"Rotated checksum-mismatch file to {rotated.name}"
                    download_file(resolved["fasta_url"], destination)
                    local_md5 = compute_md5(destination)
                    manifest_row["local_exists"] = "True"
                    manifest_row["local_bytes"] = str(destination.stat().st_size)
                    manifest_row["local_md5"] = local_md5
                    manifest_row["md5_match"] = str(local_md5 == resolved["expected_md5"])
                    if local_md5 != resolved["expected_md5"]:
                        raise ValueError(
                            f"MD5 mismatch after download for {gca_accession}: "
                            f"{local_md5} != {resolved['expected_md5']}"
                        )
                    manifest_row["status"] = "downloaded_verified"
                    n_downloaded += 1
            else:
                if local_exists and local_md5 == resolved["expected_md5"]:
                    manifest_row["status"] = "local_verified"
                    n_verified += 1
                elif local_exists:
                    manifest_row["status"] = "local_md5_mismatch"
                else:
                    manifest_row["status"] = "resolved_remote_only"
        except Exception as exc:  # noqa: BLE001
            manifest_row["status"] = "error"
            manifest_row["notes"] = str(exc)
            n_failed += 1

        manifest_rows.append(manifest_row)

    write_manifest(manifest_rows)
    print(f"Wrote manifest: {MANIFEST_PATH.relative_to(PROJECT_ROOT)}")
    print(f"Assemblies resolved: {len(manifest_rows) - n_failed}/{len(manifest_rows)}")
    print(f"Assemblies verified locally: {n_verified}")
    if args.download:
        print(f"Assemblies downloaded this run: {n_downloaded}")
    print(f"Assembly errors: {n_failed}")
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "raw"
OUTPUT_FILE = RAW_DIR / "source_file_inventory.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def infer_source_id(path: Path) -> str:
    if path.name == "AmphiBIO_v1.zip":
        return "amphibio_2017"
    if path.name == "Pyron_Beamer_conanti_fuscus_Appendix_S1.csv":
        return "pyron_beamer_2023_conanti_fuscus_appendix_s1"
    if path.parent.name == "plazi_treatments":
        return "pyron_beamer_2022_ocoee_treatmentbank"
    if path.parent.name == "nc_biodiversity_species":
        return "nc_biodiversity_amphibians"
    if path.parent.name == "web_sources":
        if "valentinei_louisianaherps" in path.name:
            return "louisianaherps_valentinei"
        if "valentinei_aquila" in path.name or "valentinei_zootaxa4263" in path.name:
            return "means_lamb_bernardo_2017_valentinei_article"
        if "pascagoula_plazi" in path.name:
            return "pyron_oconnell_lamb_beamer_2022_pascagoula_treatmentbank"
    return "unmapped_source"


def infer_source_url(path: Path) -> str:
    if path.name == "AmphiBIO_v1.zip":
        return "https://figshare.com/articles/dataset/Oliveira_et_al_AmphiBIO_v1/4644424"
    if path.name == "Pyron_Beamer_conanti_fuscus_Appendix_S1.csv":
        return "https://zenodo.org/records/7972011"
    if path.parent.name == "plazi_treatments":
        return "https://treatment.plazi.org/GgServer/xhtml/039087A5C5193F4AC4BA191F72C40EB6"
    if path.parent.name == "nc_biodiversity_species" and path.stem.startswith("Desmognathus_"):
        species = path.stem.replace("Desmognathus_", "")
        return f"https://auth1.dpr.ncparks.gov/amphibians/view.php?sciName=Desmognathus+{species}"
    if path.parent.name == "web_sources":
        if "valentinei_louisianaherps" in path.name:
            return "https://www.louisianaherps.com/southern-dusky-salamander-d.html"
        if "valentinei_aquila" in path.name:
            return "https://aquila.usm.edu/fac_pubs/17676/"
        if "valentinei_zootaxa4263" in path.name:
            return "https://www.biotaxa.org/Zootaxa/article/view/zootaxa.4263.3.3"
        if "pascagoula_plazi" in path.name:
            return "https://treatment.plazi.org/GgServer/xhtml/03CE87E1FFC1FA62FF44FB7CFC539344"
    return ""


def infer_notes(path: Path) -> str:
    if path.name == "AmphiBIO_v1.zip":
        return "Downloaded from Figshare API article 4644424 and used to generate Desmognathus trait subset"
    if path.name == "Pyron_Beamer_conanti_fuscus_Appendix_S1.csv":
        return "Downloaded directly from Zenodo and used to generate conanti-fuscus morphometric summaries"
    if path.parent.name == "plazi_treatments":
        return "XHTML snapshot of TreatmentBank rendering used for manual trait extraction"
    if path.parent.name == "nc_biodiversity_species":
        return "HTML snapshot of NC Biodiversity Project species account used for traceable trait extraction"
    if path.parent.name == "web_sources":
        return "HTML or XHTML snapshot of web source used for traceable gap-filling and validation"
    return ""


def main() -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(RAW_DIR.rglob("*")):
        if not path.is_file():
            continue
        if path.name == OUTPUT_FILE.name:
            continue
        if path.name == ".gitkeep":
            continue
        rows.append(
            {
                "source_id": infer_source_id(path),
                "local_filename": str(path.relative_to(RAW_DIR)),
                "sha256": sha256_for_file(path),
                "downloaded_on": pd.Timestamp(path.stat().st_mtime, unit="s").date().isoformat(),
                "source_url": infer_source_url(path),
                "notes": infer_notes(path),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Files: {len(out)}")


if __name__ == "__main__":
    main()

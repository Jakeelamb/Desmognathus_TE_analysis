#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import re
from html import unescape
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "raw" / "nc_biodiversity_species"
DERIVED_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"

OUTPUT_FILE = DERIVED_DIR / "nc_biodiversity_desmognathus_traits.csv"


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def standardize_species(value: str) -> str:
    return value.strip().split()[-1].lower()


def account_url(species: str) -> str:
    return f"https://auth1.dpr.ncparks.gov/amphibians/view.php?sciName=Desmognathus+{species}"


def clean_html_fragment(fragment: str | None) -> str | None:
    if fragment is None:
        return None
    cleaned = fragment
    cleaned = re.sub(r"(?i)<br\s*/?>", "\n", cleaned)
    cleaned = re.sub(r"(?i)</p>|</tr>|</td>|</li>|</fieldset>", "\n", cleaned)
    cleaned = re.sub(r"(?is)<script.*?</script>", " ", cleaned)
    cleaned = re.sub(r"(?is)<style.*?</style>", " ", cleaned)
    cleaned = re.sub(r"(?s)<[^>]+>", " ", cleaned)
    cleaned = unescape(cleaned).replace("\xa0", " ")
    lines = [re.sub(r"\s+", " ", line).strip() for line in cleaned.splitlines()]
    lines = [line for line in lines if line]
    return "\n".join(lines) if lines else None


def extract_field(html: str, label: str) -> str | None:
    pattern = re.compile(
        rf"<b>\s*{re.escape(label)}\s*:\s*</b>(.*?)(?:</td>|</tr>)",
        flags=re.IGNORECASE | re.DOTALL,
    )
    match = pattern.search(html)
    if not match:
        return None
    return clean_html_fragment(match.group(1))


def extract_page_title_common_name(html: str) -> str | None:
    match = re.search(
        r"Photo Gallery for <i>Desmognathus [^<]+</i>\s*-\s*([^<]+)</font>",
        html,
        flags=re.IGNORECASE,
    )
    if not match:
        return None
    return clean_html_fragment(match.group(1))


def parse_numeric(value: str) -> float:
    return float(value.replace(",", ""))


def parse_measurement_range(text: str | None, token: str) -> tuple[float | None, float | None]:
    if not text:
        return (None, None)
    normalized = text.replace("\u2013", "-").replace("\u2014", "-")
    patterns = {
        "svl": [
            r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
            r"(\d+(?:\.\d+)?)\s*to\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
            r"(\d+(?:\.\d+)?)\s*mm\s*SVL\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        ],
        "tl": [
            r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
            r"(\d+(?:\.\d+)?)\s*to\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
        ],
    }
    for pattern in patterns[token]:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (parse_numeric(match.group(1)), parse_numeric(match.group(2)))

    if token == "svl":
        single = re.search(r"(\d+(?:\.\d+)?)\s*mm\s*SVL", normalized, flags=re.IGNORECASE)
        if single:
            value = parse_numeric(single.group(1))
            return (None, value)
    if token == "tl":
        single = re.search(r"(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b", normalized, flags=re.IGNORECASE)
        if single:
            value = parse_numeric(single.group(1))
            return (None, value)
    return (None, None)


def parse_description_size(text: str | None) -> tuple[float | None, float | None, str | None]:
    if not text:
        return (None, None, None)
    normalized = text.replace("\u2013", "-").replace("\u2014", "-")

    mixed_stage_range_patterns = [
        r"juveniles and adults ranging from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"(?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL\s*for juveniles and adults",
        r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL\s*for transformed specimens",
    ]
    for pattern in mixed_stage_range_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (parse_numeric(match.group(1)), parse_numeric(match.group(2)), "mixed_stage_range")

    adult_range_patterns = [
        r"SVL of adults ranging from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm",
        r"snout-vent-length of\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm",
        r"adults ranging from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"adults range from (?:about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"adult size\s*\((\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"larger adult size\s*\((\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
    ]
    for pattern in adult_range_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (parse_numeric(match.group(1)), parse_numeric(match.group(2)), "adult_range")

    max_patterns = [
        r"maximum size of the adults (?:ca\.?|about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"maximum adult size (?:ca\.?|about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"relatively small \(maximum size of the adults (?:ca\.?|about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*mm\s*SVL\)",
        r"largest adults reaching (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
        r"up to\s*(\d+(?:\.\d+)?)\s*mm\s*SVL",
    ]
    for pattern in max_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (None, parse_numeric(match.group(1)), "adult_max")

    return (None, None, None)


def parse_description_total_length(text: str | None) -> tuple[float | None, float | None]:
    if not text:
        return (None, None)
    normalized = text.replace("\u2013", "-").replace("\u2014", "-")

    adult_range_patterns = [
        r"adults vary from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
        r"adults range from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
        r"adult size\s*\((\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
    ]
    for pattern in adult_range_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (parse_numeric(match.group(1)), parse_numeric(match.group(2)))

    adult_range_patterns_cm = [
        r"adults vary from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*cm\s*(?:total length|TL)\b",
        r"adults range from (?:around|about|approximately)?\s*(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)\s*cm\s*(?:total length|TL)\b",
    ]
    for pattern in adult_range_patterns_cm:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (parse_numeric(match.group(1)) * 10.0, parse_numeric(match.group(2)) * 10.0)

    adult_max_patterns = [
        r"maximum adult size (?:ca\.?|about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
        r"up to\s*(\d+(?:\.\d+)?)\s*mm\s*(?:total length|TL)\b",
    ]
    for pattern in adult_max_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (None, parse_numeric(match.group(1)))

    adult_max_patterns_cm = [
        r"maximum adult size (?:ca\.?|about|around|approximately)?\s*(\d+(?:\.\d+)?)\s*cm\s*(?:total length|TL)\b",
        r"up to\s*(\d+(?:\.\d+)?)\s*cm\s*(?:total length|TL)\b",
    ]
    for pattern in adult_max_patterns_cm:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            return (None, parse_numeric(match.group(1)) * 10.0)

    return (None, None)


def parse_elevation_range(*texts: str | None) -> tuple[float | None, float | None]:
    combined = " ".join(text for text in texts if text)
    if not combined:
        return (None, None)
    normalized = combined.replace("\u2013", "-").replace("\u2014", "-")

    range_patterns = [
        r"elevations?[^.]{0,120}?(\d{1,4}(?:,\d{3})?)\s*-\s*(\d{1,4}(?:,\d{3})?)\s*m\b",
        r"between\s+(\d{1,4}(?:,\d{3})?)\s+and\s+(\d{1,4}(?:,\d{3})?)\s*m\b",
    ]
    for pattern in range_patterns:
        match = re.search(pattern, normalized, flags=re.IGNORECASE)
        if match:
            low = parse_numeric(match.group(1))
            high = parse_numeric(match.group(2))
            if low <= high:
                return (low, high)

    above_match = re.search(r"above\s+(?:approximately\s+)?(\d{1,4}(?:,\d{3})?)\s*m\b", normalized, flags=re.IGNORECASE)
    if above_match:
        return (parse_numeric(above_match.group(1)), None)

    below_match = re.search(r"below\s+(?:approximately\s+)?(\d{1,4}(?:,\d{3})?)\s*m\b", normalized, flags=re.IGNORECASE)
    if below_match:
        return (None, parse_numeric(below_match.group(1)))

    return (None, None)


def infer_development_mode(*texts: str | None) -> tuple[str | None, str | None]:
    combined = " ".join(text.lower() for text in texts if text)
    if not combined:
        return (None, None)

    direct_keywords = [
        "direct development",
        "direct-develop",
        "hatch as miniature",
        "no aquatic larval stage",
        "no larval stage",
        "without a larval stage",
    ]
    if any(keyword in combined for keyword in direct_keywords):
        return ("direct_development", "explicit_direct_development_text")

    larval_keywords = [
        "aquatic larval stage",
        "larval period",
        "larvae transform",
        "larvae metamorphose",
        "metamorphose",
        "metamorphs",
        "hatchlings presumably move",
        "hatchlings move",
    ]
    if any(keyword in combined for keyword in larval_keywords):
        return ("aquatic_larva", "larval_stage_text")

    larval_morphology_keywords = [
        "metamorphosed specimens",
        "metamorphosed individuals",
        "paired larval spots",
        "larval spots",
    ]
    if any(keyword in combined for keyword in larval_morphology_keywords):
        return ("aquatic_larva", "larval_morphology_description")

    return (None, None)


def infer_microhabitat_class(*texts: str | None) -> tuple[str | None, float | None, str | None]:
    combined = " ".join(text.lower() for text in texts if text)
    if not combined:
        return (None, None, None)

    stream_aquatic_keywords = [
        "almost entirely aquatic",
        "almost always found in water",
        "rarely leaves the water",
        "stream proper",
        "under submerged stones",
        "shallow riffle",
        "fast current",
        "small to medium-sized streams",
    ]
    if any(keyword in combined for keyword in stream_aquatic_keywords):
        return ("stream_aquatic", 2.0, "strong_stream_aquatic_keywords")

    streamside_keywords = [
        "streamside",
        "stream habitats",
        "small streams",
        "lower-order stream",
        "seeps",
        "seepages",
        "springs",
        "headwater",
        "ravine",
        "rockfaces",
    ]
    if any(keyword in combined for keyword in streamside_keywords):
        return ("streamside", 1.0, "streamside_seep_keywords")

    swamp_keywords = ["slough", "floodplain", "bottomland", "blackwater", "peat", "mud-bottomed", "wetland"]
    if any(keyword in combined for keyword in swamp_keywords) or (
        "swamp" in combined and ("coastal plain" in combined or "bottomland" in combined)
    ):
        return ("swamp_coastal_plain", 3.0, "swamp_wetland_keywords")

    woodland_keywords = ["forest floor", "woodland", "terrestrial", "cove forest", "leaf litter"]
    if any(keyword in combined for keyword in woodland_keywords):
        return ("mountain_woodland", 0.0, "woodland_keywords")

    return (None, None, None)


def parse_species_page(path: Path) -> dict[str, object]:
    html = path.read_text(encoding="utf-8", errors="ignore")
    species = standardize_species(path.stem.replace("Desmognathus_", "Desmognathus "))
    description_text = extract_field(html, "Description")
    size_text = extract_field(html, "Adult Body Size")
    habitat_text = extract_field(html, "Habitat")
    reproductive_text = extract_field(html, "Reproductive Mode")
    aquatic_text = extract_field(html, "Aquatic Life History")
    terrestrial_text = extract_field(html, "Terrestrial Life History")
    status_comments = extract_field(html, "Status Comments")

    svl_min, svl_max = parse_measurement_range(size_text, "svl")
    tl_min, tl_max = parse_measurement_range(size_text, "tl")
    size_source_field = "adult_body_size_field" if pd.notna(svl_min) or pd.notna(svl_max) else None
    size_phrase_type = "adult_range" if pd.notna(svl_min) or pd.notna(svl_max) else None
    if pd.isna(svl_min) and pd.isna(svl_max):
        desc_svl_min, desc_svl_max, desc_size_type = parse_description_size(description_text)
        svl_min = desc_svl_min
        svl_max = desc_svl_max
        size_source_field = "description_field" if pd.notna(desc_svl_min) or pd.notna(desc_svl_max) else None
        size_phrase_type = desc_size_type
    if pd.isna(tl_min) and pd.isna(tl_max):
        desc_tl_min, desc_tl_max = parse_description_total_length(description_text)
        tl_min = desc_tl_min
        tl_max = desc_tl_max
        if pd.notna(desc_tl_min) or pd.notna(desc_tl_max):
            size_source_field = size_source_field or "description_field"
            if pd.notna(desc_tl_min):
                size_phrase_type = size_phrase_type or "adult_tl_range"
            else:
                size_phrase_type = size_phrase_type or "adult_tl_max"
    elevation_min, elevation_max = parse_elevation_range(habitat_text, status_comments)
    development_mode, development_rule = infer_development_mode(
        description_text, reproductive_text, aquatic_text, terrestrial_text
    )
    microhabitat_class, aquaticity_index, lifestyle_rule = infer_microhabitat_class(
        habitat_text, aquatic_text, terrestrial_text
    )

    return {
        "species": species,
        "nc_common_name": extract_page_title_common_name(html),
        "nc_description_text": description_text,
        "nc_body_size_text": size_text,
        "nc_habitat_text": habitat_text,
        "nc_reproductive_mode_text": reproductive_text,
        "nc_aquatic_life_history_text": aquatic_text,
        "nc_terrestrial_life_history_text": terrestrial_text,
        "nc_status_comments_text": status_comments,
        "nc_adult_svl_min_mm": svl_min,
        "nc_adult_svl_max_mm": svl_max,
        "nc_adult_svl_mid_mm": np.mean([svl_min, svl_max]) if pd.notna(svl_min) and pd.notna(svl_max) else np.nan,
        "nc_size_source_field": size_source_field,
        "nc_size_phrase_type": size_phrase_type,
        "nc_adult_tl_min_mm": tl_min,
        "nc_adult_tl_max_mm": tl_max,
        "nc_elevation_min_m": elevation_min,
        "nc_elevation_max_m": elevation_max,
        "nc_elevation_mid_m": (
            np.mean([elevation_min, elevation_max]) if pd.notna(elevation_min) and pd.notna(elevation_max) else np.nan
        ),
        "nc_development_mode": development_mode,
        "nc_development_rule": development_rule,
        "nc_microhabitat_class": microhabitat_class,
        "nc_aquaticity_index": aquaticity_index,
        "nc_lifestyle_rule": lifestyle_rule,
        "nc_source_id": "nc_biodiversity_amphibians",
        "nc_source_filename": path.name,
        "nc_source_file_sha256": sha256_for_file(path),
        "nc_source_url": account_url(species),
    }


def main() -> None:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    if not RAW_DIR.exists():
        raise FileNotFoundError(
            f"Missing raw source directory: {RAW_DIR}\n"
            "Download NC Biodiversity Project species pages into path_analysis/data/external/raw/nc_biodiversity_species/."
        )

    html_files = sorted(RAW_DIR.glob("Desmognathus_*.html"))
    if not html_files:
        raise FileNotFoundError(
            f"No HTML source files found in {RAW_DIR}\n"
            "Expected files like Desmognathus_kanawha.html."
        )

    rows = [parse_species_page(path) for path in html_files]
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    out.to_csv(OUTPUT_FILE, index=False)

    print(f"Wrote {OUTPUT_FILE}")
    print(f"Species: {len(out)}")


if __name__ == "__main__":
    main()

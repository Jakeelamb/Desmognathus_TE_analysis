from __future__ import annotations

import csv
import hashlib
import importlib.util
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlsplit

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_candidate_review", ROOT / "scripts/publication/build_candidate_review.py"
)
assert SPEC is not None and SPEC.loader is not None
review = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = review
SPEC.loader.exec_module(review)


def test_current_gallery_manifests_use_estimated_genome_size_display_wording() -> None:
    """Check reviewed captions and titles without renaming durable asset identifiers."""
    gallery_root = ROOT / "Publication/figures/candidates"
    page = Page((gallery_root / "index.html").read_text())
    candidates = set()
    for image in page.images:
        relative_path = Path(unquote(urlsplit(image["src"]).path))
        directory = (gallery_root / relative_path).parent
        candidates.add(directory.parent if directory.name == "panels" else directory)
    assert candidates, "The current review gallery must contain figure exports"

    obsolete_label = re.compile(r"\bconditional[\s-]+genome[\s-]+sizes?\b", re.IGNORECASE)
    estimated_label = re.compile(r"\bestimated[\s-]+genome[\s-]+sizes?\b", re.IGNORECASE)
    estimates_found = 0
    for directory in sorted(candidates):
        for manifest in (directory / "FIGURE_MANIFEST.csv",
                         directory / "panels/PANEL_MANIFEST.csv"):
            with manifest.open(newline="") as stream:
                for row in csv.DictReader(stream):
                    # Filenames, release statuses, input columns, and user notes
                    # are stable identifiers or provenance, not display labels.
                    for field in ("title", "legend", "alt_text"):
                        text = row.get(field, "")
                        assert not obsolete_label.search(text), (manifest, field, text)
                        estimates_found += bool(estimated_label.search(text))
    assert estimates_found > 0, "The review set should retain genome-size labels"


def test_current_gallery_omits_only_the_pgls_species_influence_panel() -> None:
    gallery_root = ROOT / "Publication/figures/candidates"
    with (gallery_root / "figure_4_pgls/panels/PANEL_MANIFEST.csv").open(newline="") as stream:
        panels = list(csv.DictReader(stream))
    assert [(row["parent_figure"], row["panel_id"]) for row in panels] == [
        ("Figure_4_PGLS", "a")
    ]
    assert panels[0]["pdf"] == "Figure_4_PGLS_a_ratio_genome_size.pdf"

    page = Page((gallery_root / "index.html").read_text())
    assert [anchor for anchor in page.sections
            if anchor in {"figure-4-pgls", "figure-4-pgls-a", "figure-4-pgls-b"}] == [
        "figure-4-pgls", "figure-4-pgls-a"
    ]
    assert "#figure-4-pgls-a" in page.links
    assert "#figure-4-pgls-b" not in page.links
    note_ids = {note["data-note-id"] for note in page.notes}
    assert "figure-4-pgls-a" in note_ids
    assert "figure-4-pgls-b" not in note_ids

    # The unrelated four-domain comparison panels also use letter b and remain
    # visible; removing the species-influence card must not hide those panels.
    for directory, anchor in (
        ("figure_4_filter_comparison", "figure-4-filter-b"),
        ("figure_4_pgls_comparison", "figure-4-pgls-filter-b"),
    ):
        with (gallery_root / directory / "panels/PANEL_MANIFEST.csv").open(newline="") as stream:
            comparison = list(csv.DictReader(stream))
        assert [row["panel_id"] for row in comparison] == ["a", "b"]
        assert "Four-domain" in comparison[1]["title"]
        assert anchor in page.sections
        assert f"#{anchor}" in page.links
        assert anchor in note_ids


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def make_candidate(tmp_path: Path, number: int, panels: str) -> Path:
    directory = tmp_path / f"candidate {number}"
    panel_dir = directory / "panels"
    panel_dir.mkdir(parents=True)

    def assets(base: Path, name: str) -> dict[str, str]:
        files = {
            kind: f"{name} & review.{ext}"
            for kind, ext in (("pdf", "pdf"), ("png", "png"), ("tiff", "tif"))
        }
        for filename in files.values():
            (base / filename).write_bytes(b"unchanged figure export")
        return files

    common = {"legend": 'A < B & "quoted".', "alt_text": 'Alt <img src="x"> & text'}
    write_csv(
        directory / "FIGURE_MANIFEST.csv",
        [
            {
                "figure_id": f"Figure_{number}",
                **common,
                **assets(directory, "collage"),
            }
        ],
    )
    write_csv(
        panel_dir / "PANEL_MANIFEST.csv",
        [
            {
                "parent_figure": f"Figure_{number}",
                "panel_id": panel,
                "title": f"Panel {panel} <title> & detail",
                **common,
                **assets(panel_dir, panel),
            }
            for panel in panels
        ],
    )
    return directory


class Page(HTMLParser):
    def __init__(self, source: str) -> None:
        super().__init__(convert_charrefs=True)
        self.sections: list[str] = []
        self.links: list[str] = []
        self.images: list[dict[str, str]] = []
        self.image_links: list[str | None] = []
        self._link_href: str | None = None
        self.tags: list[str] = []
        self.notes: list[dict[str, str]] = []
        self.labels: list[str] = []
        self.feed(source)

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        self.tags.append(tag)
        if tag == "section":
            self.sections.append(attributes["id"])
        elif tag == "a":
            self.links.append(attributes["href"])
            self._link_href = attributes["href"]
        elif tag == "img":
            self.images.append(attributes)
            self.image_links.append(self._link_href)
        elif tag == "textarea":
            self.notes.append(attributes)
        elif tag == "label":
            self.labels.append(attributes["for"])

    def handle_endtag(self, tag: str) -> None:
        if tag == "a":
            self._link_href = None


def test_two_figures_are_ordered_linked_escaped_and_leave_sources_unchanged(tmp_path: Path) -> None:
    candidates = [make_candidate(tmp_path, 1, "abcde"), make_candidate(tmp_path, 2, "abc")]
    before = {
        path: path.read_bytes()
        for directory in candidates
        for path in directory.rglob("*")
        if path.is_file()
    }
    output = tmp_path / "review" / "index.html"
    assert review.build_candidate_review(candidates, output) == output
    source = output.read_text(encoding="utf-8")
    page = Page(source)
    assert page.sections == [
        "figure-1",
        *[f"figure-1-{letter}" for letter in "abcde"],
        "figure-2",
        *[f"figure-2-{letter}" for letter in "abc"],
    ]
    assert len(page.images) == 10
    assert page.links[:10] == [f"#{anchor}" for anchor in page.sections]
    for url in page.links[10:] + [image["src"] for image in page.images]:
        assert not url.startswith(("/", "http", "file:"))
        assert (output.parent / unquote(urlsplit(url).path)).resolve().is_file()
    assert {Path(unquote(urlsplit(url).path)).suffix for url in page.links[10:]} == {
        ".pdf", ".png", ".tif"
    }
    assert len(page.links[10:]) == 40
    assert all(image["alt"] == 'Alt <img src="x"> & text' for image in page.images)
    assert "A &lt; B &amp; &quot;quoted&quot;." in source
    assert "Panel a &lt;title&gt; &amp; detail" in source
    assert page.tags.count("title") == 1
    assert "Arial,sans-serif" in source
    assert "2 draft figures · 8 isolated panels" in source
    assert "Review numbers are provisional" in source
    assert [note["data-note-id"] for note in page.notes] == page.sections
    assert page.labels == [note["id"] for note in page.notes]
    assert all(note["aria-describedby"] for note in page.notes)
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize("extension", ["pdf", "png", "tif"])
def test_asset_versions_follow_only_content_changes_and_preserve_note_keys(
    tmp_path: Path, extension: str
) -> None:
    candidate = make_candidate(tmp_path, 4, "ab")
    output = tmp_path / "review" / "index.html"

    def build() -> Page:
        return Page(review.build_candidate_review([candidate], output).read_text(encoding="utf-8"))

    def asset_versions(page: Page) -> dict[Path, str]:
        versions = {}
        for url in [*page.links, *(image["src"] for image in page.images)]:
            parsed = urlsplit(url)
            if not parsed.path:
                continue
            path = (output.parent / unquote(parsed.path)).resolve()
            expected_version = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
            assert parse_qs(parsed.query) == {"v": [expected_version]}
            versions[path] = url
        assert page.image_links == [image["src"] for image in page.images]
        assert all(page.links.count(image["src"]) == 2 for image in page.images)
        return versions

    before = build()
    original_versions = asset_versions(before)
    assert len(original_versions) == 9
    target = candidate / "panels" / f"b & review.{extension}"
    original_bytes = target.read_bytes()
    target.write_bytes(original_bytes)
    same_content = build()
    assert asset_versions(same_content) == original_versions

    target.write_bytes(b"revised panel export")
    changed = build()
    changed_versions = asset_versions(changed)
    assert changed_versions.keys() == original_versions.keys()
    assert {
        path for path in original_versions if changed_versions[path] != original_versions[path]
    } == {target.resolve()}
    assert before.sections == same_content.sections == changed.sections
    assert [note["data-note-id"] for note in before.notes] == [
        note["data-note-id"] for note in changed.notes
    ]


def test_rebuild_keeps_note_ids_and_never_replaces_author_notes(tmp_path: Path) -> None:
    candidate = make_candidate(tmp_path, 1, "ab")
    output = tmp_path / "review" / "index.html"
    notes_file = tmp_path / "review_notes.json"
    notes_file.write_text('{"author_note":"Keep the larger labels."}', encoding="utf-8")
    first = review.build_candidate_review([candidate], output).read_text(encoding="utf-8")
    manifest = candidate / "FIGURE_MANIFEST.csv"
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["title"] = "A revised descriptive title"
    write_csv(manifest, rows)
    second = review.build_candidate_review([candidate], output).read_text(encoding="utf-8")
    assert [note["data-note-id"] for note in Page(first).notes] == [
        note["data-note-id"] for note in Page(second).notes
    ]
    assert notes_file.read_text(encoding="utf-8") == '{"author_note":"Keep the larger labels."}'


@pytest.mark.parametrize(
    "failure",
    [
        "missing_asset",
        "wrong_parent",
        "duplicate_panel",
        "empty_legend",
        "missing_column",
        "multiple_figures",
    ],
)
def test_invalid_candidate_leaves_existing_page_unchanged(tmp_path: Path, failure: str) -> None:
    candidate = make_candidate(tmp_path, 1, "ab")
    output = tmp_path / "index.html"
    output.write_text("retain current review", encoding="utf-8")
    panel_path = candidate / "panels/PANEL_MANIFEST.csv"
    with panel_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if failure == "missing_asset":
        (panel_path.parent / rows[1]["tiff"]).unlink()
    elif failure == "wrong_parent":
        rows[1]["parent_figure"] = "Figure_2"
    elif failure == "duplicate_panel":
        rows[1]["panel_id"] = "a"
    elif failure == "empty_legend":
        rows[1]["legend"] = " "
    elif failure == "missing_column":
        for row in rows:
            del row["alt_text"]
    elif failure == "multiple_figures":
        figure_path = candidate / "FIGURE_MANIFEST.csv"
        with figure_path.open(newline="", encoding="utf-8") as handle:
            figure = next(iter(csv.DictReader(handle)))
        write_csv(figure_path, [figure, figure])
    write_csv(panel_path, rows)
    with pytest.raises(ValueError):
        review.build_candidate_review([candidate], output)
    assert output.read_text(encoding="utf-8") == "retain current review"


def test_duplicate_figure_ids_are_rejected(tmp_path: Path) -> None:
    candidate = make_candidate(tmp_path, 1, "a")
    with pytest.raises(ValueError, match="Duplicate figure IDs"):
        review.build_candidate_review([candidate, candidate], tmp_path / "index.html")


@pytest.mark.parametrize("target", ["FIGURE_MANIFEST.csv", "collage & review.pdf"])
def test_output_cannot_replace_source_input(tmp_path: Path, target: str) -> None:
    candidate = make_candidate(tmp_path, 1, "a")
    path = candidate / target
    original = path.read_bytes()
    with pytest.raises(ValueError, match="must not replace"):
        review.build_candidate_review([candidate], path)
    assert path.read_bytes() == original


def test_cli_accepts_repeated_candidates(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    candidates = [make_candidate(tmp_path, 1, "a"), make_candidate(tmp_path, 2, "a")]
    output = tmp_path / "index.html"
    assert (
        review.main(
            [
                "--candidate",
                str(candidates[0]),
                "--candidate",
                str(candidates[1]),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    assert capsys.readouterr().out.strip() == str(output)
    assert Page(output.read_text(encoding="utf-8")).sections == [
        "figure-1",
        "figure-1-a",
        "figure-2",
        "figure-2-a",
    ]


def test_optional_figure_title_is_used_without_inferring_from_basename(tmp_path: Path) -> None:
    candidate = make_candidate(tmp_path, 2, "a")
    manifest = candidate / "FIGURE_MANIFEST.csv"
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["title"] = "TE composition & reproductive strategy"
    rows[0]["basename"] = "obsolete_diversity_title"
    write_csv(manifest, rows)
    output = review.build_candidate_review([candidate], tmp_path / "index.html")
    source = output.read_text(encoding="utf-8")
    assert "Figure 2 — TE composition &amp; reproductive strategy — full collage" in source
    assert "obsolete_diversity_title" not in source


def test_relettered_panels_keep_original_note_keys(tmp_path: Path) -> None:
    candidate = make_candidate(tmp_path, 10, "acde")
    manifest = candidate / "panels/PANEL_MANIFEST.csv"
    with manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row, display_id in zip(rows, "abcd", strict=True):
        row["display_panel_id"] = display_id
    write_csv(manifest, rows)
    output = tmp_path / "index.html"
    source = review.build_candidate_review([candidate], output).read_text()
    page = Page(source)
    assert page.sections == ["figure-10", "figure-10-a", "figure-10-c", "figure-10-d", "figure-10-e"]
    assert [note["data-note-id"] for note in page.notes] == page.sections
    assert "Figure 10b — Panel c" in source
    assert "Figure 10c — Panel d" in source
    assert "Figure 10d — Panel e" in source
    assert "Figure 10e —" not in source


@pytest.mark.parametrize("display_id", ['<script>', "a"])
def test_invalid_or_duplicate_display_panel_ids_are_rejected(tmp_path: Path, display_id: str) -> None:
    candidate = make_candidate(tmp_path, 10, "ac")
    manifest = candidate / "panels/PANEL_MANIFEST.csv"
    with manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows[0]["display_panel_id"] = "a"
    rows[1]["display_panel_id"] = display_id
    write_csv(manifest, rows)
    with pytest.raises(ValueError, match="display_panel_id"):
        review.build_candidate_review([candidate], tmp_path / "index.html")

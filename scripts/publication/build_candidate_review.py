"""Assemble existing candidate exports into one continuous figure-review page."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote


@dataclass(frozen=True)
class Card:
    anchor: str
    title: str
    legend: str
    alt_text: str
    assets: dict[str, Path]


@dataclass(frozen=True)
class Candidate:
    figure_id: str
    title: str
    cards: list[Card]
    manifests: tuple[Path, Path]


def read_manifest(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"{path}: duplicate column names")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: manifest has no rows")
    for index, row in enumerate(rows, start=2):
        if None in row or any(not (row.get(field) or "").strip() for field in required):
            raise ValueError(f"{path}: incomplete or malformed row {index}")
    return rows


def read_candidate(directory: Path) -> Candidate:
    directory = directory.resolve()
    figure_manifest = directory / "FIGURE_MANIFEST.csv"
    panel_manifest = directory / "panels/PANEL_MANIFEST.csv"
    shared_fields = {"legend", "alt_text", "pdf", "png", "tiff"}
    figures = read_manifest(figure_manifest, shared_fields | {"figure_id"})
    if len(figures) != 1:
        raise ValueError(f"{figure_manifest}: expected exactly one figure row")
    figure = figures[0]
    figure_id = figure["figure_id"]
    if not re.fullmatch(r"[A-Za-z0-9]+(?:[_-][A-Za-z0-9]+)*", figure_id):
        raise ValueError(f"{figure_manifest}: invalid figure_id {figure_id!r}")
    anchor = figure_id.lower().replace("_", "-")
    label = figure_id.replace("_", " ")
    title = (figure.get("title") or "").strip()
    group_title = f"{label} — {title}" if title else label

    def card(row: dict[str, str], base: Path, card_anchor: str, title: str) -> Card:
        assets = {}
        for kind in ("pdf", "png", "tiff"):
            relative = Path(row[kind])
            if relative.is_absolute():
                raise ValueError(f"{base}: asset must be a relative path: {relative}")
            asset = (base / relative).resolve()
            if not asset.is_file():
                raise ValueError(f"Missing {kind} asset: {asset}")
            assets[kind] = asset
        return Card(card_anchor, title, row["legend"], row["alt_text"], assets)

    cards = [card(figure, directory, anchor, f"{group_title} — full collage")]
    panels = read_manifest(panel_manifest, shared_fields | {"parent_figure", "panel_id", "title"})
    panel_ids: set[str] = set()
    display_ids: set[str] = set()
    for panel in panels:
        if panel["parent_figure"] != figure_id:
            raise ValueError(f"{panel_manifest}: parent_figure must match {figure_id}")
        panel_id = panel["panel_id"]
        if not re.fullmatch(r"[a-z][a-z0-9]*", panel_id):
            raise ValueError(f"{panel_manifest}: invalid panel_id {panel_id!r}")
        if panel_id in panel_ids:
            raise ValueError(f"{panel_manifest}: duplicate panel_id {panel_id!r}")
        panel_ids.add(panel_id)
        display_id = (panel.get("display_panel_id") or panel_id).strip()
        if not re.fullmatch(r"[a-z][a-z0-9]*", display_id):
            raise ValueError(f"{panel_manifest}: invalid display_panel_id {display_id!r}")
        if display_id in display_ids:
            raise ValueError(f"{panel_manifest}: duplicate display_panel_id {display_id!r}")
        display_ids.add(display_id)
        cards.append(
            card(
                panel,
                panel_manifest.parent,
                f"{anchor}-{panel_id}",
                f"{label}{display_id} — {panel['title']}",
            )
        )
    return Candidate(figure_id, group_title, cards, (figure_manifest, panel_manifest))


def build_candidate_review(candidate_dirs: Sequence[Path], output: Path) -> Path:
    """Write a review HTML; validate every source first and leave assets untouched."""
    if not candidate_dirs:
        raise ValueError("At least one candidate directory is required")
    candidates = [read_candidate(Path(directory)) for directory in candidate_dirs]
    anchors = [card.anchor for candidate in candidates for card in candidate.cards]
    if len(anchors) != len(set(anchors)):
        raise ValueError("Duplicate figure IDs or colliding figure/panel anchors")
    output = Path(output).resolve()
    inputs = {
        path.resolve()
        for candidate in candidates
        for path in (
            *candidate.manifests,
            *(asset for card in candidate.cards for asset in card.assets.values()),
        )
    }
    if output in inputs or (output.exists() and any(output.samefile(path) for path in inputs)):
        raise ValueError("Output must not replace a source asset or manifest")

    def asset_url(path: Path) -> str:
        relative = Path(os.path.relpath(path, output.parent)).as_posix()
        with path.open("rb") as handle:
            revision = hashlib.file_digest(handle, "sha256").hexdigest()[:12]
        return html.escape(f"{quote(relative, safe='/')}?v={revision}", quote=True)

    navigation = []
    sections = []
    for candidate in candidates:
        navigation.append(
            '<div class="nav-group"><strong>'
            + html.escape(candidate.title)
            + "</strong> "
            + " ".join(
                f'<a href="#{card.anchor}">{html.escape(card.title)}</a>'
                for card in candidate.cards
            )
            + "</div>"
        )
        for card in candidate.cards:
            urls = {kind: asset_url(path) for kind, path in card.assets.items()}
            note_label = html.escape(f"Notes for {card.title.split(' — ')[0]}")
            sections.append(
                f'<section id="{card.anchor}" aria-labelledby="{card.anchor}-title">\n'
                f'<h2 id="{card.anchor}-title">{html.escape(card.title)}</h2>\n'
                '<p class="downloads">'
                f'<a href="{urls["pdf"]}">PDF</a> · '
                f'<a href="{urls["png"]}">Full-size PNG</a> · '
                f'<a href="{urls["tiff"]}">CMYK TIFF</a> '
                '<button type="button" data-note-action="focus">Add / edit note</button></p>\n'
                f'<a href="{urls["png"]}"><img src="{urls["png"]}" '
                f'alt="{html.escape(card.alt_text, quote=True)}" loading="lazy"></a>\n'
                '<div class="figure-notes">\n'
                f'<label for="note-{card.anchor}">{note_label}</label>\n'
                f'<textarea id="note-{card.anchor}" class="note-input" '
                f'data-note-id="{card.anchor}" rows="4" maxlength="50000" '
                f'aria-describedby="note-status-{card.anchor}" disabled '
                'placeholder="What should change? Include panel letters or labels where useful.">'
                '</textarea>\n'
                f'<p class="note-status" id="note-status-{card.anchor}" role="status">'
                'Loading saved notes…</p>\n'
                '<div class="note-conflict" hidden><p>A newer note was saved elsewhere. '
                'Your draft is still in the box above. The saved version is:</p>'
                '<pre class="remote-note"></pre><div class="note-actions">'
                '<button type="button" data-note-action="use-saved">Use saved version</button> '
                '<button type="button" data-note-action="keep-draft">Save my draft instead</button>'
                '</div></div>\n</div>\n'
                f'<p class="caption">{html.escape(card.legend)}</p>\n</section>'
            )
    panel_count = sum(len(candidate.cards) - 1 for candidate in candidates)
    draft_summary = (
        f"<p><strong>{len(candidates)} draft figures · {panel_count} isolated panels.</strong> "
        "Review numbers are provisional; final selection and refinement follow this draft pass.</p>\n"
    )
    notes_script = Path(__file__).with_name("candidate_review_notes.js").read_text(encoding="utf-8")
    document = (
        """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Figure review — collages and isolated panels</title>
<style>
body{margin:0;background:#f3f4f5;color:#222;font:16px/1.5 Arial,sans-serif}
header,main{max-width:1100px;margin:auto;padding:24px}
h1{font-size:26px;margin:0 0 8px}a{color:#075b8b}
.nav-group{display:flex;flex-wrap:wrap;gap:6px 20px;margin:12px 0}
.nav-group strong{flex-basis:100%}
section{background:white;padding:24px;margin:0 0 28px;scroll-margin-top:88px}
h2{font-size:21px;margin:0}
img{display:block;width:100%;height:auto;margin:16px auto;max-width:1000px}
.caption{max-width:950px;white-space:pre-line}.downloads{margin:6px 0}
a:focus-visible{outline:3px solid #0072b2}
button,textarea{font:inherit}
button{border:1px solid #9aadb9;border-radius:5px;background:white;color:#075b8b;
padding:7px 12px;cursor:pointer;font-size:14px}
button:hover{background:#eaf2f7}button:disabled{opacity:.55;cursor:default}
button:focus-visible,textarea:focus-visible{outline:3px solid #0072b2;outline-offset:2px}
.downloads button{margin-left:12px}
.review-toolbar{position:sticky;top:0;z-index:2;background:#fff;border-block:1px solid #d6dfe5}
.review-toolbar-inner{max-width:1100px;margin:auto;padding:10px 24px;display:flex;
align-items:center;gap:10px;flex-wrap:wrap}
#notes-summary{flex:1;min-width:180px;font-size:14px}
.figure-notes{background:#f4f8fa;border:1px solid #d6e3eb;border-radius:6px;padding:16px;margin:18px 0}
.figure-notes label{display:block;font-weight:bold;margin-bottom:8px}
.note-input{display:block;box-sizing:border-box;width:100%;padding:12px;border:1px solid #92a9b8;
border-radius:4px;background:white;color:#222;resize:vertical;line-height:1.5;scroll-margin-top:110px}
.note-status{font-size:13px;margin:7px 0 0;color:#485d6c}
.note-conflict{border-top:1px solid #b88228;margin-top:12px;padding-top:8px}
.remote-note{white-space:pre-wrap;overflow-wrap:anywhere;font:inherit;background:white;padding:12px}
[hidden]{display:none!important}
@media(max-width:600px){header,main,section{padding:14px}}
@media print{body{background:white}nav,.downloads,.review-toolbar,.figure-notes{display:none}
section{break-before:page}}
</style></head><body><header>
<h1>Figure review: collages and isolated panels</h1>
"""
        + draft_summary
        + """
<p>Each full collage is followed by its individual panels. Scroll to inspect every figure,
or select an image for its full-size PNG.</p>
<p>Use <strong>Add / edit note</strong> on any figure or panel. Notes save automatically
to this project and remain available after refreshing or rebuilding the gallery.</p>
<noscript><p>Enable JavaScript to write and save figure notes.</p></noscript>
<nav aria-label="Figures and panels">
"""
        + "\n".join(navigation)
        + "\n</nav></header>\n"
        + '<div class="review-toolbar" aria-label="Figure notes controls">'
        '<div class="review-toolbar-inner">'
        '<span id="notes-summary" role="status">Loading project notes…</span>'
        '<button type="button" id="save-notes" disabled>Save now</button>'
        '<button type="button" id="export-notes" disabled>Download notes</button>'
        '</div></div><main>\n'
        + "\n".join(sections)
        + """
</main><script>
"""
        + notes_script
        + "\n</script></body></html>\n"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    return output


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        output = build_candidate_review(args.candidate, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

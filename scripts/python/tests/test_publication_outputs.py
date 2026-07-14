from __future__ import annotations

import csv
import re
import subprocess
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
FIGURES = ROOT / "Publication" / "figures"
METHODS = ROOT / "Publication" / "methods"


class PublicationOutputTests(unittest.TestCase):
    def test_figure_manifest_and_derivatives(self) -> None:
        with (FIGURES / "FIGURE_MANIFEST.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 7)
        self.assertEqual({row["figure_id"] for row in rows}, {f"Figure_{i}" for i in range(1, 8)})
        for row in rows:
            self.assertGreaterEqual(float(row["minimum_target_text_pt"]), 7)
            for column in ("pdf", "png", "tiff"):
                path = FIGURES / row[column]
                self.assertTrue(path.exists(), path)
                self.assertGreater(path.stat().st_size, 10_000, path)

    def test_pdf_fonts_are_embedded_and_tiffs_are_cmyk_300dpi(self) -> None:
        for pdf in FIGURES.glob("Figure_*.pdf"):
            output = subprocess.run(
                ["pdffonts", str(pdf)], check=True, capture_output=True, text=True
            ).stdout.splitlines()[2:]
            self.assertTrue(output, pdf)
            for line in output:
                fields = line.split()
                self.assertGreaterEqual(len(fields), 7, line)
                self.assertEqual(fields[5], "yes", line)
        for tiff in FIGURES.glob("Figure_*.tif"):
            output = subprocess.run(
                [
                    "magick",
                    "identify",
                    "-format",
                    "%[colorspace] %x %y",
                    str(tiff),
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
            self.assertRegex(output, r"^CMYK 300(?:\.0+)? 300(?:\.0+)?$")

    def test_methods_docx_contains_figures_line_numbers_and_required_sections(self) -> None:
        path = METHODS / "Materials_and_Methods_GBE_draft.docx"
        self.assertTrue(path.exists())
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            media = [name for name in names if name.startswith("word/media/")]
            self.assertEqual(len(media), 7)
            document_xml = archive.read("word/document.xml").decode("utf-8")
            self.assertIn('w:lnNumType w:countBy="1"', document_xml)
            self.assertIn('w:restart="continuous"', document_xml)
        plain = subprocess.run(
            ["pandoc", str(path), "-t", "plain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        for heading in (
            "Genomic data acquisition and processing",
            "De novo TE library construction and dnaPipeTE quantification",
            "Characterization of the repeat landscape",
            "Cell and nucleus segmentation, linkage, and quality review",
            "Nuclear integrated optical density and relative DNA-content sensitivity",
            "Exploratory phylogenetic path analysis",
            "Data availability",
        ):
            self.assertIn(heading, plain)
        self.assertGreaterEqual(len(re.findall("AUTHOR INPUT REQUIRED", plain)), 10)
        self.assertIn("not an independently validated absolute genome size", plain)


if __name__ == "__main__":
    unittest.main()

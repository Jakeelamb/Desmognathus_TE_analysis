"""Package the eight main collages and their captions after export validation."""
import csv
import hashlib
import json
import zipfile
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
verification = json.loads((RUN / 'qa/verification.json').read_text())
assert verification['status'] == 'PASS'
files = [RUN / name for name in ('FINDINGS_AND_CAPTIONS.md', 'README.md')]
for number in range(1, 9):
    folder = RUN / f'figures/figure_{number}'
    row = next(csv.DictReader((folder / 'FIGURE_MANIFEST.csv').open()))
    files.extend(folder / row[format_] for format_ in ('pdf', 'png', 'tiff'))
manifest = {
    'description': 'Eight polished review collages; frozen numerical results; separate panels remain in the gallery.',
    'formats': {'pdf': 'Vector; embedded Arial', 'png': '300 dpi RGB', 'tif': '600 dpi CMYK'},
    'files': [{'path': str(path.relative_to(RUN)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in files],
}
package = RUN / 'Manuscript_figures_1-8.zip'
with zipfile.ZipFile(package, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    for path in files:
        archive.write(path, str(path.relative_to(RUN)))
    archive.writestr('PACKAGE_MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
with zipfile.ZipFile(package) as archive:
    assert archive.testzip() is None
print(f'Packaged {len(files)} files plus manifest: {package} ({package.stat().st_size:,} bytes)')

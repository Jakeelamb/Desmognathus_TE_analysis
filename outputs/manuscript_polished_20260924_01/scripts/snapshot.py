"""Freeze existing numerical exports and record the previous review package."""
import hashlib
import json
import shutil
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
ROOT = RUN.parents[1]
SOURCE = ROOT / 'outputs/manuscript_story_20260924_01'
dest = RUN / 'inputs/source_run'
records = []
for path in sorted(SOURCE.rglob('*')):
    rel = path.relative_to(SOURCE)
    if not path.is_file() or any(p in {'node_modules', 'revisions', '.git'} for p in rel.parts):
        continue
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    copied = path.suffix.lower() in {'.csv', '.json', '.md', '.r', '.py', '.js', '.mjs', '.tre', '.nwk', '.txt', '.tsv'} and rel.parts[0] != 'qa'
    records.append({'source': str(path), 'relative': str(rel), 'sha256': digest, 'copied': copied})
    if copied:
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
notes = ROOT / 'Publication/figures/review_notes.json'
shutil.copy2(notes, RUN / 'inputs/earlier_review_notes.json')
(RUN / 'qa/source_baseline.json').write_text(json.dumps(records, indent=2) + '\n')
print(f'Preservation baseline: {len(records)} files; frozen text inputs: {sum(r["copied"] for r in records)}')

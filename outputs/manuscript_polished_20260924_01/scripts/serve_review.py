"""Serve only the polished run and its separate durable notes."""
import argparse
import importlib.util
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
path = RUN / 'inputs/source_run/inputs/publication_scripts/serve_candidate_review.py'
spec = importlib.util.spec_from_file_location('polished_server', path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
parser = argparse.ArgumentParser()
parser.add_argument('--port', type=int, default=8768)
args = parser.parse_args()
with module.ReviewServer(RUN, RUN / 'review_notes.json', args.port) as server:
    print(f'Polished eight-figure review: http://127.0.0.1:{server.server_address[1]}/', flush=True)
    server.serve_forever()

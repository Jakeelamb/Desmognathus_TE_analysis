"""Serve this isolated run with its own durable review-note store."""
import argparse
import importlib.util
import re
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("isolated_review_server", RUN / "inputs/publication_scripts/serve_candidate_review.py")
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
module.ANCHOR = re.compile(r"(?:figure-[1-9][0-9]*(?:-[a-z][a-z0-9]*)?|table-[1-9][0-9]*)\Z")
parser = argparse.ArgumentParser()
parser.add_argument('--port',type=int,default=8767)
args=parser.parse_args()
with module.ReviewServer(RUN, RUN / 'review_notes.json', args.port) as server:
    print(f'Isolated manuscript review: http://127.0.0.1:{server.server_address[1]}/',flush=True)
    server.serve_forever()

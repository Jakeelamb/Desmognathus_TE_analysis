from __future__ import annotations

import http.client
import importlib.util
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "serve_candidate_review", ROOT / "scripts/publication/serve_candidate_review.py"
)
assert SPEC is not None and SPEC.loader is not None
review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(review)


@contextmanager
def running_server(tmp_path):
    gallery = tmp_path / "candidates"
    gallery.mkdir(exist_ok=True)
    (gallery / "index.html").write_text("<h1>Figures</h1>", encoding="utf-8")
    (gallery / "figure.png").write_bytes(b"figure bytes")
    server = review.ReviewServer(gallery, tmp_path / "review_notes.json", port=0)
    thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01})
    thread.start()
    try:
        yield server
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=3)
        assert not thread.is_alive()


def request(server, method="GET", path="/api/notes", value=None, headers=None, body=None):
    connection = http.client.HTTPConnection(*server.server_address, timeout=5)
    actual_headers = {"Content-Type": "application/json", **(headers or {})}
    if value is not None:
        body = json.dumps(value).encode()
    try:
        connection.request(method, path, body=body, headers=actual_headers)
        response = connection.getresponse()
        data = response.read()
        if response.getheader("Content-Type", "").startswith("application/json"):
            data = json.loads(data)
        return response.status, data, dict(response.getheaders())
    finally:
        connection.close()


def note(text="Make species labels larger.", revision=0):
    return {"text": text, "title": "Figure 5 · Phylogeny", "base_revision": revision}


def test_notes_roundtrip_and_restart_persistence(tmp_path):
    with running_server(tmp_path) as server:
        assert request(server)[1] == {"schema_version": 1, "notes": {}}
        assert not server.store.path.exists()
        status, saved, _ = request(server, "PUT", "/api/notes/figure-5", note("Arial — α"))
        assert status == 200
        entry = saved["entry"]
        assert entry["text"] == "Arial — α"
        assert entry["revision"] == 1
        assert entry["updated_at"].endswith("+00:00")
        assert request(server)[1]["notes"] == {"figure-5": entry}
        assert json.loads(server.store.path.read_text())["notes"]["figure-5"] == entry
    with running_server(tmp_path) as server:
        assert request(server)[1]["notes"]["figure-5"] == entry
        status, saved, _ = request(server, "PUT", "/api/notes/figure-5", note("Next", 1))
        assert status == 200
        assert saved["entry"]["revision"] == 2


def test_concurrent_writes_conflict_per_note_and_empty_note_retains_revision(tmp_path):
    with running_server(tmp_path) as server:
        with ThreadPoolExecutor(max_workers=2) as executor:
            responses = list(
                executor.map(
                    lambda text: request(server, "PUT", "/api/notes/figure-5", note(text)),
                    ["First review", "Second review"],
                )
            )
        assert sorted(result[0] for result in responses) == [200, 409]
        winner = next(result[1]["entry"] for result in responses if result[0] == 200)
        conflict = next(result[1]["entry"] for result in responses if result[0] == 409)
        assert conflict == winner
        assert request(server, "PUT", "/api/notes/figure-5-a", note("Panel"))[0] == 200
        status, cleared, _ = request(server, "PUT", "/api/notes/figure-5", note("", 1))
        assert status == 200
        assert cleared["entry"]["revision"] == 2
        status, stale, _ = request(server, "PUT", "/api/notes/figure-5", note("Old", 1))
        assert status == 409
        assert stale["entry"]["text"] == ""
        assert set(request(server)[1]["notes"]) == {"figure-5", "figure-5-a"}


@pytest.mark.parametrize(
    "contents",
    [
        b"broken json",
        b'{"schema_version":1,"notes":[]}',
        b'{"schema_version":true,"notes":{}}',
        b'{"schema_version":1,"notes":{"figure-5":{"text":"old text"}}}',
    ],
)
def test_malformed_store_is_preserved_on_read_and_write(tmp_path, contents):
    with running_server(tmp_path) as server:
        server.store.path.write_bytes(contents)
        assert request(server)[0] == 500
        assert request(server, "PUT", "/api/notes/figure-5", note())[0] == 500
        assert server.store.path.read_bytes() == contents


def test_failed_atomic_replace_preserves_previous_notes(tmp_path, monkeypatch):
    with running_server(tmp_path) as server:
        assert request(server, "PUT", "/api/notes/figure-5", note())[0] == 200
        previous = server.store.path.read_bytes()

        def fail_replace(*args):
            raise OSError("disk failure")

        monkeypatch.setattr(review.os, "replace", fail_replace)
        status, value, _ = request(server, "PUT", "/api/notes/figure-5", note("New", 1))
        assert status == 500
        assert "could not be saved" in value["error"]
        assert server.store.path.read_bytes() == previous
        assert not list(tmp_path.glob(".review_notes.json.*.tmp"))


@pytest.mark.parametrize(
    "changes",
    [
        {"base_revision": True},
        {"base_revision": -1},
        {"base_revision": 0.0},
        {"base_revision": "0"},
        {"text": None},
        {"text": ["text"]},
        {"text": "a" * 50_001},
        {"title": "a" * 513},
        {"title": 5},
        {"extra": True},
    ],
)
def test_invalid_note_types_and_limits_do_not_create_store(tmp_path, changes):
    with running_server(tmp_path) as server:
        result = request(server, "PUT", "/api/notes/figure-5", {**note(), **changes})
        assert result[0] == 400
        assert not server.store.path.exists()


@pytest.mark.parametrize(
    "path",
    [
        "/api/notes/figure-0",
        "/api/notes/figure-05",
        "/api/notes/figure-1/A",
        "/api/notes/figure-1-../a",
        "/api/notes/figure-1-%61",
        "/api/notes/unrelated",
    ],
)
def test_invalid_anchors_are_rejected(tmp_path, path):
    with running_server(tmp_path) as server:
        assert request(server, "PUT", path, note())[0] == 400
        assert not server.store.path.exists()


def test_request_body_validation(tmp_path):
    with running_server(tmp_path) as server:
        path = "/api/notes/figure-5"
        assert request(server, "PUT", path, body=b"{bad")[0] == 400
        assert request(server, "PUT", path, body=b"\xff")[0] == 400
        assert request(server, "PUT", path, value=[])[0] == 400
        assert request(server, "PUT", path, note(), {"Content-Type": "text/plain"})[0] == 415
        assert request(server, "PUT", path, body=b"a" * (review.MAX_BODY + 1))[0] == 413
        assert not server.store.path.exists()


def test_cross_origin_and_bad_host_rejected_same_origin_accepted(tmp_path):
    with running_server(tmp_path) as server:
        port = server.server_address[1]
        path = "/api/notes/figure-5"
        for headers in [
            {"Origin": "https://unrelated.example"},
            {"Origin": "null"},
            {"Host": f"evil.example:{port}"},
            {"Sec-Fetch-Site": "cross-site"},
        ]:
            assert request(server, "PUT", path, note(), headers)[0] == 403
        assert request(server, headers={"Host": "evil.example"})[0] == 403
        assert not server.store.path.exists()
        headers = {"Origin": f"http://127.0.0.1:{port}"}
        assert request(server, "PUT", path, note(), headers)[0] == 200


def test_static_files_and_notes_outside_gallery(tmp_path):
    with running_server(tmp_path) as server:
        status, body, headers = request(server, path="/")
        assert status == 200
        assert body == b"<h1>Figures</h1>"
        assert headers["Cache-Control"] == "no-store"
        assert request(server, path="/index.html?fresh=1")[2]["Cache-Control"] == "no-store"
        assert request(server, path="/figure.png")[1] == b"figure bytes"
        assert request(server, "PUT", "/api/notes/figure-5", note())[0] == 200
        assert request(server, path="/../review_notes.json")[0] == 404
        assert request(server, path="/api/missing")[0] == 404
        assert request(server, "HEAD", "/")[0] == 200

"""Serve the local figure gallery and revision-checked, durable review notes."""

from __future__ import annotations

import argparse
import json
import os
import re
import tempfile
import threading
from datetime import UTC, datetime
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[2]
ANCHOR = re.compile(r"figure-[1-9][0-9]*(?:-[a-z][a-z0-9]*)?\Z")
MAX_BODY = 250 * 1024
MAX_TEXT = 50_000
MAX_TITLE = 512


class StoreError(Exception):
    """The existing notes cannot safely be read or replaced."""


def valid_entry(entry: object) -> bool:
    return (
        isinstance(entry, dict)
        and set(entry) == {"text", "title", "revision", "updated_at"}
        and isinstance(entry["text"], str)
        and len(entry["text"]) <= MAX_TEXT
        and isinstance(entry["title"], str)
        and len(entry["title"]) <= MAX_TITLE
        and type(entry["revision"]) is int
        and entry["revision"] > 0
        and isinstance(entry["updated_at"], str)
        and bool(entry["updated_at"])
    )


class NotesStore:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.lock = threading.Lock()

    def read(self) -> dict:
        try:
            value = json.loads(self.path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            return {"schema_version": 1, "notes": {}}
        except (OSError, UnicodeError, ValueError) as exc:
            raise StoreError("Notes could not be read; the existing file was preserved.") from exc
        if (
            not isinstance(value, dict)
            or set(value) != {"schema_version", "notes"}
            or type(value["schema_version"]) is not int
            or value["schema_version"] != 1
            or not isinstance(value["notes"], dict)
            or any(
                not ANCHOR.fullmatch(key) or not valid_entry(entry)
                for key, entry in value["notes"].items()
            )
        ):
            raise StoreError("Notes have an invalid format; the existing file was preserved.")
        return value

    def write(self, value: dict) -> None:
        temporary = None
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=self.path.parent,
                delete=False,
                prefix=f".{self.path.name}.",
                suffix=".tmp",
            ) as handle:
                temporary = Path(handle.name)
                json.dump(value, handle, ensure_ascii=False, indent=2)
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, self.path)
        except (OSError, UnicodeError) as exc:
            raise StoreError("Notes could not be saved; please keep your unsaved text.") from exc
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


class ReviewHandler(SimpleHTTPRequestHandler):
    server: ReviewServer

    def json_response(self, status: int, value: dict) -> None:
        body = json.dumps(value, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def allowed_request(self, *, mutate: bool = False) -> bool:
        hosts = self.headers.get_all("Host", [])
        port = self.server.server_address[1]
        if len(hosts) != 1 or hosts[0] not in {f"127.0.0.1:{port}", f"localhost:{port}"}:
            self.json_response(403, {"error": "Only the local review host is allowed."})
            return False
        origins = self.headers.get_all("Origin", [])
        if mutate and (
            (origins and origins != [f"http://{hosts[0]}"])
            or self.headers.get("Sec-Fetch-Site") == "cross-site"
        ):
            self.json_response(403, {"error": "Cross-origin note changes are not allowed."})
            return False
        return True

    def end_headers(self) -> None:
        if self.path == "/" or urlsplit(self.path).path.endswith(".html"):
            self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        super().end_headers()

    def do_HEAD(self) -> None:
        if self.allowed_request():
            super().do_HEAD()

    def do_GET(self) -> None:
        if not self.allowed_request():
            return
        path = urlsplit(self.path).path
        if path == "/api/notes":
            try:
                with self.server.store.lock:
                    value = self.server.store.read()
                self.json_response(200, value)
            except StoreError as exc:
                self.json_response(500, {"error": str(exc)})
        elif path.startswith("/api/"):
            self.json_response(404, {"error": "Unknown review API route."})
        else:
            super().do_GET()

    def do_PUT(self) -> None:
        if not self.allowed_request(mutate=True):
            return
        path = urlsplit(self.path).path
        anchor = path.removeprefix("/api/notes/")
        if not path.startswith("/api/notes/") or not ANCHOR.fullmatch(anchor):
            self.json_response(400, {"error": "Invalid figure or panel identifier."})
            return
        lengths = self.headers.get_all("Content-Length", [])
        if self.headers.get("Transfer-Encoding") or len(lengths) != 1:
            self.json_response(400, {"error": "A single Content-Length is required."})
            return
        if not lengths[0].isascii() or not lengths[0].isdigit():
            self.json_response(400, {"error": "Invalid Content-Length."})
            return
        length = int(lengths[0])
        if length > MAX_BODY:
            self.json_response(413, {"error": "Note request is too large."})
            return
        if self.headers.get_content_type() != "application/json":
            self.json_response(415, {"error": "Notes must use application/json."})
            return
        try:
            value = json.loads(self.rfile.read(length).decode("utf-8"))
        except (ValueError, UnicodeError):
            self.json_response(400, {"error": "Notes must contain valid UTF-8 JSON."})
            return
        if (
            not isinstance(value, dict)
            or set(value) != {"text", "title", "base_revision"}
            or not isinstance(value["text"], str)
            or len(value["text"]) > MAX_TEXT
            or not isinstance(value["title"], str)
            or len(value["title"]) > MAX_TITLE
            or type(value["base_revision"]) is not int
            or value["base_revision"] < 0
        ):
            self.json_response(400, {"error": "Invalid note text, title, or base revision."})
            return
        try:
            with self.server.store.lock:
                store = self.server.store.read()
                current = store["notes"].get(anchor)
                revision = current["revision"] if current else 0
                if value["base_revision"] != revision:
                    self.json_response(
                        409,
                        {
                            "error": "This note changed elsewhere. Review the saved version first.",
                            "entry": current,
                        },
                    )
                    return
                entry = {
                    "text": value["text"],
                    "title": value["title"],
                    "revision": revision + 1,
                    "updated_at": datetime.now(UTC).isoformat(),
                }
                store["notes"][anchor] = entry
                self.server.store.write(store)
            self.json_response(200, {"entry": entry})
        except StoreError as exc:
            self.json_response(500, {"error": str(exc)})


class ReviewServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, directory: Path, notes_file: Path, port: int = 8766) -> None:
        self.store = NotesStore(notes_file)
        super().__init__(("127.0.0.1", port), partial(ReviewHandler, directory=str(directory)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=ROOT / "Publication/figures/candidates")
    parser.add_argument(
        "--notes-file", type=Path, default=ROOT / "Publication/figures/review_notes.json"
    )
    parser.add_argument("--port", type=int, default=8766)
    args = parser.parse_args()
    if not args.directory.is_dir():
        parser.error("The gallery directory does not exist.")
    with ReviewServer(args.directory.resolve(), args.notes_file.resolve(), args.port) as server:
        print(f"Figure review: http://127.0.0.1:{server.server_address[1]}/", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()

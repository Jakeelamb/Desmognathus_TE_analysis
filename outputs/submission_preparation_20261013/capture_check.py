"""Retain a command, timestamps, status, and raw output before reporting it."""

from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

name, *command = sys.argv[1:]
if command and command[0] == "--":
    command.pop(0)
if not command or Path(name).name != name:
    raise SystemExit("usage: capture_check.py NAME -- COMMAND [ARGS]")
folder = Path(__file__).resolve().parent / "checks"
folder.mkdir(exist_ok=True)
log_path = folder / f"{name}.log"
receipt_path = folder / f"{name}.json"
if receipt_path.exists():
    raise SystemExit("Use a new check name to preserve previous execution evidence")
started_at = datetime.now(UTC).isoformat()
start = time.monotonic()
with log_path.open("wb") as log:
    process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
receipt = {
    "command": command,
    "cwd": str(Path.cwd()),
    "started_at": started_at,
    "completed_at": datetime.now(UTC).isoformat(),
    "elapsed_seconds": round(time.monotonic() - start, 3),
    "exit_code": process.returncode,
    "raw_log": str(log_path),
    "raw_log_sha256": hashlib.sha256(log_path.read_bytes()).hexdigest(),
}
receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
print(log_path.read_text(errors="replace")[-12000:])
print(json.dumps(receipt, indent=2))
raise SystemExit(process.returncode)

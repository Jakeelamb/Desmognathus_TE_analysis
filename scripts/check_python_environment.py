"""Read-only interpreter check; UV's sync --check verifies locked dependencies."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    expected_python = (ROOT / ".python-version").read_text().strip()
    observed_python = ".".join(map(str, sys.version_info[:3]))
    if observed_python != expected_python:
        print(
            f"Python environment: FAIL; expected {expected_python}, found {observed_python}. "
            "Use make setup-python in a new UV_PROJECT_ENVIRONMENT to preserve an older environment.",
            file=sys.stderr,
        )
        return 1

    print(f"Python interpreter: PASS; Python {observed_python}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

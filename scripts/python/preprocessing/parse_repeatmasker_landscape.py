#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compatibility wrapper for the canonical RepeatMasker landscape parser.
"""

from __future__ import annotations

import importlib.util
import logging
import sys
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("landscape_parser_wrapper")


def _load_canonical_module():
    script_path = (
        Path(__file__).resolve().parent.parent.parent
        / "processing"
        / "parse_repeatmasker_landscape.py"
    )
    spec = importlib.util.spec_from_file_location(
        "canonical_parse_repeatmasker_landscape",
        script_path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load canonical parser from {script_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    logger.warning(
        "scripts/python/preprocessing/parse_repeatmasker_landscape.py is a compatibility "
        "wrapper. Use scripts/processing/parse_repeatmasker_landscape.py."
    )
    module = _load_canonical_module()
    module.main()


if __name__ == "__main__":
    main()

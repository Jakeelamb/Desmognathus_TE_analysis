#!/usr/bin/env python3
"""
Compatibility wrapper for the canonical preprocessing superfamily writer.
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from python.preprocessing.generate_superfamily_proportions import (
    create_superfamily_proportions,
    main as preprocessing_main,
)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    """Run the canonical preprocessing implementation."""
    logger.warning(
        "scripts/python/analysis/generate_superfamily_proportions.py is a compatibility "
        "wrapper. Use scripts/python/preprocessing/generate_superfamily_proportions.py."
    )
    preprocessing_main()


if __name__ == "__main__":
    main()

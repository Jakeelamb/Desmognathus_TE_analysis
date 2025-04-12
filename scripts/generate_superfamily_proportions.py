#!/usr/bin/env python3
"""
Generate superfamily proportions from a superfamily breakdown CSV file.

This script is deprecated. Please use scripts/python/preprocessing/generate_superfamily_proportions.py instead.
"""

import os
import sys
import pandas as pd
import numpy as np
import logging
import warnings
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from python.utils.path_utils import resolve_path, ensure_directory

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """Main function that shows a deprecation warning and runs the new script."""
    logger.warning("This script is deprecated. Please use scripts/python/preprocessing/generate_superfamily_proportions.py instead.")
    
    # Import and run the new script
    try:
        from python.preprocessing.generate_superfamily_proportions import main as new_main
        logger.info("Redirecting to the new script...")
        new_main()
    except ImportError:
        logger.error("Could not import the new script. Please run scripts/python/preprocessing/generate_superfamily_proportions.py directly.")
        sys.exit(1)

if __name__ == "__main__":
    main() 
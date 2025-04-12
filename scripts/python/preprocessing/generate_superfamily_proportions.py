#!/usr/bin/env python3
"""
Generate superfamily proportions from a superfamily breakdown CSV file.

This script reads the superfamily breakdown CSV file and calculates the proportions
of each superfamily, then saves the results to a CSV file.
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import warnings
from pathlib import Path

# Add the parent directory to the path so we can import the path utils
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from python.utils.path_utils import resolve_path, ensure_directory

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def create_superfamily_proportions(input_file, output_file):
    """
    Create a superfamily proportions file from a superfamily breakdown CSV.
    Assumes input CSV has 'Species' as the first column (or index) and superfamilies as subsequent columns.
    
    Parameters:
    -----------
    input_file : str
        Path to the superfamily breakdown CSV file
    output_file : str
        Path to save the output CSV file
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    try:
        # Read the superfamily breakdown CSV, using the first column as index
        df = pd.read_csv(input_file, index_col=0)
        logger.info(f"Read superfamily breakdown from: {input_file}")

        # Identify numerical columns (superfamilies) - exclude 'Total' if present
        numerical_cols = df.select_dtypes(include=np.number).columns
        if 'Total' in numerical_cols:
             numerical_cols = numerical_cols.drop('Total')
        
        # Calculate the total bases for each species (row sums)
        # Use only numerical columns for sum, handle potential non-numeric Totals robustly
        row_totals = df[numerical_cols].sum(axis=1)

        # Create a copy for proportions to avoid modifying original df if needed elsewhere
        prop_df = df[numerical_cols].copy()

        # Calculate proportions row-wise
        # Use .div() for safe division, axis=0 aligns row_totals with rows of prop_df
        # Fill NaN results from 0/0 division with 0
        prop_df = prop_df.div(row_totals, axis=0).fillna(0)

        # Ensure proportions sum to 1 (or 100 if multiplying by 100 later) per species
        # Add verification step (optional but recommended)
        # logger.debug(f"Proportion sums per species:\\n{prop_df.sum(axis=1)}")

        # Save the proportions to a CSV file, including the index ('Species')
        ensure_directory(os.path.dirname(output_file)) # Use ensure_directory from utils
        prop_df.to_csv(output_file, index=True) # index=True saves the 'Species' column
        logger.info(f"Saved superfamily proportions to: {output_file}")
        
        # Free memory
        del df, prop_df
        
        return True
    
    except Exception as e:
        logger.error(f"Error creating superfamily proportions: {e}")
        return False


def parse_args():
    """Parse command-line arguments."""
    import argparse
    parser = argparse.ArgumentParser(description="Generate superfamily proportions from breakdown CSV")
    
    parser.add_argument("--input", "-i", type=str,
                        help="Path to superfamily breakdown CSV")
    parser.add_argument("--output", "-o", type=str,
                        help="Path to save the output superfamily proportions CSV")
    parser.add_argument("--diversity-copy", action="store_true",
                        help="Copy the output file to the diversity directory for PCA analysis")
    parser.add_argument("--config", "-c", type=str, default="config/paths.yaml",
                        help="Path to configuration file")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose logging")
    
    return parser.parse_args()


def main():
    """Run the main function."""
    args = parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Determine paths
    if args.input:
        input_path = args.input
    else:
        input_path = os.path.join(resolve_path('data.interim.pivot_tables'), 
                                  'superfamily_breakdown.csv')
    
    if args.output:
        output_path = args.output
    else:
        output_dir = resolve_path('data.processed.te_superfamily')
        output_path = os.path.join(output_dir, 'superfamily_proportions.csv')
    
    # Check if input file exists
    if not os.path.exists(input_path):
        logger.error(f"Input file does not exist: {input_path}")
        sys.exit(1)
    
    # Create superfamily proportions
    if create_superfamily_proportions(input_path, output_path):
        logger.info("Superfamily proportions generated successfully")
        
        # Copy to diversity directory if requested
        if args.diversity_copy:
            diversity_dir = resolve_path('data.processed.diversity')
            ensure_directory(diversity_dir)
            diversity_path = os.path.join(diversity_dir, 'superfamily_proportions.csv')
            
            try:
                import shutil
                shutil.copy2(output_path, diversity_path)
                logger.info(f"Copied superfamily proportions to: {diversity_path} for PCA analysis")
            except Exception as e:
                logger.error(f"Error copying to diversity directory: {e}")
    else:
        logger.error("Failed to generate superfamily proportions")
        sys.exit(1)


if __name__ == "__main__":
    main() 
#!/usr/bin/env python3
"""
Generate canonical species-by-superfamily proportions from a breakdown table.
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
from pathlib import Path

# Add the parent directory to the path so we can import the path utils
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from python.utils.path_utils import resolve_path, ensure_directory
from python.utils.te_table_utils import coerce_feature_table_to_species_matrix

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def create_superfamily_proportions(input_file, output_file):
    """
    Create a canonical species-by-superfamily proportions file.
    
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
        df = pd.read_csv(input_file)
        logger.info(f"Read superfamily breakdown from: {input_file}")
        matrix, orientation = coerce_feature_table_to_species_matrix(
            df,
            feature_axis_labels=("superfamily",),
        )
        logger.info(f"Detected input orientation: {orientation}")

        row_totals = matrix.sum(axis=1)
        prop_df = matrix.div(row_totals.replace(0, np.nan), axis=0).fillna(0)
        prop_df.index.name = "Species"

        ensure_directory(os.path.dirname(output_file))
        prop_df.to_csv(output_file, index=True)
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
    parser.add_argument("--config", "-c", type=str, default="paths.yaml",
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
        input_path = os.path.join(resolve_path('data.interim.pivot_tables', args.config), 
                                  'superfamily_breakdown.csv')
    
    if args.output:
        output_path = args.output
    else:
        output_dir = resolve_path('data.processed.te_superfamily', args.config)
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
            diversity_dir = resolve_path('data.processed.diversity', args.config)
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

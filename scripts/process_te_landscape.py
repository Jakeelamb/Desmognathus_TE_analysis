#!/usr/bin/env python3
"""
Main workflow script for TE landscape analysis.
"""

import sys
import argparse
import subprocess
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('TE_Landscape')

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process TE landscape data')
    parser.add_argument('--skip-proportions', action='store_true', 
                      help='Skip generation of superfamily proportions')
    parser.add_argument('--skip-diversity', action='store_true', 
                      help='Skip diversity analysis')
    parser.add_argument('--skip-pca', action='store_true', 
                      help='Skip PCA analysis')
    parser.add_argument('--test', action='store_true',
                      help='Use test data instead of full dataset')
    parser.add_argument('--min-species-presence', type=int, default=3,
                      help='Minimum number of species a TE must be present in (for PCA analysis)')
    
    return parser.parse_args()

def run_command(cmd):
    """Run a command and log its output."""
    logger.info(f"Running command: {' '.join(cmd)}")
    process = subprocess.run(cmd, capture_output=True, text=True)
    
    if process.returncode != 0:
        logger.error(f"Command failed with return code {process.returncode}")
        logger.error(f"STDERR: {process.stderr}")
        sys.exit(1)
    
    logger.info(f"Command completed successfully")
    return process.stdout

def main():
    """Run the TE landscape analysis workflow."""
    args = parse_args()
    
    # Define script paths
    scripts_dir = Path(__file__).parent
    proportions_script = scripts_dir / "python" / "analysis" / "generate_superfamily_proportions.py"
    diversity_script = scripts_dir / "python" / "analysis" / "run_diversity_analysis.py"
    pca_script = scripts_dir / "R" / "analysis" / "te_pca_analysis.R"
    
    # Run superfamily proportions generation
    if not args.skip_proportions:
        logger.info("Generating superfamily proportions...")
        cmd = ["python", str(proportions_script)]
        if args.test:
            cmd.extend(["--data-dir", "data/test_data/raw"])
        run_command(cmd)
    
    # Run diversity analysis
    if not args.skip_diversity:
        logger.info("Running diversity analysis...")
        cmd = ["python", str(diversity_script)]
        if args.test:
            cmd.append("--test")
        run_command(cmd)
    
    # Run PCA analysis
    if not args.skip_pca:
        logger.info("Running PCA analysis...")
        cmd = ["Rscript", str(pca_script)]
        if args.min_species_presence:
            cmd.extend(["--min-species-presence", str(args.min_species_presence)])
        if args.test:
            cmd.append("--test")
        run_command(cmd)
    
    logger.info("TE landscape processing completed!")

if __name__ == "__main__":
    main()



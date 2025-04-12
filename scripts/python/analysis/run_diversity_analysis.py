#!/usr/bin/env python3
"""
Run diversity analysis on TE superfamily proportions.

This script calculates diversity metrics (Shannon, Simpson, etc.) for TE superfamily
proportions data and generates visualizations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
from pathlib import Path

# Add the parent directory to the path so we can import the path utils
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from python.utils.path_utils import resolve_path, ensure_directory

import logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def calculate_diversity_metrics(df, output_dir=None):
    """
    Calculate diversity metrics for TE superfamily proportions.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing TE superfamily proportions data
    output_dir : str, optional
        Directory to save output files
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with diversity metrics
    """
    # Get species columns (all except the first column)
    species_cols = df.columns[1:]
    
    # Create a dictionary to store diversity metrics
    metrics = {
        'species': [],
        'shannon_entropy': [],
        'simpson_diversity': [],
        'inverse_simpson': [],
        'richness': [],
        'pielou_evenness': []
    }
    
    # Calculate diversity metrics for each species
    for col in species_cols:
        # Get the proportions for the current species
        props = df[col].values
        
        # Remove NaN values
        props = props[~np.isnan(props)]
        
        # Skip if all values are zero
        if np.sum(props) == 0:
            logger.warning(f"Skipping {col} - all proportions are zero")
            continue
        
        # Calculate Shannon entropy (H)
        # H = -sum(p_i * log(p_i))
        # Only consider non-zero proportions for log
        non_zero_props = props[props > 0]
        shannon = -np.sum(non_zero_props * np.log(non_zero_props))
        
        # Calculate Simpson diversity (D)
        # D = 1 - sum(p_i^2)
        simpson = 1 - np.sum(props ** 2)
        
        # Calculate inverse Simpson diversity
        # 1/D = 1 / sum(p_i^2)
        inv_simpson = 1 / np.sum(props ** 2) if np.sum(props ** 2) > 0 else np.nan
        
        # Calculate richness (number of non-zero proportions)
        richness = np.sum(props > 0)
        
        # Calculate Pielou's evenness (J)
        # J = H / log(S), where S is richness
        pielou = shannon / np.log(richness) if richness > 0 else np.nan
        
        # Add metrics to the dictionary
        metrics['species'].append(col)
        metrics['shannon_entropy'].append(shannon)
        metrics['simpson_diversity'].append(simpson)
        metrics['inverse_simpson'].append(inv_simpson)
        metrics['richness'].append(richness)
        metrics['pielou_evenness'].append(pielou)
    
    # Create a DataFrame from the metrics dictionary
    metrics_df = pd.DataFrame(metrics)
    
    # Save metrics to a CSV file if output_dir is provided
    if output_dir:
        ensure_directory(output_dir)
        metrics_df.to_csv(os.path.join(output_dir, 'diversity_metrics.csv'), index=False)
        logger.info(f"Saved diversity metrics to {os.path.join(output_dir, 'diversity_metrics.csv')}")
    
    return metrics_df


def plot_diversity_metrics(metrics_df, output_dir=None):
    """
    Generate plots for diversity metrics.
    
    Parameters:
    -----------
    metrics_df : pandas.DataFrame
        DataFrame with diversity metrics
    output_dir : str, optional
        Directory to save output files
        
    Returns:
    --------
    None
    """
    # Ensure the output directory exists
    if output_dir:
        ensure_directory(output_dir)
    
    # Set the style for plotting
    sns.set(style='whitegrid', context='paper', palette='colorblind')
    
    # Plot Shannon entropy
    plt.figure(figsize=(12, 6))
    sns.barplot(x='species', y='shannon_entropy', data=metrics_df)
    plt.title('Shannon Entropy by Species')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'shannon_entropy.png'), dpi=300)
        logger.info(f"Saved Shannon entropy plot to {os.path.join(output_dir, 'shannon_entropy.png')}")
    
    plt.close()
    
    # Plot Simpson diversity
    plt.figure(figsize=(12, 6))
    sns.barplot(x='species', y='simpson_diversity', data=metrics_df)
    plt.title('Simpson Diversity by Species')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'simpson_diversity.png'), dpi=300)
        logger.info(f"Saved Simpson diversity plot to {os.path.join(output_dir, 'simpson_diversity.png')}")
    
    plt.close()
    
    # Plot richness
    plt.figure(figsize=(12, 6))
    sns.barplot(x='species', y='richness', data=metrics_df)
    plt.title('TE Superfamily Richness by Species')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'richness.png'), dpi=300)
        logger.info(f"Saved richness plot to {os.path.join(output_dir, 'richness.png')}")
    
    plt.close()
    
    # Plot Pielou's evenness
    plt.figure(figsize=(12, 6))
    sns.barplot(x='species', y='pielou_evenness', data=metrics_df)
    plt.title('Pielou Evenness by Species')
    plt.xticks(rotation=90)
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'pielou_evenness.png'), dpi=300)
        logger.info(f"Saved Pielou evenness plot to {os.path.join(output_dir, 'pielou_evenness.png')}")
    
    plt.close()
    
    # Plot a heatmap of all metrics
    plt.figure(figsize=(12, 8))
    metrics_heatmap = metrics_df.set_index('species')
    sns.heatmap(metrics_heatmap, annot=True, cmap='YlGnBu', linewidths=0.5)
    plt.title('Diversity Metrics Heatmap')
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'diversity_heatmap.png'), dpi=300)
        logger.info(f"Saved diversity heatmap to {os.path.join(output_dir, 'diversity_heatmap.png')}")
    
    plt.close()


def parse_args():
    """Parse command-line arguments."""
    import argparse
    parser = argparse.ArgumentParser(description="Calculate diversity metrics for TE superfamily proportions")
    
    parser.add_argument("--input", "-i", type=str,
                        help="Path to superfamily proportions CSV file")
    parser.add_argument("--output-dir", "-o", type=str,
                        help="Directory to save output files")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip generating plots")
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
        input_path = os.path.join(resolve_path('data.processed.diversity'), 
                                  'superfamily_proportions.csv')
    
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = resolve_path('results.tables.diversity')
    
    figures_dir = resolve_path('results.figures.diversity')
    
    # Check if input file exists
    if not os.path.exists(input_path):
        logger.error(f"Input file does not exist: {input_path}")
        sys.exit(1)
    
    # Read the superfamily proportions CSV
    try:
        df = pd.read_csv(input_path)
        logger.info(f"Read superfamily proportions from: {input_path}")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Calculate diversity metrics
    metrics_df = calculate_diversity_metrics(df, output_dir)
    
    # Generate plots if not disabled
    if not args.no_plots:
        plot_diversity_metrics(metrics_df, figures_dir)
    
    logger.info("Diversity analysis completed successfully")


if __name__ == "__main__":
    main() 
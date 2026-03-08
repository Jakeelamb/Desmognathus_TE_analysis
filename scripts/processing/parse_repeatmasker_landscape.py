#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_repeatmasker_landscape.py

Description:
    Parses RepeatMasker .align files and joins with classification metadata
    to generate a dataset for TE repeat landscape analysis.

Usage:
    python parse_repeatmasker_landscape.py [SRX_ID]

    Example:
    python parse_repeatmasker_landscape.py SRX19953421
"""

import os
import re
import sys
import logging
import pandas as pd
import numpy as np
from pathlib import Path

# Import centralized configuration
sys.path.insert(0, str(Path(__file__).parent.parent))  # Add scripts/ to path
from config import paths, PROJECT_ROOT

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('parser.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger('landscape_parser')

# Define constants and paths from centralized config
ALIGN_DIR = paths.input_data.repeatmasker
CLASS_DIR = PROJECT_ROOT / "interim"
MERGED_CLASSIFICATION_FILE = paths.results.data / "dnaPipeTE_merged_classifications.csv"
OUTPUT_DIR = paths.results.landscapes
KIMURA_PATTERN = r"Kimura \(with divCpGMod\) = (\d+\.\d+)"

_MERGED_CLASSIFICATION_CACHE = None

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_align_file(filepath):
    """
    Parse a RepeatMasker .align file to extract contig names, aligned base pairs, and Kimura distances.
    
    Args:
        filepath (str): Path to the .align file
        
    Returns:
        list: List of tuples (contig_name, aligned_bp, kimura)
    """
    logger.info(f"Parsing .align file: {filepath}")
    
    results = []
    current_contig = None
    current_aligned_bp = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # If line starts with a digit, it's a header line for a TE hit
            if re.match(r'^\d+', line):
                parts = line.split()
                if len(parts) >= 8:  # Ensure we have enough columns
                    current_contig = parts[4]
                    # Calculate aligned_bp from start and end positions
                    try:
                        start = int(parts[5])
                        end = int(parts[6])
                        current_aligned_bp = abs(end - start)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Invalid alignment coordinates in line: {line}")
                        logger.warning(f"Error: {str(e)}")
                        logger.warning(f"Parts: {parts}")
                        current_contig = None
                        current_aligned_bp = None
            
            # Look for Kimura distance line
            elif current_contig and "Kimura (with divCpGMod) =" in line:
                kimura_match = re.search(KIMURA_PATTERN, line)
                if kimura_match:
                    kimura = float(kimura_match.group(1))
                    results.append((current_contig, current_aligned_bp, kimura))
                    current_contig = None
                    current_aligned_bp = None
    
    logger.info(f"Extracted {len(results)} entries from {filepath}")
    return results


def bin_kimura(kimura):
    """
    Bin Kimura distances into integer percentage ranges.
    
    Args:
        kimura (float): Kimura distance value
        
    Returns:
        str: String representation of the bin (e.g., "0–1%")
    """
    bin_start = int(kimura)
    return f"{bin_start}–{bin_start+1}%"


def load_classification_file(filepath):
    """
    Load the classification file containing TE annotations.
    
    Args:
        filepath (str): Path to the classification file
        
    Returns:
        pd.DataFrame: DataFrame with classification data
    """
    logger.info(f"Loading classification file: {filepath}")
    
    try:
        # Read the file line by line to handle variable number of fields
        data = []
        with open(filepath, 'r') as f:
            # Read header line
            header = f.readline().strip()
            if header.startswith('#'):
                header = header[1:]
            
            # Process each line
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # Split on whitespace while preserving quoted strings
                parts = []
                current = []
                in_quotes = False
                for char in line:
                    if char == '"':
                        in_quotes = not in_quotes
                    elif char.isspace() and not in_quotes:
                        if current:
                            parts.append(''.join(current))
                            current = []
                    else:
                        current.append(char)
                if current:
                    parts.append(''.join(current))
                
                # Extract fields
                if len(parts) >= 11:  # We expect at least 11 fields
                    entry = {
                        'dnaPipeTE_contig_name': parts[2],  # comp_TRINITY_* field
                        'Species': parts[7].replace('"', ''),  # D.species field
                        'Class': parts[8].replace('"', ''),  # Class field
                        'Order': parts[9].replace('"', ''),  # Order field
                        'Superfamily': parts[10].replace('"', '') if len(parts) > 10 else 'Unknown'  # Superfamily field
                    }
                    data.append(entry)
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Clean up the data
        df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        
        # Log some statistics
        logger.info(f"Loaded {len(df)} entries from classification file")
        logger.info(f"Unique Classes: {df['Class'].nunique()}")
        logger.info(f"Unique Orders: {df['Order'].nunique()}")
        logger.info(f"Unique Superfamilies: {df['Superfamily'].nunique()}")
        
        return df
        
    except Exception as e:
        logger.error(f"Error loading classification file: {str(e)}")
        logger.warning("Creating empty DataFrame with required columns")
        return pd.DataFrame(columns=['dnaPipeTE_contig_name', 'Species', 'Class', 'Order', 'Superfamily'])


def load_classification_for_sample(srx_id):
    """
    Load classification data for a single sample.

    Prefer the canonical merged classification table in results/data and
    fall back to the older per-sample interim file only if necessary.
    """
    global _MERGED_CLASSIFICATION_CACHE

    if MERGED_CLASSIFICATION_FILE.exists():
        if _MERGED_CLASSIFICATION_CACHE is None:
            logger.info(f"Loading canonical merged classifications: {MERGED_CLASSIFICATION_FILE}")
            _MERGED_CLASSIFICATION_CACHE = pd.read_csv(MERGED_CLASSIFICATION_FILE)

        merged_df = _MERGED_CLASSIFICATION_CACHE
        if 'Source' in merged_df.columns:
            sample_df = merged_df[merged_df['Source'].astype(str) == srx_id].copy()
            if not sample_df.empty:
                logger.info(f"Loaded {len(sample_df)} canonical classification rows for {srx_id}")
                return sample_df
            logger.warning(f"No rows for {srx_id} found in canonical merged classifications")

    legacy_file = CLASS_DIR / f"{srx_id}_reads_per_component_and_annotation_processed"
    if legacy_file.exists():
        logger.info(f"Falling back to legacy interim classification file: {legacy_file}")
        return load_classification_file(legacy_file)

    logger.error(
        "No classification source found for %s. Checked %s and %s",
        srx_id,
        MERGED_CLASSIFICATION_FILE,
        legacy_file,
    )
    return pd.DataFrame(columns=['dnaPipeTE_contig_name', 'Species', 'Class', 'Order', 'Superfamily'])


def merge_align_classification(align_data, class_df):
    """
    Merge .align data with classification metadata.
    
    Args:
        align_data (list): List of tuples (contig_name, aligned_bp, kimura)
        class_df (pd.DataFrame): DataFrame with classification data
        
    Returns:
        pd.DataFrame: Merged DataFrame
    """
    # Convert align data to DataFrame
    align_df = pd.DataFrame(align_data, columns=['contig_name', 'aligned_bp', 'kimura'])
    
    # Add Kimura bins
    align_df['Kimura_bin'] = align_df['kimura'].apply(bin_kimura)
    
    # Log sample of contig names from both sources
    logger.info("Sample RepeatMasker contig names:")
    logger.info(align_df['contig_name'].head().tolist())
    logger.info("\nSample dnaPipeTE contig names:")
    logger.info(class_df['dnaPipeTE_contig_name'].head().tolist())
    
    # Clean up contig names for matching
    align_df['clean_contig'] = align_df['contig_name'].apply(lambda x: 
        re.sub(r'comp_TRINITY_', '', x.split('_i')[0])  # Remove comp_TRINITY_ prefix and everything after _i
    )
    
    class_df['clean_contig'] = class_df['dnaPipeTE_contig_name'].apply(lambda x:
        re.sub(r'comp_TRINITY_', '', x.split('_i')[0])  # Same cleaning for dnaPipeTE names
    )
    
    # Try matching on cleaned names
    merged_df = pd.merge(
        align_df,
        class_df[['clean_contig', 'Species', 'Class', 'Order', 'Superfamily']],
        on='clean_contig',
        how='left'
    )
    
    # Check match rate
    match_rate = (1 - merged_df['Species'].isna().mean()) * 100
    logger.info(f"\nMatch rate after cleaning: {match_rate:.2f}%")
    
    # For unmatched entries, try to extract classification from RM_annotation
    if 'RM_annotation' in class_df.columns:
        annotation_map = {}
        for _, row in class_df.iterrows():
            if pd.notna(row['RM_annotation']):
                parts = row['RM_annotation'].split('/')
                if len(parts) >= 2:
                    order = parts[0]
                    superfamily = parts[1] if len(parts) > 1 else 'Unknown'
                    annotation_map[row['clean_contig']] = {
                        'Order': order,
                        'Superfamily': superfamily,
                        'Class': row['Class']
                    }
        
        # Apply mapping to unmatched rows
        unmatched = merged_df['Species'].isna()
        if unmatched.any():
            for idx in merged_df[unmatched].index:
                contig = merged_df.loc[idx, 'clean_contig']
                if contig in annotation_map:
                    merged_df.loc[idx, 'Order'] = annotation_map[contig]['Order']
                    merged_df.loc[idx, 'Superfamily'] = annotation_map[contig]['Superfamily']
                    merged_df.loc[idx, 'Class'] = annotation_map[contig]['Class']
    
    # Ensure we have all required columns with valid values
    if 'Species' not in merged_df.columns or merged_df['Species'].isna().all():
        # Use species from classification data if available
        if not class_df.empty and 'Species' in class_df.columns:
            species = class_df['Species'].iloc[0]
        else:
            species = "Desmognathus_sp"
        merged_df['Species'] = species
    
    # Fill in missing classifications
    merged_df['Order'] = merged_df['Order'].fillna('Unknown')
    merged_df['Superfamily'] = merged_df['Superfamily'].fillna('Unknown')
    merged_df['Class'] = merged_df['Class'].fillna('Unknown')
    
    # Keep only necessary columns and drop any duplicates
    result_df = merged_df[['Species', 'Kimura_bin', 'Class', 'Order', 'Superfamily', 'aligned_bp']].drop_duplicates()
    
    # Log final statistics
    logger.info("\nFinal dataset statistics:")
    logger.info(f"Total entries: {len(result_df)}")
    logger.info(f"Unique Classes: {result_df['Class'].nunique()}")
    logger.info(f"Unique Orders: {result_df['Order'].nunique()}")
    logger.info(f"Unique Superfamilies: {result_df['Superfamily'].nunique()}")
    logger.info("\nClass distribution:")
    logger.info(result_df['Class'].value_counts())
    logger.info("\nOrder distribution:")
    logger.info(result_df['Order'].value_counts())
    logger.info("\nSuperfamily distribution:")
    logger.info(result_df['Superfamily'].value_counts())
    
    return result_df


def aggregate_landscape(df):
    """
    Aggregate the landscape data by species, Kimura bin, class, order, and superfamily.
    
    Args:
        df (pd.DataFrame): DataFrame with merged align and classification data
        
    Returns:
        pd.DataFrame: Aggregated DataFrame
    """
    # Handle empty dataframe case
    if df.empty:
        logger.warning("Empty DataFrame passed to aggregate_landscape")
        return pd.DataFrame(columns=['Species', 'Kimura_bin', 'Class', 'Order', 'Superfamily', 'aligned_bp'])
    
    try:
        # Group by and sum aligned base pairs
        agg_df = df.groupby(['Species', 'Kimura_bin', 'Class', 'Order', 'Superfamily'])['aligned_bp'].sum().reset_index()
        
        # Remove quotes from any string values (some files have quoted values)
        for col in ['Species', 'Class', 'Order', 'Superfamily']:
            if df[col].dtype == 'object':
                agg_df[col] = agg_df[col].str.replace('"', '')
        
        return agg_df
    except Exception as e:
        logger.error(f"Error in aggregate_landscape: {str(e)}")
        return pd.DataFrame(columns=['Species', 'Kimura_bin', 'Class', 'Order', 'Superfamily', 'aligned_bp'])


def process_single_sample(srx_id):
    """
    Process a single SRX sample.
    
    Args:
        srx_id (str): SRX identifier of the sample
        
    Returns:
        pd.DataFrame: Processed landscape data for the sample
    """
    logger.info(f"Processing sample {srx_id}")
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Check if align file exists
    align_file = os.path.join(ALIGN_DIR, f"{srx_id}_Trinity.align")
    if not os.path.exists(align_file):
        logger.error(f"Align file not found: {align_file}")
        return None
    
    # Parse .align file
    align_data = parse_align_file(align_file)
    if not align_data:
        logger.error(f"No data extracted from {align_file}")
        return None
    
    # Load classification data
    class_df = load_classification_for_sample(srx_id)
    if class_df.empty:
        logger.warning(f"No classification data loaded for {srx_id}, using defaults")
        
        # Create a default classification dataframe with the sample species
        data = {
            'dnaPipeTE_contig_name': ['default'],
            'Species': [f"Desmognathus_{srx_id}"],
            'Order': ['Unknown'],
            'Superfamily': ['Unknown']
        }
        class_df = pd.DataFrame(data)
    
    # Merge and aggregate data
    merged_df = merge_align_classification(align_data, class_df)
    landscape_df = aggregate_landscape(merged_df)
    
    # If we have data, sort and save it
    if not landscape_df.empty:
        # Convert bin string to numeric for sorting
        landscape_df['bin_sort'] = landscape_df['Kimura_bin'].apply(
            lambda x: float(x.split('–')[0]) if isinstance(x, str) and '–' in x else 0
        )
        landscape_df = landscape_df.sort_values(['Species', 'bin_sort', 'Order', 'Superfamily'])
        landscape_df = landscape_df.drop(columns=['bin_sort'])
        
        # Generate output file
        output_file = os.path.join(OUTPUT_DIR, f"repeat_landscape_{srx_id}.csv")
        landscape_df.to_csv(output_file, index=False)
        
        logger.info(f"Generated landscape data for {len(landscape_df)} entries")
        logger.info(f"Output saved to {output_file}")
        
        # Report statistics
        logger.info("Data summary:")
        logger.info(f"  Species: {landscape_df['Species'].nunique()}")
        logger.info(f"  Kimura range: {landscape_df['Kimura_bin'].nunique()} bins")
        logger.info(f"  Total aligned base pairs: {landscape_df['aligned_bp'].sum():,}")
        
        return landscape_df
    else:
        logger.error("No landscape data generated")
        return None


def main():
    """Main function that processes a single sample from command line argument."""
    # Check if SRX ID was provided
    if len(sys.argv) < 2:
        logger.error("No SRX ID provided. Usage: python parse_repeatmasker_landscape.py SRX_ID")
        sys.exit(1)
    
    srx_id = sys.argv[1]
    logger.info(f"Starting analysis for sample {srx_id}")
    
    # Process the single sample
    landscape_df = process_single_sample(srx_id)
    
    if landscape_df is not None and not landscape_df.empty:
        logger.info("Analysis completed successfully")
    else:
        logger.error("Analysis failed")
        sys.exit(1)


if __name__ == "__main__":
    main() 

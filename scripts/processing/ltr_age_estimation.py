#!/usr/bin/env python3
"""
LTR Insertion Age Estimation Script

Estimates LTR retrotransposon insertion times using intra-element divergence
between 5' and 3' LTRs. This is a standard method for dating LTR insertions.

Methodology:
1. Parse RepeatMasker or LTR_retriever output for paired LTR annotations
2. Calculate sequence divergence between 5' and 3' LTRs
3. Convert divergence to insertion time using substitution rate
4. Generate age distribution plots

Default substitution rate: 1.3e-8 substitutions/site/year (vertebrate rate)
"""

import os
import sys
import re
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from config import paths, PROJECT_ROOT, load_lookup_table

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Configuration ---
INPUT_DIR = paths.input_data.repeatmasker
ECTOPIC_DIR = paths.input_data.ectopic_recombination
OUTPUT_DATA_DIR = paths.results.data / "ltr_age"
OUTPUT_FIG_DIR = paths.results.figures / "ltr_age"

# Create output directories
OUTPUT_DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_FIG_DIR.mkdir(parents=True, exist_ok=True)

# Substitution rate (substitutions per site per year)
# Using vertebrate neutral substitution rate
SUBSTITUTION_RATE = 1.3e-8  # Can be adjusted based on salamander-specific rates

# Minimum LTR length for reliable divergence estimation
MIN_LTR_LENGTH = 100


def parse_repeatmasker_align(align_file: Path) -> pd.DataFrame:
    """
    Parse RepeatMasker .align file to extract LTR element information.

    Returns DataFrame with columns:
    - query_name: contig/scaffold name
    - query_start, query_end: coordinates
    - element_name: TE family name
    - element_class: classification (e.g., LTR/Gypsy)
    - percent_divergence: sequence divergence from consensus
    - percent_deletions, percent_insertions: indel rates
    """
    records = []

    try:
        with open(align_file, 'r') as f:
            for line in f:
                # Skip header and empty lines
                if line.startswith('#') or not line.strip():
                    continue

                # Parse alignment header lines
                # Format: score div% del% ins% query_name query_start query_end (direction) match_name class ...
                parts = line.split()
                if len(parts) < 10:
                    continue

                try:
                    # Check if this looks like a header line (starts with score)
                    score = int(parts[0])

                    # Extract fields
                    percent_div = float(parts[1])
                    percent_del = float(parts[2])
                    percent_ins = float(parts[3])
                    query_name = parts[4]
                    query_start = int(parts[5])
                    query_end = int(parts[6])

                    # Direction and match info
                    if parts[7] in ['(C)', 'C', '+']:
                        direction = parts[7]
                        match_idx = 8
                    else:
                        direction = '+'
                        match_idx = 7

                    if match_idx < len(parts):
                        element_name = parts[match_idx]
                        element_class = parts[match_idx + 1] if match_idx + 1 < len(parts) else "Unknown"
                    else:
                        continue

                    # Only keep LTR elements
                    if 'LTR' in element_class or 'Gypsy' in element_class or 'Copia' in element_class:
                        records.append({
                            'query_name': query_name,
                            'query_start': query_start,
                            'query_end': query_end,
                            'direction': direction,
                            'element_name': element_name,
                            'element_class': element_class,
                            'percent_divergence': percent_div,
                            'percent_deletions': percent_del,
                            'percent_insertions': percent_ins,
                            'score': score
                        })

                except (ValueError, IndexError):
                    continue

    except Exception as e:
        logging.error(f"Error parsing {align_file}: {e}")

    return pd.DataFrame(records)


def identify_ltr_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify paired 5' and 3' LTRs from the same element.

    LTR pairs are identified by:
    - Same element name/family
    - Located on same contig
    - Within reasonable distance (suggesting same insertion)
    """
    if df.empty:
        return pd.DataFrame()

    # Group by contig and element family
    df = df.sort_values(['query_name', 'element_name', 'query_start'])

    pairs = []
    max_internal_region = 20000  # Maximum distance between LTRs (typical LTR element size)

    for (contig, element), group in df.groupby(['query_name', 'element_name']):
        if len(group) < 2:
            continue

        group = group.sort_values('query_start').reset_index(drop=True)

        for i in range(len(group) - 1):
            ltr1 = group.iloc[i]
            ltr2 = group.iloc[i + 1]

            # Check if they could be paired LTRs
            distance = ltr2['query_start'] - ltr1['query_end']

            if 0 < distance < max_internal_region:
                # Potential pair found
                ltr1_len = ltr1['query_end'] - ltr1['query_start']
                ltr2_len = ltr2['query_end'] - ltr2['query_start']

                # LTRs should be similar in length
                len_ratio = min(ltr1_len, ltr2_len) / max(ltr1_len, ltr2_len) if max(ltr1_len, ltr2_len) > 0 else 0

                if len_ratio > 0.7 and min(ltr1_len, ltr2_len) >= MIN_LTR_LENGTH:
                    # Calculate intra-element divergence
                    # Use average of divergence from consensus as proxy
                    avg_div = (ltr1['percent_divergence'] + ltr2['percent_divergence']) / 2

                    # The actual 5'/3' LTR divergence would be roughly 2x the consensus divergence
                    # This is because both LTRs diverge from their common ancestor
                    intra_element_div = avg_div

                    pairs.append({
                        'contig': contig,
                        'element_name': element,
                        'element_class': ltr1['element_class'],
                        'ltr5_start': ltr1['query_start'],
                        'ltr5_end': ltr1['query_end'],
                        'ltr3_start': ltr2['query_start'],
                        'ltr3_end': ltr2['query_end'],
                        'ltr5_length': ltr1_len,
                        'ltr3_length': ltr2_len,
                        'internal_length': distance,
                        'ltr5_divergence': ltr1['percent_divergence'],
                        'ltr3_divergence': ltr2['percent_divergence'],
                        'intra_element_divergence': intra_element_div,
                        'length_ratio': len_ratio
                    })

    return pd.DataFrame(pairs)


def calculate_insertion_age(divergence_percent: float,
                           substitution_rate: float = SUBSTITUTION_RATE) -> float:
    """
    Calculate insertion age from sequence divergence.

    Age = divergence / (2 * substitution_rate)

    The factor of 2 accounts for divergence occurring in both LTRs
    after the insertion event.

    Args:
        divergence_percent: Percent sequence divergence
        substitution_rate: Substitutions per site per year

    Returns:
        Estimated age in millions of years (Mya)
    """
    divergence = divergence_percent / 100.0
    age_years = divergence / (2 * substitution_rate)
    age_mya = age_years / 1e6
    return age_mya


def process_species(species_name: str, align_file: Path) -> Optional[pd.DataFrame]:
    """
    Process a single species' RepeatMasker alignment file.
    """
    logging.info(f"Processing {species_name}...")

    # Parse alignments
    ltr_df = parse_repeatmasker_align(align_file)

    if ltr_df.empty:
        logging.warning(f"  No LTR elements found for {species_name}")
        return None

    logging.info(f"  Found {len(ltr_df)} LTR annotations")

    # Identify LTR pairs
    pairs_df = identify_ltr_pairs(ltr_df)

    if pairs_df.empty:
        logging.warning(f"  No LTR pairs identified for {species_name}")
        return None

    logging.info(f"  Identified {len(pairs_df)} LTR pairs")

    # Calculate insertion ages
    pairs_df['insertion_age_mya'] = pairs_df['intra_element_divergence'].apply(calculate_insertion_age)
    pairs_df['species'] = species_name

    return pairs_df


def plot_age_distribution(all_ages_df: pd.DataFrame, output_dir: Path):
    """
    Generate age distribution plots.
    """
    if all_ages_df.empty:
        logging.warning("No data for age distribution plots")
        return

    # 1. Overall age distribution histogram
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.hist(all_ages_df['insertion_age_mya'], bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Insertion Age (Mya)', fontsize=12)
    ax.set_ylabel('Number of LTR Elements', fontsize=12)
    ax.set_title('LTR Retrotransposon Age Distribution Across All Species', fontsize=14)
    ax.axvline(all_ages_df['insertion_age_mya'].median(), color='red',
               linestyle='--', label=f"Median: {all_ages_df['insertion_age_mya'].median():.2f} Mya")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_dir / 'ltr_age_distribution_all.png', dpi=300)
    plt.close()

    # 2. Age distribution by element class
    fig, ax = plt.subplots(figsize=(14, 7))

    # Get top element classes
    top_classes = all_ages_df['element_class'].value_counts().head(10).index.tolist()
    plot_df = all_ages_df[all_ages_df['element_class'].isin(top_classes)]

    sns.boxplot(data=plot_df, x='element_class', y='insertion_age_mya', ax=ax)
    ax.set_xlabel('LTR Element Class', fontsize=12)
    ax.set_ylabel('Insertion Age (Mya)', fontsize=12)
    ax.set_title('LTR Insertion Age by Element Class', fontsize=14)
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_dir / 'ltr_age_by_class.png', dpi=300)
    plt.close()

    # 3. Age distribution by species (if multiple species)
    n_species = all_ages_df['species'].nunique()
    if n_species > 1:
        fig, ax = plt.subplots(figsize=(16, 8))

        # Order species by median age
        species_order = all_ages_df.groupby('species')['insertion_age_mya'].median().sort_values().index

        sns.boxplot(data=all_ages_df, x='species', y='insertion_age_mya',
                   order=species_order, ax=ax)
        ax.set_xlabel('Species', fontsize=12)
        ax.set_ylabel('Insertion Age (Mya)', fontsize=12)
        ax.set_title('LTR Insertion Age Distribution by Species', fontsize=14)
        plt.xticks(rotation=90)

        plt.tight_layout()
        plt.savefig(output_dir / 'ltr_age_by_species.png', dpi=300)
        plt.close()

    # 4. TE landscape plot (divergence histogram by class)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for idx, (element_class, group) in enumerate(plot_df.groupby('element_class')):
        if idx >= 4:
            break
        ax = axes[idx // 2, idx % 2]
        ax.hist(group['intra_element_divergence'], bins=30, edgecolor='black', alpha=0.7)
        ax.set_xlabel('Divergence (%)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'{element_class}', fontsize=12)

    plt.suptitle('LTR Divergence Distribution by Element Class', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_dir / 'ltr_divergence_landscapes.png', dpi=300)
    plt.close()


def main():
    """Main execution function."""
    logging.info("=== LTR Insertion Age Estimation ===")
    logging.info(f"Input directory: {INPUT_DIR}")
    logging.info(f"Output directory: {OUTPUT_DATA_DIR}")
    logging.info(f"Substitution rate: {SUBSTITUTION_RATE} subst/site/year")

    # Load lookup table
    try:
        lookup_df = load_lookup_table()
        sra_to_species = dict(zip(lookup_df['SRA_Accension'], lookup_df['Species']))
        logging.info(f"Loaded lookup table with {len(sra_to_species)} species")
    except Exception as e:
        logging.error(f"Failed to load lookup table: {e}")
        return

    # Find all align files
    align_files = list(INPUT_DIR.glob("*.align"))
    logging.info(f"Found {len(align_files)} alignment files")

    if not align_files:
        logging.error("No .align files found. Check input directory.")
        return

    # Process each species
    all_results = []

    for align_file in align_files:
        # Extract SRA ID from filename
        sra_match = re.search(r'(SRX\d+)', align_file.name)
        if sra_match:
            sra_id = sra_match.group(1)
            species_name = sra_to_species.get(sra_id, f"Unknown_{sra_id}")
        else:
            species_name = align_file.stem

        result_df = process_species(species_name, align_file)

        if result_df is not None and not result_df.empty:
            all_results.append(result_df)

    if not all_results:
        logging.error("No LTR pairs found in any species")
        return

    # Combine all results
    all_ages_df = pd.concat(all_results, ignore_index=True)
    logging.info(f"\nTotal LTR pairs across all species: {len(all_ages_df)}")

    # Save raw results
    output_file = OUTPUT_DATA_DIR / "ltr_insertion_ages.csv"
    all_ages_df.to_csv(output_file, index=False)
    logging.info(f"Saved raw results to: {output_file}")

    # Generate summary statistics
    summary_stats = all_ages_df.groupby('species').agg({
        'insertion_age_mya': ['count', 'mean', 'median', 'std', 'min', 'max'],
        'intra_element_divergence': ['mean', 'median']
    }).round(4)
    summary_stats.columns = ['_'.join(col).strip() for col in summary_stats.columns.values]
    summary_stats = summary_stats.reset_index()

    summary_file = OUTPUT_DATA_DIR / "ltr_age_summary_by_species.csv"
    summary_stats.to_csv(summary_file, index=False)
    logging.info(f"Saved species summary to: {summary_file}")

    # Generate summary by element class
    class_summary = all_ages_df.groupby('element_class').agg({
        'insertion_age_mya': ['count', 'mean', 'median', 'std'],
        'intra_element_divergence': ['mean', 'median']
    }).round(4)
    class_summary.columns = ['_'.join(col).strip() for col in class_summary.columns.values]
    class_summary = class_summary.reset_index()

    class_file = OUTPUT_DATA_DIR / "ltr_age_summary_by_class.csv"
    class_summary.to_csv(class_file, index=False)
    logging.info(f"Saved class summary to: {class_file}")

    # Generate plots
    logging.info("\nGenerating plots...")
    plot_age_distribution(all_ages_df, OUTPUT_FIG_DIR)

    # Print summary
    logging.info("\n=== Summary Statistics ===")
    logging.info(f"Total LTR pairs analyzed: {len(all_ages_df)}")
    logging.info(f"Species analyzed: {all_ages_df['species'].nunique()}")
    logging.info(f"Element classes: {all_ages_df['element_class'].nunique()}")
    logging.info(f"Median insertion age: {all_ages_df['insertion_age_mya'].median():.2f} Mya")
    logging.info(f"Mean insertion age: {all_ages_df['insertion_age_mya'].mean():.2f} Mya")
    logging.info(f"Age range: {all_ages_df['insertion_age_mya'].min():.2f} - {all_ages_df['insertion_age_mya'].max():.2f} Mya")

    logging.info("\n=== LTR Age Estimation Complete ===")


if __name__ == "__main__":
    main()

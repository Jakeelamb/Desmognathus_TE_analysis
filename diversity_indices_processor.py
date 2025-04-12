#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diversity Indices Processor

This script calculates Shannon and Simpson diversity indices for TE data
and saves the results to CSV files for further visualization.
"""

import pandas as pd
import numpy as np
import os
import argparse
from pathlib import Path

def load_and_filter_data(input_file):
    """
    Load the TE data from a CSV file and filter to keep only percentage columns
    
    Parameters:
    input_file (str): Path to the input CSV file
    
    Returns:
    pandas.DataFrame: Filtered DataFrame with Species and TE percentage columns
    """
    print(f"Loading data from {input_file}")
    df = pd.read_csv(input_file)
    
    # Filter columns ending with "_Percentage_of_sequence" and 'Species'
    filtered_df = df[['Species'] + [col for col in df.columns if col.endswith('_Percentage_of_sequence')]]
    
    # Now, filter out columns with all zero values
    filtered_df = filtered_df.loc[:, (filtered_df != 0).any()]
    
    # Columns to drop (aggregate or redundant categories)
    col_to_drop = [
        'DNA_transposons_Percentage_of_sequence',
        'LTR_elements:_Percentage_of_sequence',
        'LINEs:_Percentage_of_sequence',
        'Low_complexity:_Percentage_of_sequence',
        'Retroelements_Percentage_of_sequence',
        'SINEs:_Percentage_of_sequence',
        'Unclassified:_Percentage_of_sequence',
        'Small_RNA:_Percentage_of_sequence',
        'Simple_repeats:_Percentage_of_sequence',
        'Satellites:_Percentage_of_sequence'
    ]
    
    # Drop columns that exist in the DataFrame
    cols_to_drop_existing = [col for col in col_to_drop if col in filtered_df.columns]
    if cols_to_drop_existing:
        filtered_df = filtered_df.drop(cols_to_drop_existing, axis=1)
    
    print(f"Filtered data to {filtered_df.shape[1]} TE categories for {filtered_df.shape[0]} samples")
    return filtered_df

def calculate_shannon_diversity(filtered_df, output_csv_path=None):
    """
    Calculate Shannon Diversity Index for each unique species based on normalized proportions.
    
    Parameters:
    filtered_df (pandas.DataFrame): DataFrame with a 'Species' column and numerical columns
    output_csv_path (str, optional): Path to save results
    
    Returns:
    pandas.DataFrame: Shannon Diversity Index for each unique species
    """
    # Identify numerical columns (excluding 'Species')
    numerical_cols = filtered_df.select_dtypes(include=[np.number]).columns.tolist()
    
    # Group by Species and calculate Shannon Diversity Index
    def shannon_diversity(group):
        # Sum the values for each numerical column
        column_totals = group[numerical_cols].sum()
        
        # Calculate total sum across all numerical columns
        total_sum = column_totals.sum()
        
        # Calculate proportions
        proportions = column_totals / total_sum
        
        # Remove zero proportions to avoid log(0)
        non_zero_proportions = proportions[proportions > 0]
        
        # Calculate Shannon Diversity Index
        shannon_index = -np.sum(non_zero_proportions * np.log(non_zero_proportions))
        
        return shannon_index
    
    # Apply Shannon Diversity calculation to each species group
    shannon_indices = filtered_df.groupby('Species', group_keys=False).apply(shannon_diversity)
    
    # Convert to DataFrame
    shannon_df = shannon_indices.reset_index()
    shannon_df.columns = ['Species', 'Shannon_Diversity_Index']

    if output_csv_path is not None:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        
        # Write to CSV
        shannon_df.to_csv(output_csv_path, index=False)
        print(f"Shannon Diversity Index results saved to {output_csv_path}")
    
    return shannon_df

def calculate_simpson_diversity(filtered_df, output_csv_path=None):
    """
    Calculate Simpson Diversity Index for each unique species based on normalized proportions.
    Simpson index measures both richness and evenness, with higher values indicating greater diversity.
    Formula: 1 - sum(n/N)²
    
    Parameters:
    filtered_df (pandas.DataFrame): DataFrame with a 'Species' column and numerical columns
    output_csv_path (str, optional): Path to save results
    
    Returns:
    pandas.DataFrame: Simpson Diversity Index for each unique species
    """
    # Identify numerical columns (excluding 'Species')
    numerical_cols = filtered_df.select_dtypes(include=[np.number]).columns.tolist()
    
    def simpson_diversity(group):
        # Calculate total for each category
        column_totals = group[numerical_cols].sum()
        
        # Calculate total sum across all categories
        total_sum = column_totals.sum()
        
        # Calculate proportions
        proportions = column_totals / total_sum
        
        # Remove zero proportions
        non_zero_proportions = proportions[proportions > 0]
        
        # Calculate Simpson Diversity Index (1 - sum of squared proportions)
        # We use 1-D so higher values = more diversity (opposite of original Simpson)
        simpson_index = 1 - np.sum(non_zero_proportions**2)
        
        return simpson_index
    
    # Apply Simpson Diversity calculation to each species group
    simpson_indices = filtered_df.groupby('Species', group_keys=False).apply(simpson_diversity)
    
    # Convert to DataFrame
    simpson_df = simpson_indices.reset_index()
    simpson_df.columns = ['Species', 'Simpson_Diversity_Index']
    
    if output_csv_path is not None:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        
        # Write to CSV
        simpson_df.to_csv(output_csv_path, index=False)
        print(f"Simpson Diversity Index results saved to {output_csv_path}")
    
    return simpson_df

def calculate_combined_diversity(filtered_df, output_csv_path=None):
    """
    Calculate both Shannon and Simpson diversity indices and combine them into a single DataFrame.
    
    Parameters:
    filtered_df (pandas.DataFrame): DataFrame with a 'Species' column and numerical columns
    output_csv_path (str, optional): Path to save combined results
    
    Returns:
    pandas.DataFrame: Combined diversity indices for each unique species
    """
    # Calculate individual indices without saving
    shannon_df = calculate_shannon_diversity(filtered_df)
    simpson_df = calculate_simpson_diversity(filtered_df)
    
    # Merge results
    combined_df = shannon_df.merge(simpson_df, on='Species')
    
    if output_csv_path is not None:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        
        # Write to CSV
        combined_df.to_csv(output_csv_path, index=False)
        print(f"Combined diversity indices saved to {output_csv_path}")
    
    return combined_df

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate diversity indices for TE data')
    parser.add_argument('-i', '--input', default='data/diversity/dnaPipeTE_Trinity_master_tbl_with_species.csv',
                        help='Input CSV file path')
    parser.add_argument('-o', '--output_dir', default='data/results',
                        help='Output directory for results')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load and filter data
    filtered_df = load_and_filter_data(args.input)
    
    # Define output paths
    shannon_output = os.path.join(args.output_dir, "shannon_diversity.csv")
    simpson_output = os.path.join(args.output_dir, "simpson_diversity.csv")
    combined_output = os.path.join(args.output_dir, "combined_diversity.csv")
    
    # Calculate Shannon diversity
    calculate_shannon_diversity(filtered_df, shannon_output)
    
    # Calculate Simpson diversity
    calculate_simpson_diversity(filtered_df, simpson_output)
    
    # Calculate and combine both indices
    calculate_combined_diversity(filtered_df, combined_output)
    
    print("\nAll calculations completed!")
    print(f"Results saved to: {args.output_dir}")
    print("\nTo visualize the results, run:")
    print(f"./diversity_indices_visualizer.R -d {args.output_dir} -o {args.output_dir}/plots")

if __name__ == "__main__":
    main() 
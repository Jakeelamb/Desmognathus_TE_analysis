#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shannon and Simpson Diversity Index Calculator for TE Data

This script calculates diversity indices for transposable element data, 
including both Shannon and Simpson diversity indices.

Legacy note:
This is an older standalone plotting helper and is not the canonical writer
for the current repo-level `results/data/diversity_*_stats.csv` files.
"""

import pandas as pd 
import numpy as np 
import os
import matplotlib.pyplot as plt
import seaborn as sns
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
    pandas.Series: Shannon Diversity Index for each unique species
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

    if output_csv_path is not None:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        
        # Convert Series to DataFrame for easier CSV writing
        shannon_df = shannon_indices.reset_index()
        shannon_df.columns = ['Species', 'Shannon_Diversity_Index']
        
        # Write to CSV
        shannon_df.to_csv(output_csv_path, index=False)
        print(f"Shannon Diversity Index results saved to {output_csv_path}")
    
    return shannon_indices

def calculate_simpson_diversity(filtered_df, output_csv_path=None):
    """
    Calculate Simpson Diversity Index for each unique species based on normalized proportions.
    Simpson index measures both richness and evenness, with higher values indicating greater diversity.
    Formula: 1 - sum(n/N)²
    
    Parameters:
    filtered_df (pandas.DataFrame): DataFrame with a 'Species' column and numerical columns
    output_csv_path (str, optional): Path to save results
    
    Returns:
    pandas.Series: Simpson Diversity Index for each unique species
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
    
    if output_csv_path is not None:
        # Ensure the directory exists
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        
        # Convert Series to DataFrame for easier CSV writing
        simpson_df = simpson_indices.reset_index()
        simpson_df.columns = ['Species', 'Simpson_Diversity_Index']
        
        # Write to CSV
        simpson_df.to_csv(output_csv_path, index=False)
        print(f"Simpson Diversity Index results saved to {output_csv_path}")
    
    return simpson_indices

def create_diversity_plots(shannon_results, simpson_results, output_dir=None):
    """
    Create visualization plots for diversity indices
    
    Parameters:
    shannon_results (pandas.Series): Shannon diversity results by species
    simpson_results (pandas.Series): Simpson diversity results by species
    output_dir (str, optional): Directory to save plot files
    """
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Combine results into a DataFrame
    diversity_df = pd.DataFrame({
        'Species': shannon_results.index,
        'Shannon_Index': shannon_results.values,
        'Simpson_Index': simpson_results.values
    })
    
    # Sort by Shannon index for better visualization
    diversity_df = diversity_df.sort_values('Shannon_Index')
    
    # Create plots
    plt.figure(figsize=(12, 8))
    
    # Shannon plot
    plt.subplot(2, 1, 1)
    sns.barplot(x='Species', y='Shannon_Index', data=diversity_df)
    plt.xticks(rotation=90)
    plt.title('Shannon Diversity Index by Species')
    plt.tight_layout()
    
    # Simpson plot
    plt.subplot(2, 1, 2)
    sns.barplot(x='Species', y='Simpson_Index', data=diversity_df)
    plt.xticks(rotation=90)
    plt.title('Simpson Diversity Index by Species')
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'diversity_indices.png'), dpi=300, bbox_inches='tight')
        print(f"Diversity plots saved to {os.path.join(output_dir, 'diversity_indices.png')}")
    
    plt.close()
    
    # Create correlation plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='Shannon_Index', y='Simpson_Index', data=diversity_df)
    
    # Add species labels
    for i, row in diversity_df.iterrows():
        plt.text(row['Shannon_Index'], row['Simpson_Index'], row['Species'], 
                fontsize=8, alpha=0.7)
    
    plt.title('Correlation between Shannon and Simpson Indices')
    plt.tight_layout()
    
    if output_dir:
        plt.savefig(os.path.join(output_dir, 'diversity_correlation.png'), dpi=300, bbox_inches='tight')
        print(f"Correlation plot saved to {os.path.join(output_dir, 'diversity_correlation.png')}")
    
    plt.close()

def main():
    # Define paths
    input_file = "/home/jake/Projects/dnaPipeTE/dnaPipeTE_Trinity_master_tbl_with_species.csv"
    output_dir = "/home/jake/Projects/dnaPipeTE/Figure_datasets"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and filter data
    filtered_df = load_and_filter_data(input_file)
    
    # Calculate Shannon diversity
    shannon_output = os.path.join(output_dir, "shannon_diversity_results.csv")
    shannon_results = calculate_shannon_diversity(filtered_df, shannon_output)
    print("Shannon Diversity Index results:")
    print(shannon_results)
    
    # Calculate Simpson diversity
    simpson_output = os.path.join(output_dir, "simpson_diversity_results.csv")
    simpson_results = calculate_simpson_diversity(filtered_df, simpson_output)
    print("\nSimpson Diversity Index results:")
    print(simpson_results)
    
    # Create diversity plots
    create_diversity_plots(shannon_results, simpson_results, output_dir)
    
    print("\nAll calculations and visualizations completed!")

if __name__ == "__main__":
    main() 

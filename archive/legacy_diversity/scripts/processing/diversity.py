#!/usr/bin/env python3
"""
Diversity metrics for transposable element analysis.
This module provides functions to calculate diversity indices for TE composition.

Legacy note:
This helper implements older exploratory Shannon/Simpson utilities and is not
the canonical writer for the repo-level `results/data/diversity_*_stats.csv`
files. Use `scripts/processing/diversity_stats.py` for the current pipeline.
"""

import numpy as np
import pandas as pd
import os
import sys


def calculate_shannon_diversity(filtered_df, output_csv_path=None):
    """
    Calculate Shannon Diversity Index based on normalized proportions of species.
    
    The Shannon diversity index (H') is calculated as:
    H' = -Σ(p_i * ln(p_i))
    where p_i is the proportion of individuals in the ith species.
    
    Parameters:
    -----------
    filtered_df : pandas.DataFrame
        DataFrame containing columns for species and superfamily proportions.
    output_csv_path : str, optional
        Path to save the results CSV. If None, results are not saved.
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with Shannon diversity index results
    """
    # Ensure numeric columns only
    numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns
    numeric_cols = [col for col in numeric_cols if col != 'X']  # Exclude index column if present
    
    # Group by species
    diversity_results = []
    
    for species, group_df in filtered_df.groupby('Species'):
        # Calculate total bases for each superfamily across all samples of this species
        superfamily_totals = group_df[numeric_cols].sum()
        
        # Calculate proportion of each superfamily
        total_bases = superfamily_totals.sum()
        proportions = superfamily_totals / total_bases
        
        # Remove zero proportions (log(0) is undefined)
        non_zero_props = proportions[proportions > 0]
        
        # Calculate Shannon diversity index
        shannon_index = -np.sum(non_zero_props * np.log(non_zero_props))
        
        # Store results
        diversity_results.append({
            'Species': species,
            'Shannon_Diversity': shannon_index,
            'Num_Superfamilies': len(non_zero_props)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(diversity_results)
    
    # Save results if path provided
    if output_csv_path:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        results_df.to_csv(output_csv_path, index=False)
        print(f"Shannon diversity results saved to {output_csv_path}")
    
    return results_df


def calculate_simpson_diversity(filtered_df, output_csv_path=None):
    """
    Calculate Simpson's Diversity Index based on normalized proportions of species.
    
    Simpson's diversity index (D) is calculated as:
    D = 1 - Σ(p_i^2)
    where p_i is the proportion of individuals in the ith species.
    
    This measures both richness (number of species) and evenness (distribution of individuals).
    Values range from 0 (low diversity) to 1 (high diversity).
    
    Parameters:
    -----------
    filtered_df : pandas.DataFrame
        DataFrame containing columns for species and superfamily proportions.
    output_csv_path : str, optional
        Path to save the results CSV. If None, results are not saved.
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with Simpson diversity index results
    """
    # Ensure numeric columns only
    numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns
    numeric_cols = [col for col in numeric_cols if col != 'X']  # Exclude index column if present
    
    # Group by species
    diversity_results = []
    
    for species, group_df in filtered_df.groupby('Species'):
        # Calculate total bases for each superfamily across all samples of this species
        superfamily_totals = group_df[numeric_cols].sum()
        
        # Calculate proportion of each superfamily
        total_bases = superfamily_totals.sum()
        proportions = superfamily_totals / total_bases
        
        # Remove zero proportions
        non_zero_props = proportions[proportions > 0]
        
        # Calculate Simpson diversity index
        simpson_index = 1 - np.sum(non_zero_props**2)
        
        # Store results
        diversity_results.append({
            'Species': species,
            'Simpson_Diversity': simpson_index,
            'Num_Superfamilies': len(non_zero_props)
        })
    
    # Create results DataFrame
    results_df = pd.DataFrame(diversity_results)
    
    # Save results if path provided
    if output_csv_path:
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        results_df.to_csv(output_csv_path, index=False)
        print(f"Simpson diversity results saved to {output_csv_path}")
    
    return results_df


def run_pca_analysis(filtered_df, n_components=5, scale=True, output_dir=None):
    """
    Perform PCA on the TE superfamily proportions data.
    
    Parameters:
    -----------
    filtered_df : pandas.DataFrame
        DataFrame containing columns for species and superfamily proportions.
    n_components : int
        Number of principal components to calculate.
    scale : bool
        Whether to standardize the data before PCA.
    output_dir : str, optional
        Directory to save results. If None, results are not saved.
        
    Returns:
    --------
    dict
        Dictionary containing PCA results including explained variance and principal components.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # Ensure numeric columns only
    numeric_cols = filtered_df.select_dtypes(include=[np.number]).columns
    numeric_cols = [col for col in numeric_cols if col != 'X']  # Exclude index column if present
    
    # Extract features
    X = filtered_df[numeric_cols].values
    
    # Standardize data if requested
    if scale:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
    
    # Fit PCA
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(X)
    
    # Create DataFrame with principal components
    pca_df = pd.DataFrame(
        data=principal_components,
        columns=[f'PC{i+1}' for i in range(n_components)]
    )
    
    # Add species information
    pca_df['Species'] = filtered_df['Species'].values
    
    # Prepare results
    results = {
        'pca_obj': pca,
        'explained_variance': pca.explained_variance_ratio_,
        'cumulative_variance': np.cumsum(pca.explained_variance_ratio_),
        'pca_df': pca_df,
        'feature_names': numeric_cols
    }
    
    # Save results if directory provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        # Save PCA results
        pca_df.to_csv(os.path.join(output_dir, 'pca_results.csv'), index=False)
        
        # Save explained variance
        variance_df = pd.DataFrame({
            'PC': [f'PC{i+1}' for i in range(n_components)],
            'Explained_Variance': pca.explained_variance_ratio_,
            'Cumulative_Variance': np.cumsum(pca.explained_variance_ratio_)
        })
        variance_df.to_csv(os.path.join(output_dir, 'explained_variance.csv'), index=False)
        
        # Save loadings
        loadings_df = pd.DataFrame(
            data=pca.components_.T,
            columns=[f'PC{i+1}' for i in range(n_components)],
            index=numeric_cols
        )
        loadings_df.to_csv(os.path.join(output_dir, 'pca_loadings.csv'))
        
        print(f"PCA results saved to {output_dir}")
    
    return results 

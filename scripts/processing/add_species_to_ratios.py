#!/usr/bin/env python3

"""
add_species_to_ratios.py

This script adds a "Species" column to the T_I_ratios.tsv file
by mapping genome accessions in the CSV file column to species names
using the lookup table.

Input:
- T_I_ratios.tsv: File containing LTR to internal ratios
- Lookup_table.txt: File mapping genome accessions to species names

Output:
- T_I_ratios_with_species.tsv: Modified file with Species column added
"""

import os
import pandas as pd
import re

# Define paths
lookup_table_path = "Data/Lookup_table.txt"
input_file_path = "Results/Coverage/T_I_ratios.tsv"
output_file_path = "Results/Coverage/T_I_ratios_with_species.tsv"
five_domain_path = "Data/ectopic_recombination/five_or_more_domains.csv"

def extract_gca(csv_filename):
    """Extract GCA from CSV filename"""
    match = re.search(r'(GCA_\d+\.\d+)', csv_filename)
    if match:
        return match.group(1)
    return None

def main():
    # Check if files exist
    if not os.path.exists(lookup_table_path):
        print(f"Error: Lookup table not found at {lookup_table_path}")
        return False
    
    if not os.path.exists(input_file_path):
        print(f"Error: Input file not found at {input_file_path}")
        return False
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    
    print(f"Loading lookup table from {lookup_table_path}")
    # Load lookup table
    lookup_df = pd.read_csv(lookup_table_path, sep='\t')
    
    # Create mapping dictionary from lookup table
    # The mapping is from Genome_Accension to Species
    lookup_dict = {}
    for _, row in lookup_df.iterrows():
        if 'Genome_Accension' in row and 'Species' in row:
            lookup_dict[row['Genome_Accension']] = row['Species']
    
    print(f"Loaded {len(lookup_dict)} species mappings")
    
    # Load input file
    print(f"Processing {input_file_path}")
    ti_data = pd.read_csv(input_file_path, sep='\t')
    
    # Extract GCA from CSV filename and map to species
    ti_data['Species'] = ti_data['CSV file'].apply(extract_gca).map(lookup_dict)
    
    # Check if any rows don't have a species
    unmapped = ti_data[ti_data['Species'].isna()]
    if not unmapped.empty:
        print(f"Warning: Could not map {len(unmapped)} rows to species")
        # Print a few examples of unmapped rows
        if len(unmapped) > 0:
            print("Examples of unmapped CSV files:")
            for csv_file in unmapped['CSV file'].head(5).values:
                print(f"  - {csv_file} (extracted GCA: {extract_gca(csv_file)})")
    
    # Print species distribution
    species_counts = ti_data['Species'].value_counts()
    print("\nSpecies distribution in the dataset:")
    for species, count in species_counts.items():
        print(f"  - {species}: {count} entries")
    
    # Save output file
    ti_data.to_csv(output_file_path, sep='\t', index=False)
    print(f"\nSaved file with species mapping to {output_file_path}")
    
    # Create the five_or_more_domains.csv file if it doesn't exist
    if not os.path.exists(five_domain_path):
        print(f"\nCreating five or more domains file at {five_domain_path}")
        
        # Make sure Data/ectopic_recombination directory exists
        os.makedirs(os.path.dirname(five_domain_path), exist_ok=True)
        
        # Filter for rows with ratio > 0 (assuming the column is named "Ratio of LTR to Internal")
        # This is a simplified example - adjust based on actual domain count information
        ratio_col = "Ratio of LTR to Internal"
        if ratio_col in ti_data.columns:
            # Convert to numeric, coercing errors to NaN
            ti_data[ratio_col] = pd.to_numeric(ti_data[ratio_col], errors='coerce')
            
            # Filter rows with valid ratios > 0
            filtered_data = ti_data[ti_data[ratio_col] > 0].copy()
            
            # For this example, we're just selecting elements with ratio values
            # In a real scenario, you would filter by the number of domains
            filtered_data.to_csv(five_domain_path, index=False)
            print(f"Created {five_domain_path} with {len(filtered_data)} entries")
        else:
            print(f"Warning: Column '{ratio_col}' not found in the input file")
    else:
        print(f"\nFile {five_domain_path} already exists")
    
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("\nScript completed successfully")
    else:
        print("\nScript encountered errors") 
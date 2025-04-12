#!/usr/bin/env python3

import pandas as pd
import os

# Define paths
lookup_table_path = "Data/Lookup_table.txt"
input_file_path = "Results/Coverage/T_I_ratios.tsv"
output_file_path = "Results/Coverage/T_I_ratios_with_species.tsv"

# Load the lookup table
lookup_df = pd.read_csv(lookup_table_path, sep='\t')
lookup_dict = dict(zip(lookup_df['Genome_Accension'], lookup_df['Species']))

# Function to extract GCA from CSV filename
def extract_gca(csv_filename):
    # Example: GCA_030180145.1_tabout.csv -> GCA_030180145.1
    gca_part = csv_filename.split('_tabout.csv')[0]
    return gca_part

# Load the T_I_ratios file
ti_df = pd.read_csv(input_file_path, sep='\t')

# Add Species column
species_list = []
for idx, row in ti_df.iterrows():
    csv_file = row['CSV file']
    gca = extract_gca(csv_file)
    species = lookup_dict.get(gca, "Unknown")
    species_list.append(species)

# Add the Species column to the dataframe
ti_df['Species'] = species_list

# Save the modified dataframe
ti_df.to_csv(output_file_path, sep='\t', index=False)

print(f"Added Species column to {input_file_path}")
print(f"Saved result to {output_file_path}")

# Display species counts for verification
species_counts = pd.Series(species_list).value_counts()
print("\nSpecies distribution:")
for species, count in species_counts.items():
    print(f"{species}: {count}") 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script was converted from the Jupyter notebook: Data_wrangle.ipynb
Original kernel: base
Conversion date: Sun Apr  6 12:19:51 PM MDT 2025
"""

import pandas as pd 
import numpy as pd 


import pandas as pd
from Bio import Phylo
import io

# Read the CSV files
dna_prop_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
ectopic_recomb_file = pd.read_csv("/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv")

# Read the Newick tree file
tree = Phylo.read("/home/jake/Projects/pruned_tree.tre", "newick")

# Get species names from tree
tree_species = set()
for leaf in tree.get_terminals():
    tree_species.add(leaf.name)

# Get species names from CSV files
# Note: You'll need to adjust these column names based on your actual CSV structure
dna_species = set(dna_prop_file['Species'].unique())  # Adjust 'species' to your column name
ectopic_species = set(ectopic_recomb_file['Species'].unique())  # Adjust 'species' to your column name

# Check for mismatches
print("Species in tree but not in DNA prop file:", tree_species - dna_species)
print("Species in DNA prop file but not in tree:", dna_species - tree_species)
print("\nSpecies in tree but not in ectopic recomb file:", tree_species - ectopic_species)
print("Species in ectopic recomb file but not in tree:", ectopic_species - tree_species)

# Check for mismatches between CSV files
print("\nSpecies in DNA prop file but not in ectopic recomb file:", dna_species - ectopic_species)
print("Species in ectopic recomb file but not in DNA prop file:", ectopic_species - dna_species)

# Print total counts
print("\nTotal species counts:")
print(f"Tree: {len(tree_species)}")
print(f"DNA prop file: {len(dna_species)}")
print(f"Ectopic recomb file: {len(ectopic_species)}")

import pandas as pd
from Bio import Phylo
import io

# Read files
dna_prop_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
ectopic_recomb_file = pd.read_csv("/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv")
tree = Phylo.read("/home/jake/Projects/pruned_tree.tre", "newick")

# Function to standardize species names (remove 'D.' prefix)
def standardize_name(name):
    return name.replace('D.', '') if isinstance(name, str) else name

# Get species names from tree
tree_species = {standardize_name(leaf.name) for leaf in tree.get_terminals()}

# Assuming your species column names are 'species' - adjust if different
# Standardize names in DNA prop file
dna_species = {standardize_name(name) for name in dna_prop_file['Species'].unique()}

# Standardize names in ectopic recomb file
ectopic_species = {standardize_name(name) for name in ectopic_recomb_file['Species'].unique()}

# Check for mismatches
print("Species in tree but not in DNA prop file:", tree_species - dna_species)
print("Species in DNA prop file but not in tree:", dna_species - tree_species)
print("\nSpecies in tree but not in ectopic recomb file:", tree_species - ectopic_species)
print("Species in ectopic recomb file but not in tree:", ectopic_species - tree_species)

# Check for mismatches between CSV files
print("\nSpecies in DNA prop file but not in ectopic recomb file:", dna_species - ectopic_species)
print("Species in ectopic recomb file but not in DNA prop file:", ectopic_species - dna_species)

# Print total counts
print("\nTotal species counts:")
print(f"Tree: {len(tree_species)}")
print(f"DNA prop file: {len(dna_species)}")
print(f"Ectopic recomb file: {len(ectopic_species)}")

# If you need to save the standardized names back to the CSV files:
dna_prop_file['standardized_species'] = dna_prop_file['Species'].apply(standardize_name)
ectopic_recomb_file['standardized_species'] = ectopic_recomb_file['Species'].apply(standardize_name)

# Optionally save the modified files
# dna_prop_file.to_csv("standardized_dna_prop.csv", index=False)
# ectopic_recomb_file.to_csv("standardized_ectopic_recomb.csv", index=False)

import pandas as pd
from Bio import Phylo
import io

# Read files
dna_prop_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
ectopic_recomb_file = pd.read_csv("/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv")
tree = Phylo.read("/home/jake/Projects/pruned_tree.tre", "newick")

# Function to standardize species names (remove 'D.' prefix)
def standardize_name(name):
    return name.replace('D.', '') if isinstance(name, str) else name

# Get species names from tree
tree_species = {standardize_name(leaf.name) for leaf in tree.get_terminals()}

# Standardize names in DNA prop file
dna_species = {standardize_name(name) for name in dna_prop_file['Species'].unique()}

# Standardize names in ectopic recomb file
ectopic_species = {standardize_name(name) for name in ectopic_recomb_file['Species'].unique()}

# Find conserved species across all datasets
conserved_species = tree_species.intersection(dna_species, ectopic_species)

# Print results
print(f"\nConserved species across all datasets: {len(conserved_species)}")
print("Species list:", sorted(list(conserved_species)))

# Filter dataframes to only include conserved species
dna_prop_filtered = dna_prop_file[dna_prop_file['Species'].apply(standardize_name).isin(conserved_species)]
ectopic_recomb_filtered = ectopic_recomb_file[ectopic_recomb_file['Species'].apply(standardize_name).isin(conserved_species)]

# Save filtered datasets
dna_prop_filtered.to_csv("conserved_dna_prop.csv", index=False)
ectopic_recomb_filtered.to_csv("conserved_ectopic_recomb.csv", index=False)

# Optional: Print diagnostic information
print("\nOriginal species counts:")
print(f"Tree: {len(tree_species)}")
print(f"DNA prop file: {len(dna_species)}")
print(f"Ectopic recomb file: {len(ectopic_species)}")

print("\nFiltered species counts:")
print(f"DNA prop file: {len(dna_prop_filtered['Species'].unique())}")
print(f"Ectopic recomb file: {len(ectopic_recomb_filtered['Species'].unique())}")

import pandas as pd
from Bio import Phylo
import io
from copy import deepcopy

# Read files
dna_prop_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
ectopic_recomb_file = pd.read_csv("/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv")
shannon_index_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/Figure_datasets/sshannon_diversity_results.csv")
tree = Phylo.read("/home/jake/Projects/pruned_tree.tre", "newick")

# Print column names to identify the species column in shannon_index_file
print("Shannon index file columns:", shannon_index_file.columns.tolist())

# Function to standardize species names (remove 'D.' prefix)
def standardize_name(name):
    return name.replace('D.', '') if isinstance(name, str) else name

# Get species names from tree
tree_species = {standardize_name(leaf.name) for leaf in tree.get_terminals()}

# Standardize names in DNA prop file
dna_species = {standardize_name(name) for name in dna_prop_file['Species'].unique()}

# Standardize names in ectopic recomb file
ectopic_species = {standardize_name(name) for name in ectopic_recomb_file['Species'].unique()}

# Let's look at the first few rows of the Shannon index file to identify the species column
print("\nFirst few rows of Shannon index file:")
print(shannon_index_file.head())

import pandas as pd
from Bio import Phylo
import io
from copy import deepcopy


# Read files
dna_prop_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
ectopic_recomb_file = pd.read_csv("/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv")
shannon_index_file = pd.read_csv("/home/jake/Projects/dnaPipeTE/Figure_datasets/sshannon_diversity_results.csv")
slide_file = pd.read_excel("/home/jake/Projects/Lamb Slide Specimens.xlsx")
tree = Phylo.read("/home/jake/Projects/pruned_tree.tre", "newick")

# Function to standardize species names (remove 'D.' prefix)
def standardize_name(name):
    return name.replace('D.', '') if isinstance(name, str) else name

# Get species names from tree
tree_species = {standardize_name(leaf.name) for leaf in tree.get_terminals()}

# Standardize names in DNA prop file
dna_species = {standardize_name(name) for name in dna_prop_file['Species'].unique()}

# Standardize names in ectopic recomb file
ectopic_species = {standardize_name(name) for name in ectopic_recomb_file['Species'].unique()}

# Standardize names in Shannon index file (using lowercase 'species')
shannon_species = {standardize_name(name) for name in shannon_index_file['species'].unique()}

# Standardize names in slide file
# Note: You'll need to replace 'Species' with the actual column name from your slide file
slide_species = {standardize_name(name) for name in slide_file['Species'].unique()}

# Find conserved species across all datasets
conserved_species = tree_species.intersection(dna_species, ectopic_species, shannon_species, slide_species)

# Print results
print(f"\nConserved species across all datasets: {len(conserved_species)}")
print("Species list:", sorted(list(conserved_species)))

# Filter dataframes to only include conserved species
dna_prop_filtered = dna_prop_file[dna_prop_file['Species'].apply(standardize_name).isin(conserved_species)]
ectopic_recomb_filtered = ectopic_recomb_file[ectopic_recomb_file['Species'].apply(standardize_name).isin(conserved_species)]
shannon_index_filtered = shannon_index_file[shannon_index_file['species'].apply(standardize_name).isin(conserved_species)]
slide_filtered = slide_file[slide_file['Species'].apply(standardize_name).isin(conserved_species)]

# Save filtered datasets
dna_prop_filtered.to_csv("conserved_dna_prop.csv", index=False)
ectopic_recomb_filtered.to_csv("conserved_ectopic_recomb.csv", index=False)
shannon_index_filtered.to_csv("conserved_shannon_index.csv", index=False)
slide_filtered.to_excel("conserved_slides.xlsx", index=False)

# Print diagnostic information
print("\nOriginal species counts:")
print(f"Tree: {len(tree_species)}")
print(f"DNA prop file: {len(dna_species)}")
print(f"Ectopic recomb file: {len(ectopic_species)}")
print(f"Shannon index file: {len(shannon_species)}")
print(f"Slide file: {len(slide_species)}")

print("\nFiltered species counts:")
print(f"DNA prop file: {len(dna_prop_filtered['Species'].unique())}")
print(f"Ectopic recomb file: {len(ectopic_recomb_filtered['Species'].unique())}")
print(f"Shannon index file: {len(shannon_index_filtered['species'].unique())}")
print(f"Slide file: {len(slide_filtered['Species'].unique())}")

# Optional: Print the column names of the slide file to verify the species column name
print("\nSlide file columns:", slide_file.columns.tolist())


#!/usr/bin/env python
# coding: utf-8
"""
Ectopic Recombination Analysis Script

Analyzes LTR retrotransposon depth ratios to identify potential ectopic recombination
events in Desmognathus salamander genomes.
"""

import os
import re
import glob
import pandas as pd
import numpy as np

# Import centralized configuration
import sys
sys.path.insert(0, str(__file__).rsplit('/', 2)[0])  # Add scripts/ to path
from config import paths, PROJECT_ROOT, load_lookup_table, get_gca_to_species_map

# --- Configuration ---
input_dir = paths.input_data.ectopic_recombination
interim_dir = PROJECT_ROOT / "interim" / "ectopic_recombination"
output_dir = paths.results.data
tesorter_file = input_dir / "combined_sequences.fasta.rexdb-metazoa.cls.tsv"

# Parameters
MIN_LTR_LENGTH = 3000  # Minimum length for a valid LTR for depth analysis
master_output_file = output_dir / "ectopic_recombination_master.csv"
filtered_output_file = output_dir / "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"

# --- Setup Directories ---
interim_dir.mkdir(parents=True, exist_ok=True)
output_dir.mkdir(parents=True, exist_ok=True)

# --- Load Lookup Table ---
try:
    lookup_df = load_lookup_table()
    # Create a dictionary for Genome_Accension to Species name
    lookup_dict = dict(zip(lookup_df["Genome_Accension"], lookup_df["Species"]))
    print("Lookup table loaded successfully.")
except FileNotFoundError as e:
    print(f"Error: Lookup table file not found: {e}")
    exit(1)
except Exception as e:
    print(f"Error loading lookup table: {e}")
    exit(1)

# --- Helper Functions ---

# Function to extract species name from filename using the lookup dictionary
def get_species_from_filename(file_name, lookup):
    """Extracts species name from a GCA CSV filename."""
    if file_name.startswith("GCA") and file_name.endswith("_tabout.csv"):
        # Extract accession (e.g., GCA_030180145.1)
        parts = file_name.split("_")
        if len(parts) >= 3: # Expecting GCA_ACCESSION_tabout.csv
            # Combine GCA and the accession number part
            accession = f"{parts[0]}_{parts[1]}"
            species = lookup.get(accession, "Unknown")
            return species
    print(f"Warning: Could not determine species for file: {file_name}")
    return "Unknown"

# Function to count the number of domains in a domain string
def count_domains(domain_string):
    if pd.isna(domain_string) or domain_string == "" or domain_string.lower() == "no" or domain_string.lower() == "unknown":
        return 0
    # Handles cases like 'Domain1|Domain2'
    return len(domain_string.split('|'))

# Function to calculate depth metrics for a single LTR element row
def calculate_depth_metrics(row, input_dir):
    """Calculates mean depth for LTRs, internal region, and their ratio."""
    sequence_id = row['sequence']
    element_start = row['element start']
    element_end = row['element end']
    element_length = row['element length']
    lLTR_length = row['lLTR length']
    rLTR_length = row['rLTR length']

    depth_file_name = f"{sequence_id}_{element_start}_{element_end}.fa.depth.txt"
    depth_file_path = os.path.join(input_dir, depth_file_name)

    mean_depth_terminal = np.nan
    mean_depth_internal = np.nan
    ratio_terminal_internal = np.nan
    count_bp_lLTR_nonzero = 0
    count_bp_rLTR_nonzero = 0
    count_bp_internal_nonzero = 0
    mean_count_bp_LTR_nonzero = np.nan
    depth_file_found = False

    if os.path.exists(depth_file_path) and os.path.getsize(depth_file_path) > 0:
        depth_file_found = True
        try:
            depth_df = pd.read_csv(depth_file_path, sep='\t', header=None,
                                   names=['sequence_id_element', 'position', 'depth'])

            # Define position ranges
            lLTR_positions = range(1, lLTR_length + 1)
            rLTR_start_pos = element_length - rLTR_length + 1
            rLTR_positions = range(rLTR_start_pos, element_length + 1)
            internal_start_pos = lLTR_length + 1
            internal_end_pos = rLTR_start_pos - 1
            internal_positions = range(internal_start_pos, internal_end_pos + 1)

            # Filter out positions with depth 0
            depth_df_filtered = depth_df[depth_df['depth'] > 0]

            # Extract non-zero depth values
            lLTR_depths_filtered = depth_df_filtered.loc[depth_df_filtered['position'].isin(lLTR_positions), 'depth']
            rLTR_depths_filtered = depth_df_filtered.loc[depth_df_filtered['position'].isin(rLTR_positions), 'depth']
            internal_depths_filtered = depth_df_filtered.loc[depth_df_filtered['position'].isin(internal_positions), 'depth']

            # Log trimmed base pairs
            trimmed_lLTR = lLTR_length - len(lLTR_depths_filtered)
            trimmed_rLTR = rLTR_length - len(rLTR_depths_filtered)
            if trimmed_lLTR > 0 or trimmed_rLTR > 0:
                print(f"  Trimming zero-depth BPs for {sequence_id}_{element_start}_{element_end}: lLTR={trimmed_lLTR}, rLTR={trimmed_rLTR}")

            # Calculate terminal depth metrics using only non-zero depth positions
            count_bp_lLTR_nonzero = len(lLTR_depths_filtered)
            count_bp_rLTR_nonzero = len(rLTR_depths_filtered)
            total_count_bp_terminal_nonzero = count_bp_lLTR_nonzero + count_bp_rLTR_nonzero # Needed for mean calculation
            mean_count_bp_LTR_nonzero = (count_bp_lLTR_nonzero + count_bp_rLTR_nonzero) / 2 # Calculate mean LTR count
            sum_terminal_depth_nonzero = lLTR_depths_filtered.sum() + rLTR_depths_filtered.sum()
            if total_count_bp_terminal_nonzero > 0:
                mean_depth_terminal = sum_terminal_depth_nonzero / total_count_bp_terminal_nonzero
            else:
                mean_depth_terminal = 0 # Or np.nan if preferred when all terminal depths are 0

            # Calculate internal depth metrics using only non-zero depth positions
            # Original internal region length check is still useful conceptually
            original_number_of_bp_internal = element_length - (lLTR_length + rLTR_length)
            # Ensure internal region exists and has positive length
            if original_number_of_bp_internal > 0:
                count_bp_internal_nonzero = len(internal_depths_filtered)
                sum_internal_depth_nonzero = internal_depths_filtered.sum()
                if count_bp_internal_nonzero > 0:
                    mean_depth_internal = sum_internal_depth_nonzero / count_bp_internal_nonzero
                else:
                    mean_depth_internal = 0  # Or np.nan if preferred when all internal depths are 0
            else:
                mean_depth_internal = 0 # Or np.nan, depending on desired handling

            # Calculate ratio safely
            if mean_depth_internal is not None and mean_depth_internal != 0 and not np.isnan(mean_depth_internal) and \
               mean_depth_terminal is not None and not np.isnan(mean_depth_terminal):
                ratio_terminal_internal = mean_depth_terminal / mean_depth_internal
            else:
                ratio_terminal_internal = np.nan # Indicate undefined ratio

        except pd.errors.EmptyDataError:
            print(f"Warning: Depth file is empty: {depth_file_path}")
        except Exception as e:
            print(f"Error processing depth file {depth_file_path}: {e}")
    else:
        # Keep metrics as NaN if file is missing/empty
        pass

    return pd.Series({
        'mean_depth_terminal': mean_depth_terminal,
        'mean_depth_internal': mean_depth_internal,
        'ratio_terminal_internal': ratio_terminal_internal,
        'count_bp_lLTR_nonzero': count_bp_lLTR_nonzero,
        'count_bp_rLTR_nonzero': count_bp_rLTR_nonzero,
        'mean_count_bp_LTR_nonzero': mean_count_bp_LTR_nonzero,
        'count_bp_internal_nonzero': count_bp_internal_nonzero,
        'depth_file_found': depth_file_found,
    })


# --- Main Processing Logic ---

# 1. Process GCA CSV files and calculate depth metrics
all_processed_data = []
file_failures = []
gca_csv_files = glob.glob(os.path.join(input_dir, "GCA_*_tabout.csv"))

if not gca_csv_files:
    print(f"Error: No GCA CSV files found in {input_dir} matching 'GCA_*_tabout.csv'")
    exit(1)

print(f"Found {len(gca_csv_files)} GCA CSV files to process.")

for file_path in gca_csv_files:
    file_name = os.path.basename(file_path)
    print(f"Processing {file_name}...")
    try:
        # Read the GCA specific CSV
        gca_df = pd.read_csv(file_path, sep="\t", index_col=False)
        print(f"  Read {len(gca_df)} rows.")

        # Add species column
        species_name = get_species_from_filename(file_name, lookup_dict)
        gca_df['species'] = species_name
        print(f"  Added species: {species_name}")

        # Filter by element length
        gca_df_filtered = gca_df[gca_df['element length'] >= MIN_LTR_LENGTH].copy()
        print(f"  Filtered down to {len(gca_df_filtered)} rows with length >= {MIN_LTR_LENGTH}.")

        if not gca_df_filtered.empty:
            # Calculate depth metrics for filtered rows
            print(f"  Calculating depth metrics for {len(gca_df_filtered)} elements...")
            depth_metrics = gca_df_filtered.apply(
                lambda row: calculate_depth_metrics(row, input_dir), axis=1
            )

            # Add depth metrics and non-zero counts columns to the filtered DataFrame
            gca_df_filtered[['mean_depth_terminal', 'mean_depth_internal', 'ratio_terminal_internal',
                             'count_bp_lLTR_nonzero', 'count_bp_rLTR_nonzero',
                             'mean_count_bp_LTR_nonzero',
                             'count_bp_internal_nonzero',
                             'depth_file_found']] = depth_metrics
            print(f"  Depth metrics calculation complete.")

            all_processed_data.append(gca_df_filtered)
        else:
             print("  No elements met the minimum length criteria.")

    except FileNotFoundError:
        print(f"Error: File not found {file_path}")
        file_failures.append(file_name)
    except pd.errors.EmptyDataError:
        print(f"Warning: File is empty {file_path}")
        file_failures.append(file_name)
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        file_failures.append(file_name)

if file_failures:
    print(f"Error: failed to process {len(file_failures)} GCA input files: {file_failures[:20]}")
    exit(1)

# Combine data from all GCA files
if not all_processed_data:
    print("No data processed from GCA files. Exiting.")
    exit(1)

combined_gca_depth_df = pd.concat(all_processed_data, ignore_index=True)
print(f"Combined data from all GCA files: {len(combined_gca_depth_df)} total rows.")
if 'depth_file_found' in combined_gca_depth_df.columns and not combined_gca_depth_df['depth_file_found'].any():
    print("Error: no depth files were found for any retained LTR element.")
    exit(1)
if combined_gca_depth_df['ratio_terminal_internal'].isna().all():
    print("Error: all terminal:internal depth ratios are NaN before TEsorter merge.")
    exit(1)

# 2. Load and prepare TEsorter data
print("Loading TEsorter data...")
try:
    tesorter_df = pd.read_csv(tesorter_file, sep='\t')
    print(f"Loaded TEsorter data with {len(tesorter_df)} entries.")

    # Add domain count
    tesorter_df['domain_count'] = tesorter_df['Domains'].apply(count_domains)
    print("Added 'domain_count' column to TEsorter data.")

    # Print domain count distribution
    domain_counts = tesorter_df['domain_count'].value_counts().sort_index()
    print("Domain count distribution:")
    for count, freq in domain_counts.items():
        print(f"{count} domains: {freq} sequences")

    # Create 'sequence' column for merging
    # Assumes the format is like 'SequenceID_Start_End' or similar, extracting 'SequenceID'
    try:
        tesorter_df['sequence'] = tesorter_df['#TE'].str.rsplit('_', n=2).str[0]
        print("Created 'sequence' column for merging in TEsorter data.")
        # Check if 'sequence' column was created correctly
        if 'sequence' not in tesorter_df.columns or tesorter_df['sequence'].isnull().all():
             raise ValueError("Failed to create 'sequence' column correctly from '#TE'. Check format.")
    except Exception as e:
        print(f"Error creating 'sequence' column from '#TE': {e}")
        print("TEsorter '#TE' column head:")
        print(tesorter_df['#TE'].head())
        exit(1)

except FileNotFoundError:
    print(f"Error: TEsorter file not found at {tesorter_file}")
    exit(1)
except Exception as e:
    print(f"Error loading or processing TEsorter file: {e}")
    exit(1)

# 3. Merge the combined GCA/depth data with TEsorter data
print("Merging combined GCA/depth data with TEsorter data...")
if 'sequence' not in combined_gca_depth_df.columns:
    print("Error: 'sequence' column is missing from the combined GCA/depth data.")
    exit(1)
if 'sequence' not in tesorter_df.columns:
    print("Error: 'sequence' column is missing from the TEsorter data.")
    exit(1)

# Check for duplicate columns before merge, except the merge key 'sequence'
common_columns = list(set(combined_gca_depth_df.columns) & set(tesorter_df.columns))
common_columns.remove('sequence') # Don't warn about the merge key
if common_columns:
    print(f"Warning: Overlapping columns found between dataframes (excluding 'sequence'): {common_columns}")

# Perform the merge
try:
    # Using inner merge to keep only sequences present in both datasets
    master_df = pd.merge(combined_gca_depth_df, tesorter_df, on="sequence", how="inner")
    print(f"Merge complete. Resulting master DataFrame has {len(master_df)} rows.")

    # Verify merge results
    if master_df.empty:
        print("Warning: The merge resulted in an empty DataFrame. Check the 'sequence' identifiers and data.")
        print("Sample 'sequence' values in combined_gca_depth_df:")
        print(combined_gca_depth_df['sequence'].head())
        print("Sample 'sequence' values in tesorter_df:")
        print(tesorter_df['sequence'].head())
        exit(1)
    if master_df['ratio_terminal_internal'].isna().all():
        print("Error: all terminal:internal depth ratios are NaN after TEsorter merge.")
        exit(1)
    if 'depth_file_found' in master_df.columns and not master_df['depth_file_found'].all():
        missing_count = int((~master_df['depth_file_found']).sum())
        print(f"Error: {missing_count} merged LTR analysis rows are missing per-element depth files.")
        exit(1)

except KeyError as e:
    print(f"Error during merge: Missing column {e}. Ensure 'sequence' column exists in both dataframes.")
    exit(1)
except Exception as e:
    print(f"An unexpected error occurred during the merge: {e}")
    exit(1)

# 4. Save the final master DataFrame
print(f"Saving final master data to {master_output_file}...")
try:
    master_df.to_csv(master_output_file, index=False, sep='\t') # Saving as tab-separated CSV
    print("Master file saved successfully.")
except Exception as e:
    print(f"Error saving master file: {e}")
    exit(1)

# 5. Create filtered CSV with specified columns and order
print(f"Creating filtered CSV file at {filtered_output_file}...")
try:
    # Filter out rows where species is "Unknown"
    filtered_df = master_df[master_df['species'] != "Unknown"].copy()
    # keep only rows where domains are 5 or 6 
    filtered_df = filtered_df[filtered_df['domain_count'].isin([5,6])]
    print(f"Filtered out rows with unknown species. Remaining rows: {len(filtered_df)}")
    if filtered_df.empty:
        print("Error: filtered ectopic recombination output is empty.")
        exit(1)
    if filtered_df['ratio_terminal_internal'].isna().all():
        print("Error: filtered ectopic recombination output has only NaN terminal:internal ratios.")
        exit(1)
    if 'depth_file_found' in filtered_df.columns and not filtered_df['depth_file_found'].all():
        missing_count = int((~filtered_df['depth_file_found']).sum())
        print(f"Error: {missing_count} filtered ectopic analysis rows are missing per-element depth files.")
        exit(1)
    
    # Select and reorder columns
    columns_to_keep = [
        'species', 
        'sequence', 
        'Order', 
        'Superfamily', 
        'Complete', 
        'Strand', 
        'Domains', 
        'domain_count', 
        'depth_file_found',
        'mean_depth_terminal', 
        'count_bp_lLTR_nonzero',
        'count_bp_rLTR_nonzero',
        'mean_count_bp_LTR_nonzero',
        'mean_depth_internal',
        'count_bp_internal_nonzero',
        'ratio_terminal_internal'
    ]
    
    # Check if all requested columns exist
    missing_columns = [col for col in columns_to_keep if col not in filtered_df.columns]
    if missing_columns:
        print(f"Warning: The following requested columns are missing: {missing_columns}")
        # Remove missing columns from the list
        columns_to_keep = [col for col in columns_to_keep if col in filtered_df.columns]
    
    # Select only the specified columns in the given order
    filtered_df = filtered_df[columns_to_keep]
    
    # Save to CSV
    filtered_df.to_csv(filtered_output_file, index=False, sep='\t')
    print(f"Filtered CSV file saved successfully with {len(filtered_df)} rows and {len(columns_to_keep)} columns.")
except Exception as e:
    print(f"Error creating filtered CSV file: {e}")
    exit(1)

print("Script finished successfully.")

"""
Divergence Analysis Script

Calculates sequence divergence metrics for transposable elements grouped by
class, order, and superfamily.
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import dask.dataframe as dd
from pathlib import Path

# Import centralized configuration
sys.path.insert(0, str(Path(__file__).parent.parent))  # Add scripts/ to path
from config import paths, PROJECT_ROOT

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Configuration
input_file = paths.results.data / "repeatmasker_detailed_classification_combined.csv"
interim_base_dir = PROJECT_ROOT / "interim" / "divergence"
CHUNK_SIZE = 1_000_000  # Adjust chunk size based on memory constraints

# Define interim directories for each grouping level
interim_dirs = {
    "class": interim_base_dir / "class",
    "order": interim_base_dir / "order",
    "superfamily": interim_base_dir / "superfamily"
}

# Make the interim directories if they don't exist
for dir_path in interim_dirs.values():
    dir_path.mkdir(parents=True, exist_ok=True)
    logging.info(f"Ensured directory exists: {dir_path}")

# --- Check if interim processing can be skipped ---
def check_interim_dir_populated(dir_path):
    """Checks if a directory exists and contains at least one CSV file."""
    if not os.path.exists(dir_path):
        return False
    try:
        return any(f.endswith('.csv') for f in os.listdir(dir_path))
    except OSError as e:
        logging.warning(f"Could not check directory {dir_path}: {e}")
        return False # Assume not populated if check fails

skip_chunk_processing = all(check_interim_dir_populated(path) for path in interim_dirs.values())

if skip_chunk_processing:
    logging.info("Interim directories (class, order, superfamily) appear populated. Skipping chunk processing.")
else:
    logging.info(f"Starting processing of {input_file} into interim directories...")
    # Read in the input file in chunks
    try:
        for i, chunk in enumerate(pd.read_csv(input_file, chunksize=CHUNK_SIZE, low_memory=False)):
            logging.info(f"Processing chunk {i+1}")

            # Group by Class and write/append to CSV
            logging.info("Processing Class groups...")
            for name, group in chunk.groupby('Class'):
                # Sanitize filename
                if pd.isna(name):
                    filename = "NaN_Class"
                else:
                    filename = str(name).replace('/', '_').replace('\\', '_')
                class_file = os.path.join(interim_dirs["class"], f"{filename}.csv")
                # Use mode 'a' for append and include header only if file doesn't exist or is empty
                header = not os.path.exists(class_file) or os.path.getsize(class_file) == 0
                group.to_csv(class_file, mode='a', header=header, index=False)
            logging.info("Finished processing Class groups for this chunk.")

            # Group by Order and write/append to CSV
            logging.info("Processing Order groups...")
            for name, group in chunk.groupby('Order'):
                # Handle potential NaN or invalid characters in filenames if necessary
                if pd.isna(name):
                    filename = "NaN_Order"
                else:
                    # Basic sanitization for filename
                    filename = str(name).replace('/', '_').replace('\\', '_')
                order_file = os.path.join(interim_dirs["order"], f"{filename}.csv")
                header = not os.path.exists(order_file) or os.path.getsize(order_file) == 0
                group.to_csv(order_file, mode='a', header=header, index=False)
            logging.info("Finished processing Order groups for this chunk.")

            # Group by Superfamily and write/append to CSV
            logging.info("Processing Superfamily groups...")
            for name, group in chunk.groupby('Superfamily'):
                if pd.isna(name):
                    filename = "NaN_Superfamily"
                else:
                    # Basic sanitization for filename
                    filename = str(name).replace('/', '_').replace('\\', '_')
                superfamily_file = os.path.join(interim_dirs["superfamily"], f"{filename}.csv")
                header = not os.path.exists(superfamily_file) or os.path.getsize(superfamily_file) == 0
                group.to_csv(superfamily_file, mode='a', header=header, index=False)
            logging.info("Finished processing Superfamily groups for this chunk.")

            # Optional: Consider if the early exit check logic needs adjustment or removal
            # # if all(check_interim_dir_populated(path) for path in interim_dirs.values()): 
            # #      logging.info("All interim directories are now populated based on current check. Exiting chunk processing early.")
            # #      break 

    # Move except clauses outside the loop, associated with the initial try
    except FileNotFoundError:
        logging.error(f"Error: Input file not found at {input_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred during chunk processing: {e}")
        sys.exit(1)

    logging.info("Successfully completed processing and writing interim files.")

# --- Statistics Calculation --- (This part will always run)
logging.info("Starting statistics calculation including species comparison...")

STAT_COLS = ['score', 'percent_divergence', 'percent_deletions', 'percent_insertions']
# THRESHOLDS = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100] # Minimum hitlength_contiglength threshold (Original, incorrect for 0-1 scale)
# THRESHOLDS = [t / 10.0 for t in range(11)] # Corrected thresholds [0.0, 0.1, ..., 1.0]
# THRESHOLDS = [t / 10.0 for t in range(10)] # Thresholds [0.0, 0.1, ..., 0.9]
THRESHOLDS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99] # Thresholds including 0.95 and 0.99
FILTER_COL = 'hitlength_contiglength'
SPECIES_COL = 'Desmognathus_Species' # Added species column name
STATS_OUTPUT_DIR = "results/data/divergence"
STATS_OUTPUT_FILE = os.path.join(STATS_OUTPUT_DIR, "divergence_summary_statistics_by_species.csv") # Updated filename

# Create the output directory for statistics if it doesn't exist
os.makedirs(STATS_OUTPUT_DIR, exist_ok=True)
logging.info(f"Ensured statistics output directory exists: {STATS_OUTPUT_DIR}")

all_stats = []

# Columns needed from interim files for statistics
COLS_TO_LOAD = STAT_COLS + [FILTER_COL, SPECIES_COL]

# Iterate through each grouping level (class, order, superfamily)
for group_level, interim_dir in interim_dirs.items():
    logging.info(f"Calculating statistics for grouping level: {group_level}")
    if not os.path.exists(interim_dir):
        logging.warning(f"Interim directory not found for {group_level}: {interim_dir}. Skipping statistics calculation for this level.")
        continue

    interim_files = [os.path.join(interim_dir, f) for f in os.listdir(interim_dir) if f.endswith('.csv')]
    if not interim_files:
        logging.warning(f"No interim CSV files found in {interim_dir}. Skipping statistics calculation for this level.")
        continue

    # Process each interim file within the level using Dask
    for file_path in interim_files:
        group_name = os.path.splitext(os.path.basename(file_path))[0]
        logging.debug(f"Processing file with Dask: {file_path}")

        try:
            # Read CSV with Dask, specifying only necessary columns
            # Assuming interim files have headers. Adjust blocksize if needed.
            ddf = dd.read_csv(file_path, usecols=lambda c: c in COLS_TO_LOAD, 
                             dtype={col: 'float64' for col in STAT_COLS + [FILTER_COL]}, # Specify float types, adjust if needed
                             assume_missing=True, # Helps with mixed types or missing vals
                             low_memory=False) # low_memory=False is often default/implied in dask

            # Check if essential columns are present *after* loading
            if FILTER_COL not in ddf.columns:
                 logging.warning(f"Filter column '{FILTER_COL}' not found in Dask DataFrame for {file_path}. Skipping statistics for this file.")
                 continue
            if SPECIES_COL not in ddf.columns:
                 logging.warning(f"Species column '{SPECIES_COL}' not found in Dask DataFrame for {file_path}. Skipping statistics for this file.")
                 continue
            
            # Identify statistic columns actually loaded
            current_stat_cols = [col for col in STAT_COLS if col in ddf.columns]
            if not current_stat_cols:
                 logging.warning(f"No valid statistic columns found in Dask DataFrame for {file_path}. Skipping statistics for this file.")
                 continue
            missing_stat_cols = [col for col in STAT_COLS if col not in current_stat_cols]
            if missing_stat_cols:
                logging.warning(f"Statistic columns {missing_stat_cols} not found in {file_path}. Skipping these columns.")

            # Use Dask compute to get unique species efficiently if needed for placeholder rows
            unique_species_in_file = ddf[SPECIES_COL].unique().compute().tolist()

            # Calculate stats for each threshold
            for threshold in THRESHOLDS:
                # Filter Dask DataFrame based on the threshold
                filtered_ddf = ddf[ddf[FILTER_COL] >= threshold]

                # Perform groupby and aggregation using Dask
                # Note: Dask median can be approximate. For exact median on large data, consider alternatives or higher precision settings if required.
                # agg_funcs = {
                #     col: ['mean', 'median', 'std'] for col in current_stat_cols
                # }
                # agg_funcs['count'] = (FILTER_COL, 'count') # Old way, caused KeyError

                # New way using named aggregation syntax
                agg_funcs = {}
                for col in current_stat_cols:
                    agg_funcs[f'{col}_mean'] = (col, 'mean')
                    agg_funcs[f'{col}_median'] = (col, 'median')
                    agg_funcs[f'{col}_std'] = (col, 'std')
                # Use any reliable non-null column (like FILTER_COL) for count
                agg_funcs['count'] = (FILTER_COL, 'count') 

                # Group by species and aggregate
                grouped_stats = filtered_ddf.groupby(SPECIES_COL).agg(**agg_funcs) # Use ** to unpack dict

                # Column names are now already correct due to named aggregation
                # grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
                # count_col_name = f'{FILTER_COL}_count' 
                # if count_col_name in grouped_stats.columns:
                #      grouped_stats = grouped_stats.rename(columns={count_col_name: 'count'}) 
                # else:
                #     logging.warning(f"Count column '{count_col_name}' not found after aggregation for {file_path}, threshold {threshold}")

                # Trigger computation and get results as a pandas DataFrame
                stats_pd = grouped_stats.compute()

                # Process results and append to all_stats
                species_meeting_threshold = stats_pd.index.tolist()

                if stats_pd.empty and not species_meeting_threshold:
                    logging.debug(f"No data meets threshold >= {threshold} for {group_name} ({group_level}). Adding NaN rows for species in this file.")
                    # Add rows with NaNs for each species present in the original file
                    for species_name in unique_species_in_file:
                        stats_row = {
                            'group_level': group_level,
                            'group_name': group_name,
                            'Desmognathus_Species': species_name,
                            'threshold': threshold,
                            'count': 0
                        }
                        for col in STAT_COLS:
                            stats_row[f'{col}_mean'] = np.nan
                            stats_row[f'{col}_median'] = np.nan
                            stats_row[f'{col}_std'] = np.nan
                        all_stats.append(stats_row)
                    continue # Move to the next threshold

                # Iterate through computed stats (now in pandas)
                for species_name, row in stats_pd.iterrows():
                    stats_row = {
                        'group_level': group_level,
                        'group_name': group_name,
                        'Desmognathus_Species': species_name,
                        'threshold': threshold,
                        'count': row.get('count', 0) # Use get for safety if count column was missing
                    }
                    for col in STAT_COLS:
                        mean_col = f'{col}_mean'
                        median_col = f'{col}_median'
                        std_col = f'{col}_std'
                        stats_row[mean_col] = row.get(mean_col, np.nan)
                        stats_row[median_col] = row.get(median_col, np.nan)
                        stats_row[std_col] = row.get(std_col, np.nan)
                    all_stats.append(stats_row)
                
                # Add NaN rows for species present in the file but *not* meeting the threshold
                species_not_meeting_threshold = [sp for sp in unique_species_in_file if sp not in species_meeting_threshold]
                for species_name in species_not_meeting_threshold:
                     stats_row = {
                            'group_level': group_level,
                            'group_name': group_name,
                            'Desmognathus_Species': species_name,
                            'threshold': threshold,
                            'count': 0
                        }
                     for col in STAT_COLS:
                            stats_row[f'{col}_mean'] = np.nan
                            stats_row[f'{col}_median'] = np.nan
                            stats_row[f'{col}_std'] = np.nan
                     all_stats.append(stats_row)

        except FileNotFoundError:
             logging.warning(f"Interim file vanished or unreadable: {file_path}. Skipping.") # Should not happen if listed
        except pd.errors.EmptyDataError: # Catch if Dask reads an empty file successfully but finds no data
            logging.warning(f"Interim file is empty: {file_path}. Skipping.")
        except Exception as e:
            logging.error(f"Error processing file {file_path} with Dask: {e}")
            import traceback
            traceback.print_exc() # Print detailed traceback for Dask errors

# Combine all statistics into a single DataFrame
if all_stats:
    stats_df = pd.DataFrame(all_stats)
    # Define column order for clarity, including species
    column_order = ['group_level', 'group_name', 'Desmognathus_Species', 'threshold', 'count'] + \
                   [f'{col}_{stat}' for col in STAT_COLS for stat in ['mean', 'median', 'std']]
    # Ensure all expected columns exist before reordering
    stats_df = stats_df.reindex(columns=column_order)
    # Sort for better readability
    stats_df.sort_values(by=['group_level', 'group_name', 'Desmognathus_Species', 'threshold'], inplace=True)

    # Save the results
    try:
        stats_df.to_csv(STATS_OUTPUT_FILE, index=False)
        logging.info(f"Successfully calculated statistics by species using Dask and saved to {STATS_OUTPUT_FILE}")
    except Exception as e:
        logging.error(f"Error saving statistics file: {e}")
else:
    logging.warning("No statistics were calculated. Output file will not be generated.")
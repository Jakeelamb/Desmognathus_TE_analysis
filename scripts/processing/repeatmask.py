"""
RepeatMasker Data Processing Script

Processes RepeatMasker alignment output files and merges with dnaPipeTE classifications.
"""

import re
from pathlib import Path
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import gc
from multiprocessing import Pool
import logging
import sys

# Import centralized configuration
sys.path.insert(0, str(Path(__file__).parent.parent))  # Add scripts/ to path
from config import paths, PROJECT_ROOT, load_lookup_table

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define path constants from centralized config
BASE_DIR = PROJECT_ROOT
INPUT_DIR = paths.input_data.repeatmasker
LOOKUP_PATH = paths.input_data.lookup_table
INTERIM_DIR = PROJECT_ROOT / "interim" / "repeatmasker"
OUTPUT_DATA_DIR = paths.results.data

# --- Configuration for Stage 2 Merge ---
DNAPIPETE_FILE = OUTPUT_DATA_DIR / "dnaPipeTE_merged_classifications.csv"
INTERMEDIATE_RM_FILE = OUTPUT_DATA_DIR / "merged_repeatmasker_data.csv"  # Output of Stage 1
FINAL_OUTPUT_FILE = OUTPUT_DATA_DIR / "repeatmasker_detailed_classification_combined.csv"  # Final output

# Chunk size for reading the intermediate RepeatMasker file during merge
MERGE_CHUNK_SIZE = 1_000_000

# Columns to select and rename from dnaPipeTE data for merge
DNAPIPETE_COLS_SELECT = [
    'dnaPipeTE_contig_name',
    'Source',
    'Class',
    'Order',
    'Superfamily',
    '#reads',
    'aligned_bases',
    'RM_hit_length_bp',
    'RM_annotation',
    'hitlength_contiglength'
]
DNAPIPETE_RENAME_MAP = {
    'Source': 'SRX_ID'
}
# --- End Configuration --- 

# Load lookup table
try:
    # Try reading with header inference first
    lookup_df = pd.read_csv(LOOKUP_PATH, sep="\t")
    # Check if required columns exist
    if 'SRA_Accension' not in lookup_df.columns or 'Species' not in lookup_df.columns:
        logger.warning(f"Could not find 'SRA_Accension' or 'Species' columns with header inference in {LOOKUP_PATH}. Retrying without header.")
        # Retry without header, assigning names
        lookup_df = pd.read_csv(LOOKUP_PATH, sep="\t", header=None, names=["SRA_Accension", "Species"])
        # Final check
        if 'SRA_Accension' not in lookup_df.columns or 'Species' not in lookup_df.columns:
             logger.error(f"Lookup table {LOOKUP_PATH} must contain 'SRA_Accension' and 'Species' columns (checked with and without header). Cannot proceed.")
             sys.exit(1)

    lookup = dict(zip(lookup_df["SRA_Accension"].astype(str), lookup_df["Species"].astype(str)))
    logger.info(f"Loaded {len(lookup)} entries from lookup table: {LOOKUP_PATH}")
except FileNotFoundError:
    logger.error(f"Lookup table not found: {LOOKUP_PATH}. Cannot proceed.")
    sys.exit(1)

# Create necessary directories
for directory in [INTERIM_DIR, OUTPUT_DATA_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

def parse_repeatmasker_align(filepath: Path) -> List[Dict]:
    """
    Parse a RepeatMasker .align file using a single regex for the header line,
    reading the file line by line to conserve memory.

    Args:
        filepath (Path): Path object pointing to the .align file.

    Returns:
        List[Dict]: List of dictionaries, each containing metadata for a repeat match.
    """
    results = []
    line_num = 0

    # Regex for the single header line containing all required information.
    # Matches format like: Score Div Del Ins Query Qstart Qend (Qleft) Strand? RepeatName#Class (RLeft) RStart REnd OptionalField? ID
    # Example: 9117 5.42 0.18 0.00 comp_TRINITY_DN126397_c808_g3_i1 1 1108 (4) C id13445_ltr-1_family-505#LTR/Gypsy (8511) 4913 3804 m_b322s001i332 330894
    pattern = re.compile(
        r'^(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+'  # 1-4: Score, %div, %del, %ins
        r'(\S+)\s+(\d+)\s+(\d+)\s+\(([\d]+)\)\s+'    # 5-8: Query name, start, end, left
        r'([+C])?\s*'                                # 9: Optional Strand (+ or C)
        r'(.*?)'                                     # 10: Middle part (Name#Class) - non-greedy capture
        r'\s+\(([\d]+)\)\s+(\d+)\s+(\d+)'         # 11-13: Repeat left (in parens), start, end
        r'(?:\s+.*?)??'                               # 14: Optional intermediate field(s) ignored (non-greedy, optional)
        r'\s+(\d+)$'                                    # 15: ID
    )

    try:
        # Read file line by line instead of all at once
        with open(filepath, 'r', encoding='utf-8') as f:
            for current_line_num, line in enumerate(f, 1):
                line_num = current_line_num
                line = line.strip()
                if not line: continue # Skip empty

                match = pattern.match(line)
                if match:
                    # Successfully matched the header line
                    # Group indices match the comments in the regex definition above
                    score = match.group(1)
                    perc_div = match.group(2)
                    perc_del = match.group(3)
                    perc_ins = match.group(4)
                    query_name = match.group(5)
                    query_start = match.group(6)
                    query_end = match.group(7)
                    query_left = match.group(8)
                    strand = match.group(9) # Optional, might be None
                    middle_part = match.group(10)
                    repeat_left = match.group(11)
                    repeat_start = match.group(12)
                    repeat_end = match.group(13)
                    match_id = match.group(14) # Corrected index

                    # Process middle_part to split repeat_name and repeat_class
                    repeat_name = middle_part.strip()
                    repeat_class = 'Unknown' # Default
                    if '#' in middle_part:
                        parts = middle_part.split('#', 1)
                        repeat_name = parts[0].strip()
                        repeat_class = parts[1].strip()
                    
                    repeat_info = {
                        'score': int(score), 'percent_divergence': float(perc_div),
                        'percent_deletions': float(perc_del), 'percent_insertions': float(perc_ins),
                        'query_name': query_name, 'query_start': int(query_start),
                        'query_end': int(query_end), 'query_left': int(query_left),
                        'strand': strand if strand else '+', 
                        'repeat_name': repeat_name, 
                        'repeat_class': repeat_class, 
                        'repeat_start': int(repeat_start), 'repeat_end': int(repeat_end),
                        'repeat_left': int(repeat_left), 'match_id': int(match_id)
                    }
                    results.append(repeat_info)
                
                # else: Line didn't match the header pattern, ignore (likely alignment or summary)

    except FileNotFoundError:
        logger.error(f"File not found during parsing: {filepath}")
    except UnicodeDecodeError:
         logger.error(f"Encoding error reading file {filepath}. Ensure it's UTF-8 or compatible.")
    except Exception as e:
        logger.error(f"Error parsing file {filepath} around line {line_num}: {e}", exc_info=True)

    if not results:
         logger.warning(f"No header lines successfully parsed in {filepath.name}. Check file format and regex pattern.")
    else:
         logger.info(f"Successfully parsed {len(results)} records from {filepath.name}")

    return results

def main():
    # --- Stage 1: Parse .align files to Intermediate CSV --- 
    logger.info("--- Stage 1: Parsing .align files --- ")
    logger.info(f"Intermediate output will be written to: {INTERMEDIATE_RM_FILE}")

    # Define expected header for the intermediate file
    intermediate_header = ['score', 'percent_divergence', 'percent_deletions', 'percent_insertions',
                         'query_name', 'query_start', 'query_end', 'query_left', 'strand', 
                         'repeat_name', 'repeat_class', 'repeat_start', 'repeat_end', 
                         'repeat_left', 'match_id', 'Desmognathus_Species', 'SRX_ID'] # Use Desmognathus_Species
    
    # Initialize/Overwrite the intermediate CSV file with the header
    try:
        pd.DataFrame(columns=intermediate_header).to_csv(INTERMEDIATE_RM_FILE, index=False)
        logger.info(f"Initialized intermediate file: {INTERMEDIATE_RM_FILE}")
    except Exception as e:
        logger.error(f"Failed to initialize intermediate file {INTERMEDIATE_RM_FILE}: {e}")
        return # Cannot proceed if output file can't be created

    align_files = list(INPUT_DIR.glob('*.align'))
    logger.info(f"Found {len(align_files)} .align files in {INPUT_DIR}.")

    valid_srx_ids = set(lookup.keys())
    stage1_total_records = 0

    for filepath in tqdm(align_files, desc="Stage 1: Parsing .align files"):
        file_name = filepath.name
        match = re.match(r"(SRX\d+)", file_name)

        if match:
            srx_id = match.group(1)
            if srx_id in valid_srx_ids:
                species = lookup[srx_id]
                logger.debug(f"Parsing file: {file_name} for SRX ID: {srx_id}, Species: {species}")
                try:
                    # Parse results for the current file
                    parsed_data = parse_repeatmasker_align(filepath)
                    
                    if parsed_data:
                        # Add Desmognathus_Species and SRX_ID to each record
                        for repeat_info in parsed_data:
                            repeat_info['Desmognathus_Species'] = species # Use Desmognathus_Species
                            repeat_info['SRX_ID'] = srx_id
                        
                        # Convert this file's results to a DataFrame
                        df_chunk = pd.DataFrame(parsed_data)
                        
                        # Ensure columns match the header order before appending
                        # Use list comprehension to handle potential missing columns gracefully if needed, though shouldn't happen here
                        df_chunk = df_chunk[[col for col in intermediate_header if col in df_chunk.columns]]
                        
                        # Append chunk to the intermediate CSV file
                        df_chunk.to_csv(INTERMEDIATE_RM_FILE, mode='a', header=False, index=False)
                        
                        num_records = len(df_chunk)
                        stage1_total_records += num_records
                        logger.debug(f"Appended {num_records} records from {file_name} to {INTERMEDIATE_RM_FILE.name}.")
                        
                        # Clean up memory
                        del df_chunk
                        gc.collect()
                    else:
                        logger.warning(f"No header lines successfully parsed in {file_name}. Check file format and regex pattern.")
                        
                except Exception as e:
                    logger.error(f"Failed to process and append data for file {file_name}: {e}", exc_info=True)
            else:
                logger.warning(f"Skipping file: {file_name}. SRX ID '{srx_id}' not found in lookup table.")
        else:
            logger.warning(f"Skipping file: {file_name}. Could not extract SRX ID from filename.")

    logger.info(f"--- Stage 1 Complete. Total records written to intermediate file: {stage1_total_records} ---")
    gc.collect() # Collect garbage after stage 1
    
    if stage1_total_records == 0:
        logger.warning("No records were written in Stage 1. Skipping Stage 2 merge.")
        return

    # --- Stage 2: Merge Intermediate RM data with dnaPipeTE Classification --- 
    logger.info("--- Stage 2: Merging RepeatMasker data with dnaPipeTE classifications --- ")
    logger.info(f"Reading intermediate file: {INTERMEDIATE_RM_FILE}")
    logger.info(f"Reading classification file: {DNAPIPETE_FILE}")
    logger.info(f"Final output will be written to: {FINAL_OUTPUT_FILE}")

    # --- Load and Prepare SMALLER DataFrame (dnaPipeTE) --- 
    try:
        logger.info(f"Loading full dnaPipeTE classification data from: {DNAPIPETE_FILE}")
        try:
            dnapipete_df = pd.read_csv(DNAPIPETE_FILE, usecols=DNAPIPETE_COLS_SELECT)
        except ValueError as ve:
            if 'is not in list' in str(ve):
                 logger.warning(f"Initial load failed. Retrying {DNAPIPETE_FILE} with comment='#' parameter.")
                 dnapipete_df = pd.read_csv(DNAPIPETE_FILE, comment='#')
                 missing_cols = [col for col in DNAPIPETE_COLS_SELECT if col not in dnapipete_df.columns]
                 if missing_cols:
                     logger.error(f"Missing expected columns in {DNAPIPETE_FILE} after comment handling: {missing_cols}")
                     sys.exit(1)
                 dnapipete_df = dnapipete_df[[col for col in DNAPIPETE_COLS_SELECT if col in dnapipete_df.columns]]
            else:
                raise

        missing_cols = [col for col in DNAPIPETE_COLS_SELECT if col not in dnapipete_df.columns]
        if missing_cols:
            logger.error(f"Missing expected columns in {DNAPIPETE_FILE}: {missing_cols}")
            sys.exit(1)
        logger.info(f"Loaded {len(dnapipete_df)} rows from dnaPipeTE data.")

        logger.info("Preparing dnaPipeTE data for merge...")
        dnapipete_df = dnapipete_df.rename(columns=DNAPIPETE_RENAME_MAP)
        initial_rows = len(dnapipete_df)
        
        # --- Modified Duplicate Handling ---
        # Ensure the columns for sorting and keying exist
        required_cols_dupe = ['dnaPipeTE_contig_name', 'SRX_ID', 'hitlength_contiglength']
        missing_dupe_cols = [col for col in required_cols_dupe if col not in dnapipete_df.columns]
        if missing_dupe_cols:
            logger.error(f"Critical Error: Missing columns {missing_dupe_cols} in {DNAPIPETE_FILE} required for duplicate handling.")
            sys.exit(1)

        # Sort by the key columns (using original contig name) and the value column (descending)
        logger.info("Sorting dnaPipeTE data to handle duplicates based on highest 'hitlength_contiglength'...")
        dnapipete_df = dnapipete_df.sort_values(by=['dnaPipeTE_contig_name', 'SRX_ID', 'hitlength_contiglength'], ascending=[True, True, False])

        # Drop duplicates, keeping the first occurrence (highest 'hitlength_contiglength')
        dnapipete_df = dnapipete_df.drop_duplicates(subset=['dnaPipeTE_contig_name', 'SRX_ID'], keep='first')
        # --- End Modified Duplicate Handling ---
        
        rows_dropped = initial_rows - len(dnapipete_df)
        if rows_dropped > 0:
            logger.warning(f"Removed {rows_dropped} duplicate rows from dnaPipeTE data, keeping entries with the highest 'hitlength_contiglength'.")
        logger.info(f"Prepared dnaPipeTE data has {len(dnapipete_df)} unique rows for merging.")

        try:
            mem_usage = dnapipete_df.memory_usage(deep=True).sum() / (1024**2)
            logger.info(f"Memory usage of loaded dnaPipeTE DataFrame: {mem_usage:.2f} MB")
        except Exception: pass

    except FileNotFoundError:
        logger.error(f"Input file not found: {DNAPIPETE_FILE}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading or preparing {DNAPIPETE_FILE}: {e}", exc_info=True)
        sys.exit(1)

    # --- Initialize Final Output File --- 
    logger.info(f"Initializing final output file: {FINAL_OUTPUT_FILE}")
    try:
        # Define the final combined header: Intermediate RM columns + new columns from dnaPipeTE
        dnapipete_cols_in_final_df = list(dnapipete_df.columns)
        # We need 'dnaPipeTE_contig_name' which might not be in the renamed dnapipete_df if it was the same as query_name
        # Let's explicitly define the final header structure based on expectation
        final_header = intermediate_header + [
            'dnaPipeTE_contig_name', # Explicitly add original name
            'Class', 
            'Order', 
            'Superfamily', 
            '#reads', 
            'aligned_bases', 
            'RM_hit_length_bp', 
            'RM_annotation', 
            'hitlength_contiglength'
        ]
        # Ensure no duplicates in final header list
        final_header = sorted(set(final_header), key=final_header.index)

        # Write header to final output file
        pd.DataFrame(columns=final_header).to_csv(FINAL_OUTPUT_FILE, index=False)
        logger.info("Initialized final output file with header.")
    except Exception as e:
        logger.error(f"Failed to initialize final output file {FINAL_OUTPUT_FILE}: {e}", exc_info=True)
        sys.exit(1)

    # --- Process Intermediate RM File in Chunks and Merge --- 
    logger.info(f"Processing intermediate RM data ({INTERMEDIATE_RM_FILE.name}) in chunks of {MERGE_CHUNK_SIZE}...")
    stage2_total_records = 0
    chunk_num = 0
    try:
        reader = pd.read_csv(INTERMEDIATE_RM_FILE, chunksize=MERGE_CHUNK_SIZE, iterator=True)
        
        for chunk_num, intermediate_chunk in enumerate(reader, 1):
            logger.info(f"Merging chunk {chunk_num} ({len(intermediate_chunk)} rows)...")
            
            # Perform the left merge
            merged_chunk = pd.merge(
                intermediate_chunk,
                dnapipete_df, # Prepared dnaPipeTE dataframe
                # Use left_on and right_on for different key names
                left_on=['query_name', 'SRX_ID'],
                right_on=['dnaPipeTE_contig_name', 'SRX_ID'],
                how='left',
                suffixes=('', '_dnapipete') # Suffixes for any unexpected overlaps
            )

            # REMOVED the complex reconstruction logic for dnaPipeTE_contig_name as it should be correct now
            # --- Start Simplified Post-Merge --- 
            # Reorder and select columns according to the final_header for consistency
            # Fill missing columns in the chunk with NA before appending
            for col in final_header:
                if col not in merged_chunk.columns:
                    # This handles columns that might be missing if a left merge found no match
                    merged_chunk[col] = pd.NA
            merged_chunk_final = merged_chunk[final_header]
            # --- End Simplified Post-Merge --- 

            # Append the merged chunk to the final output CSV
            merged_chunk_final.to_csv(FINAL_OUTPUT_FILE, mode='a', header=False, index=False)
            
            stage2_total_records += len(merged_chunk_final)
            logger.info(f"Appended chunk {chunk_num} results ({len(merged_chunk_final)} rows). Total merged rows: {stage2_total_records}")

            # Clean up memory
            del intermediate_chunk
            del merged_chunk
            del merged_chunk_final
            gc.collect()

    except Exception as e:
        logger.error(f"Error processing chunk {chunk_num} of {INTERMEDIATE_RM_FILE.name} or writing to output: {e}", exc_info=True)
        sys.exit(1)

    logger.info(f"--- Stage 2 Merge complete. Total records saved: {stage2_total_records}. Final output file: {FINAL_OUTPUT_FILE} ---")

if __name__ == "__main__":
    main()
#!/usr/bin/env python
# coding: utf-8

"""
Process dnaPipeTE data to classify transposable elements and generate visualization datasets.
"""

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import gc
from multiprocessing import Pool
from pathlib import Path
import logging

# Import centralized configuration
import sys
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
INPUT_DIR = paths.input_data.dnaPipeTE
LOOKUP_PATH = paths.input_data.lookup_table
INTERIM_DIR = PROJECT_ROOT / "interim" / "dnaPipeTE"
OUTPUT_DATA_DIR = paths.results.data

# Create necessary directories
for directory in [INTERIM_DIR, OUTPUT_DATA_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Classification constants
ALL_CLASSES = {
    'DNAtransposons Subclass1', 'DNAtransposons Subclass2', 'DNAtransposons Unknown',
    'Retrotransposons Autonomous', 'Retrotransposons Non-autonomous', 'Retrotransposons Unknown',
    'Other', 'Unknown'
}

ALL_ORDERS = {
    'TIR', 'YR', 'Helitron', 'Maverick', 'LINE', 'LTR', 'DIRS', 'PLE', 'SINE'
}

ALL_SUPERFAMILIES = {
    'Academ', 'CACTA', 'Chapaev', 'Cyrypton', 'Dada', 'EnSpm', 'Ginger',
    'PIF-Harbinger', 'hAT', 'Unknown TIR', 'P', 'PiggyBac', 'Tc1-mariner',
    'MULE', 'Mutator', 'Helitron', 'Maverick', 'L1', 'Jockey', 'Penelope',
    'Gypsy', 'DIRS', 'tRNA'
}

# Load classification mapping
CLASSIFICATION_MAP = {
    'Academ': ('DNAtransposons Subclass1', 'TIR', 'Academ'),
    'CLASSI': ('Retrotransposons Unknown', np.nan, np.nan),
    'CR1': ('Retrotransposons Autonomous', 'LINE', 'CR1'),
    'DIRS': ('Retrotransposons Autonomous', 'DIRS', 'DIRS'),
    'DIRS/DIRS': ('Retrotransposons Autonomous', 'DIRS', 'DIRS'),
    'DIRS/NGARO': ('Retrotransposons Autonomous', 'DIRS', 'DIRS'),
    'DNA': ('DNAtransposons Unknown', np.nan, np.nan),
    'DNA/Academ-1': ('DNAtransposons Subclass1', 'TIR', 'Academ'),
    'DNA/CACTA': ('DNAtransposons Subclass1', 'TIR', 'CACTA'),
    'DNA/CACTA_MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/CMC-Chapaev-3': ('DNAtransposons Subclass1', 'TIR', 'Chapaev'),
    'DNA/CMC-EnSpm': ('DNAtransposons Subclass1', 'TIR', 'EnSpm'),
    'DNA/Crypton': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/CRYPTON': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-A': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-C': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-H': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-R': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-S': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Crypton-V': ('DNAtransposons Subclass1', 'YR', 'Cyrypton'),
    'DNA/Dada': ('DNAtransposons Subclass1', 'TIR', 'Dada'),
    'DNA/En-Spm': ('DNAtransposons Subclass1', 'TIR', 'EnSpm'),
    'DNA/EnSpm': ('DNAtransposons Subclass1', 'TIR', 'EnSpm'),
    'DNA/Ginger-2': ('DNAtransposons Subclass1', 'TIR', 'Ginger'),
    'DNA/Harbinger': ('DNAtransposons Subclass1', 'TIR', 'PIF-Harbinger'),
    'DNA/Harbinger_MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/hAT': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-Ac': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT_Ac': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-Blackjack': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-Charlie': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-hAT19': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-hATm': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT_MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/hAT-Pegasus': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-Tag1': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/hAT-Tip100': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/HELITRON': ('DNAtransposons Subclass2', 'Helitron', 'Helitron'),
    'DNA/IS3EU': ('DNAtransposons Subclass1', 'TIR', 'IS3EU'),
    'DNA/Kolobok-T2': ('DNAtransposons Subclass1', 'TIR', 'Kolobok'),
    'DNA/Maverick': ('DNAtransposons Subclass2', 'Maverick', 'Maverick'),
    'DNA/MAVERICK': ('DNAtransposons Subclass2', 'Maverick', 'Maverick'),
    'DNA/Merlin': ('DNAtransposons Subclass1', 'TIR', 'Merlin'),
    'DNA/MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/MuDR': ('DNAtransposons Subclass1', 'TIR', 'MULE'),
    'DNA/MULE': ('DNAtransposons Subclass1', 'TIR', 'MULE'),
    'DNA/MULE-MuDR': ('DNAtransposons Subclass1', 'TIR', 'MULE'),
    'DNA/MULE-NOF': ('DNAtransposons Subclass1', 'TIR', 'MULE'),
    'DNA/Mutator': ('DNAtransposons Subclass1', 'TIR', 'Mutator'),
    'DNA/Mutator_MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/nMITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'DNA/Novosib': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'DNA/P': ('DNAtransposons Subclass1', 'TIR', 'P'),
    'DNA/PIF-Harbinger': ('DNAtransposons Subclass1', 'TIR', 'PIF-Harbinger'),
    'DNA/PIF-ISL2EU': ('DNAtransposons Subclass1', 'TIR', 'PIF-Harbinger'),
    'DNA/PiggyBac': ('DNAtransposons Subclass1', 'TIR', 'PiggyBac'),
    'DNA/Sola-2': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'DNA/Sola-3': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'DNA/TcMar': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-ISRm11': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-Mariner': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar_MITE': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-Stowaway': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-Tc1': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-Tc2': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TcMar-Tigger': ('DNAtransposons Subclass1', 'TIR', 'Tc1-mariner'),
    'DNA/TIR': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'DNA/TIR/HAT': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'DNA/TIR/P': ('DNAtransposons Subclass1', 'TIR', 'P'),
    'DNA/TIR/PIFHARBINGER': ('DNAtransposons Subclass1', 'TIR', 'PIF-Harbinger'),
    'DNA/TIR/PIGGYBAC': ('DNAtransposons Subclass1', 'TIR', 'PiggyBac'),
    'DNA/Zator': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'DNA/Zisupton': ('DNAtransposons Subclass1', 'TIR', 'Unknown TIR'),
    'Endogenous': (np.nan, np.nan, np.nan),
    'Gypsy': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'Harbinger': ('DNAtransposons Subclass1', 'TIR', 'PIF-Harbinger'),
    'hAT': ('DNAtransposons Subclass1', 'TIR', 'hAT'),
    'L1': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'LINE': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/AP': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/CR1': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/CR1-Zenon': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/CRE': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/CRE-Ambal': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Dong-R4': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Gypsy': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/I': ('Retrotransposons Autonomous', 'LINE', 'Jockey'),
    'LINE/I-Jockey': ('Retrotransposons Autonomous', 'LINE', 'Jockey'),
    'LINE/Jockey': ('Retrotransposons Autonomous', 'LINE', 'Jockey'),
    'LINE/L1': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'LINE/L1-DRE': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'LINE/L1-Tx1': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'LINE/L1_Tx1': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'LINE/L2': ('Retrotransposons Autonomous', 'LINE', 'L2'),
    'LINE/LINE_FR': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Penelope': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Proto2': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/R1': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/R1-LOA': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/R2': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/R2-NeSL': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/retrotransposon': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Rex-Babar': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTE': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTE-BovB': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTE_BovB': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTE-RTE': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTE-X': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/RTEX': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'LINE/Tad1': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'Low_complexity': ('Other', np.nan, np.nan),
    'LTR': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/5-bp': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/BEL': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Bhikhari': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Caulimovirus': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Copia': ('Retrotransposons Autonomous', 'LTR', 'Copia'),
    'LTR/COPIA': ('Retrotransposons Autonomous', 'LTR', 'Copia'),
    'LTR/DIRS': ('Retrotransposons Autonomous', 'LTR', 'DIRS'),
    'LTR/DNA/CACTA': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/CACTA_MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/Harbinger_MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/hAT': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/hAT_MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/Mutator_MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/DNA/TcMar_MITE': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV1': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV4': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV-Foamy': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVK': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVL': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVL-MaLR': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Foamy': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Gypsy': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/GYPSY': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Gypsy-Cigr': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Gypsy_Cigr': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Gypsy-Gmr1': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Gypsy_Gmr1': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/L1': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Ngaro': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Pao': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/RT': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/RT_LTR': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/unknown': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'Nonautonomous': ('DNAtransposons Unknown', np.nan, np.nan),
    'Non-LTR': ('Retrotransposons Unknown', np.nan, np.nan),
    'P': ('DNAtransposons Subclass1', 'TIR', 'P'),
    'Penelope': ('Retrotransposons Autonomous', 'PLE', 'Penelope'),
    'PENT_PU': (np.nan, np.nan, np.nan),
    'piggyBac': ('DNAtransposons Subclass1', 'TIR', 'PiggyBac'),
    'PLE': ('Retrotransposons Autonomous', 'PLE', 'Penelope'),
    'PLE/Chlamys': ('Retrotransposons Autonomous', 'PLE', 'Penelope'),
    'PLE/Naiad': ('Retrotransposons Autonomous', 'PLE', 'Penelope'),
    'Pseudogene': (np.nan, np.nan, np.nan),
    'RC/Helitron': ('DNAtransposons Subclass2', 'Helitron', 'Helitron'),
    'Retroposon': ('Retrotransposons Unknown', np.nan, np.nan),
    'Retroposon/L1-dep': ('Retrotransposons Autonomous', 'LINE', 'L1'),
    'RP5S': ('Other', np.nan, np.nan),
    'rRNA': ('Other', np.nan, np.nan),
    'RT': ('Other', np.nan, np.nan),
    'RTE': ('Other', np.nan, np.nan),
    'RTEX': ('Other', np.nan, np.nan),
    'RT_nLTR_LINE_RTE': ('Other', np.nan, np.nan),
    'RT_nLTR_SINE_tRNA': ('Other', np.nan, np.nan),
    'SAT': ('Other', np.nan, np.nan),
    'Satellite': ('Other', np.nan, np.nan),
    'Satellite/centr': ('Other', np.nan, np.nan),
    'Simple_repeat': ('Other', np.nan, np.nan),
    'SINE': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE?': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/5S': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/5S-Deu-L2': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/Alu': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/B2': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/B4': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/Deu': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/ID': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/MIR': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-Core-RTE': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-Deu': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-Deu-L2': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA_Lys': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-Meta': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-RTE': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-V': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-V-Core-L2': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/tRNA-V-CR1': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/U': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/U-L1': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'snRNA': ('Other', np.nan, np.nan),
    'tRNA': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'unknown': ('Unknown', np.nan, np.nan),
    'Unknown': ('Unknown', np.nan, np.nan),
    'Unspecified': ('Unknown', np.nan, np.nan)
}


def _is_real_category(value):
    return isinstance(value, str) and bool(value.strip())


ALL_CLASSES = ALL_CLASSES | {row[0] for row in CLASSIFICATION_MAP.values() if _is_real_category(row[0])}
ALL_ORDERS = ALL_ORDERS | {row[1] for row in CLASSIFICATION_MAP.values() if _is_real_category(row[1])}
ALL_SUPERFAMILIES = ALL_SUPERFAMILIES | {row[2] for row in CLASSIFICATION_MAP.values() if _is_real_category(row[2])}


class DataProcessor:
    """Class for processing dnaPipeTE data files."""
    
    def __init__(self, species_dict):
        """Initialize with a species dictionary."""
        self.species_dict = species_dict
        
    def apply_classifications(self, chunk, species):
        """Apply classifications using vectorized operations."""
        # Initialize classification columns with proper length
        n_rows = len(chunk)
        chunk['Species'] = species
        chunk['Class'] = pd.Series([None] * n_rows, dtype='category')
        chunk['Order'] = pd.Series([None] * n_rows, dtype='category')
        chunk['Superfamily'] = pd.Series([None] * n_rows, dtype='category')
        
        # Set categories for each column
        chunk['Class'] = chunk['Class'].cat.set_categories(list(ALL_CLASSES))
        chunk['Order'] = chunk['Order'].cat.set_categories(list(ALL_ORDERS))
        chunk['Superfamily'] = chunk['Superfamily'].cat.set_categories(list(ALL_SUPERFAMILIES))
        
        # Apply classifications efficiently using vectorized operations
        for class_key, (class_val, order_val, superfam_val) in CLASSIFICATION_MAP.items():
            mask = chunk['RM_classification'] == class_key
            if mask.any():
                if class_val in ALL_CLASSES:
                    chunk.loc[mask, 'Class'] = class_val
                if order_val in ALL_ORDERS:
                    chunk.loc[mask, 'Order'] = order_val
                if superfam_val in ALL_SUPERFAMILIES:
                    chunk.loc[mask, 'Superfamily'] = superfam_val
        
        return chunk
    
    def process_file(self, filename):
        """Process a single file, adding classification columns."""
        try:
            file_path = INPUT_DIR / filename
            output_path = INTERIM_DIR / filename
            
            # Extract SRX value and map to species
            srx = filename.replace("_reads_per_component_and_annotation", "")
            species = self.species_dict.get(srx[:11], "Unknown")
            if species == "Unknown":
                raise ValueError(f"No species mapping found for {srx[:11]} from {filename}")
            
            # Read the file
            df = pd.read_csv(file_path, sep=" ", header=None, names=[
                "#reads", "aligned_bases", "dnaPipeTE_contig_name", "RM_hit_length_bp",
                "RM_annotation", "RM_classification", "hitlength_contiglength"
            ])
            
            # Apply classifications and add species
            df = self.apply_classifications(df, species)
            
            # Save with same format as input, but with additional columns
            df.to_csv(output_path, sep=" ", index=False, header=False)
            
            row_count = len(df)
            logger.info(f"Processed {filename}: {row_count} rows, Species: {species}")
            
            # Clean up memory
            del df
            gc.collect()
            return {"filename": filename, "ok": True, "rows": row_count}
            
        except Exception as e:
            logger.error(f"Error processing file {filename}: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return {"filename": filename, "ok": False, "error": str(e)}
            
    def process_all_files(self):
        """Process all files in parallel."""
        if not INPUT_DIR.exists():
            raise FileNotFoundError(f"dnaPipeTE input directory not found: {INPUT_DIR}")

        for old_file in INTERIM_DIR.glob("*_reads_per_component_and_annotation"):
            old_file.unlink()

        files = sorted(
            f for f in os.listdir(INPUT_DIR)
            if f.endswith("_reads_per_component_and_annotation")
        )
        if not files:
            raise RuntimeError(f"No dnaPipeTE input files found in {INPUT_DIR}")
        
        logger.info(f"Processing {len(files)} files...")
        
        # Use maximum number of cores while leaving some resources for the system
        num_processes = min(len(files), max(1, os.cpu_count() - 1))
        
        # Process files in parallel with progress bar
        with Pool(processes=num_processes) as pool:
            results = list(tqdm(pool.imap(self.process_file, files), total=len(files)))

        failed = [result for result in results if not result.get("ok")]
        if failed:
            failed_names = ", ".join(result["filename"] for result in failed)
            raise RuntimeError(f"Failed to process {len(failed)} dnaPipeTE files: {failed_names}")

        processed_files = [INTERIM_DIR / filename for filename in files]
        missing_outputs = [path.name for path in processed_files if not path.exists()]
        if missing_outputs:
            raise RuntimeError(f"Missing processed interim files after dnaPipeTE run: {missing_outputs}")
        return files

def merge_files(input_dir, output_file, expected_files=None):
    """
    Merge all processed files into a single consolidated file incrementally.
    
    Args:
        input_dir: Path to the directory containing processed files
        output_file: Path where the merged file will be saved
    """
    logger.info("Starting incremental file merge process...")
    
    # Get list of all processed files
    if expected_files is None:
        files = sorted(
            f for f in os.listdir(input_dir)
            if f.endswith("_reads_per_component_and_annotation")
        )
    else:
        files = list(expected_files)
    if not files:
        raise RuntimeError(f"No processed dnaPipeTE files found in {input_dir}")
    
    logger.info(f"Found {len(files)} files to merge")
    
    first_file = True
    total_rows = 0
    
    # Define column names for reading and writing
    column_names = [
        "#reads", "aligned_bases", "dnaPipeTE_contig_name", "RM_hit_length_bp",
        "RM_annotation", "RM_classification", "hitlength_contiglength",
        "Species", "Class", "Order", "Superfamily"
    ]
    output_column_names = column_names + ['Source']

    # Process each file
    failures = []
    for filename in tqdm(files, desc="Merging files incrementally"):
        try:
            file_path = Path(input_dir) / filename
            if not file_path.exists():
                raise FileNotFoundError(file_path)
            
            # Read the file - assume interim files have no header
            df = pd.read_csv(file_path, sep=" ", header=None, names=column_names)
            
            # Add source file information
            srx = filename.replace("_reads_per_component_and_annotation", "")
            df['Source'] = srx
            
            # Append to the output file
            # Write header only for the first file
            df.to_csv(
                output_file, 
                mode='a', 
                index=False, 
                header=first_file, 
                columns=output_column_names # Ensure consistent column order
            )
            
            total_rows += len(df)
            
            if first_file:
                first_file = False
            
            # Clean up to manage memory
            del df
            gc.collect()
            
        except Exception as e:
            logger.error(f"Error processing {filename} during merge: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            failures.append(filename)

    if failures:
        raise RuntimeError(f"Failed to merge {len(failures)} processed dnaPipeTE files: {failures}")
    if total_rows == 0:
        raise RuntimeError("dnaPipeTE merge wrote zero rows")
    
    logger.info(f"Finished merging. Total rows in merged file: {total_rows:,}")
    logger.info(f"Merged file saved to: {output_file}")
    
    # Return the path to the merged file instead of the DataFrame
    return output_file

def create_visualization_datasets(merged_data_path, output_data_dir):
    """
    Create breakdown datasets (Class, Order, Superfamily) from the merged data.
    
    Args:
        merged_data_path: Path to the merged data CSV file.
        output_data_dir: Directory to save the breakdown datasets.
    """
    logger.info("Creating breakdown datasets...")

    # Define expected columns to read from the merged file
    required_columns = ['Species', 'aligned_bases', 'Class', 'Order', 'Superfamily']
    
    for category in ['Class', 'Order', 'Superfamily']:
        logger.info(f"Processing {category} breakdown...")
        # Use a dictionary to accumulate sums, reducing memory footprint
        category_aggregator = {}  # {species: {category_value: sum_bases}}
        output_file = Path(output_data_dir) / f"dnaPipeTE_{category.lower()}_breakdown.csv"
        
        try:
            # Process in chunks to manage memory
            # Specify columns to read, reducing memory further
            chunk_iter = pd.read_csv(
                merged_data_path, 
                chunksize=1000000, # Adjust chunksize as needed
                usecols=lambda col: col in required_columns # Read only necessary columns
            ) 
            
            for chunk in tqdm(chunk_iter, desc=f"{category} Chunks"):
                # Ensure the category column exists and handle potential missing values
                if category not in chunk.columns:
                    logger.warning(f"Category column '{category}' not found in chunk. Skipping chunk for this category.")
                    continue
                
                # Fill NA in category column to avoid issues during grouping
                chunk[category] = chunk[category].fillna('Unknown').astype(str) 
                
                # Group within the chunk and sum bases
                grouped_chunk = chunk.groupby(['Species', category])['aligned_bases'].sum()
                
                # Update the main aggregator dictionary
                for index, total_bases in grouped_chunk.items():
                    species, cat_value = index
                    if species not in category_aggregator:
                        category_aggregator[species] = {}
                    category_aggregator[species][cat_value] = category_aggregator[species].get(cat_value, 0) + total_bases
                
                # Clean up chunk data
                del chunk, grouped_chunk
                gc.collect()

            # Convert the aggregated dictionary to a DataFrame for pivoting/saving
            if not category_aggregator:
                raise RuntimeError(f"No data aggregated for category {category}")

            # Create DataFrame from the dictionary
            final_df = pd.DataFrame.from_dict(category_aggregator, orient='index').fillna(0)

            # Ensure all possible categories (columns) for this level are present, filled with 0
            # Determine the set of all possible columns based on the classification level
            if category == 'Class':
                all_possible_cols = ALL_CLASSES
            elif category == 'Order':
                all_possible_cols = ALL_ORDERS
            elif category == 'Superfamily':
                all_possible_cols = ALL_SUPERFAMILIES
                # --- Superfamily specific: Exclude 'Unknown TIR' ---
                if 'Unknown TIR' in all_possible_cols:
                    logger.info("Excluding 'Unknown TIR' from Superfamily breakdown.")
                    # Make a copy to avoid modifying the original set
                    all_possible_cols = all_possible_cols.copy() 
                    all_possible_cols.remove('Unknown TIR')
                # Also remove the column from the dataframe if it exists
                if 'Unknown TIR' in final_df.columns:
                    final_df = final_df.drop(columns=['Unknown TIR'])
            else: # Should not happen with current loop, but good practice
                all_possible_cols = final_df.columns 

            # Add missing columns (relevant after potential exclusion) and fill with 0
            current_cols = set(final_df.columns)
            missing_cols = list(all_possible_cols - current_cols)
            for col in missing_cols:
                final_df[col] = 0
            
            # Ensure columns are in a consistent order (using sorted list of relevant possible columns)
            # Filter all_possible_cols to only include those actually present or added
            relevant_cols = sorted([col for col in all_possible_cols if col in final_df.columns])
            final_df = final_df[relevant_cols]

            # --- Convert absolute values to percentages ---
            # Calculate row sums (total aligned bases per species for this category level)
            row_sums = final_df.sum(axis=1)
            
            # Avoid division by zero for species with no aligned bases in this category
            # Create a percentage DataFrame, initializing with zeros
            perc_df = pd.DataFrame(0.0, index=final_df.index, columns=final_df.columns)
            
            # Calculate percentages only for rows with non-zero sums
            non_zero_rows = row_sums > 0
            if non_zero_rows.any():
                 perc_df.loc[non_zero_rows] = final_df.loc[non_zero_rows].div(row_sums[non_zero_rows], axis=0) * 100
            
            # Save the final percentage breakdown table
            perc_df.to_csv(output_file)
            logger.info(f"Created {category} percentage breakdown at {output_file}")
            
            # Clean up dataframes
            del final_df, perc_df
            gc.collect()

        except FileNotFoundError:
             logger.error(f"Error: Merged data file not found at {merged_data_path}")
             raise
        except pd.errors.EmptyDataError:
            logger.error(f"Error: Merged data file {merged_data_path} is empty.")
            raise
        except KeyError as e:
             logger.error(f"Error processing {category} breakdown: Missing column {e}. Ensure '{category}', 'Species', and 'aligned_bases' columns exist in {merged_data_path}.")
             raise
        except Exception as e:
            logger.error(f"Error processing {category} breakdown: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            raise

    logger.info("Finished creating breakdown datasets.")
    # No return value needed as the function now only writes files
    # return class_df, prop_df # Removed original return

def load_species_dict(lookup_path):
    """Load the species lookup dictionary from file."""
    logger.info("Loading species lookup table...")
    try:
        lookup_df = pd.read_csv(lookup_path, sep='\t', usecols=['SRA_Accension', 'Species'])
        species_dict = dict(zip(
            lookup_df['SRA_Accension'].astype(str).str[:11],
            lookup_df['Species']
        ))
        logger.info(f"Loaded {len(species_dict)} species mappings")
        return species_dict
    except Exception as e:
        logger.error(f"Error loading species lookup: {str(e)}")
        return {}

def main():
    """Main function to run the entire processing pipeline."""
    # Step 1: Load species dictionary
    species_dict = load_species_dict(LOOKUP_PATH)
    if not species_dict:
        raise RuntimeError(f"No species mappings loaded from {LOOKUP_PATH}")
    
    # Step 2: Process individual files
    processor = DataProcessor(species_dict)
    processed_files = processor.process_all_files()
    
    # Step 3: Merge processed files
    # Use the new output directory and filename format
    merged_data_path = OUTPUT_DATA_DIR / "dnaPipeTE_merged_classifications.csv"
    
    # Ensure the target file is removed before starting append mode
    if merged_data_path.exists():
        merged_data_path.unlink()
        logger.info(f"Removed existing merged file: {merged_data_path}")
        
    final_merged_path = merge_files(INTERIM_DIR, merged_data_path, expected_files=processed_files)
    
    # Step 4: Create visualization datasets using the path and new output dir
    create_visualization_datasets(final_merged_path, OUTPUT_DATA_DIR)
    
    logger.info("Processing complete!")

if __name__ == "__main__":
    main()

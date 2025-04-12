import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import gc
from multiprocessing import Pool
import logging
import time

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('process_dnaPipeTE')

# Define paths
base_dir = os.path.dirname(os.path.abspath(__file__))
raw_dir = os.path.join(base_dir, "data", "raw", "dnaPipeTE_outputs")
lookup_path = os.path.join(base_dir, "data", "raw", "lookup", "lookup_table.txt")
interim_dir = os.path.join(base_dir, "data", "interim")
pivot_dir = os.path.join(interim_dir, "pivot_tables")

# Create output directories
os.makedirs(interim_dir, exist_ok=True)
os.makedirs(pivot_dir, exist_ok=True)

# Get all possible categories for the classification columns
all_classes = {
    'DNAtransposons Subclass1', 'DNAtransposons Subclass2', 'DNAtransposons Unknown',
    'Retrotransposons Autonomous', 'Retrotransposons Non-autonomous', 'Retrotransposons Unknown',
    'Other', 'Unknown'
}

all_orders = {
    'TIR', 'YR', 'Helitron', 'Maverick', 'LINE', 'LTR', 'DIRS', 'PLE', 'SINE'
}

all_superfamilies = {
    'Academ', 'CACTA', 'Chapaev', 'Cyrypton', 'Dada', 'EnSpm', 'Ginger',
    'PIF-Harbinger', 'hAT', 'Unknown TIR', 'P', 'PiggyBac', 'Tc1-mariner',
    'MULE', 'Mutator', 'Helitron', 'Maverick', 'L1', 'Jockey', 'Penelope',
    'Gypsy', 'DIRS', 'tRNA'
}

def apply_classifications(chunk, species):
    """Apply classifications using vectorized operations"""
    # Initialize classification columns with proper length
    n_rows = len(chunk)
    chunk['Species'] = species
    chunk['Class'] = pd.Series([None] * n_rows, dtype='category')
    chunk['Order'] = pd.Series([None] * n_rows, dtype='category')
    chunk['Superfamily'] = pd.Series([None] * n_rows, dtype='category')
    
    # Set categories for each column
    chunk['Class'] = chunk['Class'].cat.set_categories(list(all_classes))
    chunk['Order'] = chunk['Order'].cat.set_categories(list(all_orders))
    chunk['Superfamily'] = chunk['Superfamily'].cat.set_categories(list(all_superfamilies))
    
    # Apply classifications
    for class_key, (class_val, order_val, superfam_val) in classification_map.items():
        mask = chunk['RM_classification'] == class_key
        if mask.any():
            if class_val in all_classes:
                chunk.loc[mask, 'Class'] = class_val
            if order_val in all_orders:
                chunk.loc[mask, 'Order'] = order_val
            if superfam_val in all_superfamilies:
                chunk.loc[mask, 'Superfamily'] = superfam_val
    
    return chunk

def process_file(filename):
    """Process a single file with improved error handling"""
    try:
        file_path = os.path.join(raw_dir, filename)
        output_path = os.path.join(interim_dir, filename + "_processed")
        
        # Extract SRX value and map to species
        srx = filename.replace("_reads_per_component_and_annotation", "")
        species = species_dict.get(srx[:11], "Unknown")
        
        logger.info(f"Processing {filename} for species {species}")
        
        # Read the file
        df = pd.read_csv(file_path, sep=" ", header=None, names=[
            "#reads", "aligned_bases", "dnaPipeTE_contig_name", "RM_hit_length_bp",
            "RM_annotation", "RM_classification", "hitlength_contiglength"
        ])
        
        # Apply classifications and add species
        df = apply_classifications(df, species)
        
        # Save with same format as input, but with additional columns
        df.to_csv(output_path, sep=" ", index=False)
        
        logger.info(f"Processed {filename} - Total rows: {len(df)}")
        
        # Free up memory
        del df
        gc.collect()
        
        return True
        
    except Exception as e:
        logger.error(f"Error processing file {filename}: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def create_pivot(df, category):
    """Create a pivot table for a specific category (Class, Order, Superfamily)"""
    pivot = df.pivot_table(
        values='RM_hit_length_bp',
        index='Species',
        columns=category,
        aggfunc='sum',
        fill_value=0
    )
    
    # Add total column
    pivot['Total'] = pivot.sum(axis=1)
    
    return pivot

def merge_files():
    """Merge processed files and create pivot tables"""
    logger.info("Merging processed files...")
    
    all_dfs = []
    processed_dir = interim_dir
    
    # Get all processed files
    processed_files = [f for f in os.listdir(processed_dir) 
                     if f.endswith("_processed")]
    
    logger.info(f"Found {len(processed_files)} processed files to merge")
    
    for file in tqdm(processed_files, desc="Reading files"):
        try:
            df = pd.read_csv(os.path.join(processed_dir, file), sep=" ")
            # Keep only necessary columns to save memory
            df = df[['Species', 'Class', 'Order', 'Superfamily', 'RM_hit_length_bp']]
            all_dfs.append(df)
        except Exception as e:
            logger.error(f"Error reading file {file}: {str(e)}")
    
    # Merge all dataframes
    logger.info("Concatenating dataframes...")
    merged_df = pd.concat(all_dfs, ignore_index=True)
    del all_dfs
    gc.collect()
    
    # Create pivot tables
    logger.info("Creating pivot tables...")
    
    # Class pivot
    class_pivot = create_pivot(merged_df, 'Class')
    class_pivot.to_csv(os.path.join(pivot_dir, 'class_breakdown.csv'))
    logger.info(f"Saved class breakdown to {pivot_dir}/class_breakdown.csv")
    
    # Order pivot
    order_pivot = create_pivot(merged_df, 'Order')
    order_pivot.to_csv(os.path.join(pivot_dir, 'order_breakdown.csv'))
    logger.info(f"Saved order breakdown to {pivot_dir}/order_breakdown.csv")
    
    # Superfamily pivot
    superfamily_pivot = create_pivot(merged_df, 'Superfamily')
    superfamily_pivot.to_csv(os.path.join(pivot_dir, 'superfamily_breakdown.csv'))
    logger.info(f"Saved superfamily breakdown to {pivot_dir}/superfamily_breakdown.csv")
    
    # Create DNA vs Retro breakdown
    dna_vs_retro = pd.DataFrame(index=class_pivot.index)
    dna_vs_retro['DNA'] = class_pivot[[col for col in class_pivot.columns if col.startswith('DNAtransposons')]].sum(axis=1)
    dna_vs_retro['Retrotransposons'] = class_pivot[[col for col in class_pivot.columns if col.startswith('Retrotransposons')]].sum(axis=1)
    dna_vs_retro.to_csv(os.path.join(pivot_dir, 'dna_vs_retro.csv'))
    logger.info(f"Saved DNA vs Retro breakdown to {pivot_dir}/dna_vs_retro.csv")
    
    # Create Known vs Unknown breakdown
    known_vs_unknown = pd.DataFrame(index=class_pivot.index)
    known_vs_unknown['Known'] = class_pivot[[col for col in class_pivot.columns if col not in ['Unknown', 'Total']]].sum(axis=1)
    known_vs_unknown['Unknown'] = class_pivot['Unknown']
    known_vs_unknown.to_csv(os.path.join(pivot_dir, 'known_vs_unknown.csv'))
    logger.info(f"Saved Known vs Unknown breakdown to {pivot_dir}/known_vs_unknown.csv")
    
    return merged_df

def process_files_sequentially():
    """Process files one at a time to avoid overwhelming the system"""
    files = [f for f in os.listdir(raw_dir) 
             if f.endswith("_reads_per_component_and_annotation")]
    
    logger.info(f"Processing {len(files)} files sequentially...")
    
    success_count = 0
    for filename in tqdm(files, desc="Processing files"):
        success = process_file(filename)
        if success:
            success_count += 1
        # Small delay to prevent system from being overwhelmed
        time.sleep(0.1)
    
    logger.info(f"Successfully processed {success_count} out of {len(files)} files")
    return success_count == len(files)

def main():
    """Run the main workflow"""
    logger.info("Starting dnaPipeTE raw data processing")
    
    # Load lookup table
    global species_dict, classification_map
    logger.info(f"Loading lookup table from {lookup_path}...")
    
    try:
        lookup_df = pd.read_csv(lookup_path, sep='\t')
        species_dict = dict(zip(
            lookup_df['SRA_Accension'].astype(str).str[:11],
            lookup_df['Species']
        ))
        logger.info(f"Loaded {len(species_dict)} species mappings")
    except Exception as e:
        logger.error(f"Error loading lookup table: {str(e)}")
        return False
    
    # Process files one by one
    success = process_files_sequentially()
    if not success:
        logger.error("Some files failed to process")
        return False
    
    # Merge files and create pivot tables
    merge_files()
    
    logger.info("Processing completed successfully")
    return True

# Your classification_map dictionary
classification_map = {
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
    'LTR/ERV': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV1': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV4': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERV-Foamy': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVK': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVL': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/ERVL-MaLR': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Gypsy': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/GYPSY': ('Retrotransposons Autonomous', 'LTR', 'Gypsy'),
    'LTR/Ngaro': ('Retrotransposons Autonomous', 'LTR', 'DIRS'),
    'LTR/Pao': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/parasitic': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Unknown': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'LTR/Viper': ('Retrotransposons Autonomous', 'LTR', np.nan),
    'MITE': ('DNAtransposons Subclass2', 'Transposon Derivatives', 'MITE'),
    'None': ('Unknown', np.nan, np.nan),
    'Other': ('Other', np.nan, np.nan),
    'Penelope': ('Retrotransposons Autonomous', 'PLE', 'Penelope'),
    'RC': ('DNAtransposons Subclass2', 'Helitron', 'Helitron'),
    'RC/Helitron': ('DNAtransposons Subclass2', 'Helitron', 'Helitron'),
    'Retroposon': ('Retrotransposons Non-autonomous', np.nan, np.nan),
    'RTEX': ('Retrotransposons Autonomous', 'LINE', np.nan),
    'rRNA': ('Other', np.nan, np.nan),
    'SAT': ('Other', np.nan, np.nan),
    'SAT/NA': ('Other', np.nan, np.nan),
    'Simple_repeat': ('Other', np.nan, np.nan),
    'SINE': ('Retrotransposons Non-autonomous', 'SINE', np.nan),
    'SINE/5S': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/Alu': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/B2': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/B4': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/Core': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/Core-RTE': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/ID': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/MIR': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/RTE': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA-Ala': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA-Alu': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA-Cys': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA-Lys': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'SINE/tRNA-RTE': ('Retrotransposons Non-autonomous', 'SINE', 'tRNA'),
    'snRNA': ('Other', np.nan, np.nan),
    'tRNA': ('Other', np.nan, np.nan),
    'Unknown': ('Unknown', np.nan, np.nan)
}

if __name__ == "__main__":
    main()

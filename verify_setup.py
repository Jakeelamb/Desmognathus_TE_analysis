#!/usr/bin/env python3
"""
Verification script for Desmognathus TE Analysis environment.
Checks for required packages, input data, and directory structure.
"""

import sys
import os
from pathlib import Path

# Color codes for terminal output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_header(text):
    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}{text}{RESET}")
    print(f"{BLUE}{'='*60}{RESET}")

def check_mark(status):
    return f"{GREEN}✓{RESET}" if status else f"{RED}✗{RESET}"

def print_status(item, status, note=""):
    mark = check_mark(status)
    note_str = f" {YELLOW}({note}){RESET}" if note else ""
    print(f"  {mark} {item}{note_str}")

def check_python_packages():
    """Check if required Python packages are installed."""
    print_header("Python Packages")
    
    # Map display name to actual import name
    required_packages = [
        ('pandas', 'pandas'),
        ('numpy', 'numpy'),
        ('matplotlib', 'matplotlib'),
        ('seaborn', 'seaborn'),
        ('biopython', 'Bio'),
        ('tqdm', 'tqdm'),
        ('yaml', 'yaml')
    ]
    
    all_ok = True
    for display_name, import_name in required_packages:
        try:
            __import__(import_name)
            print_status(display_name, True)
        except ImportError:
            print_status(display_name, False, "Not installed")
            all_ok = False
    
    return all_ok

def check_r_availability():
    """Check if R and Rscript are available."""
    print_header("R Environment")
    
    r_available = os.system("which Rscript > /dev/null 2>&1") == 0
    print_status("Rscript available", r_available)
    
    if r_available:
        # Check R version
        import subprocess
        try:
            result = subprocess.run(
                ["Rscript", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            version = result.stderr.split('\n')[0] if result.stderr else "Unknown"
            print(f"    Version: {version}")
        except:
            pass
    
    return r_available

def check_directory_structure():
    """Check if required directories exist."""
    print_header("Directory Structure")
    
    base_dir = Path("/home/jake/Projects/Desmognathus_TE")
    
    required_dirs = {
        "Input Data": [
            "input_data",
            "input_data/dnaPipeTE",
            "input_data/repeatmasker",
            "input_data/phylogeny",
        ],
        "Output": [
            "results",
            "results/data",
            "results/figures",
            "interim"
        ],
        "Scripts": [
            "scripts/processing",
            "scripts/visualization"
        ]
    }
    
    all_ok = True
    for category, dirs in required_dirs.items():
        print(f"\n  {category}:")
        for dir_path in dirs:
            full_path = base_dir / dir_path
            exists = full_path.exists()
            print_status(dir_path, exists)
            if not exists:
                all_ok = False
    
    return all_ok

def check_input_files():
    """Check for critical input files."""
    print_header("Critical Input Files")
    
    base_dir = Path("/home/jake/Projects/Desmognathus_TE")
    
    files_to_check = {
        "input_data/lookup_table.txt": "Species lookup table",
        "input_data/phylogeny/desmo900dated_test.tre": "Phylogenetic tree",
    }
    
    for file_path, description in files_to_check.items():
        full_path = base_dir / file_path
        exists = full_path.exists()
        print_status(f"{description} ({file_path})", exists)
    
    # Check for data in directories
    print("\n  Data Directory Contents:")
    
    data_dirs = {
        "input_data/dnaPipeTE": "dnaPipeTE files",
        "input_data/repeatmasker": ".align files"
    }
    
    for dir_path, description in data_dirs.items():
        full_path = base_dir / dir_path
        if full_path.exists():
            file_count = len(list(full_path.iterdir()))
            status = file_count > 0
            note = f"{file_count} files" if status else "Empty"
            print_status(f"{description} ({dir_path})", status, note)
        else:
            print_status(f"{description} ({dir_path})", False, "Directory not found")

def check_scripts():
    """Check that processing scripts exist and are readable."""
    print_header("Processing Scripts")
    
    base_dir = Path("/home/jake/Projects/Desmognathus_TE")
    
    scripts = [
        ("scripts/processing/dnaPipe.py", "dnaPipeTE processing"),
        ("scripts/processing/repeatmask.py", "RepeatMasker processing"),
        ("scripts/processing/ec.py", "Ectopic recombination"),
        ("scripts/processing/divergence.py", "Divergence analysis"),
        ("scripts/processing/diversity_stats.py", "Canonical diversity metrics"),
        ("scripts/processing/pca.R", "PCA analysis"),
        ("scripts/processing/clean_tree_phylo.R", "Phylogeny cleaning"),
        ("scripts/processing/phylogenetic_pca_analysis.R", "Phylogenetic PCA"),
    ]
    
    all_ok = True
    for script_path, description in scripts:
        full_path = base_dir / script_path
        exists = full_path.exists() and full_path.is_file()
        print_status(f"{description} ({script_path})", exists)
        if not exists:
            all_ok = False
    
    return all_ok

def check_existing_outputs():
    """Check what output files already exist."""
    print_header("Existing Output Files")
    
    base_dir = Path("/home/jake/Projects/Desmognathus_TE")
    results_data = base_dir / "results/data"
    
    if not results_data.exists():
        print("  Results directory doesn't exist yet")
        return
    
    csv_files = list(results_data.glob("*.csv"))
    tre_files = list(results_data.glob("*.tre"))
    
    print(f"\n  CSV files: {len(csv_files)}")
    print(f"  Tree files: {len(tre_files)}")
    
    if csv_files or tre_files:
        print(f"\n  {YELLOW}Note: Some analyses may have been run already{RESET}")

def main():
    """Run all verification checks."""
    print(f"\n{BLUE}{'='*60}")
    print("Desmognathus TE Analysis - Setup Verification")
    print(f"{'='*60}{RESET}\n")
    
    checks = [
        ("Python Packages", check_python_packages),
        ("R Environment", check_r_availability),
        ("Directory Structure", check_directory_structure),
        ("Scripts", check_scripts),
    ]
    
    results = {}
    for check_name, check_func in checks:
        try:
            results[check_name] = check_func()
        except Exception as e:
            print(f"\n{RED}Error running {check_name}: {e}{RESET}")
            results[check_name] = False
    
    # Additional checks that don't return boolean
    try:
        check_input_files()
        check_existing_outputs()
    except Exception as e:
        print(f"\n{RED}Error checking files: {e}{RESET}")
    
    # Summary
    print_header("Summary")
    
    all_passed = all(results.values())
    
    for check_name, passed in results.items():
        print_status(check_name, passed)
    
    if all_passed:
        print(f"\n{GREEN}✓ All critical checks passed!{RESET}")
        print(f"\n{BLUE}You can now run the analysis scripts.{RESET}")
        print(f"{BLUE}See README.md for detailed instructions.{RESET}\n")
        return 0
    else:
        print(f"\n{YELLOW}⚠ Some checks failed.{RESET}")
        print(f"{YELLOW}Review the output above and install missing components.{RESET}\n")
        return 1

if __name__ == "__main__":
    sys.exit(main())

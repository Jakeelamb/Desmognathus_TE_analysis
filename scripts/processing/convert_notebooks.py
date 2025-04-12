#!/usr/bin/env python3

"""
convert_notebooks.py

This script converts Jupyter notebooks (.ipynb) to Python scripts (.py)
with proper formatting and documentation.

Usage:
  python scripts/processing/convert_notebooks.py
"""

import os
import sys
import subprocess
import json
import re
from pathlib import Path

def convert_notebook(notebook_path, output_path=None):
    """
    Convert a Jupyter notebook to a Python script.
    
    Args:
        notebook_path: Path to the notebook file
        output_path: Path where the Python script should be saved (optional)
    
    Returns:
        Path to the output file
    """
    # Default output path is same as input but with .py extension
    if output_path is None:
        output_path = str(notebook_path).replace('.ipynb', '.py')
    
    print(f"Converting {notebook_path} -> {output_path}")
    
    # Read the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Extract notebook metadata
    metadata = notebook.get('metadata', {})
    kernelspec = metadata.get('kernelspec', {})
    kernel_name = kernelspec.get('display_name', 'Python')
    
    # Open output file for writing
    with open(output_path, 'w', encoding='utf-8') as f:
        # Add shebang and encoding
        f.write("#!/usr/bin/env python3\n")
        f.write("# -*- coding: utf-8 -*-\n\n")
        
        # Add header with notebook info
        notebook_name = os.path.basename(notebook_path)
        f.write(f'"""\n')
        f.write(f'This script was converted from the Jupyter notebook: {notebook_name}\n')
        f.write(f'Original kernel: {kernel_name}\n')
        f.write(f'Conversion date: {os.popen("date").read().strip()}\n')
        f.write(f'"""\n\n')
        
        # Process each cell
        for i, cell in enumerate(notebook.get('cells', [])):
            cell_type = cell.get('cell_type')
            source = cell.get('source', [])
            
            # Join source if it's a list
            if isinstance(source, list):
                source = ''.join(source)
                
            # Handle markdown cells as comments
            if cell_type == 'markdown':
                if source.strip():
                    f.write(f'# {"-"*78}\n')
                    f.write('# ' + source.replace('\n', '\n# ').rstrip('# ') + '\n')
                    f.write(f'# {"-"*78}\n\n')
            
            # Handle code cells
            elif cell_type == 'code':
                if source.strip():
                    f.write(source + '\n\n')
    
    print(f"Successfully converted to {output_path}")
    return output_path

def main():
    """
    Find and convert all Jupyter notebooks in the repository.
    """
    # Get the repository root directory (assuming the script is run from there or a subdirectory)
    repo_root = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], 
                                        universal_newlines=True).strip()
    if not repo_root:
        repo_root = "."
    
    # Find all .ipynb files
    print(f"Searching for Jupyter notebooks in {repo_root}")
    ipynb_files = []
    
    for root, dirs, files in os.walk(repo_root):
        for file in files:
            if file.endswith(".ipynb") and not file.endswith(".nbconvert.ipynb"):
                ipynb_files.append(os.path.join(root, file))
    
    print(f"Found {len(ipynb_files)} Jupyter notebooks")
    
    # Define output directories for different types of notebooks
    output_mappings = {
        os.path.join("notebooks", "Data_wrangle.ipynb"): 
            os.path.join("scripts", "processing", "data_wrangle.py"),
        os.path.join("Projects", "03_repeat_landscape", "read_per_comp_vizualization_allcode.ipynb"): 
            os.path.join("Projects", "03_repeat_landscape", "read_per_comp_visualization.py"),
    }
    
    # Convert each notebook
    for nb_path in ipynb_files:
        rel_path = os.path.relpath(nb_path, repo_root)
        output_path = output_mappings.get(rel_path)
        
        if output_path:
            # Use the predefined output path
            output_path = os.path.join(repo_root, output_path)
        else:
            # Default: replace ipynb with py in the same directory
            output_path = nb_path.replace('.ipynb', '.py')
            
        convert_notebook(nb_path, output_path)
        
        # Suggest removal of original .ipynb file
        print(f"Consider removing the original file: rm '{nb_path}'")
    
    print("\nConversion complete. All notebooks have been converted to Python scripts.")
    print("You may now want to update your .gitignore to exclude .ipynb files entirely.")

if __name__ == "__main__":
    main() 
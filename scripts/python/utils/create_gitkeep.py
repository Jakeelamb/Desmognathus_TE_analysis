#!/usr/bin/env python3
"""
Create .gitkeep files in essential directories to ensure the directory structure
is maintained in the git repository even if the directories are empty.
"""

import os
import yaml
from pathlib import Path

def main():
    # Load config file
    with open('config/paths.yaml', 'r') as file:
        config = yaml.safe_load(file)

    # Essential directory structure
    essential_dirs = [
        # Data directories
        'data/raw',
        'data/interim',
        'data/processed',
        'data/interim/pivot_tables',
        'data/raw/lookup',
        'data/processed/diversity',
        'data/processed/te_superfamily',

        # Results directories
        'results/figures/pca',
        'results/figures/diversity',
        'results/tables/pca',
        'results/tables/diversity',
    ]

    # Create .gitkeep files
    for directory in essential_dirs:
        dir_path = Path(directory)
        
        # Create directory if it doesn't exist
        if not dir_path.exists():
            os.makedirs(dir_path, exist_ok=True)
            print(f"Created directory: {dir_path}")
        
        # Create .gitkeep file
        gitkeep_path = dir_path / '.gitkeep'
        if not gitkeep_path.exists():
            with open(gitkeep_path, 'w') as f:
                pass  # Create empty file
            print(f"Created .gitkeep in: {dir_path}")

    print("Directory structure prepared for Git tracking.")

if __name__ == "__main__":
    main() 
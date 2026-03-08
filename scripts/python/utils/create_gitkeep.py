#!/usr/bin/env python3
"""Create .gitkeep files in essential directories."""

import os
from pathlib import Path

def main():
    project_root = Path(__file__).resolve().parents[3]

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
        dir_path = project_root / directory
        
        # Create directory if it doesn't exist
        if not dir_path.exists():
            os.makedirs(dir_path, exist_ok=True)
            print(f"Created directory: {dir_path.relative_to(project_root)}")
        
        # Create .gitkeep file
        gitkeep_path = dir_path / '.gitkeep'
        if not gitkeep_path.exists():
            with open(gitkeep_path, 'w') as f:
                pass  # Create empty file
            print(f"Created .gitkeep in: {dir_path.relative_to(project_root)}")

    print("Directory structure prepared for Git tracking.")

if __name__ == "__main__":
    main() 

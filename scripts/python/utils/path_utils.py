#!/usr/bin/env python3
"""Path utility functions for the TE analysis project."""

from pathlib import Path
import logging
import os
import sys

import yaml

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def get_project_root():
    """
    Get the absolute path to the project root directory.
    
    Returns:
    --------
    pathlib.Path
        Path to the project root directory
    """
    current = Path(__file__).resolve().parent
    for candidate in [current, *current.parents]:
        if (candidate / "paths.yaml").exists():
            return candidate
    return Path(__file__).resolve().parent.parent.parent.parent


def load_config(config_path=None):
    """
    Load the configuration file.
    
    Parameters:
    -----------
    config_path : str, optional
        Path to the configuration file
        
    Returns:
    --------
    dict
        Configuration dictionary
    """
    if config_path is None:
        project_root = get_project_root()
        config_path = project_root / "paths.yaml"
    else:
        config_path = Path(config_path)
        if not config_path.is_absolute():
            config_path = get_project_root() / config_path
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def resolve_path(path_str, config_path=None):
    """
    Resolve a path string using the configuration file.
    
    Parameters:
    -----------
    path_str : str
        Path string to resolve (e.g., 'input_data.dnaPipeTE')
    config_path : str, optional
        Path to the configuration file
        
    Returns:
    --------
    str
        Resolved absolute path
    """
    # Load the configuration file
    config = load_config(config_path)
    
    # Split the path string by dots
    parts = path_str.split('.')
    
    # Navigate through the configuration dictionary
    current = config
    for part in parts:
        if part in current:
            current = current[part]
        else:
            raise KeyError(f"Path not found in configuration: {path_str}")
    
    # If the result is not a string, raise an error
    if not isinstance(current, str):
        raise ValueError(f"Path resolves to a non-string value: {path_str}")
    
    # Return the resolved path
    return current


def ensure_directory(path):
    """
    Ensure that a directory exists, creating it if necessary.
    
    Parameters:
    -----------
    path : str
        Path to the directory to ensure
        
    Returns:
    --------
    str
        Path to the directory
    """
    os.makedirs(path, exist_ok=True)
    return path


def find_file(filename, search_paths, return_first=True):
    """
    Find a file in a list of search paths.
    
    Parameters:
    -----------
    filename : str
        Filename to find
    search_paths : list
        List of directories to search
    return_first : bool, default=True
        Whether to return the first match or all matches
        
    Returns:
    --------
    str or list
        Path to the file if found, or list of paths if return_first=False
    """
    matches = []
    
    for search_path in search_paths:
        full_path = os.path.join(search_path, filename)
        if os.path.exists(full_path):
            if return_first:
                return full_path
            matches.append(full_path)
    
    if matches:
        return matches
    return None


if __name__ == '__main__':
    """Run path utilities as a standalone script."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Path utilities for the Desmognathus_TE project')
    parser.add_argument('--check-paths', action='store_true',
                        help='Check if all paths in config exist')
    parser.add_argument('--config', '-c', type=str, default='paths.yaml',
                        help='Path to configuration file')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    if args.check_paths:
        config = load_config(args.config)
        print("Checking paths in config...")
        
        # Recursively find all path entries and check if they exist
        def check_paths_recursive(config_obj, prefix=''):
            results = {'exists': [], 'missing': []}
            
            if isinstance(config_obj, dict):
                for key, value in config_obj.items():
                    current_path = f"{prefix}.{key}" if prefix else key
                    
                    if isinstance(value, dict):
                        # Check if this is a leaf node with a 'root' key
                        if 'root' in value and isinstance(value['root'], str):
                            path = value['root']
                            if os.path.exists(path):
                                results['exists'].append((current_path, path))
                            else:
                                results['missing'].append((current_path, path))
                        else:
                            # Recurse
                            sub_results = check_paths_recursive(value, current_path)
                            results['exists'].extend(sub_results['exists'])
                            results['missing'].extend(sub_results['missing'])
                    elif isinstance(value, str) and (
                            key.endswith('_dir') or key.endswith('_path') or
                            key in ['root', 'figures', 'tables']):
                        # This looks like a path
                        path = value
                        if os.path.exists(path):
                            results['exists'].append((current_path, path))
                        else:
                            results['missing'].append((current_path, path))
            
            return results
        
        path_results = check_paths_recursive(config)
        
        # Print results
        print(f"\nFound {len(path_results['exists'])} existing paths and {len(path_results['missing'])} missing paths")
        
        if path_results['exists']:
            print("\nExisting paths:")
            for key, path in path_results['exists']:
                print(f"  {key}: {path}")
        
        if path_results['missing']:
            print("\nMissing paths:")
            for key, path in path_results['missing']:
                print(f"  {key}: {path}")
                
            # Ask if missing paths should be created
            create = input("\nCreate missing directories? (y/n): ").lower().strip() == 'y'
            if create:
                for key, path in path_results['missing']:
                    if os.path.splitext(path)[1] == '':  # No extension, likely a directory
                        if ensure_directory(path):
                            print(f"Created directory: {path}")
                        else:
                            print(f"Failed to create directory: {path}")
    else:
        parser.print_help() 

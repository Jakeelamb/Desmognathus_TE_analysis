#!/usr/bin/env python3
"""
Path configuration utilities for the Desmognathus_TE project.
This module loads the paths.yaml configuration and provides functions to access paths.
"""

import os
import yaml
from typing import Dict, Any, Union, Optional

# Default paths
DEFAULT_CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))), 
                                 "config", "paths.yaml")
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Global variable to store the configuration
_config = None

def get_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load and return the path configuration.
    
    Args:
        config_path: Path to the configuration file. If None, uses the default path.
        
    Returns:
        Dict containing the path configuration.
    """
    global _config
    
    if _config is None:
        if config_path is None:
            config_path = DEFAULT_CONFIG_PATH
        
        # Ensure the path is absolute
        if not os.path.isabs(config_path):
            config_path = os.path.join(PROJECT_ROOT, config_path)
        
        # Load the YAML configuration
        try:
            with open(config_path, 'r') as f:
                _config = yaml.safe_load(f)
        except Exception as e:
            raise RuntimeError(f"Error loading path configuration from {config_path}: {e}")
    
    return _config

def get_path(path_key: str, default: Optional[str] = None) -> str:
    """
    Get a path from the configuration using a dot-separated key.
    
    Args:
        path_key: Dot-separated key to locate the path in the configuration.
        default: Default value to return if the key is not found.
        
    Returns:
        The requested path as a string.
        
    Example:
        >>> get_path("data.raw")
        "data/raw"
    """
    config = get_config()
    
    # Split the key into parts
    parts = path_key.split('.')
    
    # Navigate through the config dictionary
    current = config
    for part in parts:
        if part in current:
            current = current[part]
        else:
            return default
    
    # Ensure the path is a string
    if not isinstance(current, str):
        raise ValueError(f"Path key '{path_key}' does not correspond to a string path.")
    
    # Make the path absolute
    if not os.path.isabs(current):
        current = os.path.join(PROJECT_ROOT, current)
    
    return current

def ensure_dir(path: str) -> str:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        path: Path to the directory.
        
    Returns:
        The path to the directory.
    """
    os.makedirs(path, exist_ok=True)
    return path

def get_and_ensure_dir(path_key: str, default: Optional[str] = None) -> str:
    """
    Get a path from the configuration and ensure the directory exists.
    
    Args:
        path_key: Dot-separated key to locate the path in the configuration.
        default: Default value to return if the key is not found.
        
    Returns:
        The path to the directory.
    """
    path = get_path(path_key, default)
    return ensure_dir(path)

if __name__ == "__main__":
    """Print all paths when run as a script for testing."""
    config = get_config()
    
    def print_paths(data, prefix=""):
        for key, value in data.items():
            if isinstance(value, dict):
                print_paths(value, prefix=f"{prefix}{key}.")
            else:
                print(f"{prefix}{key}: {value}")
    
    print("Desmognathus_TE Path Configuration:")
    print("==================================")
    print_paths(config) 
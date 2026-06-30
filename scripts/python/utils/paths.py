#!/usr/bin/env python3
"""Path configuration utilities for the Desmognathus_TE project."""

import os
from pathlib import Path
from typing import Any, Dict, Optional

import yaml


def find_project_root(start: Optional[Path] = None) -> Path:
    """Find the nearest project root containing ``paths.yaml``."""
    current = (start or Path(__file__).resolve()).parent
    for candidate in [current, *current.parents]:
      if (candidate / "paths.yaml").exists():
        return candidate
    raise RuntimeError("Could not find project root containing paths.yaml")


PROJECT_ROOT = find_project_root()
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "paths.yaml"

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
        config_path = Path(config_path)
        if not config_path.is_absolute():
            config_path = PROJECT_ROOT / config_path
        
        # Load the YAML configuration
        try:
            with config_path.open('r') as f:
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
        >>> get_path("input_data.root")
        "input_data"
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
        current = str(PROJECT_ROOT / current)
    
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

"""
Centralized configuration management for the Desmognathus TE analysis project.

This module provides consistent path handling across all Python scripts by:
1. Auto-detecting the project root directory
2. Loading paths from the canonical root `paths.yaml`
3. Providing Path objects for all configured directories

Usage:
    from config import paths, PROJECT_ROOT

    # Access paths
    input_dir = paths.input_data.dnaPipeTE
    output_dir = paths.results.data
    figure_dir = paths.results.figures / "custom_subdir"
    lookup_file = paths.input_data.lookup_table
"""

from pathlib import Path
from typing import Any, Union
import yaml


def find_project_root() -> Path:
    """
    Find the project root by looking for the canonical path configuration.

    Searches upward from this file's location until finding a directory
    containing `paths.yaml`.

    Returns:
        Path to the project root directory.

    Raises:
        FileNotFoundError: If no canonical path configuration can be found.
    """
    current = Path(__file__).resolve().parent

    while current != current.parent:
        if (current / "paths.yaml").exists():
            return current
        current = current.parent

    raise FileNotFoundError(
        "Could not find project root (no paths.yaml found). "
        "Ensure you're running from within the project directory."
    )


class PathConfig:
    """
    A nested configuration object that provides attribute access to paths.

    Paths are stored as Path objects relative to the project root.
    """

    def __init__(self, config_dict: dict, root: Path):
        """
        Initialize PathConfig from a dictionary.

        Args:
            config_dict: Dictionary of path configurations.
            root: Project root directory.
        """
        self._root = root
        self._config = config_dict
        self._path = None

        section_root = config_dict.get("root")
        if isinstance(section_root, str):
            self._path = root / section_root

        for key, value in config_dict.items():
            if isinstance(value, dict):
                # Nested config section
                setattr(self, key, PathConfig(value, root))
            elif isinstance(value, str):
                # Convert string path to Path object
                setattr(self, key, root / value)
            else:
                setattr(self, key, value)

    def __repr__(self) -> str:
        return f"PathConfig({self._config})"

    def __str__(self) -> str:
        return str(self._path) if self._path is not None else super().__repr__()

    def __fspath__(self) -> str:
        if self._path is None:
            raise TypeError("This PathConfig section does not define a root path.")
        return str(self._path)

    def __truediv__(self, other: Union[str, Path]) -> Path:
        if self._path is None:
            raise TypeError("This PathConfig section does not define a root path.")
        return self._path / other

    @property
    def path(self) -> Path:
        if self._path is None:
            raise AttributeError("This PathConfig section does not define a root path.")
        return self._path

    def get(self, key: str, default: Any = None) -> Any:
        """Get a path by key name."""
        return getattr(self, key, default)

    def ensure_dirs(self) -> None:
        """Create all directories in this config section if they don't exist."""
        for key, value in self._config.items():
            if isinstance(value, str):
                path = self._root / value
                if not path.suffix:  # Only create if it's a directory (no file extension)
                    path.mkdir(parents=True, exist_ok=True)


def load_config(config_path: Path = None) -> PathConfig:
    """
    Load configuration from paths.yaml.

    Args:
        config_path: Optional path to config file. If None, auto-detects.

    Returns:
        PathConfig object with all paths.
    """
    if config_path is None:
        root = find_project_root()
        config_path = root / "paths.yaml"
    else:
        config_path = Path(config_path)
        root = find_project_root()
        if not config_path.is_absolute():
            config_path = root / config_path

    with open(config_path, 'r') as f:
        config_dict = yaml.safe_load(f)

    return PathConfig(config_dict, root)


# Module-level constants for easy importing
PROJECT_ROOT = find_project_root()
paths = load_config()


def load_lookup_table():
    """
    Load the species lookup table.

    Returns:
        pandas DataFrame with columns: Species, SRA_Accension, Genome_Accension
    """
    import pandas as pd

    lookup_path = paths.input_data.lookup_table
    if not lookup_path.exists():
        raise FileNotFoundError(f"Lookup table not found at: {lookup_path}")

    return pd.read_csv(lookup_path, sep='\t')


def get_sra_to_species_map() -> dict:
    """
    Get a mapping from SRA accession to species name.

    Returns:
        Dictionary mapping SRA IDs to species names.
    """
    df = load_lookup_table()
    return dict(zip(df['SRA_Accension'], df['Species']))


def get_gca_to_species_map() -> dict:
    """
    Get a mapping from GCA genome accession to species name.

    Returns:
        Dictionary mapping GCA IDs to species names.
    """
    df = load_lookup_table()
    # Extract just the GCA ID (e.g., "GCA_030264455" from "GCA_030264455.1")
    df['GCA_ID'] = df['Genome_Accension'].str.replace(r'\.\d+$', '', regex=True)
    return dict(zip(df['GCA_ID'], df['Species']))


if __name__ == "__main__":
    # Print configuration when run directly (useful for debugging)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"\nInput data paths:")
    print(f"  Root: {paths.input_data.root}")
    print(f"  dnaPipeTE: {paths.input_data.dnaPipeTE}")
    print(f"  RepeatMasker: {paths.input_data.repeatmasker}")
    print(f"  Lookup table: {paths.input_data.lookup_table}")
    print(f"\nResults paths:")
    print(f"  Data: {paths.results.data}")
    print(f"  Figures root: {paths.results.figures.path}")
    print(f"\nInterim: {paths.interim.root}")

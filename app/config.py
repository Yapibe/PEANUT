"""
Configuration loader for PEANUT.

This module loads and provides access to the application configuration from YAML.
"""

import yaml
from pathlib import Path
from typing import Dict, Any


def load_config() -> Dict[str, Any]:
    """Load configuration from YAML file.
    
    Returns:
        Dictionary containing the configuration values.
    """
    # Load from YAML
    config_path = Path(__file__).resolve().parent.parent / "config" / "app_config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    # Add computed paths (not in YAML)
    base_dir = Path(__file__).resolve().parent
    config["base_dir"] = base_dir
    config["static_dir"] = base_dir / "static"
    config["sample_data_dir"] = base_dir / "static" / "sample_data"
    
    return config


# Load the configuration once at module import
config = load_config() 
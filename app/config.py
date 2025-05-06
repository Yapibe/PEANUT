"""
Configuration module for PEANUT.

This module loads the application configuration from YAML and provides it to other modules.
"""

import yaml
from pathlib import Path
from typing import Dict, Any

def load_config() -> Dict[str, Any]:
    """Load configuration from YAML file with environment variable overrides."""
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

# Load configuration at module import time
config = load_config()
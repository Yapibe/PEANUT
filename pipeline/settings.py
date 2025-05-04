# pipeline/settings.py
"""
Settings and configuration classes for the PEANUT pipeline.

This module defines the dataclasses used to store and manage configuration 
settings for both global pipeline parameters and condition-specific settings.
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Optional, Dict, Set, Any
import pandas as pd
import logging

logger = logging.getLogger(__name__)

# Constants for default values
DEFAULT_ALPHA = 1.0
PRECOMPUTED_ALPHAS = {0.1, 0.2}
DEFAULT_MIN_GENES = 15
DEFAULT_MAX_GENES = 500
DEFAULT_FDR = 0.05


class Species(str, Enum):
    """Valid species options."""
    H_SAPIENS = "H_sapiens"
    MOUSE = "mouse"


class NetworkType(str, Enum):
    """Valid network types."""
    ANAT = "Anat"
    HUMANNET = "HumanNet"
    CUSTOM = "custom"


class PathwayFileType(str, Enum):
    """Valid pathway file types."""
    KEGG = "kegg"
    MSIG_C2_CANONICAL = "msig_c2_canonical"
    CUSTOM = "custom"


@dataclass
class Settings:
    """General arguments used throughout the pipeline.
    
    This class stores global settings that apply to the entire pipeline execution,
    including network configuration, pathway parameters, and analysis thresholds.
    """

    experiment_name: str
    species: Species = Species.H_SAPIENS
    alpha: float = DEFAULT_ALPHA
    network: NetworkType = NetworkType.ANAT
    pathway_file: PathwayFileType = PathwayFileType.KEGG
    min_genes_per_pathway: int = DEFAULT_MIN_GENES
    max_genes_per_pathway: int = DEFAULT_MAX_GENES
    fdr_threshold: float = DEFAULT_FDR
    run_gsea: bool = False
    restrict_to_network: bool = False
    create_similarity_matrix: bool = False
    similarity_matrix_path: Optional[Path] = field(init=False, default=None)
    
    # Paths and directories
    root_folder: Path = field(init=False, default=None)
    data_dir: Path = field(init=False, default=None)
    date: str = field(init=False, default=None)
    figure_title: str = "Pathway Enrichment"
    network_file_path: Path = field(init=False, default=None)
    pathway_file_dir: Path = field(init=False, default=None)
    plot_output_path: Path = field(init=False, default=None)

    def __post_init__(self):
        """Initialize derived paths and settings after basic initialization."""
        # Initialize basic paths
        self.root_folder = Path(__file__).resolve().parent
        self.data_dir = self.root_folder / "Data" / self.species
        self.date = datetime.today().strftime("%Y-%m-%dT%H-%M-%S")

        # Handle matrix settings - alpha=1.0 is special case for GSEA
        if self.alpha == DEFAULT_ALPHA:
            self.run_gsea = True
            self.similarity_matrix_path = None
        else:
            # Determine matrix path based on network and alpha
            self._setup_matrix_path()
        
        # Set network file path
        self._setup_network_path()
        
        # Set up pathway directory and plot output paths
        self.pathway_file_dir = self.data_dir / "pathways" / self.pathway_file
        self.plot_output_path = self.root_folder / "Outputs" / "Plots"
        
        # Note: Directory creation should be handled by the pipeline runner
        # or a dedicated setup function, not within this settings class.

    def _setup_matrix_path(self) -> None:
        """Set up the similarity matrix path based on network and alpha."""
        if self._is_custom_network():
            matrix_path = self.data_dir / "matrix" / "custom" / f"{Path(self.network).stem}_{self.alpha}.npz"
            self.create_similarity_matrix = True
        else:
            # For default networks
            if self.alpha in PRECOMPUTED_ALPHAS and not self.create_similarity_matrix:
                matrix_path = self.data_dir / "matrix" / f"{self.network}_{self.alpha}.npz"
            else:
                matrix_path = self.data_dir / "matrix" / "matrix" / f"{self.network}_{self.alpha}.npz"
                self.create_similarity_matrix = True
        
        self.similarity_matrix_path = matrix_path
    
    def _is_custom_network(self) -> bool:
        """Check if network is custom."""
        return self.network == NetworkType.CUSTOM or str(self.network).startswith('custom/')
    
    def _setup_network_path(self) -> None:
        """Set up the network file path."""
        if self._is_custom_network():
            self.network_file_path = Path(self.network)
        else:
            self.network_file_path = self.data_dir / "network" / self.network


@dataclass
class ConditionSettings:
    """Holds condition-specific settings and results."""

    condition_name: str
    experiment_name: str
    scores_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    output_path: Path = field(init=False, default=None)
    date: str = field(init=False, default=None)
    filtered_genes: Set[int] = field(default_factory=set)
    pathways_statistics: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize derived fields after basic initialization."""
        self.date = datetime.today().strftime("%Y-%m-%dT%H-%M-%S")
        self.output_path = (
            Path(__file__).resolve().parent
            / "Outputs"
            / f"{self.experiment_name}_{self.condition_name}_{self.date}.xlsx"
        )
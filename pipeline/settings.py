# pipeline/settings.py
"""
Settings and configuration classes for the PEANUT pipeline.

This module defines the dataclasses used to store and manage configuration 
settings for both global pipeline parameters and condition-specific settings.
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Set, Any
import pandas as pd
import logging

logger = logging.getLogger(__name__)

@dataclass
class Settings:
    """General arguments used throughout the pipeline.
    
    This class stores global settings that apply to the entire pipeline execution,
    including network configuration, pathway parameters, and analysis thresholds.
    
    Attributes:
        experiment_name: Unique name for the experiment.
        species: Species identifier (default: "H_sapiens").
        alpha: Propagation parameter between 0-1 (default: 1).
        network: Network identifier or path (default: "Anat").
        pathway_file: Pathway database to use (default: "kegg").
        minimum_gene_per_pathway: Minimum genes needed for pathway analysis (default: 15).
        maximum_gene_per_pathway: Maximum genes allowed for pathway analysis (default: 500).
        FDR_THRESHOLD: False discovery rate threshold for significance (default: 0.05).
        run_gsea: Whether to run GSEA instead of propagation (default: False).
        restrict_to_network: Whether to restrict gene analysis to network genes only (default: False).
        create_similarity_matrix: Whether to create a new similarity matrix (default: False).
        figure_title: Title for generated figures (default: "Pathway Enrichment").
    """

    experiment_name: str
    species: str = "H_sapiens"
    alpha: float = 1.0
    network: str = "Anat"
    pathway_file: str = "kegg"
    minimum_gene_per_pathway: int = 15
    maximum_gene_per_pathway: int = 500
    FDR_THRESHOLD: float = 0.05
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
        """Initialize derived paths and settings after basic initialization.
        
        This method sets up all necessary file paths and directories based on the 
        provided configuration.
        """
        self.root_folder = Path(__file__).resolve().parent
        self.data_dir = self.root_folder / "Data" / self.species
        self.date = datetime.today().strftime("%d_%m_%Y__%H_%M_%S")

        # Set up matrix path based on network type and create_similarity_matrix flag
        if self.network.startswith('custom/'):
            # For custom networks, use the custom matrix directory
            matrix_dir = self.data_dir / "matrix" / "custom"
            matrix_dir.mkdir(parents=True, exist_ok=True)
            self.similarity_matrix_path = matrix_dir / f"{Path(self.network).stem}_{self.alpha}.npz"
        else:
            # For default networks
            if self.create_similarity_matrix:
                matrix_dir = self.data_dir / "matrix" / "matrix"
                matrix_dir.mkdir(parents=True, exist_ok=True)
                self.similarity_matrix_path = matrix_dir / f"{self.network}_{self.alpha}.npz"
            else:
                self.similarity_matrix_path = self.data_dir / "matrix" / f"{self.network}_{self.alpha}.npz"

        # Handle network file path
        if self.network.startswith('custom/'):
            # For custom networks, use the path provided in the settings
            self.network_file_path = Path(self.network)
        else:
            # For default networks, just use the network name
            self.network_file_path = self.data_dir / "network" / self.network

        # Ensure the network file exists
        if not self.network_file_path.exists():
            raise FileNotFoundError(f"Network file not found: {self.network_file_path}")

        # Set up pathway directory and plot output paths
        self.pathway_file_dir = self.data_dir / "pathways" / self.pathway_file
        self.plot_output_path = self.root_folder / "Outputs" / "Plots"


@dataclass
class ConditionSettings:
    """Holds condition-specific settings and results.
    
    This class stores settings and results specific to a single experimental condition.
    
    Attributes:
        condition_name: Name identifier for the condition.
        experiment_name: Parent experiment name.
        scores_df: DataFrame containing gene scores (default: empty DataFrame).
        output_path: Path for result output (auto-generated).
        date: Timestamp for the analysis (auto-generated).
        filtered_genes: Set of genes after filtering (default: empty set).
        pathways_statistics: Dictionary storing pathway statistics (default: empty dict).
    """

    condition_name: str
    experiment_name: str
    scores_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    output_path: Path = field(init=False, default=None)
    date: str = field(init=False, default=None)
    filtered_genes: Set[int] = field(default_factory=set)
    pathways_statistics: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize derived fields after basic initialization."""
        self.date = datetime.today().strftime("%d_%m_%Y__%H_%M_%S")
        self.output_path = (
            Path(__file__).resolve().parent
            / "Outputs"
            / f"{self.experiment_name}_{self.condition_name}_{self.date}.xlsx"
        )
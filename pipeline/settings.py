from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional


@dataclass
class Settings:
    """General arguments used throughout the pipeline."""

    experiment_name: str
    species: str = "H_sapiens"
    species: str = "H_sapiens"
    alpha: float = 1
    network: str = "Anat"
    pathway_file: str = "kegg"
    network: str = "Anat"
    pathway_file: str = "kegg"
    minimum_gene_per_pathway: int = 15
    maximum_gene_per_pathway: int = 500
    FDR_THRESHOLD: float = 0.05
    JAC_THRESHOLD: float = 0.2

    # Paths and directories
    root_folder: Path = field(init=False)
    data_dir: Path = field(init=False)
    similarity_matrix_path: Optional[Path] = field(init=False)
    date: str = field(init=False)
    figure_title: str = "Pathway Enrichment"
    network_file_path: Path = field(init=False)
    genes_names_file: str = "gene_info.json"
    genes_names_file_path: Path = field(init=False)
    pathway_file_dir: Path = field(init=False)
    run_gsea: bool = field(default=False)

    def __post_init__(self):
        """Post-initialization to set up directories, file paths, and other settings."""
        self.root_folder = Path(__file__).resolve().parent
        self.data_dir = self.root_folder / "Data" / self.species
        self.date = datetime.today().strftime("%d_%m_%Y__%H_%M_%S")

        # Adjust the similarity matrix paths based on normalization type
        self.similarity_matrix_path = (
            self.data_dir
            / "matrix"
            / f"{self.network}_{self.alpha}.npz"
        )
        self.network_file_path = (
            self.data_dir / "network" / self.network
        )
        self.genes_names_file_path = (
            self.data_dir / "gene_names" / self.genes_names_file
        )
        self.pathway_file_dir = (
            self.data_dir / "pathways" / self.pathway_file
        )


@dataclass
class ConditionSettings:
    """Holds condition-specific settings."""

    condition_name: str
    experiment_name: str
    output_path: Path = field(init=False)
    date: str = field(init=False)
    filtered_genes: set = field(default_factory=set)
    filtered_pathways: dict = field(default_factory=dict)
    ks_significant_pathways_with_genes: dict = field(
        default_factory=dict
    )

    def __post_init__(self):
        """Initialize the output path based on the experiment and condition names."""
        self.date = datetime.today().strftime("%d_%m_%Y__%H_%M_%S")
        self.output_path = (
            Path(__file__).resolve().parent
            / "Outputs"
            / f"{self.experiment_name}_{self.condition_name}_{self.date}.xlsx"
        )

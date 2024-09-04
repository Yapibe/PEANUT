from datetime import datetime
from dataclasses import dataclass, field
from typing import Optional
from os import path


@dataclass
class Settings:
    """
    General arguments used throughout the pipeline.
    """
    experiment_name: str
    species: str = 'H_sapiens'
    alpha: float = 1
    network: str = 'Anat'
    pathway_file: str = 'kegg'
    minimum_gene_per_pathway: int = 15
    maximum_gene_per_pathway: int = 500
    FDR_THRESHOLD: float = 0.05
    JAC_THRESHOLD: float = 0.2


    # Paths and directories
    root_folder: str = field(init=False)
    data_dir: str = field(init=False)
    similarity_matrix_path: Optional[str] = field(init=False)
    date: str = field(init=False)
    figure_title: str = 'Pathway Enrichment'
    output_path: str = field(init=False)
    network_file_path: str = field(init=False)
    genes_names_file: str = 'gene_info.json'
    genes_names_file_path: str = field(init=False)
    pathway_file_dir: str = field(init=False)
    run_gsea: bool = field(default=False)
    filtered_genes: set = field(default_factory=set)
    filtered_pathways: dict = field(default_factory=dict)
    ks_significant_pathways_with_genes: dict = field(default_factory=dict)
    def __post_init__(self):
        """
        Post-initialization to set up directories, file paths, and other settings.
        """
        self.root_folder = path.dirname(path.abspath(__file__))
        self.data_dir = path.join(self.root_folder, 'Data', self.species)
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.output_path = path.join(self.root_folder, 'Outputs', '{}_{}.xlsx'.format(self.experiment_name, self.date))

        # Adjust the similarity matrix paths based on normalization type
        self.similarity_matrix_path = path.join(self.data_dir, 'matrix',
                                                f'{self.network}_{self.alpha}.npz')
        self.network_file_path = path.join(self.data_dir, 'network', self.network)
        self.genes_names_file_path = path.join(self.data_dir, 'gene_names', self.genes_names_file)
        self.pathway_file_dir = path.join(self.data_dir, 'pathways', self.pathway_file)



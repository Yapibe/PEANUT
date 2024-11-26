from os import path, makedirs
from datetime import datetime
from dataclasses import dataclass, field
from typing import Callable, Optional


@dataclass
class GeneralArgs:
    """
    General arguments used throughout the pipeline.
    """
    alpha: float = 1
    method: Optional[str] = None
    run_NGSEA: bool = False
    input_type: str = 'Score'
    network: str = 'H_sapiens'
    pathway_file: str = 'kegg'
    run_propagation: bool = True
    run_simulated: bool = False
    run_gsea: bool = True
    debug: bool = False
    create_similarity_matrix: bool = False
    normalization_type: str = 'symmetric'

    # Paths and directories
    root_folder: str = field(init=False)
    data_dir: str = field(init=False)
    similarity_matrix_path: Optional[str] = field(init=False)
    FDR_threshold: float = 0.05
    minimum_gene_per_pathway: int = 15
    maximum_gene_per_pathway: int = 500
    JAC_THRESHOLD: float = 0.2
    Experiment_name: str = field(init=False)
    date: str = field(init=False)
    figure_title: str = 'Pathway Enrichment'
    output_dir: str = field(init=False)
    input_dir: str = field(init=False)
    propagation_folder: str = field(init=False)
    gsea_out: Optional[str] = field(init=False)
    network_file_path: str = field(init=False)
    genes_names_file: str = 'gene_info.json'
    genes_names_file_path: str = field(init=False)
    pathway_file_dir: str = field(init=False)

    def __post_init__(self):
        """
        Post-initialization to set up directories, file paths, and other settings.
        """
        self.root_folder = path.dirname(path.abspath(__file__))
        self.data_dir = path.join(self.root_folder, 'Data', 'H_sapiens')
        self.Experiment_name = 'Simulated' if self.run_simulated else 'GSE'
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.output_dir = path.join(self.root_folder, 'Outputs')

        self.similarity_matrix_path = path.join(self.data_dir, 'matrix',
                                                f'Anat_{self.alpha}.npz')

        self.input_dir = path.join(self.root_folder, 'Inputs', 'Simulated') if self.run_simulated else path.join(
            self.root_folder,
            'Inputs',
            'experiments_data',
            self.Experiment_name,
            'XLSX')
        self.propagation_folder = self._create_output_subdir(
            path.join('Propagation_Scores', self.method or '', self.network, self.pathway_file, f'alpha_{self.alpha}'))
        self.gsea_out = self._create_output_subdir(
            path.join('GSEA', self.method or '', self.network, self.pathway_file)) if self.run_gsea else None
        self.network_file_path = path.join(self.data_dir, 'network', self.network)
        self.genes_names_file_path = path.join(self.data_dir, 'gene_names', self.genes_names_file)
        self.pathway_file_dir = path.join(self.data_dir, 'pathways', self.pathway_file)

    def _create_output_subdir(self, subdir: str) -> str:
        """
        Create a subdirectory under the output directory if it does not exist.

        Parameters:
        - subdir (str): Subdirectory path to create.

        Returns:
        - str: Full path to the created subdirectory.
        """
        full_path = path.join(self.output_dir, subdir)
        makedirs(full_path, exist_ok=True)  # Allow existing directories without error
        return full_path



@dataclass
class PropagationTask:
    """
    A task for performing propagation on a dataset.

    Attributes:
    - general_args (GeneralArgs): General arguments and settings.
    - test_name (str): Name of the comparison test.
    - results (dict): Dictionary to store results.
    - output_folder (str): Directory for output files.
    """
    general_args: GeneralArgs
    test_name: str
    results: dict = field(default_factory=dict)
    output_folder: str = field(init=False)

    def __post_init__(self):
        """
        Post-initialization to set up the output folder based on the test name.
        """
        self.test_name = self.test_name.split(".")[0]
        self.output_folder = path.join(self.general_args.propagation_folder, self.test_name)


@dataclass
class EnrichTask:
    """
    A task for performing enrichment analysis.

    Attributes:
    - name (str): Name of the task.
    - statistic_test (Callable): Statistical test function for enrichment analysis.
    - target_field (str): Field in the data to target for enrichment.
    - create_scores (bool): Flag to determine whether to create scores (default: True).
    - propagation_file (Optional[str]): Filename of the propagated gene scores.
    - results (dict): Dictionary to store results.
    - filtered_genes (set): Set of filtered genes.
    - filtered_pathways (dict): Dictionary of filtered pathways.
    - ks_significant_pathways_with_genes (dict): Dictionary of significant pathways with genes.
    """
    name: str
    statistic_test: Callable
    target_field: str
    create_scores: bool = True
    propagation_file: Optional[str] = None
    results: dict = field(default_factory=dict)
    filtered_genes: set = field(default_factory=set)
    filtered_pathways: dict = field(default_factory=dict)
    ks_significant_pathways_with_genes: dict = field(default_factory=dict)
    ks_non_significant_pathways_with_genes: dict = field(default_factory=dict)
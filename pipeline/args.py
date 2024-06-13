from os import path, makedirs
from datetime import datetime


class GeneralArgs:
    def __init__(self, alpha: float = 1, FDR_threshold: float = 0.05, figure_title: str ='Pathway Enrichment',
                 create_similarity_matrix: bool =False, run_propagation: bool =True, run_gsea: bool =True):
        """
        Initializes general arguments used throughout the pipeline.

        Parameters:
        - alpha (float): Alpha value for similarity matrix (default: 1).
        - FDR_threshold (float): False discovery rate threshold (default: 0.05).
        - figure_title (str): Name of the experiment (default: 'Pathway Enrichment').
        - create_similarity_matrix (bool): Flag to create similarity matrix (default: False).
        - run_propagation (bool): Flag to run propagation (default: True).

        Attributes:
        - alpha (float): Alpha value.
        - FDR_threshold (float): FDR threshold.
        - minimum_gene_per_pathway (int): Minimum number of genes per pathway.
        - maximum_gene_per_pathway (int): Maximum number of genes per pathway.
        - JAC_THRESHOLD (float): Jaccard threshold.
        - run_propagation (bool): Flag to run propagation.
        - Experiment_name (str): Name of the experiment.
        - date (str): Current date and time.
        - figure_title (str): Name of the experiment to be displayed on figures.
        - run_gsea (bool): Flag to run GSEA.
        - root_folder (str): Root directory of the script.
        - data_dir (str): Directory for input data.
        - output_dir (str): Directory for output data.
        - input_dir (str): Directory for input experiment data.
        - network_file (str): Network file name.
        - network_file_path (str): Path to the network file.
        - genes_names_file (str): Gene names file name.
        - genes_names_file_path (str): Path to the gene names file.
        - bio_pathways (str): Pathway file name.
        - pathway_file_dir (str): Directory for pathway files.
        - similarity_matrix_path (str): Path to the similarity matrix file.
        - create_similarity_matrix (bool): Flag to create similarity matrix.
        - temp_output_folder (str): Directory for temporary outputs.
        - propagation_folder (str): Directory for propagation scores.
        """
        # General Parameters
        self.alpha = alpha  # for no propagation use alpha=1
        self.FDR_threshold = FDR_threshold
        self.minimum_gene_per_pathway = 20
        self.maximum_gene_per_pathway = 200
        self.JAC_THRESHOLD = 0.2
        self.run_propagation = run_propagation
        # self.Experiment_name = 'Parkinson'
        self.Experiment_name = 'Simulated'
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.figure_title = figure_title
        self.run_simulated = True
        self.run_gsea = run_gsea
        self.run_hyper = False

        # root directory
        self.root_folder = path.dirname(path.abspath(__file__))

        # directory paths
        self.data_dir = path.join(self.root_folder, 'Data', 'H_sapiens')
        self.output_dir = path.join(self.root_folder, 'Outputs')
        # self.input_dir = path.join(self.root_folder, 'Inputs', 'experiments_data', self.Experiment_name)
        self.input_dir = path.join(self.root_folder, 'Inputs', 'Simulated')
        self.pathway_file = 'bio_pathways_gmt.gmt'
        # Data directory directories
        self.network_file = 'H_sapiens.net'
        self.network_file_path = path.join(self.data_dir, 'network', self.network_file)
        self.genes_names_file = 'H_sapiens.gene_info'
        self.genes_names_file_path = path.join(self.data_dir, 'genes_names', self.genes_names_file)

        self.pathway_file_dir = path.join(self.data_dir, 'pathways', self.pathway_file)
        self.similarity_matrix_path = path.join(self.data_dir, 'matrix')
        self.create_similarity_matrix = create_similarity_matrix
        self.similarity_matrix_path = path.join(self.similarity_matrix_path,
                                                f'similarity_matrix_{self.alpha}.npz')
        # Output directory directories
        self.temp_output_folder = path.join(self.output_dir, 'Temp')
        makedirs(self.temp_output_folder, exist_ok=True)
        self.propagation_folder = path.join(self.output_dir, 'Propagation_Scores')
        self.gsea_out = path.join(self.output_dir, 'GSEA')


class PropagationTask:
    def __init__(self, general_args: GeneralArgs, test_name: str):
        """
        Initialize a PropagationTask instance.

        Parameters:
        - general_args (GeneralArgs): General arguments and settings.
        - test_name (str): Name of the test for which propagation is performed.

        Attributes:
        - general_args (GeneralArgs): General arguments and settings.
        - results (dict): Dictionary to store results.
        - test_name (str): Name of the test.
        - test_file (str): Name of the test file.
        - test_file_path (str): Path to the test file.
        - remove_self_propagation (bool): Flag to remove self propagation.
        - output_folder (str): Directory for output files.
        """
        self.general_args = general_args
        self.results = {}
        self.test_name = test_name
        self.test_file = f'{test_name}.xlsx'
        self.test_file_path = path.join(self.general_args.input_dir, self.test_file)
        #TODO: understand this flag
        self.remove_self_propagation = False
        self.output_folder = path.join(self.general_args.propagation_folder, self.test_name)


class EnrichTask:
    def __init__(self, name: str, statistic_test: callable, target_field: str, create_scores: bool = True, propagation_file: str = None):
        """
        Initialize an EnrichTask instance.

        This class configures an enrichment analysis task, including setting file paths and statistical tests.

        Parameters:
        - name (str): Name of the task.
        - statistic_test (function): Statistical test function to use for enrichment analysis.
        - target_field (str): Field in the data to target for enrichment analysis.
        - create_scores (bool): Flag to determine whether to create scores (default: True).
        - propagation_file (str): Filename of the propagated gene scores (default: None).

        Attributes:
        - name (str): Name of the task.
        - statistic_test (function): Statistical test function.
        - target_field (str): Field in the data to target.
        - results (dict): Dictionary to store results.
        - create_scores (bool): Flag to determine whether to create scores.
        - filtered_genes (set): Set of filtered genes.
        - filtered_pathways (dict): Dictionary of filtered pathways.
        - ks_significant_pathways_with_genes (dict): Dictionary of significant pathways with genes.
        - propagation_file (str): Filename of the propagated gene scores.
        - temp_output_folder (str): Directory for temporary output files.
        """
        self.name = name
        self.statistic_test = statistic_test
        self.target_field = target_field
        self.results = dict()
        self.create_scores = create_scores
        self.filtered_genes = set()
        self.filtered_pathways = dict()
        self.ks_significant_pathways_with_genes = dict()
        self.propagation_file = propagation_file
        self.temp_output_folder = path.join(path.dirname(path.realpath(__file__)), 'Outputs', 'Temp')


class PathwayResults:
    def __init__(self, p_value, direction, adj_p_value):
        self.p_value = p_value
        self.direction = direction
        self.adj_p_value = adj_p_value

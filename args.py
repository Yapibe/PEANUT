from os import path, makedirs
from datetime import datetime


class PropagationTask:
    def __init__(self, general_args, test_name):
        self.general_args = general_args
        self.results = dict()
        self.test_name = test_name
        self.test_file = f'{test_name}.xlsx'
        self.test_file_path = path.join(self.general_args.input_dir, self.test_file)
        #TODO: understand this flag
        self.remove_self_propagation = False
        self.output_folder = path.join(self.general_args.propagation_folder, self.test_name)


class EnrichTask:
    def __init__(self, name, statistic_test, target_field, create_scores=True, propagation_file=None):
        """
        Initializes an enrichment task with specified parameters.

        This class configures an enrichment analysis task, including setting file paths and statistical tests
        Parameters:
        - name (str): Name of the task.
        - propagation_file (str): Filename of the propagated gene scores.
        - propagation_folder (str): Folder path where the propagation file is located.
        - statistic_test (function): Statistical test function to use for enrichment analysis.
        - target_field (str): Field in the data to target for enrichment analysis.
        - constrain_to_experiment_genes (bool): Flag to constrain analysis to experiment genes only
        Attributes:
        - Paths and parameters for running enrichment analysis.
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


class GeneralArgs:
    def __init__(self,alpha=1, FDR_threshold=0.05, figure_title='Pathway Enrichment', create_similarity_matrix=False,
                 run_propagation=True):
        """
        Contains general arguments and settings for pathway enrichment analysis
        This class encapsulates various parameters and settings used across different stages of pathway enrichment analysis
        Parameters:
        - FDR_threshold (float): False Discovery Rate threshold for statistical significance.
        - figure_title (str): Title for the figure or output
        Attributes:
        - Configurations like minimum and maximum genes per pathway, FDR threshold, and paths for output and figures.
        """
        # General Parameters
        self.minimum_gene_per_pathway = 20
        self.maximum_gene_per_pathway = 200
        self.FDR_threshold = FDR_threshold
        self.JAC_THRESHOLD = 0.2
        # for no propagation use alpha=1
        self.alpha = alpha
        self.run_propagation = run_propagation

        self.Experiment_name = 'Parkinson'
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.figure_title = figure_title

        # root directory
        self.root_folder = path.dirname(path.abspath(__file__))

        # directory paths
        self.data_dir = path.join(self.root_folder, 'Data', 'H_sapiens')
        self.output_dir = path.join(self.root_folder, 'Outputs')
        self.input_dir = path.join(self.root_folder, 'Inputs', 'experiments_data', self.Experiment_name)

        # Data directory directories
        self.network_file = 'H_sapiens.net'
        self.network_file_path = path.join(self.data_dir, 'network', self.network_file)
        self.genes_names_file = 'H_sapiens.gene_info'
        self.genes_names_file_path = path.join(self.data_dir, 'genes_names', self.genes_names_file)
        self.pathway_file = 'pathways'
        self.pathway_file_dir = path.join(self.data_dir, 'pathways', self.pathway_file)
        self.similarity_matrix_path = path.join(self.data_dir, 'matrix')
        self.create_similarity_matrix = create_similarity_matrix
        self.similarity_matrix_path = path.join(self.similarity_matrix_path,
                                                f'similarity_matrix_{self.alpha}.npz')
        # Output directory directories
        self.temp_output_folder = path.join(self.output_dir, 'Temp')
        makedirs(self.temp_output_folder, exist_ok=True)
        self.propagation_folder = path.join(self.output_dir, 'propagation_scores')


class PathwayResults:
    def __init__(self, p_value, direction, adj_p_value=None):
        """
        Stores the results of pathway analysis
        This class is used to hold the results of a pathway analysis, including the p-value, direction of change,
        and adjusted p-value if available
        Parameters:
        - p_value (float): The p-value resulting from the statistical test.
        - direction (bool): Indicates the direction of the effect (True for positive, False for negative).
        - adj_p_value (float, optional): The adjusted p-value after correction for multiple testing
        Attributes:
        - Holds the statistical results for a specific pathway.
        """
        self.p_value = p_value
        self.direction = direction
        self.adj_p_value = adj_p_value
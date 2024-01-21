from os import path
from datetime import datetime
from utils import get_root_path


class PropagationTask:
    def __init__(self, experiment_name, alpha=1, network_file='H_sapiens.net',
                 create_similarity_matrix=False):
        """
        Initializes a task for gene score propagation.

        This class sets up the necessary parameters and paths for performing gene score propagation in a network.

        Parameters:
        - experiment_name (str): Name of the experiment.
        - alpha (float): Propagation parameter, controlling the influence of network structure on propagation.
        - network_file (str): Name of the file containing the network data.
        - create_similarity_matrix (bool): Flag to indicate whether to create a new similarity matrix.

        Attributes:
        - Various paths and parameters are set up based on the input parameters to facilitate gene score propagation.
        """

        # General Parameters
        self.experiment_name = experiment_name
        self.experiment_file = 'Inputs/experiments_data/roded_M.xlsx'
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'Data'
        self.network_file = network_file
        self.genes_names_file = 'H_sapiens.gene_info'
        self.pathway_file = 'pathway_file'
        self.propagation_input_type = 'Score'
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.create_similarity_matrix = create_similarity_matrix
        self.remove_self_propagation = False
        # Propagation Parameters
        # for no propagation use alpha=1
        self.alpha = alpha
        self.results = dict()
        # Derived Parameters (Initial placeholders)
        self.data_dir = None
        self.network_file_path = None
        self.experiment_file_path = None
        self.pathway_file_dir = None
        self.genes_names_file_path = None
        self.input_dir = None
        self.similarity_matrix_path = None
        self.output_folder = None

        # Initialize derived parameters
        self.get_derived_parameters()

    def get_derived_parameters(self):
        """
        Set derived parameters based on the initial parameters.
        """
        self.data_dir = path.join(self.root_folder, self.data_file)
        self.genes_names_file_path = path.join(self.data_dir, 'H_sapiens', 'genes_names', self.genes_names_file)
        self.network_file_path = path.join(self.data_dir, 'H_sapiens', 'network', self.network_file)
        self.pathway_file_dir = path.join(self.data_dir, 'H_sapiens', 'pathways', self.pathway_file)
        self.similarity_matrix_path = path.join(self.data_dir, 'H_sapiens', 'matrix',
                                                f'similarity_matrix_{self.alpha}.npz')
        self.input_dir = path.join(self.root_folder, 'Inputs', 'experiments_data')
        self.experiment_file_path = path.join(self.input_dir, self.experiment_file)
        self.output_folder = path.join(self.root_folder, 'Outputs', 'propagation_scores', self.experiment_name)


class EnrichTask:
    def __init__(self, name, statistic_test, target_field, alpha=0.1,
                 create_propagation_matrix=False, create_scores=True, propagation_file=None, propagation_folder=None):
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
        self.propagation_file = propagation_file
        self.propagation_scores_path = path.join(get_root_path(), propagation_folder)
        self.statistic_test = statistic_test
        self.target_field = target_field
        self.results = dict()
        self.alpha = alpha
        self.create_scores = create_scores
        self.create_similarity_matrix = create_propagation_matrix
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'Data'
        self.data_dir = path.join(self.root_folder, self.data_file)
        self.similarity_matrix_path = path.join(self.data_dir, 'H_sapiens', 'matrix',
                                                f'similarity_matrix_{self.alpha}.npz')

class NewEnrichTask:
    def __init__(self, experiment_name, statistic_test):
        """
        Initializes an enrichment task with specified parameters.

        This class configures an enrichment analysis task, including setting file paths and statistical tests
        Parameters:
        - experiment_name (str): Name of the experiment.
        - statistic_test (function): Statistical test function to use for enrichment analysis.
        - network_file (str): Name of the file containing the network data. (optional if needed)
        Attributes:
        - Paths and parameters for running enrichment analysis.
        """
        self.experiment_name = experiment_name
        self.species = 'H_sapiens'
        self.statistic_test = statistic_test
        self.results = dict()
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'Data'
        # self.network_file = network_file  # optional if needed
        self.genes_names_file = 'H_sapiens.gene_info'  # optional if needed
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.pathway_file = 'pathways'

        # Derived Parameters (Initial placeholders)
        self.data_dir = path.join(self.root_folder, self.data_file)
        self.genes_names_file_path = path.join(self.data_dir, self.species, 'genes_names', self.genes_names_file)  # optional if needed
        self.pathway_file_dir = path.join(self.data_dir, self.species, 'pathways', self.pathway_file)
        # self.network_file_path = path.join(self.data_dir, 'H_sapiens', 'network', self.network_file)  # optional if needed
        self.input_dir = None
        self.output_folder = None

        # Initialize derived parameters
        self.get_derived_parameters()

    def get_derived_parameters(self):
        """
        Set derived parameters based on the initial parameters.
        """
        self.input_dir = path.join(self.root_folder, 'Inputs', 'experiments_data')
        self.experiment_file_path = path.join(self.input_dir, f'{self.experiment_name}.xlsx')
        self.output_folder = path.join(self.root_folder, 'Outputs', 'enrichment_scores', self.experiment_name)


class GeneralArgs:
    def __init__(self, network_path, genes_names_path, pathway_members_path, FDR_threshold=0.05,
                 output_folder_name=None, figure_name=None, figure_title='Pathway Enrichment '):
        """
        Contains general arguments and settings for pathway enrichment analysis
        This class encapsulates various parameters and settings used across different stages of pathway enrichment analysis
        Parameters:
        - network_path (str): Path to the network file.
        - genes_names_path (str): Path to the genes names file.
        - pathway_members_path (str): Path to the file containing pathway members.
        - FDR_threshold (float): False Discovery Rate threshold for statistical significance.
        - output_folder_name (str, optional): Name of the output folder.
        - figure_name (str, optional): Name of the figure to be generated.
        - figure_title (str): Title for the figure or output
        Attributes:
        - Configurations like minimum and maximum genes per pathway, FDR threshold, and paths for output and figures.
        """
        self.minimum_gene_per_pathway = 20
        self.maximum_gene_per_pathway = 200
        self.network_file_path = network_path
        self.genes_names_file_path = genes_names_path
        self.pathway_databases = ['_']
        self.pathway_keywords = ['_']
        self.FDR_threshold = FDR_threshold
        if output_folder_name is None:
            output_folder_name = 'Enrichment_maps'
        self.output_path = path.join(get_root_path(), 'Outputs', output_folder_name)
        self.figure_name = figure_name if figure_name is not None else 'figure'
        self.pathway_members_path = pathway_members_path
        self.figure_title = figure_title
        self.use_gsea = False


class NewGeneralArgs:
    def __init__(self, FDR_threshold=0.05, output_folder_name=None, figure_name=None, figure_title='Pathway Enrichment'):
        self.minimum_gene_per_pathway = 20
        self.maximum_gene_per_pathway = 60
        self.FDR_threshold = FDR_threshold
        if output_folder_name is None:
            output_folder_name = 'Enrichment_maps'
        self.output_path = path.join(get_root_path(), 'Outputs', output_folder_name)
        self.figure_name = figure_name if figure_name is not None else 'figure'
        self.figure_title = figure_title
        self.use_gsea = False


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
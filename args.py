from os import path
from datetime import datetime
from utils import get_root_path


class PropagationTask:
    def __init__(self, experiment_name, alpha=0.1, network_file='H_sapiens.net',
                 create_similarity_matrix=True):
        """
        Initialize the Propagation Task with default and optional parameters.
        """

        # General Parameters
        self.experiment_name = experiment_name
        self.experiment_file = 'scores.xlsx'
        self.root_folder = path.dirname(path.realpath(__file__))
        self.data_file = 'Data'
        self.network_file = network_file
        self.genes_names_file = 'H_sapiens.gene_info'
        self.pathway_file = 'pathway_file'
        self.propagation_input_type = 'abs_Score'
        self.date = datetime.today().strftime('%d_%m_%Y__%H_%M_%S')
        self.create_similarity_matrix = create_similarity_matrix
        self.remove_self_propagation = False
        # Propagation Parameters
        self.alpha = alpha

        # Derived Parameters (Initial placeholders)
        self.data_dir = None
        self.network_file_path = None
        self.experiment_file_path = None
        self.pathway_file_dir = None
        self.genes_names_file_path = None
        self.input_dir = None
        self.similarity_matrix_path = None
        self.output_folder = None
        # self.propagation_scores_path = None
        # self.random_networks_dir = None
        # self.interesting_pathway_file_dir = None

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

        # self.random_networks_dir = path.join(self.root_folder, self.random_network_file)


class EnrichTask:
    def __init__(self, name, propagation_file, propagation_folder, statistic_test, target_field,
                 constrain_to_experiment_genes):
        self.name = name
        self.propagation_file = propagation_file
        self.propagation_scores_path = path.join(get_root_path(), propagation_folder)
        self.statistic_test = statistic_test
        self.target_field = target_field
        self.constrain_to_experiment_genes = constrain_to_experiment_genes
        self.results = dict()


class RawScoreTask:
    def __init__(self, name, experiment_file_path, score_file_path, sheet_name, statistic_test, propagation_input_type,
                 constrain_to_network_genes=True):
        self.name = name
        self.score_file_path = score_file_path
        self.sheet_name = sheet_name
        self.statistic_test = statistic_test
        self.propagation_input_type = propagation_input_type
        self.constrain_to_experiment = constrain_to_network_genes
        self.results = dict()
        self.experiment_file_path = experiment_file_path


class GeneralArgs:
    def __init__(self, network_path, genes_names_path, pathway_members_path, FDR_threshold=0.1,
                 output_folder_name=None, figure_name=None, figure_title='Pathway Enrichment '):
        self.minimum_gene_per_pathway = 10
        self.maximum_gene_per_pathway = 50
        self.display_only_significant_pathways = True
        self.network_file_path = network_path
        self.genes_names_file_path = genes_names_path
        self.pathway_databases = ['_']
        self.pathway_keywords = ['_']
        self.significant_pathway_threshold = FDR_threshold
        if output_folder_name is None:
            output_folder_name = 'Enrichment_maps'
        self.output_path = path.join(get_root_path(), 'Outputs', output_folder_name)
        self.figure_name = figure_name if figure_name is not None else 'figure'
        self.pathway_members_path = pathway_members_path
        self.figure_title = figure_title


class PathwayResults:
    def __init__(self, p_value, direction, score):
        self.p_value = p_value
        self.direction = direction
        self.adj_p_value = None
        self.score = score
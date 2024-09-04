import zlib
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from os import path, makedirs
from .settings import Settings
import time
from typing import Union
from io import BytesIO


####################################################################LOAD FUNCTIONS############################################################################

def load_pathways_genes(pathways_dir):
    """
    Loads the pathways and their associated genes from a file.

    Args:
        pathways_dir (str): Path to the file containing the pathway data.

    Returns:
        dict: A dictionary mapping pathway names to lists of genes in each pathway.
    """
    pathways = {}
    # Open the file containing pathway data
    try:
        with open(pathways_dir, 'r') as file:
            for line in file:
                # Process each line, split by tab
                parts = line.strip().split()

                # Skip lines that don't have at least 3 parts
                if len(parts) < 3:
                    continue

                # Parse pathway name
                pathway_name = parts[0]

                # Extract all parts starting from the third part as genes, converting them to integers
                try:
                    genes = [int(gene) for gene in parts[2:]]
                except ValueError:
                    print(f"Non-integer gene ID found in pathway {pathway_name}, skipping pathway.")
                    continue

                pathways[pathway_name] = genes

    except FileNotFoundError:
        print(f"File not found: {pathways_dir}")
    except Exception as e:
        print(f"An error occurred while loading pathways: {e}")

    return pathways



def load_propagation_file(file_path, decompress=True):
    """
    Loads the propagation score data from a file.

    Args:
        file_path (str): The path to the file to be loaded.
        decompress (bool): Whether to decompress the file.
    Returns:
        dict: The loaded data containing propagation scores and other information.
    """
    with open(file_path, 'rb') as file:
        decompressed_data = pickle.load(file)
    if decompress:
        try:
            decompressed_data = pickle.loads(zlib.decompress(decompressed_data))
        except zlib.error:
            print('Entered an uncompressed file but asked to decompress it')
    return decompressed_data


def load_propagation_scores(propagation_file_path):
    propagation_result_dict: dict = load_propagation_file(propagation_file_path, decompress=True)
    return propagation_result_dict


def load_pathways_and_propagation_scores(settings, scores_df):
    """
    Loads the pathways based on the provided configuration and the propagation scores.
    Returns:
        tuple: A tuple containing the network graph, a list of interesting pathways, and a dictionary mapping
               pathways to their genes.
    """
    pathways_with_many_genes = load_pathways_genes(settings.pathway_file_dir)
    sorted_scores = scores_df.sort_values(by='GeneID').reset_index(drop=True)
    scores = {gene_id: (score, pvalue) for gene_id, score, pvalue
                   in zip(sorted_scores['GeneID'], sorted_scores['Score'], sorted_scores['P-value'])}
    if settings.run_gsea:
        return pathways_with_many_genes, scores
    else:
        scores_keys = set(scores.keys())
        # Filter pathways by those having gene counts within the specified range and that intersect with scored genes
        genes_by_pathway = {pathway: set(genes).intersection(scores_keys)
                            for pathway, genes in pathways_with_many_genes.items()
                            if settings.minimum_gene_per_pathway <=
                            len(set(genes).intersection(scores_keys)) <= settings.maximum_gene_per_pathway}

    return genes_by_pathway, scores

####################################################################GET FUNCTIONS############################################################################

def set_input_type(prior_data, norm: bool = False):
    """
    Modifies the 'Score' column based on the specified input type.

    Args:
        prior_data (pandas.DataFrame): DataFrame containing all experimental data.
        norm (bool): Flag indicating whether to normalize the 'Score' column. Default is False.

    Returns:
        pandas.DataFrame: DataFrame with the modified 'Score' column.
    """
    modified_prior_data = prior_data.copy()

    if norm:
        modified_prior_data['Score'] = 1
    else:
        modified_prior_data['Score'] = prior_data['Score'].abs()

    return modified_prior_data


def get_scores(score_path):
    """
        Loads gene scores and P-values from a file and returns a dictionary mapping gene IDs to their scores and P-values.

        Args:
            score_path (str): The path to the file containing the propagation scores.

        Returns:
            dict: A dictionary with gene IDs as keys and tuples of (Score, P-value) as values.
        """
    try:
        # Load propagation results from the specified file
        propagation_results = load_propagation_scores(score_path)

        # Sort the propagation scores by GeneID
        sorted_scores = propagation_results['gene_prop_scores'].sort_values(by='GeneID').reset_index(drop=True)

        # Create a dictionary mapping gene IDs to tuples of (Score, P-value)
        gene_scores = {gene_id: (score, pvalue) for gene_id, score, pvalue
                       in zip(sorted_scores['GeneID'], sorted_scores['Score'], sorted_scores['P-value'])}

        return gene_scores

    except FileNotFoundError:
        print(f"File not found: {score_path}")
        return {}
    except Exception as error:
        print(f"An error occurred: {error}")
        return {}

####################################################################READ FUNCTIONS############################################################################

def read_temp_scores(file_name):
    """
    Read scores from a file into a dictionary.

    Parameters:
    - file_name (str): Path to the file containing the scores.

    Returns:
    dict: A dictionary mapping pathways to their scores.
    """
    try:
        scores = pd.read_csv(file_name, sep=' ', header=None, names=['Pathway', 'Score'], index_col='Pathway')
        return scores['Score'].to_dict()
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return {}
    except Exception as e:
        print(f"Error reading scores from {file_name}: {e}")
        return {}


def read_network(network_filename: str) -> nx.Graph:
    """
    Read a network from a file and return a NetworkX graph.

    Parameters:
    - network_filename (str): Path to the file containing the network data.

    Returns:
    - nx.Graph: A graph object representing the network.
    """
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


def read_prior_set(excel_input: Union[str, BytesIO], is_bytes: bool = False) -> pd.DataFrame:
    """
    Read prior data set from an Excel file and apply preprocessing.

    Parameters:
    - excel_input (str or BytesIO): Path to the Excel file or a BytesIO object containing the prior data.
    - is_bytes (bool): Flag indicating if the input is a BytesIO object. Default is False.

    Returns:
    - pd.DataFrame: DataFrame containing the preprocessed prior data.
    """
    if is_bytes:
        prior_data = pd.read_excel(excel_input, engine='openpyxl')
    else:
        prior_data = pd.read_excel(excel_input, engine='openpyxl')

    # Identify duplicate GeneID values
    duplicate_gene_ids = prior_data[prior_data.duplicated(subset='GeneID', keep=False)]
    if not duplicate_gene_ids.empty:
        print("Duplicate GeneID values found:")
        print(duplicate_gene_ids)

    # Drop duplicate GeneID values
    prior_data = prior_data.drop_duplicates(subset='GeneID')

    # Remove any row with no value in the Score column
    prior_data = prior_data[prior_data['Score'].notna()]

    # Remove any row with "?" in the Score column
    prior_data = prior_data[prior_data['Score'] != '?']

    # Filter out GeneIDs that are not purely numeric (to exclude concatenated IDs)
    prior_data = prior_data[prior_data['GeneID'].apply(lambda x: str(x).isdigit())]

    # Convert GeneID to integer
    prior_data['GeneID'] = prior_data['GeneID'].astype(int)

    # Reset the DataFrame index
    prior_data = prior_data.reset_index(drop=True)

    return prior_data

####################################################################SAVE FUNCTIONS############################################################################

def save_file(obj, save_dir=None, compress=True):
    """
    Saves an object to a file, with optional compression.
    Args:
        obj (object): The object to be saved.
        save_dir (str, optional): The directory where the file will be saved.
        compress (bool, optional): Whether to compress the file.
    Returns:
        None
    """
    obj = pickle.dumps(obj)
    if compress:
        obj = zlib.compress(obj)
    with open(save_dir, 'wb') as f:
        pickle.dump(obj, f)


def save_propagation_score(propagation_scores: pd.DataFrame, prior_set: pd.DataFrame,
                           propagation_input: dict, genes_id_to_idx: dict, task,
                           save_dir: str, general_args: Settings) -> None:
    """
    Save the propagation scores to a compressed file.

    Parameters:
    - propagation_scores (pd.DataFrame): DataFrame containing the propagated scores.
    - prior_set (pd.DataFrame): DataFrame containing the prior set of gene scores.
    - propagation_input (dict): Dictionary mapping gene IDs to their input scores.
    - genes_id_to_idx (dict): Dictionary mapping gene IDs to their indices.
    - task (PropagationTask): Propagation task containing task-specific settings.
    - save_dir (str): Directory where the output file will be saved.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    file_name = f"{general_args.alpha}_{general_args.date}"
    save_dir = save_dir
    makedirs(save_dir, exist_ok=True)
    propagation_results_path = path.join(save_dir, file_name)

    save_dict = {
        'args': task, 'prior_set': prior_set, 'propagation_input': propagation_input,
        'gene_id_to_idx': genes_id_to_idx, 'gene_prop_scores': propagation_scores,
    }

    save_file(save_dict, propagation_results_path)

####################################################################ELSE######################################################

def process_condition(condition_file, experiment_file, pathways_file, all_pathways, P_VALUE_THRESHOLD=0.05):
    """
    Processes experimental data to determine the enrichment and trends of pathways based on specified conditions.

    Parameters:
    - condition_file (str): Path to the file containing score conditions.
    - experiment_file (str): Path to the file containing experimental data.
    - pathways_file (str): Path to the file containing pathways data.

    Returns:
    dict: A dictionary containing scores, enriched pathway genes, and mean scores for pathways under the given condition.
    """
    # Read scores for the condition, mapping pathways to their scores
    enriched_pathway_dict = read_temp_scores(condition_file)

    # Load experiment data and filter out entries where the score is zero.
    condition_data_df = pd.read_excel(experiment_file)
    experiment_data_filtered_df = condition_data_df[condition_data_df['Score'] != 0]

    # Extract the condition name from the file name
    condition_name = path.basename(condition_file).split('.')[-2]

    # Load pathway data mapping pathway names to lists of gene IDs
    homo_sapien_pathway_dict = load_pathways_genes(pathways_file)

    # Dictionary to store enriched pathway genes
    enriched_pathway_genes = {}

    # Loop through each pathway
    for pathway in all_pathways:
        # Initialize a dictionary for the pathway under the current condition
        all_pathways[pathway][condition_name] = {}

        # List of genes associated with the current pathway
        pathway_genes = homo_sapien_pathway_dict[pathway]

        # Ensure pathway_genes are integers
        pathway_genes = [int(gene) for gene in pathway_genes]

        # Filter the experiment data to only include genes that are part of the current pathway
        pathway_filtered_genes = experiment_data_filtered_df[experiment_data_filtered_df['GeneID'].isin(pathway_genes)]

        # Store details of filtered genes in a dictionary
        enriched_pathway_genes[pathway] = pathway_filtered_genes.set_index('GeneID')[['Symbol', 'Score', 'P-value']].to_dict(
            orient='index')

        # Filter to find significant genes based on the P-value threshold
        significant_genes = {gene_id: gene_details for gene_id, gene_details in
                             enriched_pathway_genes[pathway].items() if
                             gene_details['P-value'] <= P_VALUE_THRESHOLD}

        # Calculate the mean score of significant genes or set to 0 if none are significant
        mean_score = np.mean(
            [gene_details['Score'] for gene_details in significant_genes.values()]) if significant_genes else 0

        # Store the mean score and significant genes for the pathway under the condition
        all_pathways[pathway][condition_name]['Mean'] = mean_score
        all_pathways[pathway][condition_name]['significant_genes'] = significant_genes

        # Check if the pathway is in the enriched pathway dictionary to assign a P-value and trend
        if pathway in enriched_pathway_dict:
            all_pathways[pathway][condition_name]['P-value'] = enriched_pathway_dict[pathway]
            all_pathways[pathway][condition_name]['Trend'] = "Up*" if mean_score > 0 else "Down*"
        else:
            all_pathways[pathway][condition_name]['P-value'] = 1  # Default P-value if not in enriched dictionary
            all_pathways[pathway][condition_name]['Trend'] = "Up" if mean_score > 0 else "Down"



def filter_network_by_prior_data(network_filename: str, prior_data: pd.DataFrame) -> nx.Graph:
    """
    Filter a network to only include nodes present in the prior_data DataFrame.

    Parameters:
    - network_filename (str): Path to the network file.
    - prior_data (pd.DataFrame): DataFrame containing gene information.

    Returns:
    - nx.Graph: A filtered graph object.
    """
    # Read the network
    network = read_network(network_filename)

    # Get the set of gene IDs from prior_data
    gene_ids_in_prior_data = set(prior_data['GeneID'])

    # Filter the network
    filtered_network = network.copy()
    for node in network.nodes:
        if node not in gene_ids_in_prior_data:
            filtered_network.remove_node(node)

    return filtered_network



def calculate_time(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)  # Store the result
        end_time = time.time()
        print(f"Execution time: {end_time - start_time:.2f} seconds")
        return result  # Return the result
    return wrapper
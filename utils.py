import pandas as pd
import networkx as nx
import numpy as np
import os
from os import path
import pickle
import zlib


# noinspection PyTypeChecker
def read_network(network_filename):
    """
    Reads a network from a file and returns a NetworkX graph.

    Args:
        network_filename (str): Path to the file containing the network data.

    Returns:
        networkx.Graph: A graph object representing the network.
    """
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


def read_prior_set(excel_dir):
    """
    Reads prior data set from an Excel file.
    Args:
        excel_dir (str): Path to the Excel file containing the prior data.
    Returns:
        pandas.DataFrame: DataFrame containing the prior data.
    """
    prior_data = pd.read_excel(excel_dir, engine='openpyxl')
    # Apply data preparation logic
    prior_data = prior_data[prior_data['GeneID'].apply(lambda x: str(x).isdigit())]
    prior_data['GeneID'] = prior_data['GeneID'].astype(int)
    prior_data = prior_data.reset_index(drop=True)
    return prior_data


def get_propagation_input(prior_gene_ids, prior_data, input_type):
    """
    Generates propagation inputs based on specified type and network.

    Args:
        prior_gene_ids (set): List of gene IDs for the prior set.
        prior_data (pandas.DataFrame): DataFrame containing all experimental data.
        input_type (str): Type of input to generate (e.g., 'ones', 'abs_Score', etc.).

    Returns:
        dict: A dictionary mapping gene IDs to their corresponding input values.
    """
    inputs = dict()
    if input_type == 'ones' or input_type is None:
        inputs = {int(x): 1 for x in prior_gene_ids}
    elif input_type == 'abs_Score':
        for id in prior_gene_ids:
            try:
                # Your original code for processing each id
                inputs[int(id)] = np.abs(float(prior_data[prior_data.GeneID == id]['Score'].values[0]))
            except TypeError:
                print(f"Error processing ID: {id}")
                print(f"Data for this ID: {prior_data[prior_data.GeneID == id]}")
                # Optionally, you can break or continue based on your needs
                break

    elif input_type == 'Score':
        inputs = {int(id): float(prior_data[prior_data.GeneID == id]['Score'].values[0]) for id in prior_gene_ids}
    else:
        assert 0, '{} is not a valid input type'.format(input_type)

    inputs = {id: np.round(input_score, 3) for id, input_score in inputs.items()}
    return inputs


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
    print('File was saved in {}'.format(save_dir))


def save_propagation_score(propagation_scores, prior_set, propagation_input, genes_id_to_idx, task, save_dir=None):
    """
    Saves the propagation scores to a file.

    Args:
        propagation_scores (dict): The propagation scores to be saved.
        prior_set (dataframe): The set of prior genes.
        propagation_input (dict): The input used for propagation.
        genes_idx_to_id (dict): Mapping from gene indices to gene IDs.
        task (PropagationTask): The propagation task object containing program arguments.
        save_dir (str, optional): Directory to save the results.

    Returns:
        dict: A dictionary containing the saved data.
    """
    file_name = f"{task.experiment_name}_{task.propagation_input_type}_{task.alpha}_{task.date}"
    save_dir = save_dir or task.propagation_scores_path
    os.makedirs(save_dir, exist_ok=True)
    propagation_results_path = path.join(save_dir, file_name)

    save_dict = {
        'args': task, 'prior_set': prior_set, 'propagation_input': propagation_input,
        'gene_id_to_idx': genes_id_to_idx, 'gene_prop_scores': propagation_scores,
    }

    save_file(save_dict, propagation_results_path)
    return save_dict


def get_root_path():
    """
    Retrieves the root path of the current script.
    Returns:
        str: The root directory path of the current script file.
    """
    return path.dirname(path.realpath(__file__))


def load_pathways_genes(pathways_dir):
    """
    Loads the pathways and their associated genes from a file.

    Args:
        pathways_dir (str): Path to the file containing the pathway data.

    Returns:
        dict: A dictionary mapping pathway names to lists of genes in each pathway.
    """
    with open(pathways_dir, 'r') as f:
        lines = [str.upper(x.strip()).split('\t') for x in f]
    pathways = {x[0]: [int(y) for y in x[2:]] for x in lines}

    return pathways


def load_file(load_dir, decompress=True):
    """
    Loads an object from a file using pickle with optional decompression.

    Args:
        load_dir (str): The directory from which to load the file.
        decompress (bool, optional): Whether to decompress the file.

    Returns:
        object: The object loaded from the file.
    """
    with open(load_dir, 'rb') as f:
        file = pickle.load(f)
    if decompress:
        try:
            file = pickle.loads(zlib.decompress(file))
        except zlib.error:
            print('Entered an uncompressed file but asked to decompress it')
    return file


def load_propagation_scores(task, propagation_file_name=None):
    """
    Loads the propagation scores from a file.

    Args:
        task (PropagationTask): The propagation task object containing program arguments.
        propagation_file_name (str, optional): Name of the file to load the results from.

    Returns:
        dict: Dictionary containing the propagation scores, input, and gene index to ID mapping.
    """
    propagation_results_path = path.join(task.propagation_scores_path, propagation_file_name)
    propagation_result_dict: dict = load_file(propagation_results_path, decompress=True)

    return propagation_result_dict

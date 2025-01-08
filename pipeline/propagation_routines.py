import os
import time
import scipy.sparse as sp
from scipy.sparse.linalg import inv
import numpy as np
import networkx as nx
from args import GeneralArgs, PropagationTask
from utils import save_propagation_score, set_input_type
from tqdm import tqdm
import pandas as pd
import logging

# Configure logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

def calculate_ngsea(network: nx.Graph, prior_data: pd.DataFrame) -> dict:
    """
    Calculates the score for each gene based on the provided equation.

    Parameters:
    - network (networkx.Graph): The network graph.
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - dict: Dictionary mapping gene IDs to their calculated scores.
    """
    # Initialize all network genes with a score of 0
    gene_scores = {gene_id: 0 for gene_id in network.nodes()}

    # Create a dictionary from prior_data
    prior_data_dict = prior_data.set_index('GeneID')['Score'].to_dict()

    # Assign scores to all genes in prior_data
    for gene_id, score in prior_data_dict.items():
        gene_scores[gene_id] = abs(score)

    # Calculate the scores for genes in the prior data and in the network based on their neighbors
    for gene_id in prior_data_dict.keys():
        if gene_id in network.nodes():
            neighbors = list(network.neighbors(gene_id))

            if neighbors:
                neighbor_scores = [abs(prior_data_dict.get(neighbor, 0)) for neighbor in neighbors]
                ni = len(neighbor_scores)
                sum_neighbor_scores = sum(neighbor_scores)
                gene_scores[gene_id] += (1 / ni) * sum_neighbor_scores

    return gene_scores


def generate_similarity_matrix(network: nx.Graph, args: GeneralArgs) -> tuple:
    """
    Generates and saves a similarity matrix for network propagation, based on the provided network graph.

    Parameters:
    - network (nx.Graph or sp.sparse matrix): Network graph or sparse matrix.
    - args (GeneralArgs): Arguments related to the general configuration of the experiment.

    Returns:
    - tuple: A tuple containing the inverse of the similarity matrix and the list of genes.
    """
    start = time.time()

    genes = sorted(network.nodes())
    gene_index = {gene: index for index, gene in enumerate(genes)}

    # Convert network to a sparse matrix if it isn't one already
    if not sp.isspmatrix(network):
        row, col, data = [], [], []
        for edge in network.edges(data=True):
            weight = edge[2].get(2, 1.0)  # Use the edge weight if provided, otherwise default to 1.0
            row.append(gene_index[edge[0]])
            col.append(gene_index[edge[1]])
            data.append(weight)
            # Add the reverse edge since the graph is undirected
            row.append(gene_index[edge[1]])
            col.append(gene_index[edge[0]])
            data.append(weight)
        matrix = sp.csr_matrix((data, (row, col)), shape=(len(genes), len(genes)))
    else:
        matrix = network

    degrees = np.array(matrix.sum(axis=1)).flatten()
    # Replace zero degrees with 1 to avoid division by zero
    degrees[degrees == 0] = 1


    # Symmetric normalization: A' = D^(-1/2) * A * D^(-1/2)
    inv_sqrt_degrees = 1.0 / np.sqrt(degrees)
    norm_matrix = sp.diags(inv_sqrt_degrees)
    matrix = norm_matrix @ matrix @ norm_matrix



    n = matrix.shape[0]
    identity_matrix = sp.eye(n)
    matrix_to_invert = identity_matrix - (1 - args.alpha) * matrix
    matrix_to_invert_csc = matrix_to_invert.tocsc()

    print("Inverting the matrix")
    inverse_matrix = inv(matrix_to_invert_csc)
    print("Matrix inverted")
    inverse_matrix = args.alpha * inverse_matrix

    print("Saving the matrix")
    if not os.path.exists(os.path.dirname(args.similarity_matrix_path)):
        os.makedirs(os.path.dirname(args.similarity_matrix_path))
    sp.save_npz(args.similarity_matrix_path, inverse_matrix)

    end = time.time()
    print(f"Time elapsed: {end - start} seconds")
    return inverse_matrix, gene_index


def read_sparse_matrix_txt(network: nx.Graph, similarity_matrix_path: str) -> tuple:
    """
    Reads a precomputed sparse similarity matrix from a file.

    Parameters:
    - network (nx.Graph): Network graph used to generate the similarity matrix.
    - similarity_matrix_path (str): Path to the file containing the sparse matrix.

    Returns:
    - tuple: A tuple containing the sparse matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = {gene: index for index, gene in enumerate(genes)}

    matrix = sp.load_npz(similarity_matrix_path)
    return matrix, gene_index


def matrix_prop(propagation_input: dict, gene_indexes: dict, matrix=None) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - matrix (sp.spmatrix): Precomputed inverse matrix for propagation.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    num_genes = len(gene_indexes)
    F_0 = np.zeros(num_genes)  # Changed to a 1D array
    for gene_id, value in propagation_input.items():
        F_0[gene_indexes[gene_id]] = value
    F = matrix @ F_0

    return F

def _handle_ngsea_case(prior_data, prop_task, general_args, network):
    """
    Handle the GSE case for alpha = 1.

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.
    - network (networkx.Graph): The network graph.

    Returns:
    - None
    """
    # Calculate gene scores for genes in both prior_data and network
    gene_scores = calculate_ngsea(network, prior_data)

    # Filter the posterior set to include only genes in both prior_data and network
    posterior_set = prior_data[prior_data['GeneID'].isin(network.nodes())].copy()
    posterior_set['Score'] = posterior_set['GeneID'].map(gene_scores)

    save_propagation_score(
        prior_set=prior_data,
        propagation_input={gene_id: score for gene_id, score in zip(posterior_set['GeneID'], posterior_set['Score'])},
        propagation_scores=posterior_set,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(posterior_set['GeneID'])},
        task=prop_task,
        save_dir=prop_task.output_folder,
        general_args=general_args
    )


def _handle_no_propagation_case(prior_data, prop_task, general_args, network):
    """
    Handle the case where no propagation is performed (alpha = 1).

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.
    - network (networkx.Graph): The network graph.

    Returns:
    - None
    """
    # Filter prior_data to include only genes in both prior_data and network
    filtered_prior_data = prior_data[prior_data['GeneID'].isin(network.nodes())].copy()

    # Sort the filtered_prior_data by GeneID
    sorted_prior_data = filtered_prior_data.sort_values(by='GeneID').reset_index(drop=True)

    # Flatten the scores for saving
    gene_scores = sorted_prior_data['Score'].values.reshape((len(sorted_prior_data), 1))
    sorted_prior_data['Score'] = gene_scores.flatten()

    save_propagation_score(
        prior_set=sorted_prior_data,
        propagation_input={gene_id: score for gene_id, score in
                           zip(sorted_prior_data['GeneID'], sorted_prior_data['Score'])},
        propagation_scores=sorted_prior_data,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(sorted_prior_data['GeneID'])},
        task=prop_task,
        save_dir=prop_task.output_folder,
        general_args=general_args
    )



def _normalize_prop_scores(similarity_matrix, gene_to_index_map, raw_prop_scores) -> pd.DataFrame:
    """
    Normalize the propagation scores.

    Parameters:
    - similarity_matrix (sp.sparse.spmatrix or np.ndarray): The similarity matrix or upper triangular matrix.
    - gene_to_index_map (dict): Mapping of gene IDs to their indices in the matrix.
    - raw_prop_scores (np.ndarray): Array of propagation scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    # Set all network nodes to a score of 1
    all_nodes = list(gene_to_index_map.keys())  # All nodes in the network
    ones_scores_dict = {node: 1 for node in all_nodes}

    # Perform propagation with ones input
    normalization_scores = matrix_prop(ones_scores_dict, gene_to_index_map, matrix=similarity_matrix)
    
    # Identify genes with zero normalization scores
    zero_norm_indices = np.nonzero(normalization_scores == 0)[0]
    zero_prop_indices = np.nonzero(raw_prop_scores == 0)[0]
    indices_to_reset = list(set(zero_norm_indices).difference(zero_prop_indices))
    normalization_scores[indices_to_reset] = 1

    # Normalize non-zero propagation scores
    nonzero_prop_indices = np.nonzero(raw_prop_scores != 0)[0]
    raw_prop_scores[nonzero_prop_indices] /= np.abs(normalization_scores[nonzero_prop_indices])

    # Create DataFrame with all network genes
    network_genes = list(gene_to_index_map.keys())
    normalized_scores = raw_prop_scores[np.array([gene_to_index_map[gene_id] for gene_id in network_genes])]

    return pd.DataFrame({
        'GeneID': network_genes,
        'Score': normalized_scores
    })



def _save_propagation_results(propagation_input_df, full_propagated_scores_df, prop_task, general_args):
    """
    Save the results of the propagation process.

    Parameters:
    - propagation_input_df (pandas.DataFrame): DataFrame containing the modified input data.
    - full_propagated_scores_df (pandas.DataFrame): DataFrame containing full propagated scores.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    save_propagation_score(
        prior_set=propagation_input_df,
        propagation_input={gene_id: score for gene_id, score in
                           zip(full_propagated_scores_df['GeneID'], full_propagated_scores_df['Score'])},
        propagation_scores=full_propagated_scores_df,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(full_propagated_scores_df['GeneID'])},
        task=prop_task,
        save_dir=prop_task.output_folder,
        general_args=general_args
    )

def filter_network_genes(propagation_input_df, network, mode="both"):
    """
    Filters or combines genes based on the specified imputation mode.

    Args:
        propagation_input_df (pd.DataFrame): Input genes with scores.
        network (nx.Graph): Network graph.
        mode (str): Imputation mode: "both", "dataset", "network", "all".

    Returns:
        pd.DataFrame: Filtered genes with scores.
    """
    network_genes_df = propagation_input_df[propagation_input_df['GeneID'].isin(network.nodes())].copy()
    non_network_genes_df = propagation_input_df[~propagation_input_df['GeneID'].isin(network.nodes())].copy()

    if mode == "both":
        return network_genes_df
    elif mode == "dataset":
        return pd.concat([network_genes_df, non_network_genes_df])
    elif mode == "network":
        # Impute missing scores for genes in network but not in dataset
        missing_genes = set(network.nodes()) - set(propagation_input_df['GeneID'])
        imputed_genes = pd.DataFrame({'GeneID': list(missing_genes), 'Score': 0})  # Default score = 0
        return pd.concat([network_genes_df, imputed_genes])
    elif mode == "all":
        # Include all genes (network + dataset) with propagated/imputed scores
        missing_genes = set(network.nodes()) - set(propagation_input_df['GeneID'])
        imputed_genes = pd.DataFrame({'GeneID': list(missing_genes), 'Score': 0})  # Default score = 0
        return pd.concat([network_genes_df, non_network_genes_df, imputed_genes])

def get_similarity_matrix(network, general_args):
    if general_args.create_similarity_matrix:
        return generate_similarity_matrix(network, general_args)
    else:
        return read_sparse_matrix_txt(network, general_args.similarity_matrix_path)


def handle_no_propagation_cases(prior_data, prop_task, general_args, network):
    if general_args.run_NGSEA:
        _handle_ngsea_case(prior_data, prop_task, general_args, network)
    else:
        _handle_no_propagation_case(prior_data, prop_task, general_args, network)


def perform_propagation(test_name: str, general_args, network=None, prior_data=None):
    """
    Performs the propagation of gene scores through the network.

    Parameters:
    - test_name (str): Name of the test for which propagation is performed.
    - general_args: General arguments and settings.
    - network (networkx.Graph): The network graph.
    - prior_data (pandas.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - None
    """
    prop_task = PropagationTask(general_args=general_args, test_name=test_name)
    
    # Modify prior_data based on the input type
    if general_args.input_type == 'Abs_Score':
        propagation_input_df = set_input_type(prior_data, general_args.input_type)
    else:
        propagation_input_df = prior_data

    if general_args.alpha == 1:
        handle_no_propagation_cases(propagation_input_df, prop_task, general_args, network)
        return

    # Generate or load the similarity matrix
    matrix, network_gene_index = get_similarity_matrix(network, general_args)

    # Filter genes in the network (used for propagation input)
    filtered_propagation_input = {
        gene_id: score
        for gene_id, score in zip(propagation_input_df['GeneID'], propagation_input_df['Score'])
        if gene_id in network.nodes()
    }

    # Perform propagation
    propagation_score = matrix_prop(filtered_propagation_input, network_gene_index, matrix=matrix)

    # Normalize scores
    normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score)

    # Handle post-normalization based on imputation mode
    if general_args.imputation_mode == "both":
        # Keep only genes in both dataset and network
        final_scores_df = normalized_df.merge(
            prior_data[['GeneID']], on='GeneID', how='inner'
        )
    elif general_args.imputation_mode == "dataset":
        # Combine propagated scores for "both" and original scores for dataset-only genes
        dataset_only_genes = prior_data[~prior_data['GeneID'].isin(network.nodes())].copy()
        final_scores_df = pd.concat([normalized_df, dataset_only_genes])
    elif general_args.imputation_mode == "network":
        # Combine propagated scores for "both" and imputed scores for network-only genes
        network_only_genes = normalized_df[~normalized_df['GeneID'].isin(prior_data['GeneID'])].copy()
        final_scores_df = pd.concat([normalized_df, network_only_genes])
    elif general_args.imputation_mode == "all":
        # Keep all genes (both propagated and original scores)
        dataset_only_genes = prior_data[~prior_data['GeneID'].isin(network.nodes())].copy()
        network_only_genes = normalized_df[~normalized_df['GeneID'].isin(prior_data['GeneID'])].copy()
        final_scores_df = pd.concat([normalized_df, dataset_only_genes, network_only_genes])

    # Save the results
    _save_propagation_results(propagation_input_df, final_scores_df, prop_task, general_args)

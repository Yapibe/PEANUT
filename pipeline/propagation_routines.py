import os
import time
import scipy.sparse as sp
from scipy.sparse.linalg import inv
import json
import numpy as np
import networkx as nx
from args import GeneralArgs, PropagationTask
from utils import save_propagation_score, set_input_type
from tqdm import tqdm
import pandas as pd


def calculate_gene_scores(network: nx.Graph, prior_data: pd.DataFrame) -> dict:
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
    for gene_id in tqdm(prior_data_dict.keys(), desc="Calculating gene scores"):
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
    inv_sqrt_degrees = 1.0 / np.sqrt(degrees)
    norm_matrix = sp.diags(inv_sqrt_degrees)
    matrix = norm_matrix @ matrix @ norm_matrix

    n = matrix.shape[0]
    identity_matrix = sp.eye(n)
    matrix_to_invert = identity_matrix - (1 - args.alpha) * matrix
    matrix_to_invert_csc = matrix_to_invert.tocsc()

    print("Inverting the matrix")
    inverse_matrix = inv(matrix_to_invert_csc)
    inverse_matrix = args.alpha * inverse_matrix

    # Extract the upper triangular part of the inverse matrix
    upper_tri_indices = np.triu_indices(n)
    upper_tri_inverse_matrix = inverse_matrix[upper_tri_indices]

    # Print densities for debugging
    original_density = matrix.nnz / (matrix.shape[0] * matrix.shape[1])
    inverse_density = inverse_matrix.nnz / (inverse_matrix.shape[0] * inverse_matrix.shape[1])
    print(f"Original matrix density: {original_density}")
    print(f"Inverse matrix density: {inverse_density}")

    print("Saving the matrix")
    if not os.path.exists(os.path.dirname(args.similarity_matrix_path)):
        os.makedirs(os.path.dirname(args.similarity_matrix_path))
    sp.save_npz(args.similarity_matrix_path, inverse_matrix)
    save_upper_path = args.tri_similarity_matrix_path
    np.save(save_upper_path, upper_tri_inverse_matrix)

    end = time.time()
    print(f"Time elapsed: {end - start} seconds")
    return inverse_matrix, gene_index


def read_sparse_matrix_txt(network: nx.Graph, similarity_matrix_path: str, tri_matrix_path: str, debug: bool) -> tuple:
    """
    Reads a precomputed sparse similarity matrix from a file.

    Parameters:
    - network (nx.Graph): Network graph used to generate the similarity matrix.
    - similarity_matrix_path (str): Path to the file containing the sparse matrix.
    - tri_matrix_path (str): Path to the file containing the upper triangular part of the sparse matrix.
    - debug (bool): Flag to indicate whether to use the full matrix or the upper triangular part.

    Returns:
    - tuple: A tuple containing the sparse matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = {gene: index for index, gene in enumerate(genes)}

    if not os.path.exists(similarity_matrix_path) and not os.path.exists(tri_matrix_path):
        raise FileNotFoundError(f"Neither specified file {similarity_matrix_path} nor {tri_matrix_path} exists.")

    if debug:
        upper_tri_inverse_matrix = np.load(tri_matrix_path)
        return upper_tri_inverse_matrix, gene_index
    else:
        matrix = sp.load_npz(similarity_matrix_path)
        return matrix, gene_index

def symmetric_matrix_vector_multiply(upper_tri_matrix: np.ndarray, F_0: np.ndarray) -> np.ndarray:
    """
    Performs matrix-vector multiplication using only the upper triangular part of a symmetric matrix.

    Parameters:
    - upper_tri_matrix (np.ndarray): The upper triangular part of the symmetric matrix stored as a 2D array with dimensions (1, number of cells in the upper triangle).
    - F_0 (np.ndarray): The vector to be multiplied.

    Returns:
    - np.ndarray: The result of the matrix-vector multiplication.
    """
    num_genes = F_0.size
    result = np.zeros(num_genes)
    index = 0

    upper_tri_matrix = upper_tri_matrix.flatten()  # Ensure it is a 1D array

    for i in range(num_genes):
        for j in range(i, num_genes):
            result[i] += upper_tri_matrix[index] * F_0[j]
            if i != j:
                result[j] += upper_tri_matrix[index] * F_0[i]
            index += 1

    return result

def matrix_prop(propagation_input: dict, gene_indexes: dict, debug: bool, inverse_matrix=None) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - inverse_matrix (sp.spmatrix): Precomputed inverse matrix for propagation.
    - upper_tri_inverse_matrix (np.ndarray): Upper triangular part of the inverse matrix.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - debug (bool): Flag to indicate whether to use the full matrix or the upper triangular part.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    num_genes = len(gene_indexes)
    F_0 = np.zeros(num_genes)  # Changed to a 1D array
    for gene_id, value in propagation_input.items():
        F_0[gene_indexes[gene_id]] = value

    if debug:
        F = symmetric_matrix_vector_multiply(inverse_matrix, F_0)
    else:
        F = inverse_matrix @ F_0
    # # Compare F_full and F_upper_tri
    # difference = np.abs(F_full - F_upper_tri)
    # threshold = 1e-5  # Define a small threshold for floating-point comparison
    # num_differences = np.sum(difference > threshold)
    #
    # print(f"Number of different elements: {num_differences}")
    # print(f"Differences: {difference[difference > threshold]}")

    return F



def _handle_ngsea_case(prior_data, prop_task, general_args, network):
    """
    Handle the GSE case for alpha = 1.

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    gene_scores = calculate_gene_scores(network, prior_data)
    posterior_set = prior_data.copy()
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


def _handle_no_propagation_case(prior_data, prop_task, general_args):
    sorted_prior_data = prior_data.sort_values(by='GeneID').reset_index(drop=True)
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


def _normalize_prop_scores(matrix, network_gene_index, propagation_score, filtered_prior_data, debug: bool) -> pd.DataFrame:
    """
    Normalize the propagation scores.

    Parameters:
    - matrix (sp.sparse.spmatrix or np.ndarray): The similarity matrix or upper triangular matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices in the matrix.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.
    - debug (bool): Flag to indicate whether to use the full matrix or the upper triangular part.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    # Set input type to 'ones' and convert to dictionary
    ones_input_df = set_input_type(filtered_prior_data, 'ones')
    ones_input = ones_input_df.set_index('GeneID')['Score'].to_dict()

    # Perform propagation with ones input
    ones_gene_scores_inverse = matrix_prop(ones_input, network_gene_index, debug, inverse_matrix=matrix)

    # Identify genes with zero normalization scores
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))
    ones_gene_scores_inverse[genes_to_delete] = 1

    # Adjust the normalization scores
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] /= np.abs(ones_gene_scores_inverse[non_zero_indices])

    # Create a DataFrame for all network genes
    all_network_genes = list(network_gene_index.keys())
    all_normalized_scores = propagation_score[np.array([network_gene_index[gene_id] for gene_id in all_network_genes])]

    # Create DataFrame for normalized scores
    normalized_df = pd.DataFrame({
        'GeneID': all_network_genes,
        'Score': all_normalized_scores
    })

    return normalized_df



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


def merge_with_prior_data(final_propagation_results, prior_data, gene_name_dict):
    full_propagated_scores_df = final_propagation_results.merge(prior_data[['GeneID', 'Symbol', 'P-value']], on='GeneID', how='left')
    full_propagated_scores_df['Symbol'] = full_propagated_scores_df['Symbol'].fillna(full_propagated_scores_df['GeneID'].map(gene_name_dict).fillna(""))
    full_propagated_scores_df['P-value'] = full_propagated_scores_df['P-value'].fillna(0)
    return full_propagated_scores_df


def filter_network_genes(propagation_input_df, network):
    network_genes_df = propagation_input_df[propagation_input_df['GeneID'].isin(network.nodes())].copy()
    non_network_genes = propagation_input_df[~propagation_input_df['GeneID'].isin(network.nodes())].copy()
    filtered_propagation_input = {gene_id: score for gene_id, score in zip(network_genes_df['GeneID'], network_genes_df['Score'])}
    return network_genes_df, non_network_genes, filtered_propagation_input


def get_similarity_matrix(network, general_args):
    if general_args.create_similarity_matrix:
        return generate_similarity_matrix(network, general_args)
    else:
        return read_sparse_matrix_txt(network, general_args.similarity_matrix_path, general_args.tri_similarity_matrix_path, general_args.debug)


def handle_no_propagation_cases(prior_data, prop_task, general_args, network):
    if general_args.run_NGSEA:
        print("Running GSE")
        _handle_ngsea_case(prior_data, prop_task, general_args, network)
    else:
        print("Running GSEA")
        _handle_no_propagation_case(prior_data, prop_task, general_args)


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

    if general_args.alpha == 1:
        handle_no_propagation_cases(prior_data, prop_task, general_args, network)
        return

    print(f"Running propagation with '{general_args.input_type}'")

    matrix, network_gene_index = get_similarity_matrix(network, general_args)

    # Modify prior_data based on the input type
    propagation_input_df = set_input_type(prior_data, general_args.input_type)

    # Filter genes in the network
    network_genes_df, non_network_genes, filtered_propagation_input = filter_network_genes(propagation_input_df, network)

    # Perform network propagation
    propagation_score = matrix_prop(filtered_propagation_input, network_gene_index, general_args.debug,
                                    inverse_matrix=matrix)

    # Normalize the propagation scores and create DataFrame within the function
    normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score, network_genes_df,
                                           general_args.debug)

    # Combine network and non-network genes
    final_propagation_results = pd.concat([
        normalized_df[['GeneID', 'Score']],
        non_network_genes[['GeneID', 'Score']]
    ], ignore_index=True)

    # Load the gene_info.json file
    with open(general_args.genes_names_file_path, 'r') as f:
        gene_name_dict = json.load(f)
    full_propagated_scores_df = merge_with_prior_data(final_propagation_results, prior_data, gene_name_dict)

    # Save the results
    _save_propagation_results(propagation_input_df, full_propagated_scores_df, prop_task, general_args)
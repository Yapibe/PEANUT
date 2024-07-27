import os
import time
import json
import scipy.sparse as sp
from scipy.sparse.linalg import inv
import pandas as pd
import numpy as np
import networkx as nx
from args import GeneralArgs, PropagationTask
from utils import save_propagation_score, read_prior_set, read_network, set_input_type


def matrix_prop(propagation_input: dict, inverse_matrix: sp.spmatrix, gene_indexes: dict) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.

    Parameters:
    - seeds (list): List of seed gene IDs.
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - inverse_matrix (sp.spmatrix): Precomputed inverse matrix for propagation.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - num_genes (int): Total number of genes in the network.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    seeds = list(propagation_input.keys())
    num_genes = len(gene_indexes)
    F_0 = np.zeros(num_genes)  # Changed to a 1D array
    for seed in seeds:
        F_0[gene_indexes[seed]] = propagation_input[seed]

    F = inverse_matrix @ F_0
    return F


def generate_similarity_matrix(network: nx.Graph, args: GeneralArgs) -> tuple:
    """
    Generates and saves a similarity matrix for network propagation, based on the provided network graph.

    Parameters:
    - network (nx.Graph or sp.sparse matrix): Network graph or sparse matrix.
    - args (GeneralArgs): Arguments related to the general configuration of the experiment.

    Returns:
    - tuple: A tuple containing the inverse of the similarity matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = {gene: index for index, gene in enumerate(genes)}

    if not sp.issparse(network):
        matrix = nx.to_scipy_sparse_array(network, nodelist=genes, weight=2)
    else:
        matrix = network

    # use timer to measure time
    start = time.time()
    print("Weight normalization")
    # Normalize the matrix
    norm_matrix = sp.diags(1 / np.sqrt(matrix.sum(0).ravel()), format="csr")
    matrix = norm_matrix * matrix * norm_matrix

    print("Calculating the inverse")
    # First, let's get the shape of the matrix W
    n = matrix.shape[0]

    # Create an identity matrix of the same shape as W
    Identity = sp.eye(n)

    # Calculate (I - (1-alpha)*W)
    matrix_to_invert = Identity - (1 - args.alpha) * matrix

    print("Inverting the matrix")
    # Use scipy's sparse linear solver to find the inverse
    inverse_matrix = inv(matrix_to_invert)


    # calculate alpha * (I - (1-alpha)*W)^-1
    inverse_matrix = args.alpha * inverse_matrix

    print("Converting to CSR format")
    # Convert to CSR format before saving
    matrix_inverse_csr = sp.csr_matrix(inverse_matrix)


    print("Saving the matrix")
    # Save the matrix in .npz format
    # check if path exists, if not create it
    if not os.path.exists(os.path.dirname(args.similarity_matrix_path)):
        os.makedirs(os.path.dirname(args.similarity_matrix_path))
    sp.save_npz(args.similarity_matrix_path, matrix_inverse_csr)

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

    if not os.path.exists(similarity_matrix_path):
        raise FileNotFoundError(f"The specified file {similarity_matrix_path} does not exist.")

    matrix = sp.load_npz(similarity_matrix_path)

    return matrix, gene_index


def handle_no_propagation_cases(prior_data, prop_task, general_args):
    # Sort prior_data by GeneID
    sorted_prior_data = prior_data.sort_values(by='GeneID').reset_index(drop=True)
    gene_scores = sorted_prior_data['Score'].values.reshape((len(sorted_prior_data), 1))
    sorted_prior_data['Score'] = gene_scores.flatten()
    save_propagation_score(prior_set=sorted_prior_data, propagation_input={gene_id: score for gene_id, score in
                                                                           zip(sorted_prior_data['GeneID'],
                                                                               sorted_prior_data['Score'])},
                           propagation_scores=sorted_prior_data,
                           genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(sorted_prior_data['GeneID'])},
                           task=prop_task, save_dir=prop_task.output_folder, general_args=general_args)


def get_similarity_matrix(network, general_args):
    if general_args.create_similarity_matrix:
        return generate_similarity_matrix(network, general_args)
    else:
        return read_sparse_matrix_txt(network, general_args.similarity_matrix_path)


def filter_network_genes(propagation_input_df, network):
    network_genes_df = propagation_input_df[propagation_input_df['GeneID'].isin(network.nodes())]
    filtered_propagation_input = {gene_id: score for gene_id, score in zip(network_genes_df['GeneID'], network_genes_df['Score'])}
    return network_genes_df, filtered_propagation_input


def _normalize_prop_scores(matrix, network_gene_index, propagation_score, filtered_prior_data):
    """
    Normalize the propagation scores.

    Parameters:
    - matrix (sp.sparse.spmatrix): The similarity matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices in the matrix.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    # Set input type to 'ones'
    ones_input_df = set_input_type(filtered_prior_data, 'ones')

    # Convert DataFrame to dictionary format
    ones_input = {gene_id: score for gene_id, score in zip(ones_input_df['GeneID'], ones_input_df['Score'])}

    # Perform propagation with ones input
    ones_gene_scores_inverse = matrix_prop(ones_input, matrix, network_gene_index)

    # Identify zero normalization and zero propagation genes
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))

    # Adjust the normalization scores
    ones_gene_scores_inverse[genes_to_delete] = 1
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] = propagation_score[non_zero_indices] / np.abs(
        ones_gene_scores_inverse[non_zero_indices])

    # Create DataFrame for normalized scores
    normalized_df = pd.DataFrame({
        'GeneID': filtered_prior_data['GeneID'],
        'Score': [propagation_score[network_gene_index[gene_id]] for gene_id in filtered_prior_data['GeneID']]
    })

    return normalized_df


def merge_with_prior_data(final_propagation_results, prior_data, gene_name_dict):
    full_propagated_scores_df = final_propagation_results.merge(prior_data[['GeneID', 'Symbol', 'P-value']], on='GeneID', how='left')
    full_propagated_scores_df['Symbol'] = full_propagated_scores_df['Symbol'].fillna(full_propagated_scores_df['GeneID'].map(gene_name_dict).fillna(""))
    full_propagated_scores_df['P-value'] = full_propagated_scores_df['P-value'].fillna(0)
    return full_propagated_scores_df



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
        propagation_scores=full_propagated_scores_df,
        prior_set=propagation_input_df,
        propagation_input={gene_id: score for gene_id, score in zip(full_propagated_scores_df['GeneID'], full_propagated_scores_df['Score'])},
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(full_propagated_scores_df['GeneID'])},
        task=prop_task,
        save_dir=prop_task.output_folder,
        general_args=general_args
    )

def perform_propagation(test_name: str, general_args: GeneralArgs):
    """
    Performs the propagation of gene scores through the network.

    Parameters:
    - test_name (str): Name of the test for which propagation is performed.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    prop_task = PropagationTask(general_args=general_args, test_name=test_name)

    # Load and prepare prior set
    prior_data = read_prior_set(prop_task.test_file_path)

    if general_args.alpha == 1:
        handle_no_propagation_cases(prior_data, prop_task, general_args)
        return

    # Read the network graph from a file
    network = read_network(general_args.network_file_path)

    # create or upload similarity matrix
    matrix, network_gene_index = get_similarity_matrix(network, general_args)


    # Propagate network
    print("propagating network")
    propagation_input_df = set_input_type(prior_data, general_args.input_type)

    # Filter genes in the network
    network_genes_df, filtered_propagation_input = filter_network_genes(propagation_input_df, network)

    # Perform network propagation
    propagation_score = matrix_prop(filtered_propagation_input, matrix, network_gene_index)

    # Normalize the propagation scores and create DataFrame within the function
    normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score, network_genes_df)

    # Handle genes not in the network
    non_network_genes = propagation_input_df[~propagation_input_df['GeneID'].isin(network.nodes())].copy()

    # Combine network and non-network genes
    final_propagation_results = pd.concat([
        normalized_df[['GeneID', 'Score']],
        non_network_genes[['GeneID', 'Score']]
    ], ignore_index=True)

    # Load the gene_info.json file
    with open(general_args.genes_names_file_path, 'r') as f:
        gene_name_dict = json.load(f)
    full_propagated_scores_df = merge_with_prior_data(final_propagation_results, prior_data, gene_name_dict)

    # Print the number of scores different from the original scores
    different_scores = propagation_input_df[propagation_input_df['Score'] != final_propagation_results['Score']]
    print(f"Number of scores different from the original scores: {len(different_scores)}")

    # Save the results
    _save_propagation_results(propagation_input_df, full_propagated_scores_df, prop_task, general_args)

    del matrix

# def propagate(seeds, propagation_input, matrix, gene_indexes, num_genes, inverted_scores, alpha=0.1, n_max_iterations=100000,
#               convergence_th=1e-6):
#     # Initialize the vector F_t with zeros, with a shape of (num_genes,)
#     F_t = np.zeros(num_genes)
#
#     # If no propagation_input is given, set it to 1 for all seed genes
#     if not propagation_input:
#         propagation_input = {x: 1 for x in seeds}
#
#     # Set the initial propagation values
#     for seed in seeds:
#         if seed in gene_indexes:
#             F_t[gene_indexes[seed]] = propagation_input[seed]
#
#     # Calculate Y = alpha * F_t (initial values scaled by alpha)
#     Y = alpha * F_t
#
#     # Initialize F_t with Y
#     F_t = Y.copy()
#
#     iterations_to_similarity = -1
#
#     # Iterate to propagate values through the network
#     for iteration in range(n_max_iterations):
#         # Save the previous F_t for convergence check
#         F_t_1 = F_t.copy()
#
#         # Update F_t using the inverse matrix
#         F_t = matrix.dot(F_t_1) + Y
#
#         # Check for similarity with inverted_scores
#         if np.linalg.norm(F_t - inverted_scores) < convergence_th:
#             print(f"Converged after {iteration} iterations")
#             break
#
#     return F_t
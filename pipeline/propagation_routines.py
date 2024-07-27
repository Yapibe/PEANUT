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
        for edge in network.edges():
            row.append(gene_index[edge[0]])
            col.append(gene_index[edge[1]])
            data.append(1.0)  # Unweighted graph, so weight is 1.0 for both directions
            # Add the reverse edge since the graph is unweighted and assumed undirected
            row.append(gene_index[edge[1]])
            col.append(gene_index[edge[0]])
            data.append(1.0)
        matrix = sp.csr_matrix((data, (row, col)), shape=(len(genes), len(genes)))
    else:
        matrix = network

    print("Weight normalization")
    degrees = np.array(matrix.sum(axis=1)).flatten()
    # Replace zero degrees with 1 to avoid division by zero
    degrees[degrees == 0] = 1
    inv_sqrt_degrees = 1.0 / np.sqrt(degrees)
    norm_matrix = sp.diags(inv_sqrt_degrees)
    matrix = norm_matrix @ matrix @ norm_matrix

    print("Calculating the inverse")
    n = matrix.shape[0]
    identity_matrix = sp.eye(n)

    matrix_to_invert = identity_matrix - (1 - args.alpha) * matrix
    print(f"Matrix to invert non-zero elements: {matrix_to_invert.nnz}")

    # Convert to CSC format for efficient inversion
    matrix_to_invert_csc = matrix_to_invert.tocsc()
    # save the matrix to disk
    save_pre_inverse_path = os.path.join(args.data_dir, 'matrix', 'pre_inverse_matrix.npz')
    sp.save_npz(save_pre_inverse_path, matrix_to_invert_csc)

    print("Inverting the matrix")
    inverse_matrix = inv(matrix_to_invert_csc)
    inverse_matrix = args.alpha * inverse_matrix
    print("Inverse matrix non-zero elements:", inverse_matrix.nnz)

    # Print densities for debugging
    original_density = matrix.nnz / (matrix.shape[0] * matrix.shape[1])
    inverse_density = inverse_matrix.nnz / (inverse_matrix.shape[0] * inverse_matrix.shape[1])
    print(f"Original matrix density: {original_density}")
    print(f"Inverse matrix density: {inverse_density}")

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

    if not os.path.exists(similarity_matrix_path):
        raise FileNotFoundError(f"The specified file {similarity_matrix_path} does not exist.")

    matrix = sp.load_npz(similarity_matrix_path)

    return matrix, gene_index

def matrix_prop(propagation_input: dict, inverse_matrix: sp.spmatrix, gene_indexes: dict, num_genes: int) -> np.ndarray:
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
    F_0 = np.zeros(num_genes)  # Changed to a 1D array
    for seed in seeds:
        F_0[gene_indexes[seed]] = propagation_input[seed]

    F = inverse_matrix.dot(F_0)
    return F


def _handle_ngsea_case(prior_data, prop_task, general_args, network):
    """
    Handle the NGSEA case for alpha = 1.

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
    ones_input = set_input_type(filtered_prior_data, 'ones')
    ones_gene_scores_inverse = matrix_prop(ones_input, matrix, network_gene_index, len(network_gene_index))

    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))

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



def _save_propagation_results(prior_data, propagation_results, prop_task, general_args):
    """
    Save the results of the propagation process.

    Parameters:
    - prior_data (pd.DataFrame): The original prior data.
    - propagation_results (pd.DataFrame): DataFrame containing propagation results.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    save_propagation_score(
        prior_set=prior_data,
        propagation_input={gene_id: score for gene_id, score in
                           zip(propagation_results['GeneID'], propagation_results['PropagatedScore'])},
        propagation_scores=propagation_results,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(propagation_results['GeneID'])},
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
    network_genes_df = propagation_input_df[propagation_input_df['GeneID'].isin(network.nodes())]
    filtered_propagation_input = {gene_id: score for gene_id, score in zip(network_genes_df['GeneID'], network_genes_df['Score'])}
    return network_genes_df, filtered_propagation_input


def get_similarity_matrix(network, general_args):
    if general_args.create_similarity_matrix:
        return generate_similarity_matrix(network, general_args)
    else:
        return read_sparse_matrix_txt(network, general_args.similarity_matrix_path)


def handle_no_propagation_cases(prior_data, prop_task, general_args, network):
    if general_args.run_NGSEA:
        print("Running NGSEA")
        _handle_ngsea_case(prior_data, prop_task, general_args, network)
    else:
        print("Running GSEA")
        _handle_no_propagation_case(prior_data, prop_task, general_args)


def perform_propagation(test_name: str, general_args, network, prior_data):
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

    print(f"Running propagation with scores '{general_args.input_type}'")

    matrix, network_gene_index = get_similarity_matrix(network, general_args)
    # Modify prior_data based on the input type
    propagation_input_df = set_input_type(prior_data, general_args.input_type)

    # Filter genes in the network
    network_genes_df, filtered_propagation_input = filter_network_genes(propagation_input_df, network)

    # Perform network propagation
    propagation_score = matrix_prop(filtered_propagation_input, matrix, network_gene_index, len(network_gene_index))

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
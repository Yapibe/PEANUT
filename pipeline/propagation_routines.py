import os
import time
import scipy as sp, scipy.linalg
import numpy as np
import networkx as nx
from args import GeneralArgs, PropagationTask
from utils import save_propagation_score, read_prior_set, get_propagation_input, read_network


def propagate_with_inverse(seeds: list, propagation_input: dict, inverse_matrix: sp.matrix, gene_indexes: dict,
                           num_genes: int) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.

    Parameters:
    - seeds (list): List of seed gene IDs.
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - inverse_matrix (sp.matrix): Precomputed inverse matrix for propagation.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - num_genes (int): Total number of genes in the network.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    F_0 = np.zeros((num_genes, 1))
    for seed in seeds:
        F_0[gene_indexes[seed]] = propagation_input[seed]

    F = inverse_matrix.dot(F_0)
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

    if not sp.sparse.issparse(network):
        matrix = nx.to_scipy_sparse_array(network, nodelist=genes, weight=2)
    else:
        matrix = network

    # use timer to measure time
    start = time.time()
    print("Weight normalization")
    # Normalize the matrix
    norm_matrix = sp.sparse.diags(1 / sp.sqrt(matrix.sum(0).ravel()), format="csr")
    matrix = norm_matrix * matrix * norm_matrix

    print("Calculating the inverse")
    # First, let's get the shape of the matrix W
    n = matrix.shape[0]

    # Create an identity matrix of the same shape as W
    Identity = sp.sparse.eye(n)

    # Calculate (I - (1-alpha)*W)
    matrix_to_invert = Identity - (1 - args.alpha) * matrix

    print("Inverting the matrix")
    # Use scipy's sparse linear solver to find the inverse
    inverse_matrix = sp.sparse.linalg.inv(matrix_to_invert)


    # calculate alpha * (I - (1-alpha)*W)^-1
    inverse_matrix = args.alpha * inverse_matrix

    print("Converting to CSR format")
    # Convert to CSR format before saving
    matrix_inverse_csr = sp.sparse.csr_matrix(inverse_matrix)


    print("Saving the matrix")
    # Save the matrix in .npz format
    # check if path exists, if not create it
    if not os.path.exists(os.path.dirname(args.similarity_matrix_path)):
        os.makedirs(os.path.dirname(args.similarity_matrix_path))
    sp.sparse.save_npz(args.similarity_matrix_path, matrix_inverse_csr)

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

    matrix = sp.sparse.load_npz(similarity_matrix_path)

    return matrix, gene_index


def propagate_network(propagation_input: dict, matrix: sp.sparse.spmatrix, gene_index: dict) -> tuple:
    """
    Performs network propagation using a given input and similarity matrix.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - matrix (sp.sparse.spmatrix): Propagation matrix.
    - gene_index (dict): Mapping of gene IDs to their indices in the matrix.

    Returns:
    - tuple: A tuple containing a dictionary of gene indexes, the array of inverted gene scores,
             and a dictionary of gene indexes to scores.
    """
    inverted_gene_scores = propagate_with_inverse(
        list(propagation_input.keys()), propagation_input, matrix, gene_index, len(gene_index)
    )


    gene_indexes_scores = {
        gene_index[gene]: inverted_gene_scores[gene_index[gene]]
        for gene in propagation_input.keys() if gene in gene_index
    }

    return inverted_gene_scores, gene_indexes_scores


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
        # Sort prior_data by GeneID
        sorted_prior_data = prior_data.sort_values(by='GeneID').reset_index(drop=True)

        # Create experiment_gene_index based on sorted GeneID
        genes = sorted_prior_data['GeneID']
        experiment_gene_index = {gene_id: idx for idx, gene_id in enumerate(genes)}

        # Create propagation_input as a dictionary
        propagation_input = {gene_id: score for gene_id, score in
                             zip(sorted_prior_data['GeneID'], sorted_prior_data['Score'])}

        # Create gene_scores as a 2 dimensional ndarray of scores
        gene_scores = sorted_prior_data['Score'].values.reshape((len(sorted_prior_data), 1))
        posterior_set = sorted_prior_data.copy()
        posterior_set['Score'] = gene_scores.flatten()

        save_propagation_score(propagation_scores=posterior_set, prior_set=sorted_prior_data,
                               propagation_input=propagation_input, genes_id_to_idx=experiment_gene_index,
                               task=prop_task, save_dir=prop_task.output_folder, general_args=general_args)

        return

    # Read the network graph from a file
    network = read_network(general_args.network_file_path)
    all_genes_ids = set(network.nodes())
    # Filter prior_data to include only genes in the network
    filtered_prior_data = prior_data[prior_data['GeneID'].isin(all_genes_ids)]
    filtered_prior_gene_ids = set(filtered_prior_data['GeneID'])

    # create or upload similarity matrix
    if general_args.create_similarity_matrix:
        print("generating similarity matrix")
        matrix, network_gene_index = generate_similarity_matrix(network, general_args)
    else:
        matrix, network_gene_index = read_sparse_matrix_txt(network, general_args.similarity_matrix_path)

    # Propagate network
    print("propagating network")
    propagation_input = get_propagation_input(filtered_prior_gene_ids, filtered_prior_data, general_args.input_type)
    propagation_score, gene_score_dict = propagate_network(propagation_input, matrix, network_gene_index)

    ones_input = get_propagation_input(filtered_prior_gene_ids, filtered_prior_data, 'ones')
    ones_gene_scores_inverse, ones_gene_score_dict = propagate_network(ones_input, matrix, network_gene_index)

    del matrix

    # Identify genes with zero normalization score but non-zero propagation score
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))

    # Set the normalization score of these genes to 1
    ones_gene_scores_inverse[genes_to_delete] = 1

    # Perform the normalization
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] = propagation_score[non_zero_indices] / np.abs(
        ones_gene_scores_inverse[non_zero_indices])

    # Filter the propagation score to include only genes in the prior data
    filtered_propagation_scores = {gene_id: propagation_score[network_gene_index[gene_id]] for gene_id in filtered_prior_gene_ids
                                   if gene_id in network_gene_index}

    # Create posterior set DataFrame
    posterior_set = filtered_prior_data.copy()
    posterior_set['Score'] = posterior_set['GeneID'].map(filtered_propagation_scores)

    # Save propagation score
    save_propagation_score(propagation_scores=posterior_set, prior_set=prior_data, propagation_input=propagation_input,
                           genes_id_to_idx=network_gene_index, task=prop_task, save_dir=prop_task.output_folder,
                           general_args=general_args)



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
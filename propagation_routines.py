import os
import time
import scipy as sp
import numpy as np
import networkx as nx
from args import PropagationTask
from utils import save_propagation_score, read_prior_set, get_propagation_input, filter_network_by_prior_data


def propagate_with_inverse(seeds, propagation_input, inverse_matrix, gene_indexes, num_genes):
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.
    Args:
        seeds (list): List of seed gene IDs.
        propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
        inverse_matrix (numpy.ndarray): Precomputed inverse matrix for propagation.
        gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
        num_genes (int): Total number of genes in the network.
    Returns:
        numpy.ndarray: Array containing the final propagated values for each gene.
    """
    F_0 = np.zeros((num_genes, 1))
    for seed in seeds:
        F_0[gene_indexes[seed]] = propagation_input[seed]
    # change F_0 to a 2D array
    F_0 = F_0.reshape((F_0.shape[0], 1))

    F = inverse_matrix.dot(F_0)
    return F


def generate_similarity_matrix(network, args):
    """
    Generates and saves a similarity matrix for network propagation, based on the provided network graph.
    Args:
        network (networkx.Graph or scipy.sparse matrix): Network graph or sparse matrix.
    Returns:
        tuple: A tuple containing the inverse of the similarity matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = {gene: index for index, gene in enumerate(genes)}

    if not sp.sparse.issparse(network):
        matrix = nx.to_scipy_sparse_array(network, nodelist=genes, weight=2)
    else:
        matrix = network

    # use timer to measure time
    start = time.time()
    print("Normalizing the matrix")
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
    matrix_inverse = np.linalg.inv(matrix_to_invert.toarray())

    # calculate alpha * (I - (1-alpha)*W)^-1
    matrix_inverse = args.alpha * matrix_inverse

    print("Converting to CSR format")
    # Convert to CSR format before saving
    matrix_inverse_csr = sp.sparse.csr_matrix(matrix_inverse)

    print("Saving the matrix")
    # Save the matrix in .npz format
    # check if path exists, if not create it
    if not os.path.exists(os.path.dirname(args.similarity_matrix_path)):
        os.makedirs(os.path.dirname(args.similarity_matrix_path))
    sp.sparse.save_npz(args.similarity_matrix_path, matrix_inverse_csr)

    end = time.time()
    print(f"Time elapsed: {end - start} seconds")
    return matrix_inverse, gene_index


def read_sparse_matrix_txt(network, similarity_matrix_path):
    """
    Reads a precomputed sparse similarity matrix from a file.
    Args:
        network (networkx.Graph): Network graph used to generate the similarity matrix.
        similarity_matrix_path (str): Path to the file containing the sparse matrix.
    Returns:
        tuple: A tuple containing the sparse matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = dict([(gene, index) for (index, gene) in enumerate(genes)])
    if not os.path.exists(similarity_matrix_path):
        raise FileNotFoundError(f"The specified file {similarity_matrix_path} does not exist.")

    # Load the matrix in .npz format
    start = time.time()
    matrix = sp.sparse.load_npz(similarity_matrix_path)
    end = time.time()
    print(f"Time elapsed: {end - start} seconds")
    return matrix, gene_index


def propagate_network(propagation_input, matrix, gene_index):
    """
    Performs network propagation using a given input and similarity matrix.
    Args:
        propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
        matrix (numpy.ndarray or scipy.sparse matrix): Propagation matrix.
        gene_index (dict): Mapping of gene IDs to their indices in the matrix.
    Returns:
        tuple: A tuple containing a dictionary of gene indexes, the array of inverted gene scores,
               and a dictionary of gene indexes to scores.
    """
    inverted_gene_scores = propagate_with_inverse(list(propagation_input.keys()), propagation_input, matrix,
                                                  gene_index, len(gene_index))
    # return dictionary of gene indexes and inverted gene scores
    gene_indexes_scores = dict([(gene_index[gene], inverted_gene_scores[gene_index[gene]])
                                for gene in propagation_input.keys() if gene in gene_index])
    return inverted_gene_scores, gene_indexes_scores


def perform_propagation(test_name, general_args):
    prop_task = PropagationTask(general_args=general_args, test_name=test_name)

    # Load and prepare prior set
    prior_data = read_prior_set(prop_task.test_file_path)
    print("loaded prior data")

    if general_args.alpha == 1:
        print("Skipping propagation, saving sorted scores directly")

        # Sort prior_data by GeneID
        sorted_prior_data = prior_data.sort_values(by='GeneID').reset_index(drop=True)

        # Create experiment_gene_index based on sorted GeneID
        genes = sorted_prior_data['GeneID']
        experiment_gene_index = {gene_id: idx for idx, gene_id in enumerate(genes)}

        # Create propagation_input as a dictionary
        propagation_input = {gene_id: score for gene_id, score in
                             zip(sorted_prior_data['GeneID'], sorted_prior_data['Score'])}

        # Create gene_scores as a ndarray of scores
        gene_scores = sorted_prior_data['Score'].values
        # make gene_scores a 2 dimensional array
        gene_scores = gene_scores.reshape((len(gene_scores), 1))

        posterior_set = sorted_prior_data.copy()
        posterior_set['Score'] = gene_scores.flatten()

        save_propagation_score(propagation_scores=posterior_set, prior_set=sorted_prior_data,
                               propagation_input=propagation_input, genes_id_to_idx=experiment_gene_index,
                               task=prop_task, save_dir=prop_task.output_folder, general_args=general_args)

        return

    # reads the network graph from a file
    filtered_network = filter_network_by_prior_data(general_args.network_file_path, prior_data)
    # Intersection with network nodes
    all_genes_ids = filtered_network.nodes()

    # Filter prior_data to include only genes in the filtered network
    filtered_prior_data = prior_data[prior_data['GeneID'].isin(all_genes_ids)]

    # create or upload similarity matrix
    if general_args.create_similarity_matrix:
        print("generating similarity matrix")
        matrix, network_gene_index = generate_similarity_matrix(filtered_network, general_args)
    else:
        print("reading similarity matrix")
        matrix, network_gene_index = read_sparse_matrix_txt(filtered_network, general_args.similarity_matrix_path)
        print("uploaded similarity matrix")

    # Propagate network
    print("propagating network")
    propagation_input = get_propagation_input(all_genes_ids, prior_data)
    propagation_score, gene_score_dict = propagate_network(propagation_input, matrix, network_gene_index)

    # print("getting ones input")
    ones_input = get_propagation_input(all_genes_ids, prior_data, 'ones')
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
    propagation_score[non_zero_indices] = propagation_score[non_zero_indices] / np.abs(ones_gene_scores_inverse[non_zero_indices])

    # Create posterior set DataFrame
    posterior_set = filtered_prior_data.copy()
    posterior_set['Score'] = propagation_score

    # save propagation score
    print("saving propagation score")
    save_propagation_score(propagation_scores=posterior_set, prior_set=prior_data, propagation_input=propagation_input,
                           genes_id_to_idx=network_gene_index, task=prop_task, save_dir=prop_task.output_folder, general_args=general_args)


# def propagate(seeds, propagation_input, matrix, gene_indexes, num_genes, task: PropagationTask):
#     """
#     Propagates the influence of seed genes through the network using the specified propagation matrix.
#     Args:
#         seeds (list): List of seed gene IDs.
#         propagation_input (dict): Mapping of gene IDs to their initial propagation values.
#         matrix (numpy.ndarray or scipy.sparse matrix): Propagation matrix.
#         gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
#         num_genes (int): Total number of genes in the network.
#         task (PropagationTask): Propagation task object containing propagation parameters.
#     Returns:
#         numpy.ndarray: Array containing the final propagated values for each gene.
#     """
#     F_t = np.zeros(num_genes)
#     if not propagation_input:
#         propagation_input = {x: 1 for x in seeds}
#     for seed in seeds:
#         if seed in gene_indexes:
#             F_t[gene_indexes[seed]] = propagation_input[seed]
#     Y = task.alpha * F_t
#     matrix = (1 - task.alpha) * matrix
#     for _ in range(task.n_max_iterations):
#         F_t_1 = F_t
#         F_t = matrix.dot(F_t_1) + Y
#
#         if scipy.linalg.norm(F_t_1 - F_t) < task.convergence_th:
#             break
#     return F_t
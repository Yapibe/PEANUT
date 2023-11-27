import scipy as sp
import scipy.sparse
import scipy.linalg
import numpy as np
import pandas as pd
from tqdm import tqdm
from args import PropagationTask
import networkx as nx
import os
import time


def propagate(seeds, propagation_input, matrix, gene_indexes, num_genes, task: PropagationTask):
    """
    Propagates the influence of seed genes through the network using the specified propagation matrix.
    Args:
        seeds (list): List of seed gene IDs.
        propagation_input (dict): Mapping of gene IDs to their initial propagation values.
        matrix (numpy.ndarray or scipy.sparse matrix): Propagation matrix.
        gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
        num_genes (int): Total number of genes in the network.
        task (PropagationTask): Propagation task object containing propagation parameters.
    Returns:
        numpy.ndarray: Array containing the final propagated values for each gene.
    """
    F_t = np.zeros(num_genes)
    if not propagation_input:
        propagation_input = {x: 1 for x in seeds}
    for seed in seeds:
        if seed in gene_indexes:
            F_t[gene_indexes[seed]] = propagation_input[seed]
    Y = task.alpha * F_t
    matrix = (1 - task.alpha) * matrix
    for _ in range(task.n_max_iterations):
        F_t_1 = F_t
        F_t = matrix.dot(F_t_1) + Y

        if scipy.linalg.norm(F_t_1 - F_t) < task.convergence_th:
            break
    return F_t


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


def generate_similarity_matrix(network, similarity_matrix_path, alpha):
    """
    Generates and saves a similarity matrix for network propagation, based on the provided network graph.
    Args:
        network (networkx.Graph or scipy.sparse matrix): Network graph or sparse matrix.
        similarity_matrix_path (str): Path to save the computed similarity matrix.
        alpha (float): Propagation parameter, controls the influence of the network structure.
    Returns:
        tuple: A tuple containing the inverse of the similarity matrix and the list of genes.
    """
    genes = sorted(network.nodes())
    gene_index = dict([(gene, index) for (index, gene) in enumerate(genes)])

    if not sp.sparse.issparse(network):
        matrix = nx.to_scipy_sparse_matrix(network, genes, weight=2)
    else:
        matrix = network

    # use timer to measure time
    start = time.time()
    print("Normalizing the matrix")
    # Normalize the matrix
    norm_matrix = sp.sparse.diags(1 / sp.sqrt(matrix.sum(0).A1), format="csr")
    matrix = norm_matrix * matrix * norm_matrix

    print("Calculating the inverse")
    # First, let's get the shape of the matrix W
    n = matrix.shape[0]

    # Create an identity matrix of the same shape as W
    Identity = sp.sparse.eye(n)

    # Calculate (I - (1-alpha)*W)
    matrix_to_invert = Identity - (1 - alpha) * matrix

    print("Inverting the matrix")
    # Use scipy's sparse linear solver to find the inverse
    matrix_inverse = np.linalg.inv(matrix_to_invert.toarray())

    # calculate alpha * (I - (1-alpha)*W)^-1
    matrix_inverse = alpha * matrix_inverse

    print("Converting to CSR format")
    # Convert to CSR format before saving
    matrix_inverse_csr = sp.sparse.csr_matrix(matrix_inverse)

    print("Saving the matrix")
    # Save the matrix in .npz format
    sp.sparse.save_npz(similarity_matrix_path, matrix_inverse_csr)

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


def column_propagation(column_data, matrix, gene_index):
    """
    Optimized function for column-wise network propagation.

    Args:
        column_data (pd.DataFrame): DataFrame containing 'GeneID' and the score column.
        matrix (numpy.ndarray or scipy.sparse matrix): Propagation matrix.
        gene_index (dict): Mapping of gene IDs to their indices in the matrix.

    Returns:
        np.ndarray: Array of propagated scores.
    """
    # Assuming 'GeneID' is the first column and the score is the second
    score_column = column_data.columns[1]

    # Directly work with Pandas series
    gene_id_to_score = column_data.set_index('GeneID')[score_column].to_dict()

    # Propagate network scores
    score_gene_scores_inverse, gene_score_dict = propagate_network(gene_id_to_score, matrix, gene_index)

    return score_gene_scores_inverse.flatten()


def generate_propagation_scores(task1, updated_prior_data, network_graph, num_shuffles):
    if task1.create_similarity_matrix:
        print("generating similarity matrix")
        matrix, network_gene_index = generate_similarity_matrix(network_graph, task1.similarity_matrix_path, task1.alpha)
    else:
        print("reading similarity matrix")
        matrix, network_gene_index = read_sparse_matrix_txt(network_graph, task1.similarity_matrix_path)
        print("uploaded similarity matrix")
    column_scores = {}

    # Iterate over shuffled columns with progress bar
    print("propagating")
    for column in tqdm(updated_prior_data.columns, desc="Propagating columns"):
        if column != 'GeneID':
            # Create sub dataframe of GeneID and current column
            column_data = updated_prior_data[['GeneID', column]]
            # Perform propagation
            column_scores[column] = column_propagation(column_data, matrix, network_gene_index)
    del matrix

    # Convert column scores to dataframe and process
    column_scores_df = pd.DataFrame.from_dict(column_scores)
    column_scores_df['GeneID'] = updated_prior_data['GeneID']
    column_scores_df = column_scores_df[['GeneID'] + [col for col in column_scores_df.columns if col != 'GeneID']]

    # save to csv
    column_scores_df.to_csv(f'Outputs\\propagation_matrix\\alpha:{task1.alpha}_{num_shuffles}.csv', index=False)
    return column_scores_df


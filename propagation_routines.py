import scipy as sp
import scipy.sparse
import scipy.linalg
import numpy as np
from args import PropagationTask
import networkx as nx
import os
import time


def propagate(seeds, propagation_input, matrix, gene_indexes, num_genes, task: PropagationTask):
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
    F_0 = np.zeros(num_genes)
    for seed in seeds:
        F_0[gene_indexes[seed]] = propagation_input[seed]
    # change F_0 to a 2D array
    F_0 = F_0.reshape((F_0.shape[0], 1))

    F = inverse_matrix.dot(F_0)
    return F


def generate_similarity_matrix(network, similarity_matrix_path, alpha):
    genes = sorted(network.nodes())

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
    return matrix_inverse, genes


def read_sparse_matrix_txt(network, similarity_matrix_path):
    genes = sorted(network.nodes())
    if not os.path.exists(similarity_matrix_path):
        raise FileNotFoundError(f"The specified file {similarity_matrix_path} does not exist.")

    # rows, cols, data = [], [], []
    # with open(similarity_matrix_path, 'r') as f:
    #     for line in f:
    #         # Extract row, col, and data from the line
    #         i, j, v = line.split()
    #         i = int(i[1:-1])  # Remove parentheses and convert to int
    #         j = int(j[:-1])  # Remove trailing comma and convert to int
    #         v = float(v)  # Convert to float
    #         rows.append(i)
    #         cols.append(j)
    #         data.append(v)
    # coo_mat = sp.sparse.coo_matrix((data, (rows, cols)))
    # csr_mat = coo_mat.tocsr()

    # Load the matrix in .npz format
    start = time.time()
    matrix = sp.sparse.load_npz(similarity_matrix_path)
    end = time.time()
    print(f"Time elapsed: {end - start} seconds")
    return matrix, genes


def propagate_network(propagation_input, matrix, genes):
    num_genes = len(genes)
    gene_indexes = dict([(gene, index) for (index, gene) in enumerate(genes)])

    time_start = time.time()
    inverted_gene_scores = propagate_with_inverse([x for x in propagation_input.keys()], propagation_input, matrix,
                                                  gene_indexes, num_genes)

    time_end = time.time()
    print(f"Time to propagate: {time_end - time_start} seconds")
    # return dictionary of gene indexes and inverted gene scores
    gene_indexes_scores = dict([(gene_indexes[gene], inverted_gene_scores[gene_indexes[gene]]) for gene in propagation_input.keys() if gene in gene_indexes])
    return gene_indexes, inverted_gene_scores, gene_indexes_scores

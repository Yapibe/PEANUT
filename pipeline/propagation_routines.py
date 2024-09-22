import os
import time
import scipy.sparse as sp
from scipy.sparse.linalg import inv
import numpy as np
import networkx as nx
from .settings import Settings
from .utils import set_input_type
import pandas as pd
from .pathway_enrichment import perform_enrichment


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
    prior_data_dict = prior_data.set_index("GeneID")["Score"].to_dict()

    # Assign scores to all genes in prior_data
    for gene_id, score in prior_data_dict.items():
        gene_scores[gene_id] = abs(score)

    # Calculate the scores for genes in the prior data and in the network based on their neighbors
    for gene_id in prior_data_dict.keys():
        if gene_id in network.nodes():
            neighbors = list(network.neighbors(gene_id))

            if neighbors:
                neighbor_scores = [
                    abs(prior_data_dict.get(neighbor, 0)) for neighbor in neighbors
                ]
                ni = len(neighbor_scores)
                sum_neighbor_scores = sum(neighbor_scores)
                gene_scores[gene_id] += (1 / ni) * sum_neighbor_scores

    return gene_scores


def generate_similarity_matrix(network: nx.Graph, settings: Settings) -> tuple:
    """
    Generates and saves a similarity matrix for network propagation, based on the provided network graph.

    Parameters:
    - network (nx.Graph or sp.sparse matrix): Network graph or sparse matrix.
    - settings (Settings): Settings object containing the experiment settings.

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
            weight = edge[2].get(
                2, 1.0
            )  # Use the edge weight if provided, otherwise default to 1.0
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
    matrix_to_invert = identity_matrix - (1 - settings.alpha) * matrix
    matrix_to_invert_csc = matrix_to_invert.tocsc()

    print("Inverting the matrix")
    inverse_matrix = inv(matrix_to_invert_csc)
    print("Matrix inverted")
    inverse_matrix = settings.alpha * inverse_matrix

    print("Saving the matrix")
    if not os.path.exists(os.path.dirname(settings.similarity_matrix_path)):
        os.makedirs(os.path.dirname(settings.similarity_matrix_path))
    sp.save_npz(settings.similarity_matrix_path, inverse_matrix)

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
        raise FileNotFoundError(
            f"Specified file {similarity_matrix_path} doesn't  exists."
        )

    matrix = sp.load_npz(similarity_matrix_path)
    return matrix, gene_index


def matrix_prop(
    propagation_input: dict, gene_indexes: dict, inverse_matrix: sp.spmatrix = None
) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix for faster calculation.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values for propagation.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - inverse_matrix (sp.spmatrix): Precomputed inverse matrix for propagation.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    num_genes = len(gene_indexes)
    F_0 = np.zeros(num_genes)
    indices = [gene_indexes[gene_id] for gene_id in propagation_input.keys()]
    F_0[indices] = list(propagation_input.values())

    F = inverse_matrix @ F_0

    return F


def _handle_no_propagation_case(prior_data: pd.DataFrame, network: nx.Graph) -> pd.DataFrame:
    """
    Handle the case where no propagation is performed (alpha = 1).

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - network (networkx.Graph): The network graph.

    Returns:
    - pd.DataFrame: The filtered and sorted prior data.
    """
    # Filter prior_data to include only genes in both prior_data and network
    filtered_prior_data = prior_data[prior_data["GeneID"].isin(network.nodes())].copy()

    # Sort the filtered_prior_data by GeneID
    sorted_prior_data = filtered_prior_data.sort_values(by="GeneID").reset_index(
        drop=True
    )

    # Flatten the scores for saving
    gene_scores = sorted_prior_data["Score"].values.reshape((len(sorted_prior_data), 1))
    sorted_prior_data["Score"] = gene_scores.flatten()

    return sorted_prior_data


def _normalize_prop_scores(
    matrix: sp.spmatrix,
    network_gene_index: dict,
    propagation_score: np.ndarray,
    filtered_prior_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Normalize the propagation scores.

    Parameters:
    - matrix (sp.spmatrix): The similarity matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices in the matrix.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    # Set input type to 'ones' and convert to dictionary
    ones_input_df = set_input_type(filtered_prior_data, norm=True)
    ones_input = ones_input_df.set_index("GeneID")["Score"].to_dict()

    # Perform propagation with ones input
    ones_gene_scores_inverse = matrix_prop(
        ones_input, network_gene_index, inverse_matrix=matrix
    )

    # Identify genes with zero normalization scores
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(
        set(zero_normalization_genes).difference(zero_propagation_genes)
    )
    ones_gene_scores_inverse[genes_to_delete] = 1

    # Adjust the normalization scores
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] /= np.abs(
        ones_gene_scores_inverse[non_zero_indices]
    )

    # Create a DataFrame for all network genes
    all_network_genes = list(network_gene_index.keys())
    all_normalized_scores = propagation_score[
        np.array([network_gene_index[gene_id] for gene_id in all_network_genes])
    ]

    # Create DataFrame for normalized scores
    normalized_df = pd.DataFrame(
        {"GeneID": all_network_genes, "Score": all_normalized_scores}
    )

    return normalized_df


def filter_network_genes(
    propagation_input_df: pd.DataFrame, network: nx.Graph
) -> tuple:
    """
    Filter genes in the network.

    Parameters:
    - propagation_input_df (pd.DataFrame): DataFrame containing propagation input data.
    - network (nx.Graph): The network graph.

    Returns:
    - tuple: A tuple containing the filtered network genes DataFrame, non-network genes DataFrame,
             and filtered propagation input dictionary.
    """
    network_genes_df = propagation_input_df[
        propagation_input_df["GeneID"].isin(network.nodes())
    ].copy()
    non_network_genes = propagation_input_df[
        ~propagation_input_df["GeneID"].isin(network.nodes())
    ].copy()
    filtered_propagation_input = {
        gene_id: score
        for gene_id, score in zip(network_genes_df["GeneID"], network_genes_df["Score"])
    }
    return network_genes_df, non_network_genes, filtered_propagation_input


def perform_propagation(settings: Settings, network: nx.Graph = None, prior_data: pd.DataFrame = None) -> None:
    """
    Performs the propagation of gene scores through the network.

    Parameters:
    - settings (Settings): Settings object containing the experiment settings.
    - network (nx.Graph): The network graph.
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - None
    """

    if settings.alpha == 1:
        print("Running GSEA")
        full_propagated_scores_df = _handle_no_propagation_case(prior_data, network)
        perform_enrichment(settings, full_propagated_scores_df)
        return

    matrix, network_gene_index = read_sparse_matrix_txt(
        network, settings.similarity_matrix_path
    )

    # Modify prior_data based on the input type
    propagation_input_df = set_input_type(prior_data, norm=False)

    # Filter genes in the network
    network_genes_df, non_network_genes, filtered_propagation_input = (
        filter_network_genes(propagation_input_df, network)
    )

    # Perform network propagation
    propagation_score = matrix_prop(
        filtered_propagation_input, network_gene_index, inverse_matrix=matrix
    )

    # Normalize the propagation scores and create DataFrame within the function
    normalized_df = _normalize_prop_scores(
        matrix, network_gene_index, propagation_score, network_genes_df
    )

    # Merge the normalized scores directly with prior_data
    full_propagated_scores_df = normalized_df.join(
        prior_data.set_index('GeneID')[['Symbol', 'P-value']],
        on='GeneID',
        how='inner'
    ).reset_index()

    perform_enrichment(settings, full_propagated_scores_df)

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


def calculate_gene_scores(network: nx.Graph, prior_data: pd.DataFrame) -> dict:
    """
    Calculates the score for each gene based on its own absolute value 
    and the absolute mean score of its neighbors' original scores.

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

    # Assign scores to all genes in the network
    for gene_id in network.nodes():
        # Own absolute score (default to 0 if not in prior_data)
        own_score = abs(prior_data_dict.get(gene_id, 0))

        # Neighbors' absolute mean score (use prior_data values, default to 0 for missing genes)
        neighbors = list(network.neighbors(gene_id))
        if neighbors:
            neighbor_scores = [abs(prior_data_dict.get(neighbor, 0)) for neighbor in neighbors]
            mean_neighbor_score = sum(neighbor_scores) / len(neighbor_scores)
        else:
            mean_neighbor_score = 0

        # Final score: own score + mean neighbor score
        gene_scores[gene_id] = own_score + mean_neighbor_score

    return gene_scores



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
    gene_scores = calculate_gene_scores(network, prior_data)

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

    if args.normalization_type == 'symmetric':
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

def matrix_prop(propagation_input: dict, gene_indexes: dict, inverse_matrix: sp.spmatrix) -> np.ndarray:
    """
    Propagate scores across the network using a similarity matrix.
    
    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial scores.
    - gene_indexes (dict): Mapping of gene IDs to their matrix indices.
    - inverse_matrix (sp.spmatrix): Precomputed inverse similarity matrix.

    Returns:
    - np.ndarray: Propagated scores for all nodes in the network.
    """
    num_genes = len(gene_indexes)
    
    # Initialize scores for all genes to 0
    F_0 = np.zeros(num_genes)
    
    # Assign initial scores to genes with known values
    for gene_id, value in propagation_input.items():
        if gene_id in gene_indexes:
            F_0[gene_indexes[gene_id]] = value
    
    
    # Perform propagation using the similarity matrix
    F = inverse_matrix @ F_0

    return F



def _handle_ngsea_case(prior_data, prop_task, general_args, network):
    """
    Handle the NGSEA case for alpha = 1, incorporating all network genes.

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - prop_task (PropagationTask): The propagation task object.
    - general_args (GeneralArgs): General arguments and settings.
    - network (networkx.Graph): The network graph.

    Returns:
    - None
    """
    logger.info("Starting NGSEA propagation for alpha=1.")

    # Calculate gene scores for all genes in the network
    gene_scores = calculate_gene_scores(network, prior_data)

    # Create a DataFrame with all network genes and their scores
    all_network_genes = pd.DataFrame({
        'GeneID': list(network.nodes()),
        'Score': [gene_scores[gene] for gene in network.nodes()]
    })

    logger.info("Merging prior data with all network genes.")
    # Merge prior_data with all network genes
    full_propagated_scores_df = all_network_genes.merge(
        prior_data[['GeneID', 'Symbol', 'P-value']],
        on='GeneID',
        how='left'
    )

    # Fill in default values for genes not in the dataset
    logger.info("Assigning default values for genes not in the experiment dataset.")
    full_propagated_scores_df['Symbol'] = full_propagated_scores_df['Symbol'].fillna(full_propagated_scores_df['GeneID'].astype(str))
    full_propagated_scores_df['P-value'] = full_propagated_scores_df['P-value'].fillna(0)


    # Reorder columns for consistency
    full_propagated_scores_df = full_propagated_scores_df[['GeneID', 'Score', 'Symbol', 'P-value']]
    logger.info(f"Final NGSEA propagated DataFrame contains {len(full_propagated_scores_df)} genes.")

    # Save the results
    logger.info("Saving NGSEA propagation results.")
    save_propagation_score(
        prior_set=prior_data,
        propagation_input={gene_id: score for gene_id, score in zip(full_propagated_scores_df['GeneID'], full_propagated_scores_df['Score'])},
        propagation_scores=full_propagated_scores_df,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(full_propagated_scores_df['GeneID'])},
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
    # Filter prior_data to include only genes present in the network
    filtered_prior_data = prior_data[prior_data['GeneID'].isin(network.nodes())].copy()

    # Sort the filtered_prior_data by GeneID
    sorted_prior_data = filtered_prior_data.sort_values(by='GeneID').reset_index(drop=True)

    save_propagation_score(
        prior_set=sorted_prior_data,
        propagation_input={gene_id: score for gene_id, score in zip(sorted_prior_data['GeneID'], sorted_prior_data['Score'])},
        propagation_scores=sorted_prior_data,
        genes_id_to_idx={gene_id: idx for idx, gene_id in enumerate(sorted_prior_data['GeneID'])},
        task=prop_task,
        save_dir=prop_task.output_folder,
        general_args=general_args
    )



def _normalize_prop_scores(matrix, network_gene_index, propagation_score, prior_data) -> pd.DataFrame:
    """
    Normalize the propagation scores and include all network nodes in the output.

    Parameters:
    - matrix (sp.sparse.spmatrix or np.ndarray): The similarity matrix or upper triangular matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices in the matrix.
    - propagation_score (np.ndarray): Array of propagation scores.
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores for all genes in the network.
    """    
    # Step 1: Create a "ones input" for all genes in the network
    ones_input_df = pd.DataFrame({
        'GeneID': list(network_gene_index.keys()),
        'Score': 1  # Assign a score of 1 for all genes
    })

    # Perform propagation with ones input to calculate normalization vector
    ones_input = ones_input_df.set_index('GeneID')['Score'].to_dict()
    ones_gene_scores_inverse = matrix_prop(ones_input, network_gene_index, inverse_matrix=matrix)

    # Step 2: Handle zero normalization scores to avoid division by zero
    zero_norm_count = np.sum(ones_gene_scores_inverse == 0)
    ones_gene_scores_inverse[ones_gene_scores_inverse == 0] = 1

    # Step 3: Normalize the propagation scores
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] /= np.abs(ones_gene_scores_inverse[non_zero_indices])

    # Step 4: Include all network genes in the normalized scores DataFrame
    all_network_genes = list(network_gene_index.keys())
    all_normalized_scores = propagation_score[np.array([network_gene_index[gene_id] for gene_id in all_network_genes])]

    # Step 5: Merge with prior_data to retain Symbol and P-value, set defaults for missing genes
    normalized_df = pd.DataFrame({
        'GeneID': all_network_genes,
        'Score': all_normalized_scores
    })
    normalized_df = normalized_df.merge(
        prior_data[['GeneID', 'Symbol', 'P-value']], on='GeneID', how='left'
    )
    normalized_df['Symbol'].fillna(normalized_df['GeneID'], inplace=True)  # Default Symbol to GeneID
    normalized_df['P-value'].fillna(0, inplace=True)  # Default P-value to 0 for genes not in the dataset

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


def filter_network_genes(propagation_input_df, network):
    network_genes_df = propagation_input_df[propagation_input_df['GeneID'].isin(network.nodes())].copy()
    non_network_genes = propagation_input_df[~propagation_input_df['GeneID'].isin(network.nodes())].copy()
    filtered_propagation_input = {gene_id: score for gene_id, score in zip(network_genes_df['GeneID'], network_genes_df['Score'])}
    return network_genes_df, non_network_genes, filtered_propagation_input


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
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - None
    """
    prop_task = PropagationTask(general_args=general_args, test_name=test_name)

    # Modify prior_data based on the input type
    if general_args.input_type == 'abs_Score':
        propagation_input_df = set_input_type(prior_data, general_args.input_type)
    else:
        propagation_input_df = prior_data

    if general_args.alpha == 1:
        handle_no_propagation_cases(propagation_input_df, prop_task, general_args, network)
        return

    # Step 1: Load or generate similarity matrix
    matrix, network_gene_index = get_similarity_matrix(network, general_args)

    # Step 2: Filter genes in the network
    network_genes_df, non_network_genes, filtered_propagation_input = filter_network_genes(propagation_input_df, network)

    # Step 3: Perform network propagation
    propagation_score = matrix_prop(filtered_propagation_input, network_gene_index, inverse_matrix=matrix)

    # Step 4: Normalize the propagation scores and create DataFrame
    normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score, prior_data)

    # Step 5: Merge the normalized scores with prior_data and handle missing genes
    full_propagated_scores_df = normalized_df.merge(
        prior_data[['GeneID', 'Symbol', 'P-value']], 
        on='GeneID', 
        how='left'
    )

    # Resolve 'Symbol' column
    full_propagated_scores_df['Symbol'] = full_propagated_scores_df['Symbol_y'].fillna(
        full_propagated_scores_df['Symbol_x'].fillna(
            full_propagated_scores_df['GeneID'].astype(str)  # Default to GeneID
        )
    )

    # Resolve 'P-value' column
    full_propagated_scores_df['P-value'] = full_propagated_scores_df['P-value_y'].fillna(
        full_propagated_scores_df['P-value_x'].fillna(0)  # Default to 0
    )

    # Drop intermediate columns with suffixes
    full_propagated_scores_df.drop(columns=['Symbol_x', 'Symbol_y', 'P-value_x', 'P-value_y'], inplace=True)

    # Reorder columns for consistency
    full_propagated_scores_df = full_propagated_scores_df[['GeneID', 'Score', 'Symbol', 'P-value']]
    # Step 6: Save the results
    _save_propagation_results(propagation_input_df, full_propagated_scores_df, prop_task, general_args)


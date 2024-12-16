import logging
import time
from pathlib import Path
from scipy.sparse.linalg import cg
import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import inv

from .pathway_enrichment import perform_enrichment
from .settings import Settings

logger = logging.getLogger(__name__)

def generate_similarity_matrix(network: nx.Graph, settings: Settings) -> tuple:
    """
    Generate and save a similarity matrix for network propagation.

    Parameters:
    - network (nx.Graph): Network graph.
    - settings (Settings): Settings containing configuration.

    Returns:
    - tuple: A tuple containing the similarity matrix and the gene index mapping.
    """
    try:
        start = time.time()

        # Map genes to indices
        genes = sorted(network.nodes())
        gene_index = {gene: index for index, gene in enumerate(genes)}

        # Create adjacency matrix
        edges = list(network.edges(data=True))
        row = []
        col = []
        data = []

        for u, v, edge_data in edges:
            weight = edge_data.get(2, 1.0)
            row.extend([gene_index[u], gene_index[v]])
            col.extend([gene_index[v], gene_index[u]])
            data.extend([weight, weight])

        matrix = sp.csr_matrix((data, (row, col)), shape=(len(genes), len(genes)))

        # Symmetric normalization
        degrees = np.array(matrix.sum(axis=1)).flatten()
        degrees[degrees == 0] = 1  # Avoid division by zero
        inv_sqrt_degrees = 1.0 / np.sqrt(degrees)
        norm_matrix = sp.diags(inv_sqrt_degrees)
        matrix = norm_matrix @ matrix @ norm_matrix

        # Prepare for Conjugate Gradient (CG)
        n = matrix.shape[0]
        identity_matrix = sp.eye(n, format="csr")
        matrix_to_invert = sp.eye(n) - (1 - settings.alpha) * matrix

        # Compute the inverse using CG
        logger.info("Inverting the matrix using Conjugate Gradient method")
        inverse_matrix = []
        for i in range(n):
            # Extract column `i` as a dense vector
            e_i = identity_matrix[:, i].toarray().ravel()
            x, info = cg(matrix_to_invert, e_i)
            if info != 0:
                raise ValueError(f"CG did not converge for column {i}")
            inverse_matrix.append(x)

        # Convert to a dense matrix
        inverse_matrix = np.array(inverse_matrix).T
        inverse_matrix = settings.alpha * sp.csr_matrix(inverse_matrix)

        logger.info("Saving the matrix")
        settings.similarity_matrix_path.parent.mkdir(parents=True, exist_ok=True)
        sp.save_npz(settings.similarity_matrix_path, inverse_matrix)

        end = time.time()
        logger.info(f"Time elapsed: {end - start:.2f} seconds")
        return inverse_matrix, gene_index
    except Exception as e:
        logger.error(f"An error occurred in generate_similarity_matrix: {e}")
        raise




def read_sparse_matrix(network: nx.Graph, similarity_matrix_path: Path) -> tuple:
    """
    Read a precomputed sparse similarity matrix from a file.

    Parameters:
    - network (nx.Graph): Network graph used to generate the similarity matrix.
    - similarity_matrix_path (Path): Path to the file containing the sparse matrix.

    Returns:
    - tuple: A tuple containing the sparse matrix and the gene index mapping.
    """
    try:
        genes = sorted(network.nodes())
        gene_index = {gene: index for index, gene in enumerate(genes)}

        if not similarity_matrix_path.exists():
            raise FileNotFoundError(
                f"File {similarity_matrix_path} does not exist."
            )

        matrix = sp.load_npz(similarity_matrix_path)
        return matrix, gene_index
    except Exception as e:
        logger.error(f"An error occurred in read_sparse_matrix: {e}")
        return None, {}


def matrix_prop(propagation_input: dict, gene_indexes: dict, matrix: sp.spmatrix) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed matrix.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - matrix (sp.spmatrix): Precomputed matrix for propagation.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    try:
        num_genes = len(gene_indexes)
        F_0 = np.zeros(num_genes, dtype=np.float32)  # Initialize with zeros
        for gene_id, value in propagation_input.items():
            F_0[gene_indexes[gene_id]] = value  # Assign values directly

        # Propagate using the matrix
        F = matrix @ F_0
        return F
    except Exception as e:
        logger.error(f"An error occurred in matrix_prop: {e}")
        return np.array([])

def _handle_gsea_case(prior_data: pd.DataFrame, network: nx.Graph, restrict_to_network: bool) -> pd.DataFrame:
    """
    Handle the case where no propagation is performed (alpha = 1).

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - network (networkx.Graph): The network graph.
    - restrict_to_network (bool): Whether to restrict genes to those in the network.

    Returns:
    - pd.DataFrame: DataFrame containing the gene scores.
    """
    try:
        if restrict_to_network:
            # Filter prior_data to include only network genes
            filtered_prior_data = prior_data[prior_data["GeneID"].isin(map(str, network.nodes()))].copy()
            logger.info("Restricting to network genes.")
        else:
            # No filtering, include all genes
            filtered_prior_data = prior_data.copy()
            logger.info("Processing all genes without network restriction.")

        # Sort the filtered_prior_data by GeneID
        filtered_prior_data["GeneID"] = filtered_prior_data["GeneID"].astype(str)
        sorted_prior_data = filtered_prior_data.sort_values(by="GeneID").reset_index(drop=True)

        return sorted_prior_data
    except Exception as e:
        logger.error(f"An error occurred in _handle_gsea_case: {e}")
        return pd.DataFrame()


def _normalize_prop_scores(matrix: sp.spmatrix, network_gene_index: dict, propagation_score: np.ndarray, filtered_prior_data: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the propagation scores.

    Parameters:
    - matrix (sp.spmatrix): The similarity matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    try:
        # Set input type to 'ones' and convert to dictionary
        ones_input_df = pd.DataFrame({'GeneID': filtered_prior_data['GeneID'], 'Score': 1.0})
        ones_input = ones_input_df.set_index('GeneID')['Score'].to_dict()

        # Perform propagation with ones input
        ones_gene_scores = matrix_prop(
            ones_input, network_gene_index, matrix=matrix
        )

        # Adjust zero normalization scores
        ones_gene_scores[ones_gene_scores == 0] = 1

        # Normalize propagation scores
        propagation_score /= ones_gene_scores

        # Create DataFrame for normalized scores
        all_network_genes = list(network_gene_index.keys())
        normalized_df = pd.DataFrame(
            {"GeneID": all_network_genes, "Score": propagation_score}
        )

        # Convert GeneID to string for consistency
        normalized_df["GeneID"] = normalized_df["GeneID"].astype(str)

        return normalized_df
    except Exception as e:
        logger.error(
            f"An error occurred in _normalize_prop_scores: {e}"
        )
        return pd.DataFrame()


def filter_network_genes(propagation_input_df, network):
    network_genes_df = propagation_input_df[
        propagation_input_df["GeneID"].isin(network.nodes())
    ].copy()
    filtered_propagation_input = {
        gene_id: score
        for gene_id, score in zip(
            network_genes_df["GeneID"], network_genes_df["Score"]
        )
    }
    return network_genes_df, filtered_propagation_input


def get_similarity_matrix(network, settings):
    if settings.create_similarity_matrix:
        return generate_similarity_matrix(network, settings)
    else:
        return read_sparse_matrix(network, settings.similarity_matrix_path)
    

async def perform_propagation(settings: Settings, condition_settings, network: nx.Graph, condition_data: pd.DataFrame):
    """
    Perform the propagation of gene scores through the network for a specific condition.
    Can also run GSEA analysis in parallel if enabled.

    Parameters:
    - settings (Settings): General settings for the run.
    - condition_settings (ConditionSettings): Settings specific to the condition.
    - network (networkx.Graph): The network graph.
    - condition_data (pd.DataFrame): DataFrame containing condition gene scores.
    """
    try:
        if settings.run_gsea:
            logger.info("Running GSEA without propagation")
            full_propagated_scores_df = _handle_gsea_case(condition_data, network, settings.restrict_to_network)
            perform_enrichment(settings, condition_settings, full_propagated_scores_df)
            return

        matrix, network_gene_index = get_similarity_matrix(network, settings)

        # Take absolute value of scores
        propagation_input_df = condition_data.copy()
        propagation_input_df.columns = ["GeneID", "Score"]
        propagation_input_df["Score"] = propagation_input_df["Score"].abs()

        # Filter genes in the network
        network_genes_df, filtered_propagation_input = filter_network_genes(
            propagation_input_df, network
        )

        # Perform network propagation
        propagation_score = matrix_prop(filtered_propagation_input, network_gene_index, matrix=matrix)

        # Normalize the propagation scores
        normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score, network_genes_df)

        # Perform enrichment analysis for this condition
        perform_enrichment(settings, condition_settings, normalized_df)
    except Exception as e:
        logger.error(f"An error occurred in perform_propagation: {e}")
import logging
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import inv

from .pathway_enrichment import perform_enrichment
from .settings import Settings
from .utils import set_input_type

logger = logging.getLogger(__name__)


def generate_similarity_matrix(
    network: nx.Graph, settings: Settings
) -> tuple:
    """
    Generate and save a similarity matrix for network propagation.

    Parameters:
    - network (nx.Graph): Network graph.
    - settings (Settings): Settings containing configuration.

    Returns:
    - tuple: A tuple containing the inverse of the similarity matrix and the gene index mapping.
    """
    try:
        start = time.time()

        genes = sorted(network.nodes())
        gene_index = {gene: index for index, gene in enumerate(genes)}

        # Create adjacency matrix
        row, col, data = [], [], []
        for edge in network.edges(data=True):
            weight = edge[2].get(2, 1.0)  # Default weight is 1.0
            row.append(gene_index[edge[0]])
            col.append(gene_index[edge[1]])
            data.append(weight)
            # Add the reverse edge since the graph is undirected
            row.append(gene_index[edge[1]])
            col.append(gene_index[edge[0]])
            data.append(weight)
        matrix = sp.csr_matrix(
            (data, (row, col)), shape=(len(genes), len(genes))
        )

        degrees = np.array(matrix.sum(axis=1)).flatten()
        degrees[degrees == 0] = 1  # Avoid division by zero

        # Symmetric normalization
        inv_sqrt_degrees = 1.0 / np.sqrt(degrees)
        norm_matrix = sp.diags(inv_sqrt_degrees)
        matrix = norm_matrix @ matrix @ norm_matrix

        n = matrix.shape[0]
        identity_matrix = sp.eye(n)
        matrix_to_invert = (
            identity_matrix - (1 - settings.alpha) * matrix
        )
        matrix_to_invert_csc = matrix_to_invert.tocsc()

        logger.info("Inverting the matrix")
        inverse_matrix = inv(matrix_to_invert_csc)
        logger.info("Matrix inverted")
        inverse_matrix = settings.alpha * inverse_matrix

        logger.info("Saving the matrix")
        settings.similarity_matrix_path.parent.mkdir(
            parents=True, exist_ok=True
        )
        sp.save_npz(settings.similarity_matrix_path, inverse_matrix)

        end = time.time()
        logger.info(f"Time elapsed: {end - start:.2f} seconds")
        return inverse_matrix, gene_index
    except Exception as e:
        logger.error(
            f"An error occurred in generate_similarity_matrix: {e}"
        )
        return None, {}


def read_sparse_matrix(
    network: nx.Graph, similarity_matrix_path: Path
) -> tuple:
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


def matrix_prop(
    propagation_input: dict,
    gene_indexes: dict,
    matrix: sp.spmatrix,
) -> np.ndarray:
    """
    Propagates seed gene values through a precomputed inverse matrix.

    Parameters:
    - propagation_input (dict): Mapping of gene IDs to their initial values.
    - gene_indexes (dict): Mapping of gene IDs to their indices in the matrix.
    - inverse_matrix (sp.spmatrix): Precomputed inverse matrix for propagation.

    Returns:
    - np.ndarray: Array containing the final propagated values for each gene.
    """
    try:
        num_genes = len(gene_indexes)
        F_0 = np.zeros(num_genes, dtype=np.float32)
        for gene_id, value in propagation_input.items():
            F_0[gene_indexes[gene_id]] = value

        F = matrix @ F_0
        return F
    except Exception as e:
        logger.error(f"An error occurred in matrix_prop: {e}")
        return np.array([])


def _handle_gsea_case(
    prior_data: pd.DataFrame, network: nx.Graph, restrict_to_network: bool = False
) -> pd.DataFrame:
    """
    Handle the case where no propagation is performed (alpha = 1).

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - network (networkx.Graph): The network graph.

    Returns:
    - pd.DataFrame: DataFrame containing the gene scores.
    """
    try:
        if restrict_to_network:
            # Filter prior_data to include only network genes
            filtered_prior_data = prior_data[
                prior_data["GeneID"].isin(network.nodes())
            ].copy()
            logger.info("Restricting to network genes.")
        else:
            # No filtering, include all genes
            filtered_prior_data = prior_data.copy()
            logger.info("Processing all genes without network restriction.")


        # Sort the filtered_prior_data by GeneID
        sorted_prior_data = filtered_prior_data.sort_values(
            by="GeneID"
        ).reset_index(drop=True)

        # Flatten the scores for saving
        sorted_prior_data["Score"] = sorted_prior_data[
            "Score"
        ].astype(np.float32)

        return sorted_prior_data
    except Exception as e:
        logger.error(
            f"An error occurred in _handle_gsea_case: {e}"
        )
        return pd.DataFrame()


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
    - network_gene_index (dict): Mapping of gene IDs to their indices.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores.
    """
    try:
        # Set input type to 'ones' and convert to dictionary
        ones_input_df = set_input_type(filtered_prior_data, norm=True)
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


async def perform_propagation(
    settings: Settings,
    condition_settings,
    network: nx.Graph,
    condition_data: pd.DataFrame,
):
    """
    Perform the propagation of gene scores through the network for a specific condition.

    Parameters:
    - settings (Settings): General settings for the run.
    - condition_settings (ConditionSettings): Settings specific to the condition.
    - network (networkx.Graph): The network graph.
    - condition_data (pd.DataFrame): DataFrame containing condition gene scores.

    Returns:
    - None
    """
    try:
        if settings.alpha == 1:
            logger.info("Running GSEA without propagation")
            full_propagated_scores_df = _handle_gsea_case(
                condition_data, network
            )
            perform_enrichment(
                settings,
                condition_settings,
                full_propagated_scores_df,
            )
            return

        # Read the similarity matrix
        matrix, network_gene_index = read_sparse_matrix(
            network, settings.similarity_matrix_path
        )

        if matrix is None:
            # Generate the similarity matrix if it doesn't exist
            matrix, network_gene_index = generate_similarity_matrix(
                network, settings
            )

        # Modify condition_data based on the input type
        propagation_input_df = set_input_type(
            condition_data, norm=False
        )

        # Filter genes in the network
        network_genes_df, filtered_propagation_input = filter_network_genes(
            propagation_input_df, network
        )

        # Perform network propagation
        propagation_score = matrix_prop(
            filtered_propagation_input,
            network_gene_index,
            matrix=matrix,
        )

        # Normalize the propagation scores
        normalized_df = _normalize_prop_scores(
            matrix,
            network_gene_index,
            propagation_score,
            network_genes_df,
        )

        # Merge the normalized scores with condition_data
        full_propagated_scores_df = normalized_df.merge(
            condition_data[["GeneID", "Symbol", "P-value"]],
            on="GeneID",
            how="inner",
        )[["GeneID", "Score", "Symbol", "P-value"]]

        # Optimize data types
        full_propagated_scores_df["GeneID"] = (
            full_propagated_scores_df["GeneID"].astype(int)
        )
        full_propagated_scores_df["Score"] = (
            full_propagated_scores_df["Score"].astype(np.float32)
        )
        full_propagated_scores_df["P-value"] = (
            full_propagated_scores_df["P-value"].astype(np.float32)
        )
        full_propagated_scores_df["Symbol"] = (
            full_propagated_scores_df["Symbol"].astype("category")
        )
        condition_settings.scores_df = full_propagated_scores_df
        # Perform enrichment analysis for this condition
        perform_enrichment(
            settings, condition_settings, full_propagated_scores_df
        )
    except Exception as e:
        logger.error(f"An error occurred in perform_propagation: {e}")

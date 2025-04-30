"""
Network propagation module.

This module implements network propagation algorithms for spreading gene scores
through a biological interaction network. It includes functions for generating
similarity matrices, performing propagation, and normalizing results.
"""

import logging
import time
from pathlib import Path
from typing import Dict, Tuple
import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import cg
from .pathway_enrichment import perform_enrichment
from .settings import Settings, ConditionSettings

logger = logging.getLogger(__name__)

def generate_similarity_matrix(network: nx.Graph, settings: Settings) -> Tuple[sp.csr_matrix, Dict[str, int]]:
    """
    Generate and save a similarity matrix for network propagation.

    Parameters:
    - network (nx.Graph): Network graph.
    - settings (Settings): Settings containing configuration.

    Returns:
    - tuple: A tuple containing the similarity matrix and the gene index mapping.
    
    Raises:
    - ValueError: If the conjugate gradient method fails to converge.
    """
    try:
        start = time.time()
        logger.info("Starting similarity matrix generation...")

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

        # Convert to a sparse matrix
        inverse_matrix = np.array(inverse_matrix).T
        inverse_matrix = settings.alpha * sp.csr_matrix(inverse_matrix)

        logger.info(f"Saving similarity matrix to {settings.similarity_matrix_path}")
        settings.similarity_matrix_path.parent.mkdir(parents=True, exist_ok=True)
        sp.save_npz(settings.similarity_matrix_path, inverse_matrix)

        end = time.time()
        logger.info(f"Matrix generation completed in {end - start:.2f} seconds")
        return inverse_matrix, gene_index
    except Exception as e:
        logger.error(f"An error occurred in generate_similarity_matrix: {e}", exc_info=True)
        raise




def read_sparse_matrix(network: nx.Graph, similarity_matrix_path: Path) -> Tuple[sp.csr_matrix, Dict[str, int]]:
    """
    Read a precomputed sparse similarity matrix from a file.

    Parameters:
    - network (nx.Graph): Network graph used to generate the similarity matrix.
    - similarity_matrix_path (Path): Path to the file containing the sparse matrix.

    Returns:
    - tuple: A tuple containing the sparse matrix and the gene index mapping.
    
    Raises:
    - FileNotFoundError: If the similarity matrix file does not exist.
    """
    try:
        logger.info(f"Reading similarity matrix from {similarity_matrix_path}")
        
        genes = sorted(network.nodes())
        gene_index = {gene: index for index, gene in enumerate(genes)}

        if not similarity_matrix_path.exists():
            raise FileNotFoundError(
                f"Similarity matrix file {similarity_matrix_path} does not exist."
            )

        matrix = sp.load_npz(similarity_matrix_path)
        
        # Verify matrix dimensions match network size
        if matrix.shape[0] != len(genes) or matrix.shape[1] != len(genes):
            logger.warning(
                f"Matrix dimensions ({matrix.shape}) do not match network size ({len(genes)}). "
                "This may cause issues during propagation."
            )
            
        logger.info(f"Successfully loaded similarity matrix of shape {matrix.shape}")
        return matrix, gene_index
    except Exception as e:
        logger.error(f"An error occurred in read_sparse_matrix: {e}", exc_info=True)
        raise


def matrix_prop(propagation_input: Dict[str, float], gene_indexes: Dict[str, int], matrix: sp.spmatrix) -> np.ndarray:
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
        
        # Assign values to the corresponding indices
        for gene_id, value in propagation_input.items():
            if gene_id in gene_indexes:
                F_0[gene_indexes[gene_id]] = value
            else:
                # Skip genes not in the network
                continue

        # Propagate using the matrix
        F = matrix @ F_0
        return F
    except Exception as e:
        logger.error(f"An error occurred in matrix_prop: {e}", exc_info=True)
        raise

def _handle_gsea_case(prior_data: pd.DataFrame, network: nx.Graph, restrict_to_network: bool) -> pd.DataFrame:
    """
    Handle the case where no propagation is performed (alpha = 1).

    Parameters:
    - prior_data (pd.DataFrame): The prior data.
    - network (nx.Graph): The network graph.
    - restrict_to_network (bool): Whether to restrict genes to those in the network.

    Returns:
    - pd.DataFrame: DataFrame containing the gene scores.
    """
    try:
        if restrict_to_network:
            # Filter prior_data to include only network genes
            network_nodes = set(map(str, network.nodes()))
            filtered_prior_data = prior_data[prior_data["GeneID"].astype(str).isin(network_nodes)].copy()
            logger.info(f"Restricted to {len(filtered_prior_data)} network genes out of {len(prior_data)} total genes.")
        else:
            # No filtering, include all genes
            filtered_prior_data = prior_data.copy()
            logger.info(f"Processing all {len(prior_data)} genes without network restriction.")

        # Sort the filtered_prior_data by GeneID
        filtered_prior_data["GeneID"] = filtered_prior_data["GeneID"].astype(str)
        sorted_prior_data = filtered_prior_data.sort_values(by="GeneID").reset_index(drop=True)

        return sorted_prior_data
    except Exception as e:
        logger.error(f"An error occurred in _handle_gsea_case: {e}", exc_info=True)
        raise


def _normalize_prop_scores(
    matrix: sp.spmatrix, 
    network_gene_index: Dict[str, int], 
    propagation_score: np.ndarray, 
    filtered_prior_data: pd.DataFrame
) -> pd.DataFrame:
    """
    Normalize the propagation scores and return scores only for genes in both the network and the dataset.

    Parameters:
    - matrix (sp.spmatrix): The similarity matrix.
    - network_gene_index (dict): Mapping of gene IDs to their indices.
    - propagation_score (np.ndarray): Array of propagation scores.
    - filtered_prior_data (pd.DataFrame): DataFrame containing prior gene scores.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID and normalized Scores for relevant genes.
    """
    try:
        # Set input type to 'ones' and convert to dictionary
        ones_input_df = pd.DataFrame({'GeneID': filtered_prior_data['GeneID'], 'Score': 1.0})
        ones_input = ones_input_df.set_index('GeneID')['Score'].to_dict()

        # Perform propagation with ones input
        ones_gene_scores = matrix_prop(ones_input, network_gene_index, matrix=matrix)

        # Adjust zero normalization scores
        ones_gene_scores[ones_gene_scores == 0] = 1

        # Normalize propagation scores
        propagation_score /= ones_gene_scores

        # Create DataFrame for normalized scores
        all_network_genes = list(network_gene_index.keys())
        normalized_df = pd.DataFrame({"GeneID": all_network_genes, "Score": propagation_score})

        # Ensure consistent data type for merging
        normalized_df["GeneID"] = normalized_df["GeneID"].astype(str)
        filtered_prior_data["GeneID"] = filtered_prior_data["GeneID"].astype(str)

        # Return only genes that are in the original dataset
        relevant_scores_df = normalized_df.merge(
            filtered_prior_data[['GeneID']], on='GeneID', how='inner'
        )
        logger.info(f"Normalized propagation scores for {len(relevant_scores_df)} genes.")
        return relevant_scores_df

    except Exception as e:
        logger.error(f"An error occurred in _normalize_prop_scores: {e}", exc_info=True)
        raise




def filter_network_genes(propagation_input_df: pd.DataFrame, network: nx.Graph) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Filter input genes to only include those present in the network.
    
    Parameters:
    - propagation_input_df (pd.DataFrame): DataFrame with gene IDs and scores.
    - network (nx.Graph): Network graph.
    
    Returns:
    - Tuple[pd.DataFrame, Dict[str, float]]: Filtered DataFrame and dictionary mapping gene IDs to scores.
    """
    try:
        # Convert network nodes to strings for consistent comparison
        network_nodes = set(map(str, network.nodes()))
        
        # Filter the DataFrame to include only network genes
        network_genes_df = propagation_input_df[
            propagation_input_df["GeneID"].astype(str).isin(network_nodes)
        ].copy()
        
        # Create a dictionary mapping gene IDs to scores
        filtered_propagation_input = {
            str(gene_id): score
            for gene_id, score in zip(
                network_genes_df["GeneID"], network_genes_df["Score"]
            )
        }
        
        logger.info(f"Filtered from {len(propagation_input_df)} to {len(network_genes_df)} genes in network.")
        return network_genes_df, filtered_propagation_input
    except Exception as e:
        logger.error(f"An error occurred in filter_network_genes: {e}", exc_info=True)
        raise


def get_similarity_matrix(network: nx.Graph, settings: Settings) -> Tuple[sp.csr_matrix, Dict[str, int]]:
    """
    Get the similarity matrix, either by generating a new one or loading an existing one.
    
    Parameters:
    - network (nx.Graph): Network graph.
    - settings (Settings): Settings object.
    
    Returns:
    - Tuple[sp.csr_matrix, Dict[str, int]]: Similarity matrix and gene index mapping.
    """
    try:
        if settings.create_similarity_matrix:
            logger.info("Generating new similarity matrix...")
            return generate_similarity_matrix(network, settings)
        else:
            logger.info("Loading existing similarity matrix...")
            return read_sparse_matrix(network, settings.similarity_matrix_path)
    except Exception as e:
        logger.error(f"An error occurred in get_similarity_matrix: {e}", exc_info=True)
        raise
    

async def perform_propagation(
    settings: Settings, 
    condition_settings: ConditionSettings, 
    network: nx.Graph, 
    condition_data: pd.DataFrame
) -> None:
    """
    Perform the propagation of gene scores through the network for a specific condition.
    Can also run GSEA analysis in parallel if enabled.

    Parameters:
    - settings (Settings): General settings for the run.
    - condition_settings (ConditionSettings): Settings specific to the condition.
    - network (nx.Graph): The network graph.
    - condition_data (pd.DataFrame): DataFrame containing condition gene scores.
    """
    try:
        logger.info(f"Starting propagation for condition: {condition_settings.condition_name}")
        
        # For GSEA, skip propagation
        if settings.run_gsea:
            logger.info("Running GSEA without propagation")
            full_propagated_scores_df = _handle_gsea_case(condition_data, network, settings.restrict_to_network)
            perform_enrichment(settings, condition_settings, full_propagated_scores_df)
            return

        # Get the similarity matrix
        matrix, network_gene_index = get_similarity_matrix(network, settings)

        # Ensure DataFrame has the correct columns
        propagation_input_df = condition_data.copy()
        propagation_input_df.columns = ["GeneID", "Score"]
        propagation_input_df["Score"] = propagation_input_df["Score"].abs()

        # Filter genes in the network
        network_genes_df, filtered_propagation_input = filter_network_genes(
            propagation_input_df, network
        )
        
        if not filtered_propagation_input:
            logger.warning("No genes found in the network. Propagation cannot be performed.")
            return

        # Perform network propagation
        propagation_score = matrix_prop(filtered_propagation_input, network_gene_index, matrix=matrix)

        # Normalize the propagation scores
        normalized_df = _normalize_prop_scores(matrix, network_gene_index, propagation_score, network_genes_df)

        # Perform enrichment analysis for this condition
        perform_enrichment(settings, condition_settings, normalized_df)
        
        logger.info(f"Completed propagation for condition: {condition_settings.condition_name}")
        
    except Exception as e:
        logger.error(f"An error occurred in perform_propagation: {e}", exc_info=True)
        raise
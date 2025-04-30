"""
Utility functions for the PEANUT pipeline.

This module contains various utility functions for loading, processing, and
managing data files used throughout the pipeline.
"""

import logging
import pickle
import time
import zlib
from io import BytesIO
from pathlib import Path
from typing import Union, Dict, List, Tuple, Set, Any, Optional
import numpy as np
import networkx as nx
import pandas as pd
import yaml

from .settings import Settings, ConditionSettings

logger = logging.getLogger(__name__)

# ###################################################################
# LOAD FUNCTIONS
# ###################################################################

def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """Load configuration from YAML file.
    
    Args:
        config_path: Optional path to config file. If None, uses default location.
        
    Returns:
        Dictionary containing configuration settings
    """
    if config_path is None:
        config_path = Path(__file__).resolve().parent.parent / "config" / "pipeline_config.yaml"
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        logger.info(f"Loaded configuration from {config_path}")
        return config
    except Exception as e:
        logger.error(f"Error loading configuration from {config_path}: {str(e)}")
        raise

    
def load_pathways_genes(pathways_dir: Path) -> Dict[str, List[str]]:
    """
    Load the pathways and their associated genes from a file.

    Args:
        pathways_dir (Path): Path to the file containing the pathway data.

    Returns:
        Dict[str, List[str]]: A dictionary mapping pathway names to lists of genes in each pathway.
        
    Raises:
        FileNotFoundError: If the pathway file is not found.
    """
    pathways = {}
    try:
        with pathways_dir.open("r") as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                pathway_name = parts[0]
                try:
                    genes = [str(gene) for gene in parts[2:]]
                    pathways[pathway_name] = genes
                except ValueError:
                    logger.warning(
                        f"Non-integer gene ID found in pathway {pathway_name}, skipping pathway."
                    )
                    continue
        logger.info(f"Loaded {len(pathways)} pathways from {pathways_dir}")
        return pathways
    except FileNotFoundError:
        logger.error(f"File not found: {pathways_dir}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while loading pathways: {e}", exc_info=True)
        raise


def load_propagation_file(file_path: Path, decompress: bool = True) -> Dict[str, Any]:
    """
    Load the propagation score data from a file.

    Args:
        file_path (Path): The path to the file to be loaded.
        decompress (bool): Whether to decompress the file.

    Returns:
        Dict[str, Any]: The loaded data containing propagation scores and other information.
        
    Raises:
        FileNotFoundError: If the propagation file is not found.
    """
    try:
        with file_path.open("rb") as file:
            decompressed_data = pickle.load(file)
        if decompress:
            try:
                decompressed_data = pickle.loads(
                    zlib.decompress(decompressed_data)
                )
            except zlib.error:
                logger.warning(
                    "Entered an uncompressed file but asked to decompress it"
                )
        return decompressed_data
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except Exception as e:
        logger.error(
            f"An error occurred while loading propagation file: {e}", exc_info=True
        )
        raise


def load_propagation_scores(propagation_file_path: Path) -> Dict[str, Any]:
    """
    Load the propagation scores from a file.

    Args:
        propagation_file_path (Path): Path to the propagation scores file.

    Returns:
        Dict[str, Any]: The propagation scores.
    """
    return load_propagation_file(
        propagation_file_path, decompress=True
    )


def load_pathways_and_propagation_scores(settings: Settings, scores_df: pd.DataFrame) -> Tuple[Dict[str, Set[str]], Dict[str, float]]:
    """
    Load the pathways based on the provided configuration and the propagation scores.

    Args:
        settings (Settings): General settings for the pipeline.
        scores_df (pd.DataFrame): DataFrame containing the gene scores.

    Returns:
        Tuple[Dict[str, Set[str]], Dict[str, float]]: A tuple containing:
            - A dictionary mapping pathway names to sets of gene IDs
            - A dictionary mapping gene IDs to their scores
    """
    try:
        # Load pathways and filter based on gene count criteria
        pathways_with_many_genes = load_pathways_genes(settings.pathway_file_dir)
        
        # Validate scores dataframe
        if scores_df.empty:
            logger.warning("Empty scores dataframe provided")
            return {}, {}
        
        # Ensure scores_df has the expected columns
        if 'GeneID' not in scores_df.columns or 'Score' not in scores_df.columns:
            raise ValueError(f"scores_df must contain 'GeneID' and 'Score' columns. Got: {scores_df.columns}")
            
        # Sort and process scores
        sorted_scores = scores_df.sort_values(by="GeneID").reset_index(drop=True)

        # Create the scores dictionary with only GeneID and Score
        scores = {
            str(gene_id): float(score)
            for gene_id, score in zip(
                sorted_scores["GeneID"], sorted_scores["Score"]
            )
        }

        scores_keys = set(scores.keys())
        
        # Filter pathways based on gene counts and overlap with scored genes
        genes_by_pathway = {}
        for pathway, genes in pathways_with_many_genes.items():
            # Convert all genes to strings for consistent comparison
            string_genes = set(map(str, genes))
            # Find overlap with scored genes
            intersection = string_genes.intersection(scores_keys)
            
            # Check if pathway meets size criteria
            if (settings.minimum_gene_per_pathway <= len(intersection) <= 
                settings.maximum_gene_per_pathway):
                genes_by_pathway[pathway] = intersection

        logger.info(f"Filtered to {len(genes_by_pathway)} pathways with size between "
                   f"{settings.minimum_gene_per_pathway} and {settings.maximum_gene_per_pathway} genes")

        return genes_by_pathway, scores
    except Exception as e:
        logger.error(f"Error in load_pathways_and_propagation_scores: {e}", exc_info=True)
        raise


# ###################################################################
# GET FUNCTIONS
# ###################################################################
def get_scores(score_path: Path) -> Dict[int, Tuple[float, float]]:
    """
    Load gene scores and P-values from a file.

    Args:
        score_path (Path): The path to the file containing the propagation scores.

    Returns:
        Dict[int, Tuple[float, float]]: A dictionary with gene IDs as keys and tuples of (Score, P-value) as values.
    """
    try:
        propagation_results = load_propagation_scores(score_path)

        # Sort and optimize data types
        sorted_scores = (
            propagation_results["gene_prop_scores"]
            .sort_values(by="GeneID")
            .reset_index(drop=True)
        )
        sorted_scores["GeneID"] = sorted_scores["GeneID"].astype(int)
        sorted_scores["Score"] = pd.to_numeric(
            sorted_scores["Score"], downcast="float"
        )
        sorted_scores["P-value"] = pd.to_numeric(
            sorted_scores["P-value"], downcast="float"
        )

        # Create dictionary mapping gene IDs to (Score, P-value)
        gene_scores = {
            gene_id: (score, pvalue)
            for gene_id, score, pvalue in zip(
                sorted_scores["GeneID"],
                sorted_scores["Score"],
                sorted_scores["P-value"],
            )
        }
        return gene_scores
    except FileNotFoundError:
        logger.error(f"File not found: {score_path}", exc_info=True)
        raise
    except Exception as e:
        logger.error(f"An error occurred in get_scores: {e}", exc_info=True)
        raise


# ###################################################################
# READ FUNCTIONS
# ###################################################################

def read_network(network_filename: Path) -> nx.Graph:
    """
    Read a network from a file and return a NetworkX graph.

    Args:
        network_filename (Path): Path to the file containing the network data.

    Returns:
        nx.Graph: A graph object representing the network.
        
    Raises:
        FileNotFoundError: If the network file is not found.
    """
    try:
        logger.info(f"Reading network from: {network_filename}")
        if not network_filename.exists():
            logger.error(f"Network file not found: {network_filename}")
            raise FileNotFoundError(f"Network file not found: {network_filename}")

        # Read network file and create graph
        network = pd.read_table(
            network_filename, header=None, usecols=[0, 1, 2],
            names=['source', 'target', 'weight']
        )
        graph = nx.from_pandas_edgelist(network, 'source', 'target', 'weight')
        
        logger.info(f"Successfully loaded network with {len(graph.nodes())} nodes and {len(graph.edges())} edges")
        return graph
    except FileNotFoundError:
        logger.error(f"File not found: {network_filename}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while reading network: {e}", exc_info=True)
        raise


def read_prior_set(
    excel_input: Union[str, BytesIO, Path]) -> pd.DataFrame:
    """
    Read prior data set from an Excel file and apply preprocessing.

    Args:
        excel_input (str or BytesIO or Path): Path to the Excel file or a BytesIO object containing the prior data.

    Returns:
        pd.DataFrame: DataFrame containing the preprocessed prior data.
    """
    try:
        logger.info(f"Reading prior data from: {excel_input}")
        prior_data = pd.read_excel(excel_input, engine="openpyxl")

        # Data preprocessing
        initial_rows = len(prior_data)
        prior_data = prior_data.drop_duplicates(subset="GeneID")
        prior_data = prior_data[prior_data["Score"].notna()]
        prior_data = prior_data[prior_data["Score"] != "?"]
        prior_data = prior_data[
            prior_data["GeneID"].apply(lambda x: str(x).isdigit())
        ]
        prior_data["GeneID"] = prior_data["GeneID"].astype(int)
        prior_data = prior_data.reset_index(drop=True)
        
        logger.info(f"Processed prior data: {len(prior_data)} valid rows out of {initial_rows} total")
        return prior_data
    except FileNotFoundError:
        logger.error(f"File not found: {excel_input}")
        raise
    except Exception as e:
        logger.error(
            f"An error occurred while reading prior set: {e}", exc_info=True
        )
        raise


def read_enriched_pathways(condition_output_path: Path) -> Dict[str, float]:
    """
    Read the enriched pathways from the condition's output file.

    Args:
        condition_output_path (Path): Path to the condition's output Excel file.

    Returns:
        Dict[str, float]: Dictionary of pathways to FDR q-values.
    """
    try:
        logger.info(f"Reading enriched pathways from: {condition_output_path}")
        enriched_pathway_df = pd.read_excel(condition_output_path)
        
        # Check if the required columns exist
        if "Term" not in enriched_pathway_df.columns or "FDR q-val" not in enriched_pathway_df.columns:
            logger.error(f"Required columns 'Term' or 'FDR q-val' not found in {condition_output_path}")
            return {}
            
        enriched_pathway_dict = dict(
            zip(
                enriched_pathway_df["Term"],
                enriched_pathway_df["FDR q-val"],
            )
        )
        logger.info(f"Loaded {len(enriched_pathway_dict)} enriched pathways")
        return enriched_pathway_dict
    except FileNotFoundError:
        logger.error(f"File not found: {condition_output_path}")
        return {}
    except Exception as e:
        logger.error(
            f"An error occurred while reading enriched pathways: {e}", exc_info=True
        )
        return {}


# ###################################################################
# SAVE FUNCTIONS
# ###################################################################


def save_file(obj: Any, save_dir: Path, compress: bool = True) -> None:
    """
    Save an object to a file, with optional compression.

    Args:
        obj (Any): The object to be saved.
        save_dir (Path): The directory where the file will be saved.
        compress (bool): Whether to compress the file.

    Returns:
        None
    """
    try:
        # Serialize the object
        logger.info(f"Saving file to {save_dir}")
        obj_bytes = pickle.dumps(obj)
        if compress:
            obj_bytes = zlib.compress(obj_bytes)
            
        # Ensure the parent directory exists
        save_dir.parent.mkdir(parents=True, exist_ok=True)
            
        # Save to file
        with save_dir.open("wb") as f:
            pickle.dump(obj_bytes, f)
        logger.info(f"Successfully saved file to {save_dir}")
    except Exception as e:
        logger.error(f"An error occurred while saving file: {e}", exc_info=True)
        raise


def save_propagation_score(
    propagation_scores: pd.DataFrame,
    prior_set: pd.DataFrame,
    propagation_input: Dict[str, float],
    genes_id_to_idx: Dict[str, int],
    task: Any,
    save_dir: Path,
    general_args: Settings,
) -> None:
    """
    Save the propagation scores to a compressed file.

    Args:
        propagation_scores (pd.DataFrame): DataFrame containing the propagated scores.
        prior_set (pd.DataFrame): DataFrame containing the prior set of gene scores.
        propagation_input (Dict[str, float]): Mapping of gene IDs to their input scores.
        genes_id_to_idx (Dict[str, int]): Mapping of gene IDs to their indices.
        task: Propagation task containing task-specific settings.
        save_dir (Path): Directory where the output file will be saved.
        general_args (Settings): General arguments and settings.

    Returns:
        None
    """
    try:
        file_name = f"{general_args.alpha}_{general_args.date}"
        save_dir.mkdir(parents=True, exist_ok=True)
        propagation_results_path = save_dir / file_name

        # Create dictionary with all results to save
        save_dict = {
            "args": task,
            "prior_set": prior_set,
            "propagation_input": propagation_input,
            "gene_id_to_idx": genes_id_to_idx,
            "gene_prop_scores": propagation_scores,
        }

        # Save the dictionary to a compressed file
        save_file(save_dict, propagation_results_path)
        logger.info(f"Saved propagation scores for {len(propagation_scores)} genes to {propagation_results_path}")
    except Exception as e:
        logger.error(
            f"An error occurred while saving propagation score: {e}", exc_info=True
        )
        raise


# ###################################################################
# MISCELLANEOUS FUNCTIONS
# ###################################################################


def filter_network_by_prior_data(
    network_filename: Path, prior_data: pd.DataFrame
) -> nx.Graph:
    """
    Filter a network to only include nodes present in the prior_data DataFrame.

    Args:
        network_filename (Path): Path to the network file.
        prior_data (pd.DataFrame): DataFrame containing gene information.

    Returns:
        nx.Graph: A filtered graph object.
    """
    try:
        # Load the network
        network = read_network(network_filename)
        # Extract gene IDs from prior data and convert to set for faster lookup
        gene_ids_in_prior_data = set(prior_data["GeneID"].astype(str))
        # Create a subgraph containing only the nodes present in prior_data
        filtered_network = network.subgraph(
            gene_id for gene_id in network.nodes() 
            if str(gene_id) in gene_ids_in_prior_data
        ).copy()
        
        logger.info(f"Filtered network from {len(network.nodes())} to {len(filtered_network.nodes())} nodes")
        return filtered_network
    except Exception as e:
        logger.error(
            f"An error occurred while filtering network: {e}", exc_info=True
        )
        raise


def calculate_time(func):
    """
    Calculate the execution time of a function.

    Args:
        func: The function to be timed.

    Returns:
        function: The wrapped function with timing.
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logger.info(
            f"Execution time for {func.__name__}: {end_time - start_time:.2f} seconds"
        )
        return result

    return wrapper


def process_single_pathway(
    pathway: str,
    pathway_genes: Set[str],
    scores: Dict[str, float],
    pathway_pvalues: Dict[str, float]
) -> Dict[str, float]:
    """
    Process a single pathway for a given condition.

    Parameters:
    - pathway (str): Name of the pathway.
    - pathway_genes (Set[str]): Set of gene IDs in the pathway.
    - scores (Dict[str, float]): Mapping of gene IDs to scores.
    - pathway_pvalues (Dict[str, float]): Dictionary of pathway p-values.

    Returns:
    - Dict[str, float]: Dictionary containing mean score and p-value.
    """
    try:
        # Collect scores for genes in the pathway
        pathway_gene_scores = [
            scores[gene_id] for gene_id in pathway_genes 
            if gene_id in scores
        ]
        
        # Calculate mean score
        mean_score = np.mean(pathway_gene_scores) if pathway_gene_scores else 0
        
        # Get p-value for the pathway, or default to 1.0 if not found
        p_value = pathway_pvalues.get(pathway, 1.0)
        
        return {"Mean": mean_score, "P-value": p_value}
    except Exception as e:
        logger.error(f"Error processing pathway {pathway}: {e}", exc_info=True)
        return {"Mean": 0, "P-value": 1.0}


def process_condition(
    condition_settings: ConditionSettings,
    settings: Settings,
    condition_data_df: pd.DataFrame,
    all_pathways: Dict[str, Dict[str, Any]],
    pathway_sizes: Dict[str, int],
    genes_by_pathway: Dict[str, Set[str]],
    scores: Dict[str, float],
) -> None:
    """
    Processes experimental data for a condition and aggregates pathway data.

    Parameters:
    - condition_settings: Condition-specific settings.
    - settings: General pipeline settings.
    - condition_data_df: DataFrame of condition data.
    - all_pathways: Dictionary to store significant pathways by condition.
    - pathway_sizes: Dictionary to store pathway sizes.
    - genes_by_pathway: Dictionary to store gene sets for pathways.
    - scores: Dictionary to store gene scores.

    Returns:
    - None
    """
    try:
        logger.info(f"Processing condition: {condition_settings.condition_name}")
        
        # Get p-values for pathways from condition_settings
        pathway_pvalues = {
            pathway: data.get("FDR q-val", 1.0) 
            for pathway, data in condition_settings.pathways_statistics.items()
        }
        
        # Load pathways and update global dictionary with current condition's pathways
        pathways_with_many_genes = load_pathways_genes(settings.pathway_file_dir)

        # Ensure data types are consistent
        condition_data_df["GeneID"] = condition_data_df["GeneID"].astype(str)
        condition_data_df["Score"] = condition_data_df["Score"].astype(np.float32)

        # Update global scores dictionary with current condition's scores
        scores.update({
            gene_id: score 
            for gene_id, score in zip(condition_data_df["GeneID"], condition_data_df["Score"])
        })

        # Set of all gene IDs with scores
        scores_keys = set(scores.keys())
        
        # Create local mapping of pathways to their genes, filtered by score availability
        local_genes_by_pathway = {}
        for pathway, genes in pathways_with_many_genes.items():
            # Convert genes to strings and find intersection with scored genes
            genes_as_strings = set(map(str, genes))
            intersection = genes_as_strings.intersection(scores_keys)
            
            # Apply size filters
            if (settings.minimum_gene_per_pathway <= len(intersection) <= 
                settings.maximum_gene_per_pathway):
                local_genes_by_pathway[pathway] = intersection
        
        # Update global genes_by_pathway dictionary
        genes_by_pathway.update(local_genes_by_pathway)

        # Update pathway sizes
        for pathway, genes in local_genes_by_pathway.items():
            pathway_sizes[pathway] = len(genes)

        # Process each pathway and update all_pathways for significant ones
        pathways_processed = 0
        for pathway, pathway_genes in local_genes_by_pathway.items():
            # Initialize the pathway entry if it doesn't exist
            all_pathways.setdefault(pathway, {})
            
            # Process the pathway
            result = process_single_pathway(pathway, pathway_genes, scores, pathway_pvalues)
            
            # Only include significant pathways in the final dictionary
            if result["P-value"] < settings.FDR_THRESHOLD:
                all_pathways[pathway][condition_settings.condition_name] = result
                pathways_processed += 1
        
        logger.info(f"Processed {len(local_genes_by_pathway)} pathways, {pathways_processed} were significant")
        
    except Exception as e:
        logger.error(f"Error processing condition {condition_settings.condition_name}: {e}", exc_info=True)
        raise
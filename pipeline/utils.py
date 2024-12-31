import logging
import pickle
import time
import zlib
from io import BytesIO
from pathlib import Path
from typing import Union, Dict, List, Tuple, Set
import numpy as np
import networkx as nx
import pandas as pd

from .settings import Settings, ConditionSettings

logger = logging.getLogger(__name__)

# ###################################################################LOAD FUNCTIONS############################################################################


def load_pathways_genes(pathways_dir: Path) -> dict:
    """
    Load the pathways and their associated genes from a file.

    Args:
        pathways_dir (Path): Path to the file containing the pathway data.

    Returns:
    - dict: A dictionary mapping pathway names to lists of genes in each pathway.
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
                except ValueError:
                    logger.warning(
                        f"Non-integer gene ID found in pathway {pathway_name}, skipping pathway."
                    )
                    continue
                pathways[pathway_name] = genes
    except FileNotFoundError:
        logger.error(f"File not found: {pathways_dir}")
    except Exception as e:
        logger.error(f"An error occurred while loading pathways: {e}")
    return pathways


def load_propagation_file(
    file_path: Path, decompress: bool = True
) -> dict:
    """
    Load the propagation score data from a file.

    Args:
        file_path (Path): The path to the file to be loaded.
        decompress (bool): Whether to decompress the file.

    Returns:
    - dict: The loaded data containing propagation scores and other information.
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
        return {}
    except Exception as e:
        logger.error(
            f"An error occurred while loading propagation file: {e}"
        )
        return {}


def load_propagation_scores(propagation_file_path: Path) -> dict:
    """
    Load the propagation scores from a file.

    Args:
        propagation_file_path (Path): Path to the propagation scores file.

    Returns:
        dict: The propagation scores.
    """
    return load_propagation_file(
        propagation_file_path, decompress=True
    )


def load_pathways_and_propagation_scores(
    settings: Settings, scores_df: pd.DataFrame
) -> Tuple[Dict[str, Set[int]], Dict[int, float]]:
    """
    Load the pathways based on the provided configuration and the propagation scores.

    Args:
        settings (Settings): General settings for the pipeline.
        scores_df (pd.DataFrame): DataFrame containing the gene scores.

    Returns:
        Tuple: A tuple containing the pathways and a dictionary mapping gene IDs to their scores.
    """
    pathways_with_many_genes = load_pathways_genes(settings.pathway_file_dir)
    sorted_scores = scores_df.sort_values(by="GeneID").reset_index(drop=True)

    # Create the scores dictionary with only GeneID and Score
    scores = {
        gene_id: score
        for gene_id, score in zip(
            sorted_scores["GeneID"], sorted_scores["Score"]
        )
    }

    scores_keys = set(scores.keys())

    # Filter pathways based on gene counts within the specified range and intersection with scored genes
    genes_by_pathway = {
        pathway: set(map(str, genes)).intersection(scores_keys)
        for pathway, genes in pathways_with_many_genes.items()
        if settings.minimum_gene_per_pathway
        <= len(set(map(str, genes)).intersection(scores_keys))
        <= settings.maximum_gene_per_pathway
}


    return genes_by_pathway, scores




# ###################################################################GET FUNCTIONS############################################################################
def get_scores(score_path: Path) -> dict:
    """
    Load gene scores and P-values from a file.

    Args:
        score_path (Path): The path to the file containing the propagation scores.

    Returns:
        dict: A dictionary with gene IDs as keys and tuples of (Score, P-value) as values.
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
        logger.error(f"File not found: {score_path}")
        return {}
    except Exception as error:
        logger.error(f"An error occurred: {error}")
        return {}


# ###################################################################READ FUNCTIONS############################################################################

def read_network(network_filename: Path) -> nx.Graph:
    """
    Read a network from a file and return a NetworkX graph.

    Args:
        network_filename (Path): Path to the file containing the network data.

    Returns:
        nx.Graph: A graph object representing the network.
    """
    try:
        logger.info(f"Reading network from: {network_filename}")
        if not network_filename.exists():
            logger.error(f"Network file not found: {network_filename}")
            raise FileNotFoundError(f"Network file not found: {network_filename}")

        network = pd.read_table(
            network_filename, header=None, usecols=[0, 1, 2]
        )
        graph = nx.from_pandas_edgelist(network, 0, 1, 2)
        logger.info(f"Successfully loaded network with {len(graph.nodes())} nodes and {len(graph.edges())} edges")
        return graph
    except FileNotFoundError:
        logger.error(f"File not found: {network_filename}")
        raise
    except Exception as e:
        logger.error(f"An error occurred while reading network: {e}")
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
        prior_data = pd.read_excel(excel_input, engine="openpyxl")

        # Data preprocessing
        prior_data = prior_data.drop_duplicates(subset="GeneID")
        prior_data = prior_data[prior_data["Score"].notna()]
        prior_data = prior_data[prior_data["Score"] != "?"]
        prior_data = prior_data[
            prior_data["GeneID"].apply(lambda x: str(x).isdigit())
        ]
        prior_data["GeneID"] = prior_data["GeneID"].astype(int)
        prior_data = prior_data.reset_index(drop=True)
        return prior_data
    except FileNotFoundError:
        logger.error(f"File not found: {excel_input}")
        return pd.DataFrame()
    except Exception as e:
        logger.error(
            f"An error occurred while reading prior set: {e}"
        )
        return pd.DataFrame()


def read_enriched_pathways(condition_output_path: Path) -> dict:
    """
    Read the enriched pathways from the condition's output file.

    Args:
        condition_output_path (Path): Path to the condition's output Excel file.

    Returns:
        dict: Dictionary of pathways to FDR q-values.
    """
    try:
        enriched_pathway_df = pd.read_excel(condition_output_path)
        enriched_pathway_dict = dict(
            zip(
                enriched_pathway_df["Term"],
                enriched_pathway_df["FDR q-val"],
            )
        )
        return enriched_pathway_dict
    except FileNotFoundError:
        logger.error(f"File not found: {condition_output_path}")
        return {}
    except Exception as e:
        logger.error(
            f"An error occurred while reading enriched pathways: {e}"
        )
        return {}


# ###################################################################SAVE FUNCTIONS############################################################################


def save_file(obj, save_dir: Path, compress: bool = True) -> None:
    """
    Save an object to a file, with optional compression.

    Args:
        obj (object): The object to be saved.
        save_dir (Path): The directory where the file will be saved.
        compress (bool): Whether to compress the file.

    Returns:
    - None
    """
    try:
        obj = pickle.dumps(obj)
        if compress:
            obj = zlib.compress(obj)
        with save_dir.open("wb") as f:
            pickle.dump(obj, f)
    except Exception as e:
        logger.error(f"An error occurred while saving file: {e}")


def save_propagation_score(
    propagation_scores: pd.DataFrame,
    prior_set: pd.DataFrame,
    propagation_input: dict,
    genes_id_to_idx: dict,
    task,
    save_dir: Path,
    general_args: Settings,
) -> None:
    """
    Save the propagation scores to a compressed file.

    Args:
        propagation_scores (pd.DataFrame): DataFrame containing the propagated scores.
        prior_set (pd.DataFrame): DataFrame containing the prior set of gene scores.
        propagation_input (dict): Mapping of gene IDs to their input scores.
        genes_id_to_idx (dict): Mapping of gene IDs to their indices.
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

        save_dict = {
            "args": task,
            "prior_set": prior_set,
            "propagation_input": propagation_input,
            "gene_id_to_idx": genes_id_to_idx,
            "gene_prop_scores": propagation_scores,
        }

        save_file(save_dict, propagation_results_path)
    except Exception as e:
        logger.error(
            f"An error occurred while saving propagation score: {e}"
        )


# ###################################################################ELSE######################################################


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
        network = read_network(network_filename)
        gene_ids_in_prior_data = set(prior_data["GeneID"])
        filtered_network = network.subgraph(
            gene_ids_in_prior_data
        ).copy()
        return filtered_network
    except Exception as e:
        logger.error(
            f"An error occurred while filtering network: {e}"
        )
        return nx.Graph()


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
    pathway_genes: Set[int],
    scores: Dict[int, float],
    pathway_pvalues: Dict[str, float]
) -> Dict:
    """
    Process a single pathway for a given condition.

    Parameters:
    - pathway (str): Name of the pathway.
    - pathway_genes (Set[int]): Set of gene IDs in the pathway.
    - scores (Dict[int, float]): Mapping of gene IDs to scores.
    - pathway_pvalues (Dict[str, float]): Dictionary of pathway p-values.

    Returns:
    - Dict: Dictionary containing mean score and p-value.
    """
    pathway_gene_scores = [scores[gene_id] for gene_id in pathway_genes if gene_id in scores]
    mean_score = np.mean(pathway_gene_scores) if pathway_gene_scores else 0
    return {"Mean": mean_score, "P-value": pathway_pvalues.get(pathway, 1.0)}



def process_condition(
    condition_settings: ConditionSettings,
    settings: Settings,
    condition_data_df: pd.DataFrame,
    all_pathways: Dict[str, Dict],
    pathway_sizes: Dict[str, int],
    genes_by_pathway: Dict[str, Set[str]],
    scores: Dict[str, float],
):
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
    pathway_pvalues = {pathway: data.get("FDR q-val", 1.0) for pathway, data in condition_settings.pathways_statistics.items()}
    pathways_with_many_genes = load_pathways_genes(settings.pathway_file_dir)

    condition_data_df["GeneID"] = condition_data_df["GeneID"].astype(str)
    condition_data_df["Score"] = condition_data_df["Score"].astype(np.float32)

    scores.update({gene_id: score for gene_id, score in zip(condition_data_df["GeneID"], condition_data_df["Score"])})

    scores_keys = set(scores.keys())
    local_genes_by_pathway = {
        pathway: set(map(str, genes)).intersection(scores_keys)
        for pathway, genes in pathways_with_many_genes.items()
        if settings.minimum_gene_per_pathway <= len(set(map(str, genes)).intersection(scores_keys)) <= settings.maximum_gene_per_pathway
    }
    genes_by_pathway.update(local_genes_by_pathway)

    for pathway, genes in local_genes_by_pathway.items():
        pathway_sizes[pathway] = len(genes)

    for pathway, pathway_genes in local_genes_by_pathway.items():
        all_pathways.setdefault(pathway, {})
        result = process_single_pathway(pathway, pathway_genes, scores, pathway_pvalues)
        if result["P-value"] < settings.FDR_THRESHOLD:
            all_pathways[pathway][condition_settings.condition_name] = result
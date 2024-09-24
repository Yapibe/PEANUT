import logging
import pickle
import time
import zlib
from io import BytesIO
from pathlib import Path
from typing import Union

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
        dict: A dictionary mapping pathway names to lists of genes in each pathway.
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
                    genes = [int(gene) for gene in parts[2:]]
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
        dict: The loaded data containing propagation scores and other information.
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
) -> tuple:
    """
    Load the pathways based on the provided configuration and the propagation scores.

    Args:
        settings (Settings): General settings for the pipeline.
        scores_df (pd.DataFrame): DataFrame containing the gene scores.

    Returns:
        tuple: A tuple containing the pathways and a dictionary mapping gene IDs to their scores and p-values.
    """
    pathways_with_many_genes = load_pathways_genes(
        settings.pathway_file_dir
    )
    sorted_scores = scores_df.sort_values(by="GeneID").reset_index(
        drop=True
    )

    # Optimize data types
    sorted_scores["GeneID"] = sorted_scores["GeneID"].astype(int)
    sorted_scores["Score"] = pd.to_numeric(
        sorted_scores["Score"], downcast="float"
    )
    sorted_scores["P-value"] = pd.to_numeric(
        sorted_scores["P-value"], downcast="float"
    )

    scores = {
        gene_id: (score, pvalue)
        for gene_id, score, pvalue in zip(
            sorted_scores["GeneID"],
            sorted_scores["Score"],
            sorted_scores["P-value"],
        )
    }

    if settings.run_gsea:
        return pathways_with_many_genes, scores
    else:
        scores_keys = set(scores.keys())
        # Filter pathways based on gene counts within specified range and intersection with scored genes
        genes_by_pathway = {
            pathway: set(genes).intersection(scores_keys)
            for pathway, genes in pathways_with_many_genes.items()
            if settings.minimum_gene_per_pathway
            <= len(set(genes).intersection(scores_keys))
            <= settings.maximum_gene_per_pathway
        }
    return genes_by_pathway, scores


# ###################################################################GET FUNCTIONS############################################################################


def set_input_type(
    prior_data: pd.DataFrame, norm: bool = False
) -> pd.DataFrame:
    """
    Modify the 'Score' column based on the specified input type.

    Args:
        prior_data (pd.DataFrame): DataFrame containing all experimental data.
        norm (bool): Flag indicating whether to normalize the 'Score' column.

    Returns:
        pd.DataFrame: DataFrame with the modified 'Score' column.
    """
    modified_prior_data = prior_data.copy()
    if norm:
        modified_prior_data["Score"] = 1
    else:
        modified_prior_data["Score"] = prior_data["Score"].abs()
    return modified_prior_data


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


def read_temp_scores(file_name: Path) -> dict:
    """
    Read scores from a file into a dictionary.

    Args:
        file_name (Path): Path to the file containing the scores.

    Returns:
        dict: A dictionary mapping pathways to their scores.
    """
    try:
        scores = pd.read_csv(
            file_name,
            sep=" ",
            header=None,
            names=["Pathway", "Score"],
            index_col="Pathway",
        )
        return scores["Score"].to_dict()
    except FileNotFoundError:
        logger.error(f"File not found: {file_name}")
        return {}
    except Exception as e:
        logger.error(f"Error reading scores from {file_name}: {e}")
        return {}


def read_network(network_filename: Path) -> nx.Graph:
    """
    Read a network from a file and return a NetworkX graph.

    Args:
        network_filename (Path): Path to the file containing the network data.

    Returns:
        nx.Graph: A graph object representing the network.
    """
    try:
        network = pd.read_table(
            network_filename, header=None, usecols=[0, 1, 2]
        )
        return nx.from_pandas_edgelist(network, 0, 1, 2)
    except FileNotFoundError:
        logger.error(f"File not found: {network_filename}")
        return nx.Graph()
    except Exception as e:
        logger.error(f"An error occurred while reading network: {e}")
        return nx.Graph()


def read_prior_set(
    excel_input: Union[str, BytesIO, Path], is_bytes: bool = False
) -> pd.DataFrame:
    """
    Read prior data set from an Excel file and apply preprocessing.

    Args:
        excel_input (str or BytesIO or Path): Path to the Excel file or a BytesIO object containing the prior data.
        is_bytes (bool): Flag indicating if the input is a BytesIO object.

    Returns:
        pd.DataFrame: DataFrame containing the preprocessed prior data.
    """
    try:
        if is_bytes:
            prior_data = pd.read_excel(excel_input, engine="openpyxl")
        else:
            prior_data = pd.read_excel(excel_input, engine="openpyxl")

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
        None
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
    pathway_genes: list,
    experiment_data_filtered_df: pd.DataFrame,
    enriched_pathway_dict: dict,
    condition_name: str,
    FDR_THRESHOLD: float,
) -> dict:
    """
    Process a single pathway for a given condition.

    Args:
        pathway (str): Name of the pathway.
        pathway_genes (list): List of gene IDs in the pathway.
        experiment_data_filtered_df (pd.DataFrame): Filtered experimental data.
        enriched_pathway_dict (dict): Dictionary of enriched pathways and their FDR q-values.
        condition_name (str): Name of the condition.
        FDR_THRESHOLD (float): Threshold for significant P-values.

    Returns:
        dict: Dictionary containing mean score, significant genes, P-value, and trend.
    """
    pathway_genes_set = set(map(int, pathway_genes))
    genes_in_pathway_mask = experiment_data_filtered_df[
        "GeneID"
    ].isin(pathway_genes_set)
    pathway_filtered_genes = experiment_data_filtered_df.loc[
        genes_in_pathway_mask
    ]

    significant_genes_mask = (
        pathway_filtered_genes["P-value"] <= FDR_THRESHOLD
    )
    significant_genes = pathway_filtered_genes.loc[
        significant_genes_mask
    ]

    mean_score = (
        significant_genes["Score"].mean()
        if not significant_genes.empty
        else 0
    )

    if pathway in enriched_pathway_dict:
        p_value = enriched_pathway_dict[pathway]
        trend = "Up*" if mean_score > 0 else "Down*"
    else:
        p_value = 1
        trend = "Up" if mean_score > 0 else "Down"

    result = {
        "Mean": mean_score,
        "significant_genes": significant_genes.set_index(
            "GeneID"
        ).to_dict("index"),
        "P-value": p_value,
        "Trend": trend,
    }

    return result


def filter_experiment_data(condition_data_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter and preprocess experimental data.

    Parameters:
    - condition_data_df: Raw condition data DataFrame.

    Returns:
    - pd.DataFrame: Filtered and preprocessed DataFrame.
    """
    experiment_data_filtered_df = condition_data_df[
        condition_data_df['Score'] != 0
    ].copy()
    experiment_data_filtered_df['GeneID'] = experiment_data_filtered_df[
        'GeneID'
    ].astype(int)
    experiment_data_filtered_df['Score'] = pd.to_numeric(
        experiment_data_filtered_df['Score'], downcast='float'
    )
    experiment_data_filtered_df['P-value'] = pd.to_numeric(
        experiment_data_filtered_df['P-value'], downcast='float'
    )
    experiment_data_filtered_df['Symbol'] = experiment_data_filtered_df[
        'Symbol'
    ].astype('category')
    return experiment_data_filtered_df


def process_condition(
    condition_settings: ConditionSettings,
    settings: Settings,
    condition_data_df: pd.DataFrame,
    all_pathways: dict,
):
    """
    Processes experimental data for a specific condition.

    Parameters:
    - condition_settings: Settings specific to the condition.
    - settings: General settings for the pipeline.
    - condition_data_df: DataFrame containing data for the condition.
    - all_pathways: Dictionary to store pathway data across conditions.

    Returns:
    - None
    """
    enriched_pathway_dict = {
        pathway: data.get("FDR q-val", 1)
        for pathway, data in condition_settings.filtered_pathways.items()
    }

    homo_sapien_pathway_dict = load_pathways_genes(
        settings.pathway_file_dir
    )

    experiment_data_filtered_df = filter_experiment_data(
        condition_data_df
    )

    for pathway, pathway_genes in homo_sapien_pathway_dict.items():
        if pathway not in all_pathways:
            all_pathways[pathway] = {}

        result = process_single_pathway(
            pathway=pathway,
            pathway_genes=pathway_genes,
            experiment_data_filtered_df=experiment_data_filtered_df,
            enriched_pathway_dict=enriched_pathway_dict,
            condition_name=condition_settings.condition_name,
            FDR_THRESHOLD=settings.FDR_THRESHOLD,
        )

        all_pathways[pathway][condition_settings.condition_name] = result

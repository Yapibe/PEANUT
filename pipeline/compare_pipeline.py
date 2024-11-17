import os
import time
import logging
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml

from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import (
    load_pathways_genes,
    read_network,
    read_prior_set,
)

# Project and Directory Configuration
def get_project_root():
    """Get the project root directory dynamically."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PROJECT_ROOT = get_project_root()

config_path = os.path.join(PROJECT_ROOT, 'pipeline', 'config.yaml')
# Load configuration from config.yaml
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

# Configure logging using config.yaml
logging.basicConfig(level=config['logging']['level'], format=config['logging']['format'])
logger = logging.getLogger(__name__)

# Retrieve List of Excel Files
def get_file_list(input_dir, limit=None):
    """List all .xlsx files in the input directory."""
    files = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
    if limit:
        files = files[:limit]
    return files

# Retrieve List of Excel Files
def get_file_list(input_dir, limit=None):
    """List all .xlsx files in the input directory."""
    files = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
    if limit:
        files = files[:limit]
    return files

# Load Networks and Pathways Data
def load_data(networks, pathway_files, pathways_dir):
    """Load network and pathway data."""
    networks_data = {
        name: read_network(os.path.join(PROJECT_ROOT, 'pipeline', 'Data', 'H_sapiens', 'network', name))
        for name in networks
    }
    pathways_data = {
        file: load_pathways_genes(os.path.join(pathways_dir, file))
        for file in pathway_files
    }
    return networks_data, pathways_data

# Convert relative paths to absolute paths based on PROJECT_ROOT
DIRECTORIES = {key: os.path.join(PROJECT_ROOT, value) for key, value in config['directories'].items()}
# Ensure necessary directories exist
os.makedirs(DIRECTORIES['summary_base'], exist_ok=True)
# Analysis Settings
NETWORKS = config['analysis_settings']['networks']
PATHWAY_FILES = config['analysis_settings']['pathway_files']
PROP_METHODS = config['analysis_settings']['prop_methods']
ALPHAS = config['analysis_settings']['alphas']
MAX_WORKERS = config['analysis_settings']['max_workers']

# Pre-Calculated Data
PRE_CALCULATED_DATA = config['pre_calculated_data']

def get_pre_calculated_data(disease_name):
    """Retrieve pre-calculated related pathways for a given disease."""
    return PRE_CALCULATED_DATA.get(disease_name, [])


# Analysis Functions
def run_analysis(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file, related_pathways):
    """
    Perform propagation and enrichment analysis.
    
    Args:
        test_name (str): Name of the test.
        prior_data: Prior data for the analysis.
        network: Network data.
        network_name (str): Name of the network.
        alpha (float): Alpha value for the analysis.
        method (str): Propagation method.
        output_path (str): Path to save the output.
        pathway_file (str): Pathway file used.
        related_pathways (list): Related pathways for the disease.
    """
    general_args = GeneralArgs(
        network=network_name,
        pathway_file=pathway_file,
        method=method,
        alpha=alpha if method in ['PEANUT', 'ABS_PROP'] else 1,
        run_propagation= False,
        input_type='abs_Score' if method in ['PEANUT', 'ABS_PROP', 'ABS_SCORE'] else 'Score'
    )
    if method == 'NGSEA':
        general_args.run_NGSEA = True

    if method == 'PEANUT':
        general_args.run_gsea = False


    if general_args.run_propagation:
        perform_propagation(test_name, general_args, network, prior_data)

    perform_enrichment(test_name, general_args, output_path, related_pathways)

def get_pathway_rank(output_path, pathway_name):
    """
    Determine the rank and FDR q-value of a specific pathway.

    Args:
        output_path (str): Path to the results file.
        pathway_name (str): Name of the pathway.

    Returns:
        tuple: (rank, fdr_q_val, results_df) if found, else (None, None, None)
    """
    try:
        results_df = pd.read_excel(output_path)
        if 'FDR q-val' not in results_df.columns:
            logger.error(f"'FDR q-val' column not found in {output_path}.")
            return None, None, None

        results_df = results_df.sort_values(by='FDR q-val', ascending=True).reset_index(drop=True)
        results_df['Rank'] = results_df.index + 1

        pathway_row = results_df[results_df['Term'] == pathway_name]
        if not pathway_row.empty:
            rank = pathway_row['Rank'].values[0]
            fdr_q_val = pathway_row['FDR q-val'].values[0]
            return rank, fdr_q_val, results_df
        else:
            logger.error(f"Pathway '{pathway_name}' not found in {output_path}.")
    except Exception as e:
        logger.error(f"Error reading {output_path}: {e}")
    return None, None, None

def process_single_file(network, pathway_file, network_name, alpha, method, file_name, loaded_pathways, related_pathways):
    """
    Process a single file: perform analysis and collect results.
    """
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    prior_data = read_prior_set(os.path.join(DIRECTORIES['input'], file_name))
    output_dir = os.path.join(
        DIRECTORIES['output_base'], method, network_name, pathway_file, f"alpha_{alpha}"
    )
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    run_analysis(file_name, prior_data, network, network_name, alpha, method, output_file_path, pathway_file, related_pathways)
    rank, fdr_q_val, higher_rank_df = get_pathway_rank(output_file_path, pathway_name)

    higher_ranked = []
    higher_ranked_ranks = []
    genes = []
    if rank is not None:
        logger.debug(f"Available columns in {output_file_path}: {higher_rank_df.columns.tolist()}")

        higher_ranked_df = higher_rank_df[higher_rank_df['Rank'] < rank]
        higher_ranked = higher_ranked_df['Term'].tolist()
        higher_ranked_ranks = higher_ranked_df['Rank'].tolist()

        for pathway, pathway_rank in zip(higher_ranked, higher_ranked_ranks):
            gene_list = loaded_pathways[pathway_file].get(pathway, [])
            genes.extend(gene_list)

    return {
        'result': {
            'Dataset': dataset_name,
            'Pathway': pathway_name,
            'Network': network_name,
            'Pathway file': pathway_file,
            'Alpha': alpha,
            'Method': method,
            'Rank': rank,
            'FDR q-val': fdr_q_val,
            'Significant': int(float(fdr_q_val) < 0.05) if fdr_q_val else 0,
        },
        'higher_ranked_pathways': higher_ranked,
        'higher_ranked_ranks': higher_ranked_ranks,
        'higher_ranked_genes': genes
    }


# Summary Generation
def generate_summary_df(filtered_df, prop_methods):
    """
    Generate a summary DataFrame with pivoted results.

    Args:
        filtered_df (pd.DataFrame): Filtered DataFrame of results.
        prop_methods (list): List of propagation methods.

    Returns:
        pd.DataFrame: Summary DataFrame.
    """
    # Create pivot table
    pivot_df = filtered_df.pivot_table(
        index=['Dataset', 'Pathway'],
        columns='Method',
        values=['Rank', 'Significant'],
        aggfunc='first'
    ).reset_index()

    # Flatten the MultiIndex columns
    pivot_df.columns = [
        f'{metric} {method}' if method else metric
        for metric, method in pivot_df.columns
    ]

    # Define the ordered columns
    ordered_columns = ['Dataset', 'Pathway']
    for method in prop_methods:
        for metric in ['Rank', 'Significant']:
            column_name = f'{metric} {method}'
            ordered_columns.append(column_name)

    # Reorder the columns
    pivot_df = pivot_df[ordered_columns]

    # Handle missing methods by filling with NaN
    for method in prop_methods:
        for metric in ['Rank', 'Significant']:
            column_name = f'{metric} {method}'
            if column_name not in pivot_df.columns:
                pivot_df[column_name] = np.nan

    # Calculate average values for numerical columns
    avg_values = pivot_df.select_dtypes(include=[np.number]).mean().to_dict()
    avg_row = {**{col: '' for col in ['Dataset', 'Pathway']},
               **avg_values}

    # Append the average row to the DataFrame
    summary_df = pd.concat([pivot_df, pd.DataFrame([avg_row])], ignore_index=True)

    # Sort the DataFrame by 'Dataset'
    summary_df.sort_values(by='Dataset', inplace=True)

    return summary_df



FILE_LIST = get_file_list(DIRECTORIES['input'])
logger.info(f"Found {len(FILE_LIST)} .xlsx files in {DIRECTORIES['input']}")

NETWORKS_DATA, PATHWAYS_DATA = load_data(NETWORKS, PATHWAY_FILES, DIRECTORIES['pathways'])

# Main Processing Function
def main():
    start_time = time.time()
    futures = []
    results = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for network_name in NETWORKS:
            for pathway_file in PATHWAY_FILES:
                for alpha in ALPHAS:
                    for file_name in FILE_LIST:
                        dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                        if pathway_name in PATHWAYS_DATA[pathway_file]:
                            related_pathways = get_pre_calculated_data(pathway_name)
                            for prop_method in PROP_METHODS:
                                futures.append(
                                    executor.submit(
                                        process_single_file,
                                        NETWORKS_DATA[network_name],
                                        pathway_file,
                                        network_name,
                                        alpha,
                                        prop_method,
                                        file_name,
                                        PATHWAYS_DATA,
                                        related_pathways
                                    )
                                )

        for future in as_completed(futures):
            result = future.result()
            results.append(result['result'])

    results_df = pd.DataFrame(results)

    # Save Summary DataFrames
    for network_name in NETWORKS:
        for pathway_file in PATHWAY_FILES:
            for alpha in ALPHAS:
                filtered_df = results_df[
                    (results_df['Network'] == network_name) &
                    (results_df['Pathway file'] == pathway_file) &
                    (results_df['Alpha'] == alpha)
                ]
                summary_df = generate_summary_df(filtered_df, PROP_METHODS)

                summary_output_dir = os.path.join(
                    DIRECTORIES['summary_base'],
                    network_name,
                    pathway_file,
                    f"alpha_{alpha}"
                )
                os.makedirs(summary_output_dir, exist_ok=True)
                rankings_output_path = os.path.join(
                    summary_output_dir,
                    f'rankings_summary_{network_name}_{pathway_file}_alpha_{alpha}.xlsx'
                )
                summary_df.to_excel(rankings_output_path, index=False)
                logger.info(f"Rankings summary saved to {rankings_output_path}")

    elapsed_time = time.time() - start_time
    logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
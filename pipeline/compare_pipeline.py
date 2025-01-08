import os
import time
import logging
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml

from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation, filter_network_genes
from utils import (
    load_pathways_genes,
    read_network,
    read_prior_set,
    set_input_type
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
IMPUTATION_MODES = config['analysis_settings']['imputation_modes']

# Analysis Functions
def run_analysis(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file, imputation_mode):
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
        imputation_mode (str): Imputation mode.
    """
    general_args = GeneralArgs(
        network=network_name,
        pathway_file=pathway_file,
        method=method,
        alpha=alpha if method == 'PEANUT' else 1,
        run_propagation= True,
        input_type='Abs_Score',
        imputation_mode=imputation_mode
    )
    if general_args.run_propagation:
        perform_propagation(test_name, general_args, network, prior_data)

    perform_enrichment(test_name, general_args, output_path)


def get_pathway_rank(output_path: str, pathway_name: str) -> tuple:
    """
    Determine the rank, FDR q-value, and group of a specific pathway.
    
    Args:
        output_path (str): Path to the results file.
        pathway_name (str): Name of the pathway.
    
    Returns:
        tuple: (rank, fdr_q_val, group, results_df) if found, else (None, None, None, None)
    """
    try:
        results_df = pd.read_excel(output_path)

        # Pathways are already grouped and sorted
        pathway_row = results_df[results_df['Name'] == pathway_name]
        if not pathway_row.empty:
            rank = pathway_row['Rank'].values[0]
            fdr_q_val = pathway_row['FDR q-val'].values[0]
            group = pathway_row['Group'].values[0] if 'Group' in pathway_row else None
            return rank, fdr_q_val, group, results_df
        else:
            logger.error(f"Pathway '{pathway_name}' not found in {output_path}.")
    except Exception as e:
        logger.error(f"Error reading {output_path}: {e}")
    return None, None, None, None


def process_single_file(
    network: dict,
    pathway_file: str,
    network_name: str,
    alpha: float,
    method: str,
    file_name: str,
    loaded_pathways: dict,
    imputation_mode: str
) -> dict:
    """
    Process a single file: perform analysis and collect results.

    Args:
        network (dict): Network data.
        pathway_file (str): Pathway file name.
        network_name (str): Name of the network.
        alpha (float): Alpha value for analysis.
        method (str): Propagation method.
        file_name (str): Name of the input file.
        loaded_pathways (dict): Loaded pathways data.
        imputation_mode (str): Imputation mode.

    Returns:
        dict: Results and higher-ranked pathways information.
    """
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    prior_data = read_prior_set(os.path.join(DIRECTORIES['input'], file_name))
    output_dir = os.path.join(
        DIRECTORIES['output_base'], method, network_name, pathway_file, f"alpha_{alpha}"
    )
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    run_analysis(
        test_name=file_name,
        prior_data=prior_data,
        network=network,
        network_name=network_name,
        alpha=alpha,
        method=method,
        output_path=output_file_path,
        pathway_file=pathway_file,
        imputation_mode=imputation_mode
    )

    # Get rank and FDR q-val of the associated pathway
    rank, fdr_q_val, group, results_df = get_pathway_rank(output_file_path, pathway_name)

    higher_ranked = []
    higher_ranked_ranks = []
    genes = []
    ks_significant = False  # Default value

    if rank is not None and results_df is not None:
        # Get pathways ranked higher than the associated pathway
        higher_ranked_df = results_df[results_df['Rank'] < rank]
        higher_ranked = higher_ranked_df['Name'].tolist()
        higher_ranked_ranks = higher_ranked_df['Rank'].tolist()

        for pathway, pathway_rank in zip(higher_ranked, higher_ranked_ranks):
            gene_list = loaded_pathways[pathway_file].get(pathway, [])
            genes.extend(gene_list)

        # Retrieve KS significance for PEANUT method only
        if method == 'PEANUT':
            pathway_row = results_df[results_df['Name'] == pathway_name]
            if not pathway_row.empty:
                ks_significant = bool(pathway_row['ks_significant'].values[0])
                significant = int(ks_significant and (fdr_q_val < 0.05)) if fdr_q_val is not None else 0
        else:
            # For other methods, only use FDR q-val
            significant = int(fdr_q_val < 0.05) if fdr_q_val is not None else 0
    else:
        significant = 0
    prior_data_scores = set_input_type(prior_data)
    filtered_prior_data = filter_network_genes(prior_data_scores, network)
    # Get genes for the associated pathway
    pathway_genes = loaded_pathways[pathway_file].get(pathway_name, [])
    
    # Get scores for pathway genes that are in the network
    pathway_gene_scores = filtered_prior_data[filtered_prior_data['GeneID'].isin(pathway_genes)]['Score']
    
    # Calculate mean score to determine regulation direction
    mean_score = pathway_gene_scores.mean() if not pathway_gene_scores.empty else 0
    is_upregulated = mean_score > 0
    return {
        'result': {
            'Dataset': dataset_name,
            'Pathway': pathway_name,
            'Group': group,
            'Network': network_name,
            'Pathway file': pathway_file,
            'Alpha': alpha,
            'Method': method,
            'Rank': rank,
            'FDR q-val': fdr_q_val,
            'Significant': significant,
            'Up-regulated': is_upregulated,
            'Imputation Mode': imputation_mode
        },
        'higher_ranked_pathways': higher_ranked,
        'higher_ranked_ranks': higher_ranked_ranks,
        'higher_ranked_genes': genes
    }



# Summary Generation
def generate_summary_df(filtered_df):
    """
    Generate a summary DataFrame with pivoted results for imputation modes.

    Args:
        filtered_df (pd.DataFrame): Filtered DataFrame of results.

    Returns:
        pd.DataFrame: Summary DataFrame.
    """
    # Create pivot table including 'Group'
    pivot_df = filtered_df.pivot_table(
        index=['Dataset', 'Pathway'],
        columns='Imputation Mode',
        values=['Rank', 'Significant', 'Group'],
        aggfunc='first'
    ).reset_index()

    # Flatten the MultiIndex columns
    pivot_df.columns = [
        f'{metric} {mode}' if mode else metric
        for metric, mode in pivot_df.columns
    ]

    # Define the ordered columns
    ordered_columns = ['Dataset', 'Pathway']
    for mode in IMPUTATION_MODES:
        for metric in ['Rank', 'Significant', 'Group']:
            column_name = f'{metric} {mode}'
            ordered_columns.append(column_name)

    # Reorder the columns
    pivot_df = pivot_df[ordered_columns]

    # Handle missing imputation modes by filling with NaN
    for mode in IMPUTATION_MODES:
        for metric in ['Rank', 'Significant', 'Group']:
            column_name = f'{metric} {mode}'
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
                                for prop_method in PROP_METHODS:
                                    for imputation_mode in IMPUTATION_MODES:
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
                                                imputation_mode
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
                # Generate the summary focused on imputation modes
                summary_df = generate_summary_df(filtered_df)

                summary_output_dir = os.path.join(
                    DIRECTORIES['summary_base'],
                    network_name,
                    pathway_file,
                    f"alpha_{alpha}",
                )
                os.makedirs(summary_output_dir, exist_ok=True)
                rankings_output_path = os.path.join(
                    summary_output_dir,
                    f'rankings_summary.xlsx'
                )
                summary_df.to_excel(rankings_output_path, index=False)
                logger.info(f"Rankings summary saved to {rankings_output_path}")

    elapsed_time = time.time() - start_time
    logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()
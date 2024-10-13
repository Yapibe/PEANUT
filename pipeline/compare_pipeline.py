import os
import time
import logging
import pandas as pd
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import (
    load_pathways_genes,
    read_network,
    read_prior_set
)

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Project and Directory Configuration
def get_project_root():
    """Get the project root directory dynamically."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

PROJECT_ROOT = get_project_root()

DIRECTORIES = {
    'input': os.path.join(PROJECT_ROOT, 'pipeline', 'Inputs', 'experiments_data', 'GSE', 'XLSX'),
    'output_base': os.path.join(PROJECT_ROOT, 'pipeline', 'Outputs', 'NGSEA'),
    'summary_base': os.path.join(PROJECT_ROOT, 'pipeline', 'Outputs', 'NGSEA', 'Summary'),
    'pathways': os.path.join(PROJECT_ROOT, 'pipeline', 'Data', 'H_sapiens', 'pathways')
}

# Analysis Settings
NETWORKS = ['H_sapiens']
PATHWAY_FILES = ['kegg']
PROP_METHODS = ['ABS_PROP']
ALPHAS = [0.1]
MAX_WORKERS = 24

# Ensure necessary directories exist
os.makedirs(DIRECTORIES['summary_base'], exist_ok=True)

# Retrieve List of Excel Files
def get_file_list(input_dir, limit=None):
    """List all .xlsx files in the input directory."""
    files = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
    if limit:
        files = files[:limit]
    return files

FILE_LIST = get_file_list(DIRECTORIES['input'])
logger.info(f"Found {len(FILE_LIST)} .xlsx files in {DIRECTORIES['input']}")

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

NETWORKS_DATA, PATHWAYS_DATA = load_data(NETWORKS, PATHWAY_FILES, DIRECTORIES['pathways'])
logger.info("Networks loaded and pathway densities calculated.")

# Pre-Calculated Data
PRE_CALCULATED_DATA = {
    'KEGG_THYROID_CANCER': (1.954022989, 29, 4),
    'KEGG_NON_SMALL_CELL_LUNG_CANCER': (2.062678063, 54, 4),
    'KEGG_ACUTE_MYELOID_LEUKEMIA': (2.039721946, 57, 5),
    'KEGG_COLORECTAL_CANCER': (2.107879429, 62, 4),
    'KEGG_GLIOMA': (2.039072039, 65, 4),
    'KEGG_RENAL_CELL_CARCINOMA': (2.11032967, 70, 5),
    'KEGG_PANCREATIC_CANCER': (2.140724947, 70, 5),
    'KEGG_PROSTATE_CANCER': (2.106543291, 89, 5),
    'KEGG_DILATED_CARDIOMYOPATHY': (2.454394693, 90, 7),
    'KEGG_PARKINSONS_DISEASE': (1.949191984, 130, 5),
    'KEGG_ALZHEIMERS_DISEASE': (2.150090558, 166, 5),
    'KEGG_HUNTINGTONS_DISEASE': (1.758788428, 182, 4)
}

def get_pre_calculated_data(pathway_name):
    """Retrieve pre-calculated data for a given pathway."""
    return PRE_CALCULATED_DATA.get(pathway_name, (None, None, None))

# Analysis Functions
def run_analysis(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file):
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
    """
    general_args = GeneralArgs(
        network=network_name,
        pathway_file=pathway_file,
        method=method,
        alpha=alpha if method == 'ABS_PROP' else 1,
        run_propagation=True
    )

    if method == 'ABS_PROP':
        general_args.input_type = 'abs_Score'
    elif method == 'NGSEA':
        general_args.run_NGSEA = True

    if general_args.run_propagation:
        perform_propagation(test_name, general_args, network, prior_data)

    perform_enrichment(test_name, general_args, output_path)

def get_pathway_rank(output_path, pathway_name):
    """
    Determine the rank and FDR q-value of a specific pathway.
    
    Args:
        output_path (str): Path to the results file.
        pathway_name (str): Name of the pathway.
    
    Returns:
        tuple: (rank, fdr_q_val) if found, else (None, None)
    """
    try:
        results_df = pd.read_excel(output_path)
        if 'FDR q-val' not in results_df.columns:
            logger.error(f"'FDR q-val' column not found in {output_path}.")
            return None, None

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

def process_single_file(network, pathway_file, network_name, alpha, method, file_name, loaded_pathways):
    """
    Process a single file: perform analysis and collect results.
    """
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    # prior_data = read_prior_set(os.path.join(DIRECTORIES['input'], file_name))
    output_dir = os.path.join(
        DIRECTORIES['output_base'], method, network_name, pathway_file, f"alpha_{alpha}"
    )
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    pathway_density, num_genes, avg_diameter = get_pre_calculated_data(pathway_name)

    # run_analysis(file_name, prior_data, network, network_name, alpha, method, output_file_path, pathway_file)
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
            'Density': pathway_density,
            'Num Genes': num_genes,
            'Avg Diameter': avg_diameter
        },
        'higher_ranked_pathways': higher_ranked,
        'higher_ranked_ranks': higher_ranked_ranks,
        'higher_ranked_genes': genes
    }

# Aggregation Functions
def aggregate_higher_ranked(pathways_dict, genes_dict, ranks_dict, summary_dir):
    """
    Aggregate and save the most common pathways and genes, along with average ranks.
    """
    pathway_counter = defaultdict(int)
    gene_counter = defaultdict(int)
    pathway_rank_accumulator = defaultdict(list)

    for dataset, pathways in pathways_dict.items():
        for pathway in pathways:
            pathway_counter[pathway] += 1

    for dataset, ranks in ranks_dict.items():
        for pathway, rank in zip(pathways_dict[dataset], ranks):
            pathway_rank_accumulator[pathway].append(rank)

    for dataset, genes in genes_dict.items():
        for gene in genes:
            gene_counter[gene] += 1

    most_common_pathways = sorted(pathway_counter.items(), key=lambda x: x[1], reverse=True)
    most_common_genes = sorted(gene_counter.items(), key=lambda x: x[1], reverse=True)

    # Compute average ranks
    average_rankings = {pathway: np.mean(ranks) for pathway, ranks in pathway_rank_accumulator.items()}

    # Create DataFrame for pathways with average ranks
    pathways_df = pd.DataFrame(most_common_pathways, columns=['Pathway', 'Count'])

    # **Filter out pathways with count less than 6**
    pathways_df = pathways_df[pathways_df['Count'] >= 6]

    pathways_df['Average Rank'] = pathways_df['Pathway'].map(average_rankings)
    pathways_df.to_csv(os.path.join(summary_dir, 'common_pathways.csv'), index=False)

    # Create DataFrame for genes
    genes_df = pd.DataFrame(most_common_genes, columns=['Gene', 'Count'])
    genes_df.to_csv(os.path.join(summary_dir, 'common_genes.csv'), index=False)

    logger.info("Aggregated common pathways and genes saved with average ranks.")

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
        index=['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter'],
        columns='Method',
        values=['Rank', 'Significant'],
        aggfunc='first'
    ).reset_index()

    # Define the new column names based on PROP_METHODS
    new_columns = ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter']
    for method in prop_methods:
        for metric in ['Rank', 'Significant']:
            new_columns.append(f'{method} {metric}')

    # Assign the new column names
    pivot_df.columns = new_columns

    # Define the ordered columns
    ordered_columns = ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter']
    for method in prop_methods:
        ordered_columns += [f'{method} Rank', f'{method} Significant']

    # Reorder the columns
    pivot_df = pivot_df[ordered_columns]

    # Handle missing methods by filling with NaN
    for method in prop_methods:
        if f'{method} Rank' not in pivot_df.columns:
            pivot_df[f'{method} Rank'] = np.nan
        if f'{method} Significant' not in pivot_df.columns:
            pivot_df[f'{method} Significant'] = np.nan

    # Calculate average values for numerical columns
    avg_values = pivot_df.select_dtypes(include=[np.number]).mean().to_dict()
    avg_row = {**{col: '' for col in ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter']},
               **avg_values}
    
    # Append the average row to the DataFrame
    summary_df = pd.concat([pivot_df, pd.DataFrame([avg_row])], ignore_index=True)

    # Sort the DataFrame by 'Dataset'
    summary_df.sort_values(by='Dataset', inplace=True)
    
    return summary_df

# Main Processing Function
def main():
    start_time = time.time()
    futures = []
    results = []
    aggregated_pathways = defaultdict(list)
    aggregated_pathway_ranks = defaultdict(list)
    aggregated_genes = defaultdict(list)

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for network_name in NETWORKS:
            for pathway_file in PATHWAY_FILES:
                for alpha in ALPHAS:
                    for file_name in FILE_LIST:
                        dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                        if pathway_name in PATHWAYS_DATA[pathway_file]:
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
                                        PATHWAYS_DATA
                                    )
                                )

        for future in as_completed(futures):
            result = future.result()
            results.append(result['result'])
            # Aggregate pathways and genes
            dataset = result['result']['Dataset']
            aggregated_pathways[dataset].extend(result['higher_ranked_pathways'])
            aggregated_genes[dataset].extend(result['higher_ranked_genes'])
            aggregated_pathway_ranks[dataset].extend(result['higher_ranked_ranks'])

    results_df = pd.DataFrame(results)
    aggregate_higher_ranked(aggregated_pathways, aggregated_genes, aggregated_pathway_ranks, DIRECTORIES['summary_base'])

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
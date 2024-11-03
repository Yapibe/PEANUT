import os
import time
import logging
import pandas as pd
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml

from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import (
    load_pathways_genes,
    read_network,
    read_prior_set,
    load_disease_to_pathway,
    extract_subgraph,
    compute_density,
    compute_average_diameter,
    save_pathway_attributes_to_file,
    load_pathway_attributes_from_file
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

# Pre-Calculated Data
PRE_CALCULATED_DATA = config['pre_calculated_data']

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
    # Include method in test_name to prevent overwriting
    test_name_with_method = f"{method}_{test_name}"

    general_args = GeneralArgs(
        network=network_name,
        pathway_file=pathway_file,
        method=method,
        alpha=alpha if method in ['PEANUT', 'ABS_PROP'] else 1,
        run_propagation= True,
        input_type='abs_Score' if method in ['PEANUT', 'ABS_PROP', 'ABS_SCORE'] else 'Score'
    )
    if method == 'NGSEA':
        general_args.run_NGSEA = True

    if method == 'PEANUT':
        general_args.run_gsea = False


    if general_args.run_propagation:
        perform_propagation(test_name_with_method, general_args, network, prior_data)

    perform_enrichment(test_name_with_method, general_args, output_path)

def build_search_terms(pathway_name, disease_to_pathways=None):
    """
    Build a set of search terms for a given pathway, including associated terms.

    Args:
        pathway_name (str): Name of the pathway.
        disease_to_pathways (dict, optional): Mapping of pathways to their associated pathway terms.

    Returns:
        set: Set of search terms including the original pathway and associated terms.
    """
    search_terms = {pathway_name}  # Always include the original pathway
    if disease_to_pathways is not None and pathway_name in disease_to_pathways:
        search_terms.update(disease_to_pathways[pathway_name])
    return search_terms

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

def process_single_file(network, pathway_file, network_name, alpha, method, file_name, loaded_pathways):
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

    pathway_density, num_genes, avg_diameter = get_pre_calculated_data(pathway_name)

    # Include method in test_name to prevent overwriting
    test_name = f"{dataset_name}_{pathway_name}"
    run_analysis(test_name, prior_data, network, network_name, alpha, method, output_file_path, pathway_file)
    rank, fdr_q_val, higher_rank_df = get_pathway_rank(output_file_path, pathway_name)

    higher_ranked = []
    higher_ranked_ranks = []
    genes = []
    if rank is not None and higher_rank_df is not None:
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
            'Significant': int(float(fdr_q_val) <= 0.05) if fdr_q_val is not None else 0,
            'Density': pathway_density,
            'Num Genes': num_genes,
            'Avg Diameter': avg_diameter
        },
        'higher_ranked_pathways': higher_ranked,
        'higher_ranked_ranks': higher_ranked_ranks,
        'higher_ranked_genes': genes
    }

def aggregate_higher_ranked_non_disease(
    pathways_dict, genes_dict, ranks_dict, summary_dir, method, network,
    search_terms_dict, pathway_attributes=None, min_appearance=6
):
    """
    Aggregate and save the most common non-disease-associated pathways and genes per method,
    along with average ranks and pathway attributes.

    Args:
        pathways_dict (dict): Dictionary of pathways per combined dataset_pathway key.
        genes_dict (dict): Dictionary of genes per combined dataset_pathway key.
        ranks_dict (dict): Dictionary of pathway ranks per combined dataset_pathway key.
        summary_dir (str): Directory to save the aggregated results.
        method (str): The propagation method.
        network (networkx.Graph): The network graph.
        search_terms_dict (dict): Pre-built search terms per pathway.
        pathway_attributes (dict, optional): Pre-loaded pathway attributes. If None, they will be computed.
        min_appearance (int): Minimum number of times a pathway must appear to be included.
    """
    pathway_counter = defaultdict(int)
    pathway_rank_accumulator = defaultdict(list)
    if pathway_attributes is None:
        pathway_attributes = {}

    for combined_key, pathways in pathways_dict.items():
        dataset, pathway = combined_key.split('_', 1)
        search_terms = search_terms_dict.get(pathway, {pathway})

        for higher_pathway in pathways:
            # Check if the pathway contains any of the search terms
            if any(term.lower() in higher_pathway.lower() for term in search_terms):
                continue  # Skip disease-related pathways

            pathway_counter[higher_pathway] += 1
            # Accumulate ranks
            try:
                rank = ranks_dict[combined_key][pathways.index(higher_pathway)]
            except (IndexError, ValueError) as e:
                logger.error(f"Error retrieving rank for pathway '{higher_pathway}' in '{combined_key}': {e}")
                continue  # Skip this pathway

            pathway_rank_accumulator[higher_pathway].append(rank)

    # Filter pathways that appear at least `min_appearance` times
    frequent_pathways = [
        pathway for pathway, count in pathway_counter.items() if count >= min_appearance
    ]

    # If pathway_attributes were not provided, compute and save them
    if not pathway_attributes:
        for pathway in frequent_pathways:
            # Collect genes for the pathway using PATHWAYS_DATA
            genes = set()
            for pathway_file, pathways_in_file in PATHWAYS_DATA.items():
                if pathway in pathways_in_file:
                    genes.update(pathways_in_file[pathway])
                    break  # Found the pathway

            # Ensure genes are in the network
            genes_in_network = [gene for gene in genes if gene in network.nodes]
            if not genes_in_network:
                continue  # Skip pathways with no genes in the network

            # Compute attributes
            subgraph = extract_subgraph(network, genes_in_network)
            density = compute_density(subgraph)
            avg_diameter = compute_average_diameter(subgraph)
            num_genes = subgraph.number_of_nodes()

            pathway_attributes[pathway] = {
                'Gene Count': num_genes,
                'Density': density,
                'Avg Diameter': avg_diameter
            }

        # Save pathway attributes to a file for future use
        attributes_file = os.path.join(summary_dir, f'pathway_attributes_{method}.json')
        save_pathway_attributes_to_file(pathway_attributes, attributes_file)
        logger.info(f"Pathway attributes saved to {attributes_file}")

    # Compute average ranks for the frequent pathways
    average_rankings = {
        pathway: np.mean(ranks) 
        for pathway, ranks in pathway_rank_accumulator.items() 
        if pathway in frequent_pathways
    }

    # Prepare data for DataFrame
    pathway_data_list = []
    for pathway in frequent_pathways:
        if pathway in pathway_attributes:
            attributes = pathway_attributes[pathway]
            pathway_data_list.append({
                'Pathway': pathway,
                'Count': pathway_counter[pathway],
                'Average Rank': average_rankings.get(pathway, np.nan),
                'Gene Count': attributes.get('Gene Count', np.nan),
                'Density': attributes.get('Density', np.nan),
                'Avg Diameter': attributes.get('Avg Diameter', np.nan)
            })

    # Create DataFrame and save
    false_positive_df = pd.DataFrame(pathway_data_list)
    output_file = os.path.join(summary_dir, f'false_positive_pathways_{method}.csv')
    false_positive_df.to_csv(output_file, index=False)
    logger.info(f"False positive pathways saved to {output_file}")

def analyze_disease_pathway_rankings(search_terms_dict, pathways_df, summary_dir, method):
    """
    For each disease and its associated pathways, collect pathway attributes for analysis.

    Args:
        search_terms_dict (dict): Pre-built search terms per pathway.
        pathways_df (pd.DataFrame): DataFrame containing pathway rankings and attributes.
        summary_dir (str): Directory to save the analysis results.
        method (str): The propagation method.
    """
    analysis_results = []

    # Create mappings for quick lookup
    pathway_rank_mapping = pathways_df.set_index('Pathway')['Average Rank'].to_dict()
    attribute_columns = [col for col in pathways_df.columns if col not in ['Pathway', 'Count', 'Average Rank']]
    pathway_attribute_mappings = {col: pathways_df.set_index('Pathway')[col].to_dict() for col in attribute_columns}

    for pathway_name, search_terms in search_terms_dict.items():
        # Find pathways in the ranking data that contain any of the search terms
        matching_pathways = [
            pathway for pathway in pathway_rank_mapping.keys()
            if any(term.lower() in pathway.lower() for term in search_terms)
        ]

        if not matching_pathways:
            logger.warning(f"No matching pathways found for '{pathway_name}' in method {method}")
            continue

        for pathway in matching_pathways:
            pathway_data = {
                'Original Pathway': pathway_name,
                'Matched Pathway': pathway,
                'Pathway Average Rank': pathway_rank_mapping[pathway],
            }
            # Add attributes
            for attr in attribute_columns:
                pathway_data[attr] = pathway_attribute_mappings[attr].get(pathway, np.nan)

            analysis_results.append(pathway_data)

    # Convert the results to a DataFrame
    analysis_df = pd.DataFrame(analysis_results)

    # Save the analysis to a CSV file
    analysis_output_path = os.path.join(summary_dir, f'disease_pathway_attribute_analysis_{method}.csv')
    analysis_df.to_csv(analysis_output_path, index=False)
    logger.info(f"Disease pathway attribute analysis saved to {analysis_output_path}")

def generate_summary_df(filtered_df, prop_methods):
    """
    Generate a summary DataFrame with pivoted results across all methods and include an average row.

    Args:
        filtered_df (pd.DataFrame): Filtered DataFrame of results across methods.
        prop_methods (list): List of propagation methods.

    Returns:
        pd.DataFrame: Summary DataFrame containing rank and significance per method, along with an average row.
    """
    # Create pivot table to separate Rank and Significant columns per method
    pivot_df = filtered_df.pivot_table(
        index=['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter'],
        columns='Method',
        values=['Rank', 'Significant'],
        aggfunc='first'
    ).reset_index()

    # Flatten the MultiIndex columns in a way that includes both the method and metric names
    pivot_df.columns = [f'{col[1]} {col[0]}' if col[1] else col[0] for col in pivot_df.columns]

    # Reorder the columns to ensure they are correctly aligned
    ordered_columns = ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter']
    for method in prop_methods:
        ordered_columns.append(f'{method} Rank')
        ordered_columns.append(f'{method} Significant')

    # Ensure all required columns are present
    pivot_df = pivot_df[ordered_columns]

    # Calculate average rank and percentage of significant pathways for each method
    avg_row = {}
    for method in prop_methods:
        rank_col = f'{method} Rank'
        sig_col = f'{method} Significant'
        
        avg_row[rank_col] = pivot_df[rank_col].mean()  # Average rank
        avg_row[sig_col] = pivot_df[sig_col].mean() * 100  # Percentage of significant pathways (in percentage)

    # Add the dataset and pathway columns as empty for the average row
    avg_row.update({
        'Dataset': 'Average',
        'Pathway': '',
        'Density': '',
        'Num Genes': '',
        'Avg Diameter': ''
    })

    # Append the average row to the DataFrame
    summary_df = pd.concat([pivot_df, pd.DataFrame([avg_row])], ignore_index=True)

    # Sort the DataFrame by 'Dataset' or any relevant column
    summary_df.sort_values(by='Dataset', inplace=True)

    return summary_df


FILE_LIST = get_file_list(DIRECTORIES['input'])
logger.info(f"Found {len(FILE_LIST)} .xlsx files in {DIRECTORIES['input']}")

NETWORKS_DATA, PATHWAYS_DATA = load_data(NETWORKS, PATHWAY_FILES, DIRECTORIES['pathways'])
logger.info("Networks loaded and pathway densities calculated.")

def main():
    start_time = time.time()
    futures = []
    results = []
    
    # Aggregators for methods, pathways, genes, and ranks
    aggregated_pathways = defaultdict(lambda: defaultdict(list))  # method -> dataset -> list of pathways
    aggregated_genes = defaultdict(lambda: defaultdict(list))     # method -> dataset -> list of genes
    aggregated_pathway_ranks = defaultdict(lambda: defaultdict(list))  # method -> dataset -> list of ranks

    # Load the network data (assuming you're using a single network)
    network_name = NETWORKS[0]
    network = NETWORKS_DATA[network_name]

    # Load Disease to Pathway Mappings
    disease_pathway_file = os.path.join(DIRECTORIES['pathways'], 'disease_to_pathway.csv')
    disease_to_pathways = load_disease_to_pathway(disease_pathway_file)
    logger.info("Loaded disease to pathway mappings.")

     # Build search terms for each pathway once
    search_terms_dict = {}
    for file_name in FILE_LIST:
        _, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
        search_terms_dict[pathway_name] = build_search_terms(pathway_name, disease_to_pathways)
    logger.info("Built search terms for all pathways.")

    # Processing the datasets using the ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for pathway_file in PATHWAY_FILES:
            for alpha in ALPHAS:
                for file_name in FILE_LIST:
                    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                    if pathway_name in PATHWAYS_DATA[pathway_file]:
                        for prop_method in PROP_METHODS:
                            futures.append(
                                executor.submit(
                                    process_single_file,
                                    network,
                                    pathway_file,
                                    network_name,
                                    alpha,
                                    prop_method,
                                    file_name,
                                    PATHWAYS_DATA
                                )
                            )

        # Collecting results and aggregating pathways, genes, and ranks per method
        for future in as_completed(futures):
            result = future.result()
            if result is None:
                continue  # Skip if processing failed
            results.append(result['result'])
            
            # Aggregate pathways and genes per method and combined dataset_pathway key
            method = result['result']['Method']
            dataset = result['result']['Dataset']
            pathway = result['result']['Pathway']
            
            # Create a combined key
            combined_key = f"{dataset}_{pathway}"
            
            # Aggregate using the combined key
            aggregated_pathways[method][combined_key].extend(result['higher_ranked_pathways'])
            aggregated_genes[method][combined_key].extend(result['higher_ranked_genes'])
            aggregated_pathway_ranks[method][combined_key].extend(result['higher_ranked_ranks'])

    # Convert the collected results into a DataFrame
    results_df = pd.DataFrame(results)

    # Now, aggregate and save the common pathways and genes per method
    for method in PROP_METHODS:
        # Check if attributes file exists
        attributes_file = os.path.join(DIRECTORIES['summary_base'], f'pathway_attributes_{method}.json')
        if os.path.exists(attributes_file):
            pathway_attributes = load_pathway_attributes_from_file(attributes_file)
            logger.info(f"Loaded pathway attributes from file for method {method}.")
        else:
            pathway_attributes = None  # Will be computed in the function

        # Aggregate pathways, genes, and ranks while excluding disease-associated pathways
        aggregate_higher_ranked_non_disease(
            aggregated_pathways[method],
            aggregated_genes[method],
            aggregated_pathway_ranks[method],
            DIRECTORIES['summary_base'],
            method,
            network,
            search_terms_dict,                    # Pass the entire search_terms_dict
            pathway_attributes,                   # Pass the loaded pathway_attributes if available
            min_appearance=6                      # Set your desired minimum appearance
        )

    # Now, for each method, load the aggregated summaries and perform further analysis
    for method in PROP_METHODS:
        # Load Existing Aggregated Summaries
        common_pathways_path = os.path.join(DIRECTORIES['summary_base'], f'common_pathways_{method}.csv')
        common_genes_path = os.path.join(DIRECTORIES['summary_base'], f'common_genes_{method}.csv')

        if os.path.exists(common_pathways_path) and os.path.exists(common_genes_path):
            pathways_df = pd.read_csv(common_pathways_path)
            logger.info(f"Loaded existing aggregated common pathways and genes for method {method}.")
        else:
            logger.error(f"Aggregated summary files not found for method {method}. Please ensure '{common_pathways_path}' and '{common_genes_path}' exist.")
            continue  # Skip to next method if files not found

        # Perform Disease Pathway Ranking Analysis
        analyze_disease_pathway_rankings(
            search_terms_dict,            # Pass the entire search_terms_dict
            pathways_df,
            DIRECTORIES['summary_base'],
            method
        )

    # Save a single combined summary DataFrame with all methods side by side
    combined_summary_df = generate_summary_df(results_df, PROP_METHODS)
    
    summary_output_dir = os.path.join(DIRECTORIES['summary_base'])
    os.makedirs(summary_output_dir, exist_ok=True)
    combined_summary_output_path = os.path.join(summary_output_dir, 'combined_rankings_summary_all_methods.xlsx')
    combined_summary_df.to_excel(combined_summary_output_path, index=False)
    logger.info(f"Combined rankings summary for all methods saved to {combined_summary_output_path}")
    
    elapsed_time = time.time() - start_time
    logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

    # Remove GSEA output directory and its contents
    gsea_output_dir = os.path.join(PROJECT_ROOT, 'pipeline', 'Outputs', 'GSEA')
    if os.path.exists(gsea_output_dir):
        import shutil
        shutil.rmtree(gsea_output_dir)
        logger.info(f"Removed GSEA output directory: {gsea_output_dir}")
    else:
        logger.info(f"GSEA output directory does not exist: {gsea_output_dir}")

if __name__ == "__main__":
    main()

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import networkx as nx
import logging
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import load_pathways_genes, read_network, read_prior_set
import time
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a logger
logger = logging.getLogger(__name__)

# Define directories
input_dir = os.path.join('Inputs', 'experiments_data', 'GSE', 'XLSX')
output_base_dir = os.path.join('Outputs', 'NGSEA')
plot_output_dir = os.path.join(output_base_dir, 'Plots')
summary_base_dir = os.path.join(output_base_dir, 'Summary')
pathways_dir = os.path.join('Data', 'Human', 'pathways')

# Ensure plot and summary output directories exist
os.makedirs(plot_output_dir, exist_ok=True)
os.makedirs(summary_base_dir, exist_ok=True)


def run_propagation_and_enrichment(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file):
    if method in ['PROP', 'ABS_PROP']:
        # Set alpha before initializing GeneralArgs for PROP and ABS_PROP
        general_args = GeneralArgs(network=network_name, pathway_file=pathway_file, method=method, alpha=alpha, normalization_type='row')
        if method == 'ABS_PROP':
            general_args.input_type = 'abs_Score'
    elif method in ['GSEA', 'NGSEA']:
        # Initialize GeneralArgs for GSEA without modifying alpha
        general_args = GeneralArgs(network=network_name, pathway_file=pathway_file, method=method)
        if method == 'NGSEA':
            general_args.run_NGSEA = True

    perform_propagation(test_name, general_args, network, prior_data)
    perform_enrichment(test_name, general_args, output_path)


# Updated get_pathway_rank function
def get_pathway_rank(gsea_output_path, pathway_name):
    try:
        results_df = pd.read_excel(gsea_output_path)
        pathway_row = results_df[results_df['Term'] == pathway_name]
        if not pathway_row.empty:
            rank = pathway_row.index[0]
            fdr_p_val = pathway_row['FDR q-val'].values[0]
            return rank, fdr_p_val
    except Exception as e:
        logger.error(f"Error reading {gsea_output_path}: {e}")
    return None, None


# Calculate pathway density
# Calculate pathway density, number of genes, and average diameter
def calculate_pathway_density(network, genes):
    # Create a subgraph with only the genes in the pathway
    subgraph = network.subgraph(genes)

    # Calculate the number of connected components
    connected_components = list(nx.connected_components(subgraph))

    # Calculate the percentage of genes that are connected
    num_connected_genes = sum(len(component) for component in connected_components if len(component) > 1)
    percent_connected = num_connected_genes / len(genes) if len(genes) > 0 else 0

    # Calculate the average shortest path length among connected components
    avg_shortest_path_lengths = []
    avg_diameters = []
    for component in connected_components:
        if len(component) > 1:  # Only consider components with more than one node
            component_subgraph = subgraph.subgraph(component)
            try:
                avg_length = nx.average_shortest_path_length(component_subgraph)
                diameter = nx.diameter(component_subgraph)
            except nx.NetworkXError:
                avg_length = np.inf
                diameter = np.inf
            avg_shortest_path_lengths.append(avg_length)
            avg_diameters.append(diameter)

    # Calculate the average distance among connected components
    if avg_shortest_path_lengths:
        avg_distance_among_connected = np.mean(avg_shortest_path_lengths)
    else:
        avg_distance_among_connected = np.inf

    # Calculate the average diameter among connected components
    if avg_diameters:
        avg_diameter_among_connected = np.mean(avg_diameters)
    else:
        avg_diameter_among_connected = np.inf

    # Combine the measures as suggested by your PI
    combined_density_measure = percent_connected * avg_distance_among_connected

    return combined_density_measure if combined_density_measure != 0 else np.inf, len(genes), avg_diameter_among_connected



# The rest of your script would remain largely the same, with the relevant adjustments:
def process_file(network, pathway_file, network_name, alpha, prop_method, file_name, pathway_density, num_genes, avg_diameter):
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    prior_data = read_prior_set(os.path.join(input_dir, file_name))
    output_dir = os.path.join(output_base_dir, prop_method, network_name, pathway_file, f"alpha_{alpha}")
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    run_propagation_and_enrichment(file_name, prior_data, network, network_name, alpha, prop_method, output_file_path, pathway_file)
    prop_rank, fdr_q_val = get_pathway_rank(output_file_path, pathway_name)
    significant = 1 if fdr_q_val is not None and fdr_q_val < 0.05 else 0

    return {
        'Dataset': dataset_name,
        'Pathway': pathway_name,
        'Network': network_name,
        'Pathway file': pathway_file,
        'Alpha': alpha,
        'Method': prop_method,
        'Rank': prop_rank,
        'FDR q-val': fdr_q_val,
        'Significant': significant,
        'Density': pathway_density,
        'Num Genes': num_genes,
        'Avg Diameter': avg_diameter
    }



# Start timing the entire process
start_time = time.time()

networks = ['H_sapiens']
pathway_files = ['kegg']
prop_methods = ['PROP', 'ABS_PROP', 'GSEA', 'NGSEA']
alphas = [0.1, 0.2]

loaded_networks = {}
loaded_pathways = {}
pathway_densities = {}
pathways_to_consider = set(file_name.replace('.xlsx', '').split('_', 1)[1] for file_name in os.listdir(input_dir) if file_name.endswith('.xlsx'))
file_list = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]

logger.info("Loading networks and calculating pathway densities...")

# Load networks and calculate pathway densities once
for network_name in networks:
    network_file = os.path.join('Data', 'Human', 'network', network_name)
    network = read_network(network_file)
    loaded_networks[network_name] = network
    for pathway_file in pathway_files:
        pathways = load_pathways_genes(os.path.join(pathways_dir, pathway_file))
        if pathway_file not in loaded_pathways:
            loaded_pathways[pathway_file] = pathways
        for pathway_name in pathways_to_consider:
            pathway_genes = pathways[pathway_name]
            if network_name not in pathway_densities:
                pathway_densities[network_name] = {}
            if pathway_file not in pathway_densities[network_name]:
                pathway_densities[network_name][pathway_file] = {}
            pathway_densities[network_name][pathway_file][pathway_name] = calculate_pathway_density(network, pathway_genes)
logger.info("Networks loaded and pathway densities calculated.")

# Process files in parallel using ProcessPoolExecutor
futures = []
with ProcessPoolExecutor(max_workers=60) as executor:  # Adjust max_workers based on your CPU capabilities
    for network_name in tqdm(networks, desc='Networks'):
        network = loaded_networks[network_name]
        for pathway_file in tqdm(pathway_files, desc='Pathway Files', leave=False):
            pathways = loaded_pathways[pathway_file]
            for alpha in tqdm(alphas, desc='Alphas', leave=False):
                for file_name in file_list:
                    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                    if pathway_name in pathways:
                        pathway_density, num_genes, avg_diameter = pathway_densities[network_name][pathway_file].get(pathway_name, (np.inf, 0, np.inf))
                        if pathway_density != np.inf:
                            logger.info(f"Parsing {dataset_name} and {pathway_name} with density {pathway_density}, num genes {num_genes}, and avg diameter {avg_diameter}")
                        for prop_method in prop_methods:
                            futures.append(executor.submit(process_file, network, pathway_file, network_name, alpha, prop_method, file_name, pathway_density, num_genes, avg_diameter))

results = []
for future in tqdm(as_completed(futures), total=len(futures), desc='Processing Files'):
    results.append(future.result())

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Generate summary DataFrames for each combination of network_name, pathway_file, and alpha
for network_name in networks:
    for pathway_file in pathway_files:
        for alpha in alphas:
            filtered_df = results_df[(results_df['Network'] == network_name) &
                                     (results_df['Pathway file'] == pathway_file) &
                                     (results_df['Alpha'] == alpha)]

            # Pivot table to get the desired format
            pivot_df = filtered_df.pivot_table(index=['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter'], columns='Method', values=['Rank', 'Significant'], aggfunc='first').reset_index()
            pivot_df.columns = ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter'] + [f'{col[1]} {col[0]}' for col in pivot_df.columns[5:]]

            # Ensure all expected columns are present
            for method in prop_methods:
                if f'{method} Rank' not in pivot_df.columns:
                    pivot_df[f'{method} Rank'] = np.nan
                if f'{method} Significant' not in pivot_df.columns:
                    pivot_df[f'{method} Significant'] = np.nan

            # Reorder columns to the desired order
            column_order = ['Dataset', 'Pathway', 'Density', 'Num Genes', 'Avg Diameter']
            for method in prop_methods:
                column_order.append(f'{method} Rank')
                column_order.append(f'{method} Significant')
            pivot_df = pivot_df[column_order]

            # Calculate average rank and percent significant for each method
            avg_ranks = filtered_df.groupby('Method')['Rank'].mean().reset_index()
            avg_ranks.columns = ['Method', 'Average Rank']
            sig_counts = filtered_df.groupby('Method')['Significant'].sum().reset_index()
            total_counts = filtered_df.groupby('Method')['Significant'].count().reset_index()
            sig_percent = pd.merge(sig_counts, total_counts, on='Method')
            sig_percent['Percentage Significant'] = (sig_percent['Significant_x'] / sig_percent['Significant_y']) * 100
            sig_percent = sig_percent[['Method', 'Percentage Significant']]

            # Create DataFrame for Average Rank and Percent Significant rows
            avg_rank_row = pd.DataFrame([['Average Rank'] + [''] * 4 + [
                avg_ranks[avg_ranks['Method'] == method]['Average Rank'].values[0] if not avg_ranks[
                    avg_ranks['Method'] == method].empty else '' for method in prop_methods for _ in range(2)]],
                                        columns=pivot_df.columns)
            percent_sig_row = pd.DataFrame([['Percent Significant'] + [''] * 4 + [
                sig_percent[sig_percent['Method'] == method]['Percentage Significant'].values[0] if not sig_percent[
                    sig_percent['Method'] == method].empty else '' for method in prop_methods for _ in range(2)]],
                                           columns=pivot_df.columns)

            # Append the summary rows to the pivot DataFrame
            summary_df = pd.concat([pivot_df, avg_rank_row, percent_sig_row], ignore_index=True)
            # sort alphabetically by pathway name
            summary_df = summary_df.sort_values(by='Pathway')
            # Save the summary DataFrame for each network, pathway_file, and alpha
            summary_output_dir = os.path.join(summary_base_dir, network_name, pathway_file, f"alpha {alpha}")
            os.makedirs(summary_output_dir, exist_ok=True)
            rankings_output_path = os.path.join(summary_output_dir,
                                                f'rankings_summary_{network_name}_{pathway_file}_alpha_{alpha}_row.xlsx')
            summary_df.to_excel(rankings_output_path, index=False)
            logger.info(f"Rankings summary saved to {rankings_output_path}")

        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

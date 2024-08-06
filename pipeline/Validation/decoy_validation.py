from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import networkx as nx
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import load_pathways_genes, read_network, read_prior_set
import time
from tqdm import tqdm

# Define directories
input_dir = os.path.join('Inputs', 'experiments_data', 'GSE', 'XLSX')
output_base_dir = os.path.join('Outputs', 'NGSEA')
plot_output_dir = os.path.join(output_base_dir, 'Plots')
summary_base_dir = os.path.join(output_base_dir, 'Summary')
pathways_dir = os.path.join('Data', 'Human', 'pathways')

# Ensure plot and summary output directories exist
os.makedirs(plot_output_dir, exist_ok=True)
os.makedirs(summary_base_dir, exist_ok=True)

# Updated run_propagation_and_enrichment function
def run_propagation_and_enrichment(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file):
    general_args = GeneralArgs(network=network_name, pathway_file=pathway_file, method=method)
    if method == 'NGSEA':
        general_args.run_NGSEA = True
    elif method == 'PROP':
        general_args.set_alpha(alpha)

    elif method == 'ABS_PROP':
        general_args.set_alpha(alpha)
        general_args.input_type = 'abs_Score'

    perform_propagation(test_name, general_args, network, prior_data)
    perform_enrichment(test_name, general_args, output_path)

# Updated get_pathway_rank function
def get_pathway_rank(gsea_output_path, pathway_name):
    results_df = pd.read_excel(gsea_output_path)
    pathway_row = results_df[results_df['Term'] == pathway_name]
    if not pathway_row.empty:
        rank = pathway_row.index[0]
        fdr_p_val = pathway_row['FDR q-val'].values[0]
        return rank, fdr_p_val
    else:
        return None, None

# Calculate pathway density
def calculate_pathway_density(network, genes):
    subgraph = network.subgraph(genes)
    if subgraph.number_of_edges() == 0:
        return np.inf
    try:
        avg_shortest_path_length = nx.average_shortest_path_length(subgraph)
    except nx.NetworkXError:
        avg_shortest_path_length = np.inf
    return avg_shortest_path_length

def process_file(network, pathway_file, network_name, alpha, prop_method, file_name, pathway_density):
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    prior_data = read_prior_set(os.path.join(input_dir, file_name))

    # Ensure the output directory is unique and exists
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
        'Density': pathway_density
    }


# Start timing the entire process
start_time = time.time()

# Hyperparameter loops and main analysis
networks = ['H_sapiens', 'HumanNet']
pathway_files = ['c2', 'kegg']
prop_methods = ['GSEA', 'NGSEA', 'PROP', 'ABS_PROP']
alphas = [0.1, 0.2]

# Dictionary to store loaded networks
loaded_networks = {}
# Dictionary to store loaded pathways
loaded_pathways = {}
# Nested dictionary to store pathway densities
pathway_densities = {}
# Set of pathways to consider based on file_list
pathways_to_consider = set(file_name.replace('.xlsx', '').split('_', 1)[1] for file_name in os.listdir(input_dir) if file_name.endswith('.xlsx'))

file_list = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]

# Load networks and calculate pathway densities once
print("Loading networks and calculating pathway densities...")
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
print("Networks loaded and pathway densities calculated.")
# Process files in parallel
futures = []
with ThreadPoolExecutor(max_workers=32) as executor:  # Adjust max_workers based on your CPU capabilities
    for network_name in tqdm(networks, desc='Networks'):
        network = loaded_networks[network_name]
        for pathway_file in tqdm(pathway_files, desc='Pathway Files', leave=False):
            pathways = loaded_pathways[pathway_file]
            for alpha in tqdm(alphas, desc='Alphas', leave=False):
                for file_name in file_list:
                    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                    if pathway_name in pathways:
                        if pathway_name in pathway_densities[network_name][pathway_file]:
                            pathway_density = pathway_densities[network_name][pathway_file][pathway_name]
                            if pathway_density != np.inf:
                                print(f"Parsing {dataset_name} and {pathway_name} with density {pathway_density}")
                        else:
                            pathway_density = np.inf

                        for prop_method in prop_methods:
                            futures.append(executor.submit(process_file, network, pathway_file, network_name, alpha, prop_method, file_name, pathway_density))

results = []
for future in tqdm(as_completed(futures), total=len(futures), desc='Processing Files'):
    results.append(future.result())

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Pivot table to get the desired format
pivot_df = results_df.pivot_table(index=['Dataset', 'Pathway', 'Density'], columns='Method', values=['Rank', 'Significant'], aggfunc='first').reset_index()
pivot_df.columns = ['Dataset', 'Pathway', 'Density'] + [f'{col[1]} {col[0]}' for col in pivot_df.columns[3:]]

# Ensure all expected columns are present
for method in prop_methods:
    if f'{method} Rank' not in pivot_df.columns:
        pivot_df[f'{method} Rank'] = np.nan
    if f'{method} Significant' not in pivot_df.columns:
        pivot_df[f'{method} Significant'] = np.nan

# Reorder columns to the desired order
column_order = ['Dataset', 'Pathway', 'Density']
for method in prop_methods:
    column_order.append(f'{method} Rank')
    column_order.append(f'{method} Significant')
pivot_df = pivot_df[column_order]

# Calculate average rank and percent significant for each method
avg_ranks = results_df.groupby('Method')['Rank'].mean().reset_index()
avg_ranks.columns = ['Method', 'Average Rank']
sig_counts = results_df.groupby('Method')['Significant'].sum().reset_index()
total_counts = results_df.groupby('Method')['Significant'].count().reset_index()
sig_percent = pd.merge(sig_counts, total_counts, on='Method')
sig_percent['Percentage Significant'] = (sig_percent['Significant_x'] / sig_percent['Significant_y']) * 100
sig_percent = sig_percent[['Method', 'Percentage Significant']]

# Create DataFrame for Average Rank and Percent Significant rows
avg_rank_row = pd.DataFrame([['Average Rank'] + [''] * 2 + [avg_ranks[avg_ranks['Method'] == method]['Average Rank'].values[0] if not avg_ranks[avg_ranks['Method'] == method].empty else '' for method in prop_methods for _ in range(2)]], columns=pivot_df.columns)
percent_sig_row = pd.DataFrame([['Percent Significant'] + [''] * 2 + [sig_percent[sig_percent['Method'] == method]['Percentage Significant'].values[0] if not sig_percent[sig_percent['Method'] == method].empty else '' for method in prop_methods for _ in range(2)]], columns=pivot_df.columns)

# Append the summary rows to the pivot DataFrame
summary_df = pd.concat([pivot_df, avg_rank_row, percent_sig_row], ignore_index=True)

# Save the summary DataFrame for each network, pathway_file, and alpha
for network_name in networks:
    for pathway_file in pathway_files:
        for alpha in alphas:
            summary_output_dir = os.path.join(summary_base_dir, network_name, pathway_file, f"alpha {alpha}")
            os.makedirs(summary_output_dir, exist_ok=True)
            rankings_output_path = os.path.join(summary_output_dir, f'rankings_summary_{network_name}_{pathway_file}_alpha_{alpha}.xlsx')
            summary_df.to_excel(rankings_output_path, index=False)
            print(f"Rankings summary saved to {rankings_output_path}")

# End timing the entire process
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total time taken: {elapsed_time:.2f} seconds")

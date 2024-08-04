import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import json
import networkx as nx
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import load_pathways_genes
from utils import read_network, read_prior_set
import time
from scipy import stats
from sklearn.metrics import precision_recall_curve, auc
from tqdm import tqdm



# Define directories
input_dir = os.path.join('Inputs', 'experiments_data', 'GSE', 'XLSX')
output_base_dir = os.path.join( 'Outputs', 'NGSEA')
plot_output_dir = os.path.join(output_base_dir, 'Plots')
summary_base_dir = os.path.join(output_base_dir, 'Summary')
pathways_dir = os.path.join('Data', 'Human', 'pathways')

# Ensure plot and summary output directories exist
os.makedirs(plot_output_dir, exist_ok=True)
os.makedirs(summary_base_dir, exist_ok=True)

# Updated run_propagation_and_enrichment function
def run_propagation_and_enrichment(test_name, prior_data, network, network_name, alpha, method, output_path, pathway_file):
    # Default arguments
    general_args = GeneralArgs(network=network_name, pathway_file=pathway_file, method=method)

    if method == 'NGSEA':
        general_args.run_NGSEA = True
    elif method == 'PROP':
        general_args.alpha = alpha
    elif method == 'ABS_PROP':
        general_args.alpha = alpha
        general_args.input_type = 'abs_Score'

    # Perform propagation
    perform_propagation(test_name, general_args, network, prior_data)

    # Perform enrichment analysis on propagated data
    perform_enrichment(test_name, general_args, output_path)

# Updated get_pathway_rank function
def get_pathway_rank(gsea_output_path, pathway_name):
    results_df = pd.read_csv(gsea_output_path)
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
        return np.inf  # Return a high value if there are no edges
    try:
        avg_shortest_path_length = nx.average_shortest_path_length(subgraph)
    except nx.NetworkXError:
        # This exception is raised if the graph is not connected
        avg_shortest_path_length = np.inf  # or some other value indicating high distance
    return avg_shortest_path_length

# Start timing the entire process
start_time = time.time()

# Hyperparameter loops and main analysis
networks = ['H_sapiens', 'HumanNet']
pathway_files = ['c2.gmt', 'kegg.gmt']
prop_methods = ['GSEA', 'NGSEA', 'PROP', 'ABS_PROP']
alphas = [0.1, 0.2]

for network_name in tqdm(networks, desc='Networks'):
    # Read the network once per network name
    network_file = os.path.join( 'Data', 'Human', 'network', network_name)
    network = read_network(network_file)

    for pathway_file in tqdm(pathway_files, desc='Pathway Files', leave=False):
        # Load pathways and their genes
        pathways = load_pathways_genes(os.path.join(pathways_dir, pathway_file))
        #drop .gmt from pathway_file
        pathway_file = pathway_file.split('.')[0]
        for alpha in tqdm(alphas, desc='Alphas', leave=False):
            rankings_df = pd.DataFrame(columns=['Dataset', 'Pathway', 'Network', 'Pathway file', 'Alpha', 'Method', 'Rank', 'FDR q-val', 'Significant', 'Density'])

            file_list = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
            for file_name in tqdm(file_list, desc='Files', leave=False):
                # Parse the dataset and pathway from the filename
                dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                prior_data = read_prior_set(os.path.join(input_dir, file_name))

                # Get genes of the pathway
                if pathway_name in pathways:
                    pathway_genes = pathways[pathway_name]
                    # Calculate pathway density once per pathway
                    pathway_density = calculate_pathway_density(network, pathway_genes)
                    if pathway_density != np.inf:
                        print(f"Parsing {dataset_name} and {pathway_name} with density {pathway_density}")
                else:
                    pathway_density = np.inf

                for prop_method in tqdm(prop_methods, desc='Methods', leave=False):
                    output_dir = os.path.join(output_base_dir, prop_method, network_name, pathway_file)
                    os.makedirs(output_dir, exist_ok=True)
                    print(f"Running analysis for {dataset_name} and {pathway_name} with network {network_name}, pathway {pathway_file}, method {prop_method}, alpha {alpha}")

                    # Run the propagation and enrichment analysis
                    run_propagation_and_enrichment(file_name, prior_data, network, network_name, alpha, prop_method, output_dir, pathway_file)

                    # Get pathway ranks and FDR p-values
                    prop_rank, fdr_q_val = get_pathway_rank(output_dir, pathway_name)
                    significant = 1 if fdr_q_val is not None and fdr_q_val < 0.05 else 0

                    # Append the ranks to the DataFrame
                    new_row = pd.DataFrame([{
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
                    }])
                    rankings_df = pd.concat([rankings_df, new_row], ignore_index=True)

            # Calculate average ranking for each method
            avg_rankings = rankings_df.groupby('Method')['Rank'].mean().reset_index()
            avg_rankings.columns = ['Method', 'Average Rank']

            # Calculate the percentage of significant FDR q-val for each method
            sig_counts = rankings_df.groupby('Method')['Significant'].sum().reset_index()
            total_counts = rankings_df.groupby('Method')['Significant'].count().reset_index()
            sig_percent = pd.merge(sig_counts, total_counts, on='Method')
            sig_percent['Percentage Significant'] = (sig_percent['Significant_x'] / sig_percent['Significant_y']) * 100
            sig_percent = sig_percent[['Method', 'Percentage Significant']]

            # Append the average rankings and percentage significant to the DataFrame
            avg_row = pd.DataFrame([{
                'Dataset': 'Average',
                'Pathway': '',
                'Network': '',
                'Pathway file': '',
                'Alpha': '',
                'Method': row['Method'],
                'Rank': row['Average Rank'],
                'FDR q-val': '',
                'Significant': '',
                'Density': '',
                'Percentage Significant': sig_percent[sig_percent['Method'] == row['Method']]['Percentage Significant'].values[0]
            } for index, row in avg_rankings.iterrows()])

            rankings_df = pd.concat([rankings_df, avg_row], ignore_index=True)

            # Save the rankings DataFrame
            summary_output_dir = os.path.join(summary_base_dir, network_name, pathway_file)
            os.makedirs(summary_output_dir, exist_ok=True)
            rankings_output_path = os.path.join(summary_output_dir, f'rankings_summary_{network_name}_{pathway_file}.xlsx')
            rankings_df.to_excel(rankings_output_path, index=False)
            print(f"Rankings summary saved to {rankings_output_path}")

# End timing the entire process
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Total time taken: {elapsed_time:.2f} seconds")

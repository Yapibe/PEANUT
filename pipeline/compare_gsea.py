import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import json
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import read_network, read_prior_set
import time
from scipy import stats
from sklearn.metrics import precision_recall_curve, auc

# Define home and pipeline directories
home_dir = r'C:\Users\pickh\PycharmProjects\Yair_propagation'
pipeline_dir = os.path.join(home_dir, 'pipeline')

# Define directories
input_dir = os.path.join(pipeline_dir, 'Inputs', 'experiments_data', 'GSE', 'XLSX')
output_base_dir = os.path.join(pipeline_dir, 'Outputs', 'NGSEA')
summary_base_dir = os.path.join(output_base_dir, 'Summary')


# Placeholder function for running propagation and enrichment analysis
def run_propagation_and_enrichment(test_name, prior_data, network, alpha, method, output_path):
    if method == 'GSEA':
        general_args = GeneralArgs(network=network)
    elif method == 'NGSEA':
        general_args = GeneralArgs(run_NGSEA=True, network=network)
    elif method == 'PROP':
        general_args = GeneralArgs(alpha=alpha, network=network)
    else:
        general_args = GeneralArgs(alpha=alpha, input_type='abs_Score', network=network)
    # Perform propagation
    perform_propagation(test_name, general_args, network, prior_data)

    # Perform enrichment analysis on propagated data
    perform_enrichment(test_name, general_args, output_path)


# Function to extract pathway ranking from the saved results
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

# Hyperparameter loops and main analysis
networks = ['H_sapiens', 'HumanNet']
pathway_files = ['c2', 'kegg']
prop_methods = ['PROP', 'ABS_PROP', 'GSEA', 'NGSEA']
alphas = [0.1, 0.2]

for network_name in networks:
    for pathway_file in pathway_files:
        network_file = os.path.join(pipeline_dir, 'Data', 'Human', 'network', network_name)
        network = read_network(network_file)

        for alpha in alphas:
            rankings_df = pd.DataFrame(columns=['Dataset', 'Pathway', 'Network', 'Pathway file', 'Alpha', 'Method', 'Rank', 'FDR q-val'])

            for file_name in os.listdir(input_dir):
                if file_name.endswith('.xlsx'):
                    # Parse the dataset and pathway from the filename
                    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
                    prior_data = read_prior_set(os.path.join(input_dir, file_name))

                    for prop_method in prop_methods:
                        output_dir = os.path.join(output_base_dir, prop_method, network_name, pathway_file)
                        os.makedirs(output_dir, exist_ok=True)
                        print(f"Running analysis for {dataset_name} and {pathway_name} with network {network_name}, pathway {pathway_file}, method {prop_method}, alpha {alpha}")

                        # Run the propagation and enrichment analysis
                        run_propagation_and_enrichment(file_name, prior_data, network, alpha, prop_method, output_dir)

                        # Get pathway ranks and FDR p-values
                        prop_rank, fdr_q_val = get_pathway_rank(output_dir, pathway_name)

                        # Append the ranks to the DataFrame
                        new_row = pd.DataFrame([{
                            'Dataset': dataset_name,
                            'Pathway': pathway_name,
                            'Network': network_name,
                            'Pathway file': pathway_file,
                            'Alpha': alpha,
                            'Method': prop_method,
                            'Rank': prop_rank,
                            'FDR q-val': fdr_q_val
                        }])
                        rankings_df = pd.concat([rankings_df, new_row], ignore_index=True)

            # Calculate average ranking for each method
            avg_rankings = rankings_df.groupby('Method')['Rank'].mean().reset_index()
            avg_rankings.columns = ['Method', 'Average Rank']

            # Save the rankings DataFrame and average rankings
            summary_output_dir = os.path.join(summary_base_dir, network_name, pathway_file)
            os.makedirs(summary_output_dir, exist_ok=True)
            rankings_output_path = os.path.join(summary_output_dir, f'rankings_summary_{network_name}_{pathway_file}.xlsx')
            avg_rankings_output_path = os.path.join(summary_output_dir, f'average_rankings_{network_name}_{pathway_file}.xlsx')

            rankings_df.to_excel(rankings_output_path, index=False)
            avg_rankings.to_excel(avg_rankings_output_path, index=False)
            print(f"Rankings summary saved to {rankings_output_path}")
            print(f"Average rankings saved to {avg_rankings_output_path}")
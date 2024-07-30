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
from utils import read_network
from utils import read_prior_set
import time


# Define directories
input_dir = 'Inputs/experiments_data/NGSEA/XLSX'
output_dir_gsea = 'Outputs/NGSEA/GSEA'
output_dir_ngsea = 'Outputs/NGSEA/NGSEA'
output_dir_prop = 'Outputs/NGSEA/PROP'
output_dir_abs_prop = 'Outputs/NGSEA/ABS_PROP'
plot_output_dir = 'Outputs/NGSEA/Plots'
network_file = 'Data/H_sapiens/network/HumanNet'

# Ensure output directories exist
os.makedirs(output_dir_gsea, exist_ok=True)
os.makedirs(output_dir_ngsea, exist_ok=True)
os.makedirs(output_dir_prop, exist_ok=True)
os.makedirs(output_dir_abs_prop, exist_ok=True)
os.makedirs(plot_output_dir, exist_ok=True)

# Function to compare score distributions and save the plot
def compare_score_distributions(scores_gsea, scores_neighbor, scores_propagation, scores_abs_propagation, title, plot_path):
    # Extract the scores from the tuples
    scores_gsea_values = [score[0] for score in scores_gsea.values()]
    scores_neighbor_values = [score[0] for score in scores_neighbor.values()]
    scores_propagation_values = [score[0] for score in scores_propagation.values()]
    scores_abs_propagation_values = [score[0] for score in scores_abs_propagation.values()]

    plt.figure(figsize=(12, 6))
    plt.hist(scores_gsea_values, bins=50, alpha=0.5, label='GSEA')
    plt.hist(scores_neighbor_values, bins=50, alpha=0.5, label='Neighbor Averaging')
    plt.hist(scores_propagation_values, bins=50, alpha=0.5, label='Network Propagation')
    plt.hist(scores_abs_propagation_values, bins=50, alpha=0.5, label='Abs Network Propagation')
    plt.xlabel('Scores')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.savefig(plot_path)
    plt.close()


# Placeholder function for running propagation and enrichment analysis
def run_propagation_and_enrichment(test_name, prior_data, alpha, run_gsea, run_NGSEA, output_path, network,  score='Score'):
    # Setup general arguments
    general_args = GeneralArgs(alpha=alpha, run_gsea=run_gsea, run_NGSEA=run_NGSEA, input_type=score)

    # Perform propagation
    perform_propagation(test_name, general_args, network, prior_data)

    # Perform enrichment analysis on propagated data
    scores = perform_enrichment(test_name, general_args, output_path)

    return scores

# Function to extract pathway ranking from the saved results
def get_pathway_rank(gsea_output_path, pathway_name):
    results_df = pd.read_csv(gsea_output_path)
    # Ensure the pathway name matches the format in your saved results
    pathway_row = results_df[results_df['Term'] == pathway_name]
    if not pathway_row.empty:
        rank = pathway_row.index[0]
        return rank
    else:
        return None

# Initialize a DataFrame to store the rankings
rankings_df = pd.DataFrame(columns=['Dataset', 'Pathway', 'GSEA', 'NGSEA', 'PROP', 'ABS-PROP'])
network = read_network(network_file)
start = time.time()
# Process each file in the input directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('.xlsx'):
        # Parse the dataset and pathway from the filename
        dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
        prior_data = read_prior_set(os.path.join(input_dir, file_name))

        ngsea_output_path = os.path.join(output_dir_ngsea, file_name)
        gsea_output_path = os.path.join(output_dir_gsea, file_name)
        prop_output_path = os.path.join(output_dir_prop, file_name)
        abs_prop_output_path = os.path.join(output_dir_abs_prop, file_name)
        file_name = file_name.replace('.xlsx', '')
        print(f"Running analysis for {dataset_name} and {pathway_name}")
        # Run network propagation GSEA
        prop_scores = run_propagation_and_enrichment(file_name, prior_data, alpha=0.1, run_gsea=True, run_NGSEA=False,
                                                     output_path=prop_output_path, network=network)
        prop_rank = get_pathway_rank(prop_output_path, pathway_name)

        # Run abs network propagation GSEA
        prop_scores_abs = run_propagation_and_enrichment(file_name, prior_data, alpha=0.1, run_gsea=True,
                                                         run_NGSEA=False,
                                                         output_path=abs_prop_output_path, score='abs_Score',
                                                         network=network)
        prop_rank_abs = get_pathway_rank(abs_prop_output_path, pathway_name)

        # Run NGSEA
        ngsea_scores = run_propagation_and_enrichment(file_name, prior_data, alpha=1, run_gsea=True, run_NGSEA=True, output_path=ngsea_output_path, network=network)
        ngsea_rank = get_pathway_rank(ngsea_output_path, pathway_name)

        # Run normal GSEA
        gsea_scores = run_propagation_and_enrichment(file_name,prior_data, alpha=1, run_gsea=True, run_NGSEA=False, output_path=gsea_output_path, network=network, )
        gsea_rank = get_pathway_rank(gsea_output_path, pathway_name)

        # Define plot path
        plot_path = os.path.join(plot_output_dir, f"{dataset_name}_score_distributions.png")

        # Compare score distributions and save the plot
        compare_score_distributions(gsea_scores, ngsea_scores, prop_scores, prop_scores_abs,
                                    title=f"Score Distributions for {dataset_name}", plot_path=plot_path)

        # Append the ranks to the DataFrame
        new_row = pd.DataFrame([{
            'Dataset': dataset_name,
            'Pathway': pathway_name,
            'GSEA': gsea_rank,
            'NGSEA': ngsea_rank,
            'PROP': prop_rank,
            'ABS-PROP': prop_rank_abs
        }])
        rankings_df = pd.concat([rankings_df, new_row], ignore_index=True)
end = time.time()

print(f"Analysis completed in {end - start:.2f} seconds")
# Save the rankings DataFrame
rankings_output_path = 'Outputs/NGSEA/rankings_summary.xlsx'
rankings_df.to_excel(rankings_output_path, index=False)
print(f"Rankings summary saved to {rankings_output_path}")

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
from scipy import stats


network_name = 'HumanNet'
network_file = f'Data/H_sapiens/network/{network_name}'

# Define directories
input_dir = 'Inputs/experiments_data/NGSEA/XLSX'
output_dir_gsea = f'Outputs/NGSEA/GSEA/{network_name}'
output_dir_ngsea = f'Outputs/NGSEA/NGSEA/{network_name}'
output_dir_prop = f'Outputs/NGSEA/PROP/{network_name}'
output_dir_abs_prop = f'Outputs/NGSEA/ABS_PROP/{network_name}'
plot_output_dir = f'Outputs/NGSEA/Plots/{network_name}'



# Ensure output directories exist
os.makedirs(output_dir_gsea, exist_ok=True)
os.makedirs(output_dir_ngsea, exist_ok=True)
os.makedirs(output_dir_prop, exist_ok=True)
os.makedirs(output_dir_abs_prop, exist_ok=True)
os.makedirs(plot_output_dir, exist_ok=True)


def calculate_score_changes(prior_data, post_propagation_data, network):
    """
    Calculate the change in scores from the prior data to the post-propagation scores.

    Parameters:
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.
    - post_propagation_data (dict): Dictionary containing post-propagation gene scores and p-values.
    - network (networkx.Graph): The network graph.

    Returns:
    - pd.DataFrame: DataFrame containing GeneID, prior scores, post-propagation scores, and score changes.
    """
    # Filter genes that are both in the experiment and the network
    common_genes = prior_data[prior_data['GeneID'].isin(network.nodes())]

    # Create a DataFrame from post_propagation_data
    post_propagation_df = pd.DataFrame.from_dict(post_propagation_data, orient='index', columns=['Score_post', 'P_value_post'])
    post_propagation_df.index.name = 'GeneID'
    post_propagation_df.reset_index(inplace=True)

    # Merge prior data and post-propagation data on GeneID
    merged_data = pd.merge(common_genes, post_propagation_df, on='GeneID')

    # Calculate the change in scores
    merged_data['Score_change'] = merged_data['Score_post'] - merged_data['Score']

    return merged_data[['GeneID', 'Score', 'Score_post', 'Score_change']]


def plot_score_changes(score_changes, title, plot_path):
    """
    Plot the distribution of score changes.

    Parameters:
    - score_changes (pd.Series): Series containing the score changes.
    - title (str): Title for the plot.
    - plot_path (str): Path to save the plot.
    """
    plt.figure(figsize=(12, 6))
    plt.hist(score_changes, bins=50, alpha=0.75, color='blue')
    plt.xlabel('Score Change')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.savefig(plot_path)
    plt.close()

def run_validation(prior_data, post_propagation_data, network, plot_path, dataset_name, pathway_name):
    """
    Run the validation by calculating and plotting score changes.

    Parameters:
    - prior_data (pd.DataFrame): DataFrame containing prior gene scores.
    - post_propagation_data (dict): Dictionary containing post-propagation gene scores and p-values.
    - network (networkx.Graph): The network graph.
    - plot_path (str): Path to save the plot.
    """
    # Calculate score changes
    score_changes_df = calculate_score_changes(prior_data, post_propagation_data, network)

    # Plot the distribution of score changes
    plot_title = 'Score Changes Distribution for {} - {}'.format(dataset_name, pathway_name)
    plot_score_changes(score_changes_df['Score_change'], plot_title, plot_path)

    # Perform a statistical Comparison (e.g., t-Comparison) to evaluate the significance of the changes
    t_stat, p_value = stats.ttest_1samp(score_changes_df['Score_change'], 0)
    print(f"T-test results: t-statistic = {t_stat}, p-value = {p_value}")



# Function to compare score distributions and save the plot
def compare_score_distributions_in_dataset(scores_in_exp_not_net, scores_in_net_not_exp, scores_in_both, title, plot_path):
    """
    Compare the distribution of scores for genes in three categories:
    1. In experiment but not in network.
    2. In network but not in experiment.
    3. In both experiment and network.

    Parameters:
    - scores_in_exp_not_net (list): Scores of genes in the experiment but not in the network.
    - scores_in_net_not_exp (list): Scores of genes in the network but not in the experiment.
    - scores_in_both (list): Scores of genes in both the experiment and the network.
    - title (str): Title for the plot.
    - plot_path (str): Path to save the plot.
    """
    # Extract the first element of each tuple in the lists
    scores_in_exp_not_net_values = [score[0] for score in scores_in_exp_not_net]
    scores_in_net_not_exp_values = [score[0] for score in scores_in_net_not_exp]
    scores_in_both_values = [score[0] for score in scores_in_both]

    plt.figure(figsize=(12, 6))
    plt.hist(scores_in_exp_not_net_values, bins=50, alpha=0.5, label='In Experiment, Not in Network')
    plt.hist(scores_in_net_not_exp_values, bins=50, alpha=0.5, label='In Network, Not in Experiment')
    plt.hist(scores_in_both_values, bins=50, alpha=0.5, label='In Both Experiment and Network')
    plt.xlabel('Scores')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.savefig(plot_path)
    plt.close()


def compare_score_distributions_across_dataset(scores_gsea, scores_neighbor, scores_propagation, scores_abs_propagation, title, plot_path):
    """
    Compare the distribution of scores for genes in three categories:
    1. In experiment but not in network.
    2. In network but not in experiment.
    3. In both experiment and network.

    Parameters:
    - scores_in_exp_not_net (list): Scores of genes in the experiment but not in the network.
    - scores_in_net_not_exp (list): Scores of genes in the network but not in the experiment.
    - scores_in_both (list): Scores of genes in both the experiment and the network.
    - title (str): Title for the plot.
    - plot_path (str): Path to save the plot.
    """
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
rankings_df = pd.DataFrame(columns=['Dataset', 'Pathway'])
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
        changes_plot_path = os.path.join(plot_output_dir,
                                         'Changes' f"{dataset_name}_{pathway_name}_score_changes_distribution.png")
        comparison_plot_path = os.path.join(plot_output_dir, 'Comparison', f"{dataset_name}_{pathway_name}_comparison.png")
        distributions_plot_path = os.path.join(plot_output_dir, 'Distribution', f"{dataset_name}_score_distributions.png")


        file_name = file_name.replace('.xlsx', '')
        print(f"Running analysis for {dataset_name} and {pathway_name}")

        # Run network propagation GSEA
        prop_scores = run_propagation_and_enrichment(file_name, prior_data, alpha=0.1, run_gsea=True, run_NGSEA=False,
                                                     output_path=prop_output_path, network=network)
        # # Run the validation
        # run_validation(prior_data, prop_scores, network, changes_plot_path, dataset_name, pathway_name)
        # # Separate the scores into three groups
        # scores_in_exp_not_net = [score for gene, score in prop_scores.items() if gene not in network.nodes()]
        # scores_in_net_not_exp = [prop_scores[gene] for gene in network.nodes() if
        #                          gene not in prior_data['GeneID'].values]
        # scores_in_both = [score for gene, score in prop_scores.items() if
        #                   gene in network.nodes() and gene in prior_data['GeneID'].values]

        # # Compare the distributions
        # plot_title = f"Score Distribution Comparison for {dataset_name} - {pathway_name}"
        # compare_score_distributions_in_dataset(scores_in_exp_not_net, scores_in_net_not_exp, scores_in_both, plot_title, comparison_plot_path)

        prop_rank = get_pathway_rank(prop_output_path, pathway_name)

        # Run abs network propagation GSEA
        prop_scores_abs = run_propagation_and_enrichment(file_name, prior_data, alpha=0.1, run_gsea=True,
                                                         run_NGSEA=False,
                                                         output_path=abs_prop_output_path, score='abs_Score',
                                                         network=network)
        prop_rank_abs = get_pathway_rank(abs_prop_output_path, pathway_name)

        # Run GSE
        ngsea_scores = run_propagation_and_enrichment(file_name, prior_data, alpha=1, run_gsea=True, run_NGSEA=True, output_path=ngsea_output_path, network=network)
        ngsea_rank = get_pathway_rank(ngsea_output_path, pathway_name)

        # Run normal GSEA
        gsea_scores = run_propagation_and_enrichment(file_name,prior_data, alpha=1, run_gsea=True, run_NGSEA=False, output_path=gsea_output_path, network=network, )
        gsea_rank = get_pathway_rank(gsea_output_path, pathway_name)

        # Compare score distributions and save the plot
        compare_score_distributions_across_dataset(gsea_scores, ngsea_scores, prop_scores, prop_scores_abs,
                                    title=f"Score Distributions for {dataset_name}", plot_path=distributions_plot_path)

        # Append the ranks to the DataFrame
        new_row = pd.DataFrame([{
            'Dataset': dataset_name,
            'Pathway': pathway_name,
            'GSEA': gsea_rank,
            'GSE': ngsea_rank,
            'PROP': prop_rank,
            'ABS-PROP': prop_rank_abs
        }])
        rankings_df = pd.concat([rankings_df, new_row], ignore_index=True)

end = time.time()

print(f"Analysis completed in {end - start:.2f} seconds")
# Save the rankings DataFrame
rankings_output_path = f'Outputs/NGSEA/Summary/rankings_summary_{network_name}.xlsx'
rankings_df.to_excel(rankings_output_path, index=False)
print(f"Rankings summary saved to {rankings_output_path}")

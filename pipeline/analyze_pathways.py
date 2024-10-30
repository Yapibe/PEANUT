import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from utils import read_network, load_pathways_genes
import yaml
import networkx as nx
import numpy as np

# Project and Directory Configuration
def get_project_root():
    """Get the project root directory dynamically."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

PROJECT_ROOT = get_project_root()

# Load pre-calculated data from config.yaml
config_path = os.path.join(PROJECT_ROOT, 'pipeline', 'config.yaml')
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

PRE_CALCULATED_DATA = config['pre_calculated_data']

def get_pre_calculated_data(pathway_name):
    """Retrieve pre-calculated data for a given pathway."""
    return PRE_CALCULATED_DATA.get(pathway_name, (None, None, None))

# Convert pre-calculated data into a DataFrame for easier comparison
def create_pre_calculated_df():
    pre_calculated_list = []
    for pathway, data in PRE_CALCULATED_DATA.items():
        pre_calculated_list.append({
            'Pathway': pathway,
            'Density': data[0],  # Map to 'Density'
            'Gene Count': data[1],  # Map to 'Gene Count' (renamed from Num Genes)
            'Avg Diameter': data[2],  # Map to 'Avg Diameter'
            'Category': 'Associated Pathway'  # Label these as associated pathways
        })
    return pd.DataFrame(pre_calculated_list)

# Create pre-calculated DataFrame for associated pathways
pre_calculated_df = create_pre_calculated_df()

# Print to verify
print("Associated Pathways DataFrame:\n", pre_calculated_df.head())

# Load false positive pathways
false_positive_df = pd.read_csv(os.path.join(PROJECT_ROOT, 'pipeline', 'Outputs', 'NGSEA', 'Summary', 'false_positive_pathways_PEANUT.csv'))
print("False Positive Pathways DataFrame:\n", false_positive_df.head())

# Add a Category column to identify false positive pathways
false_positive_df['Category'] = 'False Positive Pathway'

# Ensure both DataFrames have aligned column names
print("Both DataFrames have the Pathway column for merging.")

# Merge pre-calculated associated pathways and false positives into one DataFrame
combined_df = pd.concat([false_positive_df, pre_calculated_df], ignore_index=True)
print("Combined DataFrame after merge:\n", combined_df.head())

# Load the network and pathways data
network = read_network(os.path.join(PROJECT_ROOT, 'pipeline', 'Data', 'H_sapiens', 'network', 'H_sapiens'))
pathways_data = load_pathways_genes(os.path.join(PROJECT_ROOT, 'pipeline', 'Data', 'H_sapiens', 'pathways', 'kegg'))

# Compute advanced network metrics for all pathways
def compute_advanced_network_metrics(network, pathways_df, pathways_data):
    advanced_metrics = []
    for pathway in pathways_df['Pathway']:
        # Get the genes associated with the pathway from pathways_data
        genes = pathways_data.get(pathway, [])
        
        if not genes:
            print(f"Pathway {pathway} not found in pathways_data. Skipping.")
            continue

        # Extract the subgraph for the given pathway genes
        genes_in_network = [gene for gene in genes if gene in network.nodes]  # Only keep genes present in the network
        subgraph = network.subgraph(genes_in_network)

        # Calculate network attributes for this pathway's subgraph
        betweenness = np.mean(list(nx.betweenness_centrality(subgraph).values())) if len(subgraph) > 0 else 0
        clustering = np.mean(list(nx.clustering(subgraph).values())) if len(subgraph) > 0 else 0
        degree_distribution = np.mean([d for n, d in subgraph.degree()]) if len(subgraph) > 0 else 0

        advanced_metrics.append({
            'Pathway': pathway,
            'Betweenness': betweenness,
            'Clustering Coefficient': clustering,
            'Degree Distribution': degree_distribution
        })

    return pd.DataFrame(advanced_metrics)

# Compute advanced metrics for all pathways
advanced_metrics_df = compute_advanced_network_metrics(network, combined_df, pathways_data)

# Merge advanced metrics into the combined DataFrame
combined_df = combined_df.merge(advanced_metrics_df, on='Pathway', how='left')

# Save the combined DataFrame
combined_df.to_csv('combined_pathway_analysis_with_advanced_metrics.csv', index=False)

# Analyze the combined data
def analyze_combined_data(df):
    # Boxplot for Gene Count by Category
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Category', y='Gene Count', data=df)
    plt.title('Gene Count Comparison: False Positives vs Associated Pathways')
    plt.savefig('gene_count_comparison.png', dpi=300)
    plt.show()

    # Boxplot for Density by Category
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Category', y='Density', data=df)
    plt.title('Density Comparison: False Positives vs Associated Pathways')
    plt.savefig('density_comparison.png', dpi=300)
    plt.show()

    # Scatterplot: Betweenness vs Average Rank by Category
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Betweenness', y='Average Rank', hue='Category', data=df)
    plt.title('Betweenness vs. Average Rank for False Positives and Associated Pathways')
    plt.savefig('betweenness_vs_rank_comparison.png', dpi=300)
    plt.show()

    # Scatterplot: Clustering Coefficient vs Average Rank by Category
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Clustering Coefficient', y='Average Rank', hue='Category', data=df)
    plt.title('Clustering Coefficient vs. Average Rank for False Positives and Associated Pathways')
    plt.savefig('clustering_vs_rank_comparison.png', dpi=300)
    plt.show()

    # Scatterplot: Degree Distribution vs Average Rank by Category
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='Degree Distribution', y='Average Rank', hue='Category', data=df)
    plt.title('Degree Distribution vs. Average Rank for False Positives and Associated Pathways')
    plt.savefig('degree_vs_rank_comparison.png', dpi=300)
    plt.show()

# Analyze the combined DataFrame
analyze_combined_data(combined_df)

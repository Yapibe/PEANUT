import matplotlib.pyplot as plt
import networkx as nx
from utils import read_network, load_pathways_genes
import os
from matplotlib.lines import Line2D
import numpy as np

# Hard-coded paths for network and pathways data
NETWORK_FILE = "pipeline/Data/H_sapiens/network/H_sapiens"
PATHWAYS_FILE = "pipeline/Data/H_sapiens/pathways/kegg"
OUTPUT_DIR = "pipeline/Data/H_sapiens/visualizations/"

def filter_connected_components(subgraph, unique_genes_pathway1, unique_genes_pathway2):
    """
    Filters connected components to include only those containing genes from both pathways.

    Args:
        subgraph (nx.Graph): The subgraph of the network containing relevant genes.
        unique_genes_pathway1 (set): Genes unique to the first pathway.
        unique_genes_pathway2 (set): Genes unique to the second pathway.

    Returns:
        nx.Graph: A subgraph containing only the valid connected components.
    """
    components = list(nx.connected_components(subgraph))
    valid_components = [
        component for component in components
        if any(gene in unique_genes_pathway1 for gene in component) and
           any(gene in unique_genes_pathway2 for gene in component)
    ]
    return subgraph.subgraph(set.union(*valid_components)).copy()


def filter_relevant_nodes(subgraph, unique_genes_pathway2):
    nodes_to_keep = set(unique_genes_pathway2)
    for gene in unique_genes_pathway2:
        if gene in subgraph:
            nodes_to_keep.update(subgraph.neighbors(gene))
    return subgraph.subgraph(nodes_to_keep).copy()


def calculate_connectivity(subgraph, shared_genes, unique_genes_pathway1, unique_genes_pathway2):
    """
    Calculates and prints connectivity metrics for the subgraph.

    Args:
        subgraph (nx.Graph): The subgraph of the network containing relevant genes.
        shared_genes (set): Genes shared between the two pathways.
        unique_genes_pathway1 (set): Genes unique to the first pathway.
        unique_genes_pathway2 (set): Genes unique to the second pathway.

    Returns:
        tuple: A tuple containing the number of shared neighbors and edges between pathways.
    """
    shared_neighbors = set()
    for gene in shared_genes:
        if gene in subgraph:
            shared_neighbors.update(subgraph.neighbors(gene))
        else:
            print(f"Warning: Shared gene {gene} is not in the subgraph.")

    pathway_edges = sum(
        1 for node1 in unique_genes_pathway1 if node1 in subgraph
        for node2 in unique_genes_pathway2 if node2 in subgraph
        if subgraph.has_edge(node1, node2)
    )

    density = nx.density(subgraph)
    print(f"Shared neighbors: {len(shared_neighbors)}")
    print(f"Edges between pathways: {pathway_edges}")
    print(f"Density of the subgraph: {density:.4f}")

    return shared_neighbors, pathway_edges

def calculate_functional_metrics(genes_pathway1, genes_pathway2, shared_neighbors, pathway_edges, subgraph):
    isi = pathway_edges / len(shared_neighbors) if len(shared_neighbors) > 0 else 0

    total_possible_edges = len(genes_pathway1) * len(genes_pathway2)
    connectivity_ratio = pathway_edges / total_possible_edges if total_possible_edges > 0 else 0
    union_genes = genes_pathway1.union(genes_pathway2)
    shared_neighbors_ratio = len(shared_neighbors) / len(union_genes)

    total_edges = subgraph.number_of_edges()
    s_score = 1 - (pathway_edges / total_edges) if total_edges > 0 else 1

    print(f"Interaction Strength Index (ISI): {isi:.4f}")
    print(f"Network Separation (S) Score: {s_score:.4f}")
    print(f"Functional Connectivity Score (FCS): {connectivity_ratio + shared_neighbors_ratio:.4f}")

def filter_relevant_nodes(subgraph, unique_genes_pathway2):
    """
    Filters the subgraph to include only nodes related to the second pathway
    and their neighbors.

    Args:
        subgraph (nx.Graph): The subgraph of the network containing relevant genes.
        unique_genes_pathway2 (set): Genes unique to the second pathway.

    Returns:
        nx.Graph: A subgraph containing only the relevant nodes.
    """
    # Identify nodes to keep: those in the second pathway or neighbors of those nodes
    nodes_to_keep = set(unique_genes_pathway2)
    for gene in unique_genes_pathway2:
        if gene in subgraph:
            nodes_to_keep.update(subgraph.neighbors(gene))

    # Create a subgraph with only the relevant nodes
    return subgraph.subgraph(nodes_to_keep).copy()

def calculate_average_shortest_path_length(graph, nodes1, nodes2):
    path_lengths = []
    for node1 in nodes1:
        for node2 in nodes2:
            if node1 != node2 and nx.has_path(graph, node1, node2):
                path_lengths.append(nx.shortest_path_length(graph, node1, node2))
    return np.mean(path_lengths) if path_lengths else float('inf')

def calculate_network_separation_score(subgraph, genes_pathway1, genes_pathway2):
    d_AB = calculate_average_shortest_path_length(subgraph, genes_pathway1, genes_pathway2)
    d_AA = calculate_average_shortest_path_length(subgraph, genes_pathway1, genes_pathway1)
    d_BB = calculate_average_shortest_path_length(subgraph, genes_pathway2, genes_pathway2)
    
    s_score = d_AB - (d_AA + d_BB) / 2
    return s_score

def plot_subgraph(subgraph, pos, genes_pathway1, genes_pathway2, shared_genes, pathway_name1, pathway_name2, title, ax):
    """
    Plots the given subgraph on the specified axes.

    Args:
        subgraph (nx.Graph): The subgraph to plot.
        pos (dict): Layout positions for the nodes.
        genes_pathway1 (set): Genes in pathway 1.
        genes_pathway2 (set): Genes in pathway 2.
        shared_genes (set): Genes shared between both pathways.
        pathway_name1 (str): Name of pathway 1.
        pathway_name2 (str): Name of pathway 2.
        title (str): Title for the subplot.
        ax (matplotlib.axes.Axes): The axes on which to draw the plot.
    """
    # Define node colors
    node_colors = []
    for node in subgraph.nodes:
        if node in shared_genes:
            node_colors.append('purple')  # Shared genes
        elif node in genes_pathway1:
            node_colors.append('red')     # Pathway 1 genes
        elif node in genes_pathway2:
            node_colors.append('blue')    # Pathway 2 genes
        else:
            node_colors.append('gray')    # Fallback color

    # Define node sizes
    node_sizes = [30 if node in shared_genes else 20 for node in subgraph.nodes]

    # Draw edges
    nx.draw_networkx_edges(subgraph, pos, alpha=0.5, ax=ax)

    # Draw nodes
    nx.draw_networkx_nodes(subgraph, pos, node_color=node_colors, node_size=node_sizes, ax=ax)

    # Draw labels
    nx.draw_networkx_labels(subgraph, pos, font_size=8, ax=ax)

    # Set title and remove axis
    ax.set_title(title)
    ax.axis('off')

    # Add legend only to the first subplot
    if title == "Original Subgraph":
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=f'{pathway_name1} Genes',
                   markerfacecolor='red', markersize=10),
            Line2D([0], [0], marker='o', color='w', label=f'{pathway_name2} Genes',
                   markerfacecolor='blue', markersize=10),
            Line2D([0], [0], marker='o', color='w', label='Shared Genes',
                   markerfacecolor='purple', markersize=10)
        ]
        ax.legend(handles=legend_elements)

def visualize_two_pathways(pathway_name1: str, pathway_name2: str, output_file: str = None) -> None:
    """
    Visualizes two specified pathways within the network, highlighting their genes.
    Shared genes are distinctly colored, and genes with no connections are excluded.

    Args:
        pathway_name1 (str): The name of the first pathway to visualize.
        pathway_name2 (str): The name of the second pathway to visualize.
        output_file (str, optional): File path to save the visualization. If None, displays the plot.

    Raises:
        ValueError: If either of the specified pathways is not found in the pathways data.
    """
    # Load the network graph
    network = read_network(NETWORK_FILE)

    # Load pathways and their genes
    pathways = load_pathways_genes(PATHWAYS_FILE)

    # Validate pathway names
    missing_pathways = [p for p in [pathway_name1, pathway_name2] if p not in pathways]
    if missing_pathways:
        raise ValueError(f"Pathway(s) '{', '.join(missing_pathways)}' not found in the pathways data.")

    # Extract genes for both pathways
    genes_pathway1 = set(pathways[pathway_name1])
    genes_pathway2 = set(pathways[pathway_name2])

    # Identify shared and unique genes
    shared_genes = genes_pathway1.intersection(genes_pathway2)
    unique_genes_pathway1 = genes_pathway1 - shared_genes
    unique_genes_pathway2 = genes_pathway2 - shared_genes

    # Create subgraph with genes present in the network
    genes_in_network = (genes_pathway1 | genes_pathway2) & set(network.nodes)
    subgraph = network.subgraph(genes_in_network).copy()

    # Filter to include only connected components with both pathways
    filtered_subgraph = filter_connected_components(subgraph, unique_genes_pathway1, unique_genes_pathway2)

    # Create a figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    plt.tight_layout()

    # Calculate layout positions once
    pos = nx.spring_layout(subgraph, seed=42, scale=2)

    # Plot original subgraph
    plot_subgraph(
        subgraph, pos, genes_pathway1, genes_pathway2, shared_genes,
        pathway_name1, pathway_name2, "Original Subgraph",
        ax=axes[0]
    )

    # Plot filtered subgraph
    plot_subgraph(
        filtered_subgraph, pos, genes_pathway1, genes_pathway2, shared_genes,
        pathway_name1, pathway_name2, "Filtered Subgraph",
        ax=axes[1]
    )

    # Adjust layout and save or show the figure
    plt.suptitle(f"Comparison of {pathway_name1} and {pathway_name2}")
    plt.subplots_adjust(top=0.92)  # Adjust to make room for the suptitle
    if output_file:
        plt.savefig(output_file, format='PNG', bbox_inches='tight')
        plt.close()
        print(f"Visualization saved to {output_file}\n")
    else:
        plt.show()
        plt.close()

    # Calculate connectivity metrics
    shared_neighbors, pathway_edges = calculate_connectivity(filtered_subgraph, shared_genes, unique_genes_pathway1, unique_genes_pathway2)

    # Calculate similarity metrics
    calculate_functional_metrics(genes_pathway1, genes_pathway2, shared_neighbors, pathway_edges, filtered_subgraph)

    # Calculate Network Separation Score
    s_score = calculate_network_separation_score(filtered_subgraph, genes_pathway1, genes_pathway2)
    print(f"Network Separation (S) Score: {s_score:.4f}")

def main():
    # Start of Selection
    pathway_name1 = "KEGG_ACUTE_MYELOID_LEUKEMIA"  # Replace with your first pathway name
    pathway_names2 = [
        "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS",
    ]

    # Ensure the output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for pathway_name2 in pathway_names2:
        output_file = f"{OUTPUT_DIR}{pathway_name1}_{pathway_name2}.png"
        try:
            visualize_two_pathways(pathway_name1, pathway_name2, output_file)
        except ValueError as ve:
            print(f"Error: {ve}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
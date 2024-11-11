import os
import json
import pandas as pd
from utils import load_pathways_genes, load_disease_to_pathway, read_network
from statistical_methods import jaccard_index
import pickle
import numpy as np
import networkx as nx

def get_project_root():
    """Get the project root directory dynamically."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

PROJECT_ROOT = get_project_root()

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
    if valid_components:
        return subgraph.subgraph(set.union(*valid_components)).copy()
    else:
        return nx.Graph()  # Return an empty graph if no valid components are found

def filter_relevant_nodes(subgraph, unique_genes_pathway2):
    """
    Filters the subgraph to include only nodes related to the second pathway and their neighbors.

    Args:
        subgraph (nx.Graph): The subgraph of the network containing relevant genes.
        unique_genes_pathway2 (set): Genes unique to the second pathway.

    Returns:
        nx.Graph: A subgraph containing only the relevant nodes.
    """
    nodes_to_keep = set(unique_genes_pathway2)
    for gene in unique_genes_pathway2:
        if gene in subgraph:
            nodes_to_keep.update(subgraph.neighbors(gene))
    return subgraph.subgraph(nodes_to_keep).copy()


def define_related_pathways(genes_by_pathway, isi_values, threshold, known_related_pathways):
    """
    Define related pathways for each pathway based on the given ISI threshold.

    Args:
        genes_by_pathway (dict): All pathways and their genes.
        isi_values (dict): ISI values between pathways.
        threshold (float): ISI threshold.
        known_related_pathways (dict): Pathways of interest and their known related pathways.

    Returns:
        dict: Mapping of pathways to their related pathways.
    """
    related_pathways = {}

    for pathway in genes_by_pathway:
        related = set()

        # Include known related pathways if available
        if pathway in known_related_pathways:
            related.update(known_related_pathways[pathway])

        # Add pathways that exceed the ISI threshold
        for other_pathway in genes_by_pathway:
            if other_pathway != pathway and isi_values.get(pathway, {}).get(other_pathway, 0) >= threshold:
                related.add(other_pathway)

        # Store the related pathways for the current pathway
        related_pathways[pathway] = list(related)

    return related_pathways


def compute_isi(genes1, genes2, network):
    """
    Compute the Interaction Strength Index (ISI) between two pathways,
    incorporating filtering stages to focus on important subgraphs.

    Args:
        genes1 (set): Genes in pathway 1.
        genes2 (set): Genes in pathway 2.
        network (nx.Graph): The network graph.

    Returns:
        float: ISI value between pathway 1 and pathway 2.
    """
    # Identify shared and unique genes
    shared_genes = genes1 & genes2
    unique_genes1 = genes1 - shared_genes
    unique_genes2 = genes2 - shared_genes

    # Create subgraph with genes present in the network
    genes_in_network = (genes1 | genes2) & set(network.nodes)
    subgraph = network.subgraph(genes_in_network).copy()

    # **Filtering Step 1**: Filter connected components to include only those with genes from both pathways
    subgraph = filter_connected_components(subgraph, unique_genes1, unique_genes2)

    # **Filtering Step 2 (Optional)**: Filter relevant nodes (e.g., nodes related to unique_genes2 and their neighbors)
    subgraph = filter_relevant_nodes(subgraph, unique_genes2)

    # Recompute unique genes after filtering
    unique_genes1 = set(node for node in unique_genes1 if node in subgraph)
    unique_genes2 = set(node for node in unique_genes2 if node in subgraph)
    shared_genes = set(node for node in shared_genes if node in subgraph)

    # Compute shared_neighbors: neighbors of shared_genes within the filtered subgraph
    shared_neighbors = set()
    for gene in shared_genes:
        shared_neighbors.update(subgraph.neighbors(gene))

    # Compute pathway_edges: edges between unique_genes1 and unique_genes2 within the filtered subgraph
    pathway_edges = subgraph.subgraph(unique_genes1 | unique_genes2).number_of_edges()

    # Compute ISI
    if len(shared_neighbors) > 0:
        isi = pathway_edges / len(shared_neighbors)
    else:
        isi = 0

    return isi



def compute_pairwise_isi(genes_by_pathway, disease_to_pathways, network):
    """
    Compute the pairwise ISI between pathways in disease_to_pathways and all pathways.

    Args:
        genes_by_pathway (dict): All pathways and their genes.
        disease_to_pathways (dict): Pathways of interest and their related pathways.
        network (nx.Graph): The network graph.

    Returns:
        dict: Nested dictionary with ISI values.
    """
    isi_values = {p1: {} for p1 in disease_to_pathways.keys()}

    for p1 in disease_to_pathways:
        genes1 = set(genes_by_pathway[p1])
        for p2 in genes_by_pathway:
            if p1 != p2:
                genes2 = set(genes_by_pathway[p2])
                isi = compute_isi(genes1, genes2, network)
                isi_values[p1][p2] = isi

    return isi_values


def find_best_isi_threshold(genes_by_pathway, disease_to_pathways, isi_values):
    """
    Find the best ISI threshold by evaluating against known related pathways.

    Args:
        genes_by_pathway (dict): All pathways and their genes.
        disease_to_pathways (dict): Pathways of interest and their related pathways.
        isi_values (dict): Precomputed ISI values.

    Returns:
        float: The best ISI threshold.
    """
    # Collect all ISI values
    all_isi_values = []
    for p1 in isi_values:
        for p2 in isi_values[p1]:
            all_isi_values.append(isi_values[p1][p2])

    # Remove zero ISI values
    all_isi_values = [v for v in all_isi_values if v > 0]

    if not all_isi_values:
        print("All ISI values are zero. Cannot find a meaningful threshold.")
        return None, []

    # Define thresholds based on percentiles
    percentiles = np.linspace(0, 100, num=100)
    thresholds = np.percentile(all_isi_values, percentiles)

    best_threshold = None
    best_f1 = -1
    performance_metrics = []

    highest_threshold_with_tp = 0
    for t in thresholds:
        total_tp = total_fp = total_fn = 0
        for pathway, related_pathways in disease_to_pathways.items():
            predicted_related = set()
            for other_pathway in genes_by_pathway:
                if other_pathway != pathway and isi_values[pathway].get(other_pathway, 0) >= t:
                    predicted_related.add(other_pathway)

            actual_related = set(related_pathways)
            tp = len(predicted_related & actual_related)
            fp = len(predicted_related - actual_related)
            fn = len(actual_related - predicted_related)

            total_tp += tp
            total_fp += fp
            total_fn += fn

        precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
        recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

        performance_metrics.append({'threshold': t, 'precision': precision, 'recall': recall, 'f1': f1})

        if total_tp > 0 and t > highest_threshold_with_tp:
            highest_threshold_with_tp = t  # Update the highest threshold with at least one TP

        if f1 > best_f1:
            best_f1 = f1
            best_threshold = t

    # Print metrics for the highest threshold with at least one TP
    for metric in performance_metrics:
        if metric['threshold'] == highest_threshold_with_tp:
            print(f"Threshold: {metric['threshold']}, Precision: {metric['precision']}, Recall: {metric['recall']}, F1: {metric['f1']}")

    return best_threshold, performance_metrics



def main():
    # Paths to your data
    pathways_dir = os.path.join(PROJECT_ROOT, 'pipeline/Data/H_sapiens/pathways/kegg')
    related_pathways_file = os.path.join(PROJECT_ROOT, 'pipeline/Data/H_sapiens/pathways/disease_to_pathway.csv')
    network_file = os.path.join(PROJECT_ROOT, 'pipeline/Data/H_sapiens/network/H_sapiens')
    isi_values_file = 'isi_values.pkl'  # File to save/load ISI values

    # Load the network graph
    network = read_network(network_file)

    # Load all pathways and their genes
    genes_by_pathway = load_pathways_genes(pathways_dir)
    # Filter pathways to have between 15 and 500 genes
    genes_by_pathway = {k: v for k, v in genes_by_pathway.items() if 15 <= len(v) <= 500}

    # Load known related pathways
    disease_to_pathways = load_disease_to_pathway(related_pathways_file)

    # Check if ISI values file exists
    if os.path.exists(isi_values_file):
        print(f"Loading ISI values from {isi_values_file}")
        with open(isi_values_file, 'rb') as f:
            isi_values = pickle.load(f)
    else:
        print("Computing ISI values...")
        # Compute ISI values with filtering steps
        isi_values = compute_pairwise_isi(genes_by_pathway, disease_to_pathways, network)
        # Save ISI values to file
        with open(isi_values_file, 'wb') as f:
            pickle.dump(isi_values, f)
        print(f"ISI values saved to {isi_values_file}")

    # Find the best ISI threshold
    best_threshold, performance_metrics = find_best_isi_threshold(
        genes_by_pathway, disease_to_pathways, isi_values)
    print(f"Best ISI threshold: {best_threshold}")

    # Print performance metrics
    for metric in performance_metrics:
        print(metric)

    # Define related pathways using the best threshold
    related_pathways = define_related_pathways(genes_by_pathway, isi_values, best_threshold, disease_to_pathways)

    # Save the related pathways to a file
    output_file = 'related_pathways_isi_filtered.json'
    with open(output_file, 'w') as f:
        json.dump(related_pathways, f, indent=4)
    print(f"Related pathways saved to {output_file}")

if __name__ == "__main__":
    main()
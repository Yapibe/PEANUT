import os
import json
import pandas as pd
from utils import load_pathways_genes, load_disease_to_pathway
from statistical_methods import jaccard_index

def get_project_root():
    """Get the project root directory dynamically."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

PROJECT_ROOT = get_project_root()

def compute_pairwise_jaccard_indices(genes_by_pathway, disease_to_pathways):
    """
    Compute the pairwise Jaccard indices between the specific pathways in disease_to_pathways
    and all pathways in genes_by_pathway.

    Args:
        genes_by_pathway (dict): Mapping of all pathways to their constituent genes.
        disease_to_pathways (dict): Mapping of the 12 specific pathways to their related pathways.

    Returns:
        dict: Nested dictionary where jaccard_indices[p1][p2] = Jaccard index between p1 and p2.
    """
    jaccard_indices = {p1: {} for p1 in disease_to_pathways.keys()}

    # Compute Jaccard index for each pathway in disease_to_pathways against all pathways in genes_by_pathway
    for p1 in disease_to_pathways:
        genes1 = set(genes_by_pathway[p1])
        for p2 in genes_by_pathway:
            if p1 != p2:  # Skip self-comparison
                genes2 = set(genes_by_pathway[p2])
                ji = jaccard_index(genes1, genes2)
                jaccard_indices[p1][p2] = ji

    return jaccard_indices

def find_best_jaccard_threshold(genes_by_pathway, disease_to_pathways, jaccard_indices):
    """
    Find the best Jaccard threshold by evaluating against known related pathways.

    Args:
        genes_by_pathway (dict): Mapping of all pathways to their genes.
        disease_to_pathways (dict): Mapping of the 12 specific pathways to their related pathways.
        jaccard_indices (dict): Precomputed pairwise Jaccard indices for the 12 specific pathways.

    Returns:
        float: The best Jaccard threshold.
    """
    thresholds = [i / 100 for i in range(1, 100)]
    best_threshold = None
    best_f1 = -1
    performance_metrics = []

    highest_threshold_with_tp = 0
    for t in thresholds:
        total_tp = total_fp = total_fn = 0
        for pathway, related_pathways in disease_to_pathways.items():
            predicted_related = set()
            for other_pathway in genes_by_pathway:
                if other_pathway != pathway and jaccard_indices[pathway].get(other_pathway, 0) >= t:
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

        if total_tp > 0:
            highest_threshold_with_tp = t  # Update the highest threshold that gives at least one TP

        if f1 > best_f1:
            best_f1 = f1
            best_threshold = t

    # Find and print the F1 values for the highest threshold with at least one TP
    for metric in performance_metrics:
        if metric['threshold'] == highest_threshold_with_tp:
            print(f"Threshold: {metric['threshold']}, Precision: {metric['precision']}, Recall: {metric['recall']}, F1: {metric['f1']}")

    return best_threshold, performance_metrics


def define_related_pathways(genes_by_pathway, jaccard_indices, threshold, known_related_pathways):
    """
    Define related pathways for each pathway based on the given Jaccard threshold.

    Args:
        genes_by_pathway (dict): Mapping of pathways to their genes.
        jaccard_indices (dict): Pairwise Jaccard indices between pathways.
        threshold (float): Jaccard threshold to define related pathways.
        known_related_pathways (dict): Mapping of pathways to their known related pathways.

    Returns:
        dict: Mapping of pathways to their related pathways.
    """
    related_pathways = {}

    for pathway in genes_by_pathway:
        related = set()

        # Include known related pathways if available
        if pathway in known_related_pathways:
            related.update(known_related_pathways[pathway])

        # Add pathways that exceed the Jaccard threshold
        for other_pathway in genes_by_pathway:
            if other_pathway != pathway and jaccard_indices.get(pathway, {}).get(other_pathway, 0) >= threshold:
                related.add(other_pathway)

        # Store the related pathways for the current pathway
        related_pathways[pathway] = list(related)

    return related_pathways


def main():
    # Paths to your data
    pathways_dir = os.path.join(PROJECT_ROOT, 'pipeline/Data/H_sapiens/pathways/kegg')
    related_pathways_file = os.path.join(PROJECT_ROOT, 'pipeline/Data/H_sapiens/pathways/disease_to_pathway.csv')

    # Load all pathways and their genes
    genes_by_pathway = load_pathways_genes(pathways_dir)
    # Filter pathways to be more than 14 and less than 501
    genes_by_pathway = {k: v for k, v in genes_by_pathway.items() if 15 <= len(v) <= 500}

    # Load known related pathways from the CSV file
    disease_to_pathways = load_disease_to_pathway(related_pathways_file)

    # Compute pairwise Jaccard indices only for the keys in disease_to_pathways
    jaccard_indices = compute_pairwise_jaccard_indices(genes_by_pathway, disease_to_pathways)

    # Find the best Jaccard threshold and print F1 values for the relevant threshold
    best_threshold, performance_metrics = find_best_jaccard_threshold(
        genes_by_pathway, disease_to_pathways, jaccard_indices)
    print(f"Best Jaccard threshold: {best_threshold}")

    # Optionally, print all performance metrics
    for metric in performance_metrics:
        print(metric)

    # Define related pathways using the best threshold
    related_pathways = define_related_pathways(genes_by_pathway, jaccard_indices, best_threshold, disease_to_pathways)

    # Save the related pathways to a file
    output_file = 'related_pathways.json'
    with open(output_file, 'w') as f:
        json.dump(related_pathways, f, indent=4)
    print(f"Related pathways saved to {output_file}")


if __name__ == "__main__":
    main()

import numpy as np
import pandas as pd
from scipy.stats import norm, hypergeom, ks_2samp, rankdata


def compute_mw(experiment_ranks, control_ranks):
    """
    Compute the Mann-Whitney U test manually using rank sums to determine the statistical difference
    between two independent samples.

    Parameters:
    - experiment_ranks (list): Ranks of the experimental group.
    - control_ranks (list): Ranks of the control group.
    Returns:
    tuple: The Mann-Whitney U statistic and the corresponding p-value.
    """

    # Calculate the sum of ranks for each group
    R1 = np.sum(experiment_ranks)
    R2 = np.sum(control_ranks)

    # Number of observations in each group
    n1 = len(experiment_ranks)
    n2 = len(control_ranks)

    # Compute the Mann-Whitney U statistics for both groups
    U1 = R1 - n1 * (n1 + 1) / 2  # U statistic for the experimental group
    U2 = R2 - n2 * (n2 + 1) / 2  # U statistic for the control group

    # Use the smaller U statistic as the test statistic
    U = min(U1, U2)

    # Calculate the mean and standard deviation for the U distribution under H0 (null hypothesis)
    mean_U = n1 * n2 / 2
    std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)

    # Calculate the Z-score associated with U statistic
    Z = (U - mean_U) / std_U

    # Calculate the two-tailed p-value from the Z-score
    p_value = 2 * norm.cdf(-np.abs(Z))  # Two-tailed test

    return U, p_value


def jaccard_index(set1: set, set2: set) -> float:
    """
    Calculate the Jaccard index, a measure of similarity between two sets.

    Parameters:
    - set1 (set): First set of elements.
    - set2 (set): Second set of elements.

    Returns:
    - float: Jaccard index (intersection over union).
    """
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    if not union:
        return 0.0
    return len(intersection) / len(union)


def compute_pairwise_jaccard_indices(genes_by_pathway):
    """
    Compute the pairwise Jaccard indices between all pathways.

    Args:
        genes_by_pathway (dict): Mapping of pathways to their constituent genes.

    Returns:
        dict: Nested dictionary where jaccard_indices[p1][p2] = Jaccard index between p1 and p2.
    """
    pathways = list(genes_by_pathway.keys())
    jaccard_indices = {p1: {} for p1 in pathways}

    for i, p1 in enumerate(pathways):
        genes1 = set(genes_by_pathway[p1])
        for j in range(i + 1, len(pathways)):
            p2 = pathways[j]
            genes2 = set(genes_by_pathway[p2])
            ji = jaccard_index(genes1, genes2)
            jaccard_indices[p1][p2] = ji
            jaccard_indices[p2][p1] = ji  # Symmetric

    return jaccard_indices


def find_best_jaccard_threshold(genes_by_pathway, known_related_pathways):
    """
    Find the best Jaccard threshold by evaluating against known related pathways.

    Args:
        genes_by_pathway (dict): Mapping of pathways to their genes.
        known_related_pathways (dict): Mapping of pathways to their known related pathways.

    Returns:
        float: The best Jaccard threshold.
    """
    # Compute pairwise Jaccard indices
    jaccard_indices = compute_pairwise_jaccard_indices(genes_by_pathway)

    # Define thresholds to test
    thresholds = [i / 100 for i in range(1, 100)]  # 0.01 to 0.99

    best_threshold = None
    best_f1 = -1  # Initialize with a value lower than possible F1 score
    performance_metrics = []

    for t in thresholds:
        tp = fp = fn = 0
        for pathway in known_related_pathways:
            predicted_related = set()
            for other_pathway in genes_by_pathway:
                if other_pathway != pathway and jaccard_indices[pathway].get(other_pathway, 0) >= t:
                    predicted_related.add(other_pathway)

            actual_related = set(known_related_pathways[pathway])
            tp += len(predicted_related & actual_related)
            fp += len(predicted_related - actual_related)
            fn += len(actual_related - predicted_related)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

        performance_metrics.append({'threshold': t, 'precision': precision, 'recall': recall, 'f1': f1})

        if f1 > best_f1:
            best_f1 = f1
            best_threshold = t

    return best_threshold, performance_metrics



def run_hyper(genes_by_pathway: dict, scores_keys: set, significant_p_vals: dict) -> list:
    """
    Run the hypergeometric test to identify pathways with significant enrichment.

    Parameters:
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores_keys (set): Set of gene IDs with scores.
    - significant_p_vals (dict): Mapping of gene IDs to significant P-values.

    Returns:
    - list: List of significant pathways identified by the hypergeometric test.
    """
    # Total number of scored genes
    M = len(scores_keys)
    # Number of genes with significant P-values
    n = len(significant_p_vals)

    # Prepare lists to hold the hypergeometric P-values and corresponding pathway names
    hypergeom_p_values = []
    pathway_names = []

    # Calculate hypergeometric P-values for each pathway with enough genes
    for pathway_name, pathway_genes in genes_by_pathway.items():
        N = len(pathway_genes)  # Number of genes in the current pathway
        x = len(set(pathway_genes).intersection(significant_p_vals.keys()))  # Enriched genes in the pathway
        # Apply hypergeometric test; if fewer than 5 enriched genes, assign a P-value of 1 (non-significant)
        pval = hypergeometric_sf(x, M, N, n) if x >= 5 else 1
        hypergeom_p_values.append(pval)
        pathway_names.append(pathway_name)

    # Identify pathways with significant hypergeometric P-values
    significant_pathways = [
        pathway for i, pathway in enumerate(pathway_names) if hypergeom_p_values[i] < 0.05
    ]

    return significant_pathways


def hypergeometric_sf(x: int, M: int, N: int, n:int) -> float:
    """
    Calculate the survival function (complement of the CDF) for the hypergeometric distribution.

    Parameters:
    - x (int): The number of successful draws.
    - M (int): The total size of the population.
    - N (int): The number of success states in the population.
    - n (int): The number of draws.

    Returns:
    float: The survival function probability.
    """
    # Calculate the survival function using hypergeometric distribution
    probability = hypergeom.sf(x - 1, M, N, n)
    return probability


def global_gene_ranking(scores: dict):
    """
    Rank all genes globally based on their scores.

    Parameters:
    - scores (dict): Mapping of gene IDs to their scores.

    Returns:
    - pd.Series: A series with gene IDs as index and their global ranks as values.
    
    Raises:
    - ValueError: If scores dictionary is empty or contains invalid values
    """
    if not scores:
        raise ValueError("Scores dictionary is empty")

    try:
        # Extract scores (no longer tuples, just direct values)
        gene_ids = list(scores.keys())
        gene_scores = list(scores.values())

        # Validate scores
        if any(not isinstance(score, (int, float)) for score in gene_scores):
            raise ValueError("Invalid score values detected. All scores must be numeric.")

        # Create a DataFrame for easier handling
        df = pd.DataFrame({'GeneID': gene_ids, 'Score': gene_scores})

        # Rank the scores (higher scores get lower rank numbers)
        df['Rank'] = df['Score'].rank(ascending=False, method='average')

        # Create a Series with GeneID as index and Rank as values
        global_ranking = df.set_index('GeneID')['Rank']

        return global_ranking

    except Exception as e:
        raise ValueError(f"Error in global gene ranking: {str(e)}")




def kolmogorov_smirnov_test_with_ranking(pathway_genes, global_ranking, alternative='greater'):
    # Ensure pathway_genes is a list
    pathway_genes = list(pathway_genes)

    # Ensure all genes are in the global ranking
    pathway_genes_in_ranking = [gene for gene in pathway_genes if gene in global_ranking.index]

    if not pathway_genes_in_ranking:
        return 1.0

    # Get the ranks for pathway genes
    pathway_ranks = global_ranking[pathway_genes_in_ranking].values

    # Get the ranks for background genes
    background_genes = global_ranking.index.difference(pathway_genes_in_ranking)
    background_ranks = global_ranking[background_genes].values

    # Perform the one-sided KS test
    D, p_value = ks_2samp(
        pathway_ranks,
        background_ranks,
        alternative=alternative
    )

    return p_value



def kolmogorov_smirnov_test_with_scores(pathway_genes, scores, alternative='less'):
    """
    Perform the one-sided Kolmogorov-Smirnov test using gene scores.

    Parameters:
    - pathway_genes (iterable): Iterable of gene IDs in the pathway.
    - scores (dict): Mapping of gene IDs to their scores.
    - alternative (str): Defines the alternative hypothesis ('less', 'greater', 'two-sided').

    Returns:
    - float: The p-value from the KS test indicating statistical difference.
    """
    # Get scores for pathway genes
    pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes if gene_id in scores]

    # Get scores for background genes
    background_genes = set(scores.keys()) - set(pathway_genes)
    background_scores = [scores[gene_id][0] for gene_id in background_genes]

    # Check for empty lists
    if not pathway_scores or not background_scores:
        return 1.0  # Return a neutral p-value if no data

    # Perform the one-sided KS test
    D, p_value = ks_2samp(
        pathway_scores,
        background_scores,
        alternative=alternative
    )

    return p_value


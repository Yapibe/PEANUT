import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, hypergeom, ks_2samp


def kolmogorov_smirnov_test(experiment_scores, control_scores):
    """
    Perform the Kolmogorov-Smirnov test to compare two samples.

    Parameters:
    - experiment_scores (list): Scores from the experimental group.
    - control_scores (list): Scores from the control group.

    Returns:
    float: The P-value from the KS test indicating statistical difference.
    """
    # Convert lists to numpy arrays for efficient operations
    experiment_scores = np.sort(experiment_scores).ravel()
    control_scores = np.sort(control_scores).ravel()

    # Calculate the length of each sample
    en1 = len(experiment_scores)
    en2 = len(control_scores)

    # Combine the scores and compute cumulative distribution functions
    data_all = np.concatenate([experiment_scores, control_scores])
    cdf_experiment = np.searchsorted(experiment_scores, data_all, side='right') / en1
    cdf_control = np.searchsorted(control_scores, data_all, side='right') / en2

    # Compute the maximum distance between the two CDFs
    D = np.max(np.abs(cdf_experiment - cdf_control))

    # Calculate the KS statistic
    en = np.sqrt(en1 * en2 / (en1 + en2))
    p_value = ks((en + 0.12 + 0.11 / en) * D)

    return p_value


def ks(alam):
    """
    Compute the Kolmogorov-Smirnov probability given a lambda value.

    Parameters:
    - alam (float): Lambda value for the KS statistic.

    Returns:
    float: The probability associated with the KS statistic.
    """
    EPS1 = 1e-6  # Precision for the convergence of term's absolute value
    EPS2 = 1e-10  # Precision for the convergence of the series_sum's relative value
    a2 = -2.0 * alam ** 2  # Adjust lambda for exponential calculation
    fac = 2.0
    series_sum = 0.0
    previous_term = 0.0

    # Sum the series until convergence criteria are met
    for j in range(1, 101):
        term = fac * np.exp(a2 * j ** 2)  # Calculate term of the series
        series_sum += term  # Add to series_sum

        # Check for convergence
        if np.abs(term) <= EPS1 * previous_term or np.abs(term) <= EPS2 * series_sum:
            return series_sum

        fac = -fac  # Alternate the sign
        previous_term = np.abs(term)  # Update the term before flag

    # Return 1.0 if the series does not converge within 100 terms
    return 1.0


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
    """
    # Extract scores
    gene_ids = list(scores.keys())
    gene_scores = [score[0] for score in scores.values()]

    # Create a DataFrame for easier handling
    df = pd.DataFrame({'GeneID': gene_ids, 'Score': gene_scores})

    # Rank the scores (higher scores get lower rank numbers)
    df['Rank'] = df['Score'].rank(ascending=False, method='average')

    # Create a Series with GeneID as index and Rank as values
    global_ranking = df.set_index('GeneID')['Rank']

    return global_ranking


def kolmogorov_smirnov_test_with_ranking(pathway_genes, global_ranking):
    """
    Perform the Kolmogorov-Smirnov test using global ranking.

    Parameters:
    - pathway_genes (iterable): Iterable of gene IDs in the pathway.
    - global_ranking (pd.Series): Global ranking of all genes.

    Returns:
    - float: The P-value from the KS test indicating statistical difference.
    """
    # Convert pathway_genes to a list if it's not already
    pathway_genes = list(pathway_genes)

    # Ensure all genes are in the global ranking
    pathway_genes_in_ranking = [gene for gene in pathway_genes if gene in global_ranking.index]

    if not pathway_genes_in_ranking:
        # If no genes are found, return a p-value of 1.0
        return 1.0

    # Get the ranks for pathway genes
    pathway_ranks = global_ranking[pathway_genes_in_ranking].values

    # Get the ranks for background genes
    background_genes = global_ranking.index.difference(pathway_genes_in_ranking)
    background_ranks = global_ranking[background_genes].values

    # Perform the KS test
    D, p_value = ks_2samp(
        pathway_ranks,
        background_ranks,
        alternative='two-sided',
        mode='auto'
    )

    return p_value
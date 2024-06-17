import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm, hypergeom


def wilcoxon_rank_sums_test(experiment_scores: list, control_scores: list, alternative: str = 'two-sided') -> float:
    """
    Perform the Wilcoxon rank-sum test to compare two independent samples.

    Parameters:
    - experiment_scores (list): Scores from the experimental group.
    - control_scores (list): Scores from the control group.
    - alternative (str): Defines the alternative hypothesis ('two-sided', 'less', 'greater').

    Returns:
    - float: The P-value from the Wilcoxon rank-sum test.
    """
    from scipy.stats import ranksums
    p_vals = ranksums(experiment_scores, control_scores, alternative=alternative).pvalue
    return p_vals


def bh_correction(p_values: np.ndarray) -> np.ndarray:
    """
    Apply the Benjamini-Hochberg procedure for controlling the false discovery rate in multiple hypothesis testing.

    Parameters:
    - p_values (array_like): Array of p-values obtained from multiple statistical tests.

    Returns:
    - numpy.ndarray: Array of adjusted p-values.
    """
    p_vals_rank = rankdata(p_values, 'max') - 1
    p_vals_rank_ord = rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) / (p_vals_rank + 1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_values_sorted]

    return adj_p_vals


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


def compute_mw_python(experiment_ranks, control_ranks):
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
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

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
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - pd.Series: A series with gene IDs as index and their global ranks as values.
    """
    # Extract scores and rank them
    gene_ids = list(scores.keys())
    gene_scores = [score[0] for score in scores.values()]

    # Rank the scores (higher scores get lower rank numbers)
    ranks = rankdata(gene_scores, method='average')
    global_ranking = pd.Series(ranks, index=gene_ids)

    return global_ranking


def kolmogorov_smirnov_test_with_ranking(pathway_genes, global_ranking):
    """
    Perform the Kolmogorov-Smirnov test using global ranking.

    Parameters:
    - pathway_genes (list): List of gene IDs in the pathway.
    - global_ranking (pd.Series): Global ranking of all genes.

    Returns:
    - float: The P-value from the KS test indicating statistical difference.
    """
    # Convert pathway genes to a list if it's a set
    pathway_genes = list(pathway_genes)

    # Get the ranks for pathway genes
    pathway_ranks = global_ranking[pathway_genes].sort_values().values
    background_ranks = global_ranking.drop(pathway_genes).values

    # Calculate the KS statistic
    en1 = len(pathway_ranks)
    en2 = len(background_ranks)
    cdf_pathway = np.searchsorted(pathway_ranks, global_ranking, side='right') / en1
    cdf_background = np.searchsorted(background_ranks, global_ranking, side='right') / en2
    D = np.max(np.abs(cdf_pathway - cdf_background))

    en = np.sqrt(en1 * en2 / (en1 + en2))
    p_value = ks((en + 0.12 + 0.11 / en) * D)

    return p_value
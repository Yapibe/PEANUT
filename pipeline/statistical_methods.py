import numpy as np
import pandas as pd
from scipy.stats import hypergeom, ks_2samp, mannwhitneyu, rankdata


def wilcoxon_rank_sums_test(
    experiment_scores: np.ndarray,
    control_scores: np.ndarray,
    alternative: str = "two-sided",
) -> float:
    """
    Perform the Wilcoxon rank-sum test to compare two independent samples.

    Parameters:
    - experiment_scores (array_like): Scores from the experimental group.
    - control_scores (array_like): Scores from the control group.
    - alternative (str): Defines the alternative hypothesis
      ('two-sided', 'less', 'greater').

    Returns:
    - float: The P-value from the Wilcoxon rank-sum test.
    """
    from scipy.stats import ranksums

    p_value = ranksums(
        experiment_scores, control_scores, alternative=alternative
    ).pvalue
    return p_value


def bh_correction(p_values: np.ndarray) -> np.ndarray:
    """
    Apply the Benjamini-Hochberg procedure for controlling the false discovery rate.

    Parameters:
    - p_values (array_like): Array of p-values from multiple tests.

    Returns:
    - np.ndarray: Array of adjusted p-values.
    """
    p_values = np.asarray(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_pvalues = p_values[sorted_indices]
    adjusted_pvalues = np.empty(n, dtype=np.float64)

    cumulative_min = 1.0
    for i in reversed(range(n)):
        rank = i + 1
        adjusted_pval = sorted_pvalues[i] * n / rank
        cumulative_min = min(cumulative_min, adjusted_pval)
        adjusted_pvalues[sorted_indices[i]] = cumulative_min

    return adjusted_pvalues


def kolmogorov_smirnov_test(
    experiment_scores: np.ndarray, control_scores: np.ndarray
) -> float:
    """
    Perform the Kolmogorov-Smirnov test to compare two samples.

    Parameters:
    - experiment_scores (array_like): Scores from the experimental group.
    - control_scores (array_like): Scores from the control group.

    Returns:
    - float: The P-value from the KS test indicating statistical difference.
    """
    result = ks_2samp(
        experiment_scores,
        control_scores,
        alternative="two-sided",
        mode="auto",
    )
    return result.pvalue


def compute_mw_python(
    experiment_ranks: np.ndarray, control_ranks: np.ndarray
) -> tuple:
    """
    Compute the Mann-Whitney U test using ranks.

    Parameters:
    - experiment_ranks (array_like): Ranks of the experimental group.
    - control_ranks (array_like): Ranks of the control group.

    Returns:
    - tuple: The Mann-Whitney U statistic and the corresponding p-value.
    """
    U, p_value = mannwhitneyu(
        experiment_ranks, control_ranks, alternative="two-sided"
    )
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
    intersection_size = len(set1 & set2)
    union_size = len(set1 | set2)
    return intersection_size / union_size if union_size > 0 else 0.0


def run_hyper(
    genes_by_pathway: dict, scores_keys: set, significant_p_vals: dict
) -> list:
    """
    Run the hypergeometric test to identify pathways with significant enrichment.

    Parameters:
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores_keys (set): Set of gene IDs with scores.
    - significant_p_vals (dict): Mapping of gene IDs to significant P-values.

    Returns:
    - list: List of significant pathways identified by the hypergeometric test.
    """
    M = len(scores_keys)  # Total number of scored genes
    n = len(significant_p_vals)  # Number of significant genes

    significant_pathways = []
    for pathway_name, pathway_genes in genes_by_pathway.items():
        N = len(pathway_genes)  # Number of genes in the pathway
        x = len(
            set(pathway_genes) & significant_p_vals.keys()
        )  # Overlap
        if x >= 5:
            p_value = hypergeom.sf(x - 1, M, N, n)
            if p_value < 0.05:
                significant_pathways.append(pathway_name)
    return significant_pathways


def global_gene_ranking(scores: dict) -> pd.Series:
    """
    Rank all genes globally based on their scores.

    Parameters:
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - pd.Series: Series with gene IDs as index and their global ranks as values.
    """
    gene_ids = list(scores.keys())
    gene_scores = np.array(
        [score[0] for score in scores.values()], dtype=np.float64
    )

    # Rank scores (higher scores get higher rank numbers)
    ranks = rankdata(-gene_scores, method="average")
    global_ranking = pd.Series(ranks, index=gene_ids)

    return global_ranking


def kolmogorov_smirnov_test_with_ranking(
    pathway_genes: set, global_ranking: pd.Series
) -> float:
    """
    Perform the Kolmogorov-Smirnov test using global ranking.

    Parameters:
    - pathway_genes (set): Set of gene IDs in the pathway.
    - global_ranking (pd.Series): Global ranking of all genes.

    Returns:
    - float: The P-value from the KS test indicating statistical difference.
    """
    pathway_ranks = global_ranking.loc[
        global_ranking.index.intersection(pathway_genes)
    ].values
    background_ranks = global_ranking.loc[
        ~global_ranking.index.isin(pathway_genes)
    ].values

    # Perform KS test
    result = ks_2samp(
        pathway_ranks,
        background_ranks,
        alternative="two-sided",
        mode="auto",
    )
    return result.pvalue

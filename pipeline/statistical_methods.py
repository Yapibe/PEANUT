import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests

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


def perform_permutation_test(filtered_pathways: dict, genes_by_pathway: dict, scores: dict, n_iterations: int = 10000, pval_threshold: float = 0.05) -> dict:
    """
    Perform a permutation test to assess the significance of the observed pathway scores.
    Updates pathways with adjusted permutation p-values and adds a boolean flag for significance.
    
    Parameters:
    - filtered_pathways (dict): Dictionary of pathways with their details.
    - genes_by_pathway (dict): Mapping of pathway names to their gene sets.
    - scores (dict): Mapping of gene IDs to their scores.
    - n_iterations (int): Number of permutations to perform.
    - pval_threshold (float): Threshold for adjusted p-values to determine significance.
    
    Returns:
    - dict: Updated filtered_pathways with permutation p-values and significance flags.
    """
    
    score_values = np.array([score[0] for score in scores.values()])
    
    # Store observed means and p-values for all pathways
    observed_means = []
    raw_p_values = []
    pathways_list = []
    valid_pathways = []
    
    for pathway, data in filtered_pathways.items():
        pathway_genes = genes_by_pathway[pathway]
    
        # Find the intersection of pathway genes with available scores
        valid_genes = pathway_genes.intersection(scores.keys())
        pathway_size = len(valid_genes)
    
        if pathway_size == 0:
            data['Permutation p-val'] = np.nan
            data['permutation_significant'] = False
            continue
    
        # Calculate the observed mean score for the pathway
        observed_mean = np.mean([scores[gene_id][0] for gene_id in valid_genes])
    
        # Generate null distribution for this pathway size
        permuted_means = np.array([
            np.mean(np.random.choice(score_values, size=pathway_size, replace=False))
            for _ in range(n_iterations)
        ])
        empirical_p_val = np.mean(permuted_means >= observed_mean)
    
        # Collect data for multiple testing correction
        observed_means.append(observed_mean)
        raw_p_values.append(empirical_p_val)
        pathways_list.append(pathway)
        valid_pathways.append((pathway, data))
    
    # Adjust p-values for multiple testing using Benjamini-Hochberg
    if raw_p_values:
        adjusted_pvals = multipletests(raw_p_values, method='fdr_bh')[1]
    
        # Update pathways with adjusted p-values and significance flags
        for i, (pathway, data) in enumerate(valid_pathways):
            adjusted_pval = adjusted_pvals[i]
            data['Permutation p-val'] = adjusted_pval
            data['permutation_significant'] = adjusted_pval <= pval_threshold
    else:
        # Set permutation_significant to False for all pathways
        for data in filtered_pathways.values():
            data['permutation_significant'] = False
    
    return filtered_pathways
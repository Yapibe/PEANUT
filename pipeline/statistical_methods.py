import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, mannwhitneyu, rankdata
from statsmodels.stats.multitest import multipletests
import multiprocessing
from functools import partial
from typing import Dict, List, Tuple
import logging

logger = logging.getLogger(__name__)

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


def global_gene_ranking(scores: Dict[int, float]) -> pd.Series:
    """
    Rank all genes globally based on their scores.

    Parameters:
    - scores (Dict[int, float]): Mapping of gene IDs to their scores.

    Returns:
    - pd.Series: Series with gene IDs as index and their global ranks as values.
    """
    gene_ids = list(scores.keys())
    gene_scores = np.array(list(scores.values()), dtype=np.float64)

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


def perform_permutation_for_pathway(
    args: Tuple[str, Dict],
    genes_by_pathway: Dict[str, List[int]],
    scores: Dict[int, Tuple[float, float]],
    score_values: np.ndarray,
    n_iterations: int,
    pval_threshold: float
) -> Tuple[str, Dict]:
    """
    Perform permutation test for a single pathway.

    Returns:
    - Tuple containing pathway name and updated data.
    """
    pathway, data = args  # Unpack the tuple
    pathway_genes = set(genes_by_pathway[pathway])
    valid_genes = pathway_genes.intersection(scores.keys())
    pathway_size = len(valid_genes)

    if pathway_size == 0:
        logger.warning(f"No valid genes found for pathway '{pathway}'. Skipping permutation test.")
        data['Permutation p-val'] = np.nan
        data['permutation_significant'] = False
        return pathway, data

    # Calculate the observed mean score for the pathway
    observed_mean = np.mean([scores[gene_id] for gene_id in valid_genes])

    # Generate null distribution for this pathway size
    permuted_means = np.array([
        np.mean(np.random.choice(score_values, size=pathway_size, replace=False))
        for _ in range(n_iterations)
    ])
    empirical_p_val = np.mean(permuted_means >= observed_mean)

    # Temporary: Assign raw p-value
    data['Permutation p-val'] = empirical_p_val
    data['permutation_significant'] = empirical_p_val <= pval_threshold

    return pathway, data


def perform_permutation_test_parallel(
    filtered_pathways: Dict[str, Dict],
    genes_by_pathway: Dict[str, List[int]],
    scores: Dict[int, Tuple[float, float]],
    n_iterations: int = 10000,
    pval_threshold: float = 0.05
) -> Dict[str, Dict]:
    """
    Perform permutation tests in parallel for all pathways.

    Parameters:
    - filtered_pathways (Dict[str, Dict]): Dictionary of pathways with their details.
    - genes_by_pathway (Dict[str, List[int]]): Mapping of pathway names to their gene sets.
    - scores (Dict[int, Tuple[float, float]]): Mapping of gene IDs to their scores.
    - n_iterations (int): Number of permutations to perform.
    - pval_threshold (float): Threshold for adjusted p-values to determine significance.

    Returns:
    - Dict[str, Dict]: Updated filtered_pathways with permutation p-values and significance flags.
    """
    if not filtered_pathways:
        logger.info("No pathways to test with the permutation test.")
        return {}

    # Extract score values for permutation
    score_values = np.array(list(scores.values()))

    # Prepare arguments for parallel processing
    pathways = list(filtered_pathways.keys())
    pathway_data = [(pathway, filtered_pathways[pathway]) for pathway in pathways]
    
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    func = partial(
        perform_permutation_for_pathway,
        genes_by_pathway=genes_by_pathway,
        scores=scores,
        score_values=score_values,
        n_iterations=n_iterations,
        pval_threshold=pval_threshold
    )

    results = pool.map(func, pathway_data)
    pool.close()
    pool.join()

    # Collect results
    for pathway, data in results:
        filtered_pathways[pathway] = data

    # Adjust p-values using Benjamini-Hochberg method
    raw_p_values = [
        data['Permutation p-val'] for pathway, data in filtered_pathways.items()
        if not np.isnan(data['Permutation p-val'])
    ]

    if raw_p_values:
        adjusted_pvals = multipletests(raw_p_values, method='fdr_bh')[1]
        idx = 0
        for pathway, data in filtered_pathways.items():
            if not np.isnan(data['Permutation p-val']):
                data['Permutation p-val'] = adjusted_pvals[idx]
                data['permutation_significant'] = adjusted_pvals[idx] <= pval_threshold
                idx += 1
    else:
        logger.info("No valid p-values obtained from permutation test.")
        # Set permutation_significant to False for all pathways
        for data in filtered_pathways.values():
            data['permutation_significant'] = False

    return filtered_pathways
import decimal
import logging
from typing import Dict, List, Tuple
import gseapy as gp
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

from .settings import Settings, ConditionSettings
from .statistical_methods import (
    compute_mw_python,
    global_gene_ranking,
    kolmogorov_smirnov_test_with_ranking,
)
from .utils import load_pathways_and_propagation_scores

logger = logging.getLogger(__name__)

import multiprocessing
from functools import partial

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
    observed_mean = np.mean([scores[gene_id][0] for gene_id in valid_genes])

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

    score_values = np.array([score[0] for score in scores.values()])

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

def perform_statist_ks(
    settings: Settings,
    condition_settings: ConditionSettings,
    genes_by_pathway: Dict[str, List[int]],
    scores: Dict[int, Tuple[float, float]],
) -> None:
    """
    Perform statistical enrichment analysis using the Kolmogorov-Smirnov test.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - genes_by_pathway (Dict[str, List[int]]): Mapping of pathways to their constituent genes.
    - scores (Dict[int, Tuple[float, float]]): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    try:
        # Collect all genes from the pathways
        condition_settings.filtered_genes = {
            gene_id for genes in genes_by_pathway.values() for gene_id in genes
        }

        # Generate global ranking
        global_ranking = global_gene_ranking(scores)

        pathways = list(genes_by_pathway.keys())


        # Perform KS test for each pathway
        ks_p_values = [
            kolmogorov_smirnov_test_with_ranking(genes_by_pathway[pathway], global_ranking)
            for pathway in pathways
        ]

        
        # Adjust p-values using Benjamini-Hochberg method
        adjusted_p_values = multipletests(ks_p_values, method="fdr_bh")[1]

        # Store all pathways with their adjusted KS p-values and significance status
        condition_settings.pathways_statistics = {
            pathway: {
                'genes': genes_by_pathway[pathway],
                'ks_p_value': adj_pval,
                'ks_significant': adj_pval < settings.FDR_THRESHOLD
            }
            for pathway, adj_pval in zip(pathways, adjusted_p_values)
        }

    except Exception as e:
        logger.error(f"An error occurred during KS test: {e}", exc_info=True)
        raise



def perform_statist_mann_whitney(
    condition_settings: ConditionSettings,
    scores: Dict[int, Tuple[float, float]],
) -> None:
    """
    Perform Mann-Whitney U test on all pathways, followed by permutation testing.
    Updates pathway statistics with MW p-values, permutation p-values, and significance flags.

    Parameters:
    - condition_settings (ConditionSettings): Condition-specific settings.
    - scores (Dict[int, Tuple[float, float]]): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    try:
        pathways = list(condition_settings.pathways_statistics.keys())

        # Prepare for ranking
        filtered_scores = np.array(
            [scores[gene_id][0] for gene_id in condition_settings.filtered_genes],
            dtype=np.float32
        )

        # Compute ranks
        ranks = rankdata(filtered_scores, method="average")
        scores_rank = dict(zip(condition_settings.filtered_genes, ranks))

        mw_p_values = []
        for pathway in pathways:
            pathway_genes = condition_settings.pathways_statistics[pathway]['genes']
            pathway_genes_set = set(pathway_genes)

            background_genes = condition_settings.filtered_genes - pathway_genes_set

            pathway_ranks = np.array(
                [scores_rank[gene_id] for gene_id in pathway_genes_set if gene_id in scores_rank],
                dtype=np.float32,
            )
            background_ranks = np.array(
                [scores_rank[gene_id] for gene_id in background_genes if gene_id in scores_rank],
                dtype=np.float32,
            )

            # Perform Mann-Whitney U test
            stat, rmw_pval = compute_mw_python(
                pathway_ranks, background_ranks
            )
            mw_p_values.append(rmw_pval)

            # Update pathway information
            condition_settings.pathways_statistics[pathway].update({
                "NOM p-val": rmw_pval,
                "Tag %": (len(pathway_genes_set) / len(condition_settings.filtered_genes)) * 100,
                "Gene %": (len(pathway_genes_set) / len(scores)) * 100,
                "Lead_genes": ",".join(map(str, sorted(pathway_genes_set))),
            })

        if not mw_p_values:
            logger.info("No pathways to test. Skipping FDR correction.")
            return

        # Apply Benjamini-Hochberg correction to adjust the MW p-values
        adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]
        for i, pathway in enumerate(pathways):
            data = condition_settings.pathways_statistics[pathway]
            data['FDR q-val'] = adjusted_mw_p_values[i]
            data['mw_significant'] = adjusted_mw_p_values[i] < 0.05

        # Perform Permutation Test
        condition_settings.pathways_statistics = perform_permutation_test_parallel(
            filtered_pathways=condition_settings.pathways_statistics,
            genes_by_pathway={pathway: data['genes'] for pathway, data in condition_settings.pathways_statistics.items()},
            scores=scores,
            n_iterations=100,
            pval_threshold=0.05
        )

        # Assign ranks based on FDR q-values
        sorted_pathways = sorted(
            condition_settings.pathways_statistics.items(),
            key=lambda item: (item[1]["FDR q-val"], item[0])  # Sort by q-val then pathway name
        )

        for rank, (pathway, data) in enumerate(sorted_pathways, start=1):
            condition_settings.pathways_statistics[pathway]["Rank"] = rank
            # Example: Combine KS and MW p-values for NES (optional)
            # data["NES"] = compute_nes(data)

        logger.info(f"Completed Mann-Whitney U test and ranking for {len(pathways)} pathways.")

    except Exception as e:
        logger.error(f"An error occurred during Mann-Whitney test: {e}", exc_info=True)
        raise



def perform_enrichment(
    settings: Settings,
    condition_settings: ConditionSettings,
    scores_df: pd.DataFrame,
):
    """
    Perform pathway enrichment analysis for a given condition.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - scores_df (pd.DataFrame): DataFrame containing gene scores for the condition.

    Returns:
    - None
    """
    try:
        genes_by_pathway, scores = (
            load_pathways_and_propagation_scores(settings, scores_df)
        )

        if settings.run_gsea:
            # Prepare data for GSEA
            gene_expression_data = pd.DataFrame(
                {
                    "gene": list(scores.keys()),
                    "logFC": [score[0] for score in scores.values()],
                }
            )

            gene_expression_data["gene"] = gene_expression_data[
                "gene"
            ].astype(str)
            gene_expression_data.sort_values(
                by="logFC", ascending=False, inplace=True
            )

            # Run GSEA
            gsea_results = gp.prerank(
                rnk=gene_expression_data,
                gene_sets=genes_by_pathway,
                verbose=True,
                no_plot=True,
            )
            # Save results
            gsea_results.res2d.to_excel(
                condition_settings.output_path
            )
        else:
            # Perform Kolmogorov-Smirnov test
            perform_statist_ks(settings, condition_settings, genes_by_pathway, scores)
            
            # Perform Mann-Whitney U test followed by permutation testing
            perform_statist_mann_whitney(condition_settings, scores)
        
            # Proceed only if there are pathways to process
            if condition_settings.pathways_statistics:
                # Convert pathways_statistics dictionary to DataFrame
                pathways_df = pd.DataFrame.from_dict(
                    condition_settings.pathways_statistics, orient='index'
                ).reset_index().rename(columns={'index': 'Name'})
            
                # Ensure correct data types using a dictionary for multiple columns
                dtype_mapping = {
                    'NOM p-val': float,
                    'FDR q-val': float,
                    'ks_significant': bool,
                    'mw_significant': bool,
                    'permutation_significant': bool,
                    'Permutation p-val': float
                }
                for col, dtype in dtype_mapping.items():
                    if col in pathways_df.columns:
                        pathways_df[col] = pathways_df[col].astype(dtype)
            
                # Define grouping criteria using vectorized conditions
                conditions = [
                    (pathways_df['ks_significant']) & (pathways_df['mw_significant']) & (pathways_df['permutation_significant']),
                    (pathways_df['ks_significant']) & (pathways_df['mw_significant']) & (~pathways_df['permutation_significant']),
                    (pathways_df['ks_significant']) & (~pathways_df['mw_significant']),
                    (~pathways_df['ks_significant'])
                ]
                choices = [1, 2, 3, 4]
                pathways_df['Group'] = np.select(conditions, choices, default=4)
            
                # Sort pathways by Group and then by Mann-Whitney p-value within each group
                pathways_df.sort_values(by=['Group', 'NOM p-val'], ascending=[True, True], inplace=True)
                pathways_df.reset_index(drop=True, inplace=True)
                pathways_df['Rank'] = pathways_df.index + 1
            
                # Reorder columns to prioritize important information
                desired_columns = ['Rank', 'Name', 'Group', 'FDR q-val', 'ks_significant', 
                                   'mw_significant', 'permutation_significant', 'Permutation p-val', 
                                   'NOM p-val']
                # Add any additional columns that might exist
                additional_columns = [col for col in pathways_df.columns if col not in desired_columns]
                pathways_df = pathways_df[desired_columns + additional_columns]
            
                # Save the DataFrame to Excel
                pathways_df.to_excel(condition_settings.output_path, index=False)
            
                logger.info(f"Enrichment analysis completed. Results saved to {condition_settings.output_path}.")
            else:
                logger.info("No pathways to output.")
    except Exception as e:
        logger.error(f"An error occurred during enrichment analysis: {e}", exc_info=True)
        raise

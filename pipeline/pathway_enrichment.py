"""
Pathway enrichment analysis module.

This module contains functions for performing enrichment analysis on gene sets,
including both standard statistical tests (KS test, Mann-Whitney U test) and GSEA.
"""

import logging
from typing import Dict, List
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
    perform_permutation_test_parallel
)
from .utils import load_pathways_and_propagation_scores

logger = logging.getLogger(__name__)

def perform_statist_ks(
    settings: Settings,
    condition_settings: ConditionSettings,
    genes_by_pathway: Dict[str, List[int]],
    scores: Dict[int, float],
) -> None:
    """
    Perform statistical enrichment analysis using the Kolmogorov-Smirnov test.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - genes_by_pathway (Dict[str, List[int]]): Mapping of pathways to their constituent genes.
    - scores (Dict[int, float]): Mapping of gene IDs to their scores.

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
    scores: Dict[int, float],
) -> None:
    """
    Perform Mann-Whitney U test on all pathways, followed by permutation testing.
    Updates pathway statistics with MW p-values, permutation p-values, and significance flags.

    Parameters:
    - condition_settings (ConditionSettings): Condition-specific settings.
    - scores (Dict[int, float]): Mapping of gene IDs to their scores.

    Returns:
    - None
    """
    try:
        pathways = list(condition_settings.pathways_statistics.keys())
        if not pathways:
            logger.info("No pathways to analyze with Mann-Whitney U test.")
            return

        # Prepare for ranking
        filtered_scores = np.array(
            [scores[gene_id] for gene_id in condition_settings.filtered_genes if gene_id in scores],
            dtype=np.float32
        )

        # Compute ranks
        ranks = rankdata(filtered_scores, method="average")
        scores_rank = dict(zip(
            [gene_id for gene_id in condition_settings.filtered_genes if gene_id in scores], 
            ranks
        ))

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

            # Skip if insufficient data for the test
            if len(pathway_ranks) < 2 or len(background_ranks) < 2:
                logger.warning(f"Insufficient data for Mann-Whitney test in pathway {pathway}. Skipping.")
                mw_p_values.append(1.0)  # Use 1.0 as a non-significant p-value
                continue

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

        logger.info(f"Completed Mann-Whitney U test and ranking for {len(pathways)} pathways.")

    except Exception as e:
        logger.error(f"An error occurred during Mann-Whitney test: {e}", exc_info=True)
        raise



def perform_enrichment(
    settings: Settings,
    condition_settings: ConditionSettings,
    scores_df: pd.DataFrame,
) -> None:
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
        genes_by_pathway, scores = load_pathways_and_propagation_scores(settings, scores_df)
        
        if not genes_by_pathway:
            logger.warning("No pathways available for enrichment analysis.")
            return

        if settings.run_gsea:
            _run_gsea_analysis(settings, condition_settings, genes_by_pathway, scores)
        else:
            # Perform Kolmogorov-Smirnov test
            perform_statist_ks(settings, condition_settings, genes_by_pathway, scores)
            
            # Perform Mann-Whitney U test followed by permutation testing
            perform_statist_mann_whitney(condition_settings, scores)
        
            # Process and save results if there are pathways to analyze
            if condition_settings.pathways_statistics:
                _save_enrichment_results(settings, condition_settings)
            else:
                logger.info("No pathways to output.")
                
    except Exception as e:
        logger.error(f"An error occurred during enrichment analysis: {e}", exc_info=True)
        raise


def _run_gsea_analysis(
    settings: Settings,
    condition_settings: ConditionSettings,
    genes_by_pathway: Dict[str, List[int]],
    scores: Dict[int, float]
) -> None:
    """
    Run Gene Set Enrichment Analysis (GSEA) using the gseapy package.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - genes_by_pathway (Dict[str, List[int]]): Mapping of pathways to their genes.
    - scores (Dict[int, float]): Mapping of gene IDs to their scores.

    Returns:
    - None
    """
    try:
        # Prepare data for GSEA
        gene_expression_data = pd.DataFrame({
            "gene": list(scores.keys()),
            "logFC": list(scores.values()),  # Use scores directly as logFC
        })

        gene_expression_data["gene"] = gene_expression_data["gene"].astype(str)
        gene_expression_data.sort_values(by="logFC", ascending=False, inplace=True)

        # Run GSEA
        gsea_results = gp.prerank(
            rnk=gene_expression_data,
            gene_sets=genes_by_pathway,
            verbose=True,
            no_plot=True,
            min_size=settings.minimum_gene_per_pathway,
            max_size=settings.maximum_gene_per_pathway,
        )
    
        # Save results
        gsea_results.res2d.reset_index(inplace=True)
        gsea_results.res2d.index = gsea_results.res2d.index + 1
        gsea_results.res2d.rename(columns={'index': 'Rank'}, inplace=True)
        gsea_results.res2d.to_excel(condition_settings.output_path, index=False)

        # Add pathway statistics to condition_settings
        condition_settings.pathways_statistics = {
            row['Term']: {
                'Rank': row['Rank'],
                'NES': row['NES'],
                'FDR q-val': row['FDR q-val'],
                'genes': genes_by_pathway.get(row['Term'], []),
                'gsea_significant': row['FDR q-val'] < settings.FDR_THRESHOLD
            }
            for _, row in gsea_results.res2d.iterrows()
        }
        
        logger.info(f"GSEA analysis completed with {len(condition_settings.pathways_statistics)} pathways.")
        
    except Exception as e:
        logger.error(f"Error during GSEA analysis: {e}", exc_info=True)
        raise


def _save_enrichment_results(
    settings: Settings,
    condition_settings: ConditionSettings
) -> None:
    """
    Process and save enrichment analysis results to Excel.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.

    Returns:
    - None
    """
    try:
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
        
    except Exception as e:
        logger.error(f"Error when saving enrichment results: {e}", exc_info=True)
        raise

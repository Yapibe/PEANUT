import decimal
import logging

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


def perform_statist_ks(
    settings: Settings,
    condition_settings: ConditionSettings,
    genes_by_pathway: dict,
    scores: dict,
):
    """
    Perform statistical enrichment analysis using the Kolmogorov-Smirnov test.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    try:
        # Collect all genes from the pathways
        condition_settings.filtered_genes = set(
            gene_id
            for genes in genes_by_pathway.values()
            for gene_id in genes
        )

        # Generate global ranking
        global_ranking = global_gene_ranking(scores)

        significant_pathways = list(genes_by_pathway.keys())

        # Perform KS test for each pathway
        ks_p_values = []
        for pathway in significant_pathways:
            genes = genes_by_pathway[pathway]
            p_value = kolmogorov_smirnov_test_with_ranking(
                genes, global_ranking
            )
            ks_p_values.append(p_value)

        if not ks_p_values:
            logger.info(
                "No significant pathways found after hypergeometric test. Skipping KS test."
            )
            return

        # Adjust p-values using Benjamini-Hochberg method
        adjusted_p_values = multipletests(
            ks_p_values, method="fdr_bh"
        )[1]

        # Store significant pathways
        condition_settings.ks_significant_pathways_with_genes = {
            pathway: (genes_by_pathway[pathway], adjusted_p_values[i])
            for i, pathway in enumerate(significant_pathways)
            if adjusted_p_values[i] < settings.FDR_THRESHOLD
        }

        if not condition_settings.ks_significant_pathways_with_genes:
            logger.info(
                "No significant pathways found after KS test."
            )
    except Exception as e:
        logger.error(f"An error occurred during KS test: {e}")
        raise


def perform_statist_mann_whitney(
    settings: Settings,
    condition_settings: ConditionSettings,
    scores: dict,
):
    """
    Perform Mann-Whitney U test on pathways that passed the KS test and
    filter significant pathways.

    Parameters:
    - settings (Settings): General settings for the experiment.
    - condition_settings (ConditionSettings): Condition-specific settings.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    try:
        mw_p_values = []
        filtered_scores = [
            scores[gene_id][0]
            for gene_id in condition_settings.filtered_genes
        ]
        filtered_scores = np.array(filtered_scores, dtype=np.float32)

        # Compute ranks
        ranks = rankdata(filtered_scores, method="average")
        scores_rank = dict(
            zip(condition_settings.filtered_genes, ranks)
        )

        for pathway, (
            pathway_genes,
            _,
        ) in (
            condition_settings.ks_significant_pathways_with_genes.items()
        ):
            pathway_genes_set = set(pathway_genes)
            pathway_scores = [
                scores[gene_id][0] for gene_id in pathway_genes_set
            ]

            background_genes = (
                condition_settings.filtered_genes - pathway_genes_set
            )

            pathway_ranks = np.array(
                [
                    scores_rank[gene_id]
                    for gene_id in pathway_genes_set
                ],
                dtype=np.float32,
            )
            background_ranks = np.array(
                [
                    scores_rank[gene_id]
                    for gene_id in background_genes
                ],
                dtype=np.float32,
            )

            # Perform Mann-Whitney U test
            _, rmw_pval = compute_mw_python(
                pathway_ranks, background_ranks
            )
            mw_p_values.append(rmw_pval)

            condition_settings.filtered_pathways[pathway] = {
                "Rank": len(condition_settings.filtered_pathways) + 1,
                "Term": pathway,
                "ES": (
                    np.mean(pathway_scores) if pathway_scores else 0
                ),
                "NES": None,
                "NOM p-val": rmw_pval,
                "FDR q-val": None,
                "FWER p-val": None,
                "Tag %": len(pathway_genes_set)
                / len(condition_settings.filtered_genes)
                * 100,
                "Gene %": len(pathway_genes_set) / len(scores) * 100,
                "Lead_genes": ",".join(map(str, pathway_genes_set)),
            }

        # Adjust p-values
        adjusted_mw_p_values = multipletests(
            mw_p_values, method="fdr_bh"
        )[1]

        # Update FDR q-values
        for i, pathway in enumerate(
            condition_settings.filtered_pathways.keys()
        ):
            condition_settings.filtered_pathways[pathway][
                "FDR q-val"
            ] = decimal.Decimal(adjusted_mw_p_values[i])
    except Exception as e:
        logger.error(
            f"An error occurred during Mann-Whitney test: {e}"
        )
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
            # Perform statistical tests
            perform_statist_ks(
                settings, condition_settings, genes_by_pathway, scores
            )
            if condition_settings.ks_significant_pathways_with_genes:
                perform_statist_mann_whitney(
                    settings, condition_settings, scores
                )
            else:
                logger.info("Skipping Mann-Whitney test.")

            if condition_settings.filtered_pathways:
                pathways_df = pd.DataFrame.from_dict(
                    condition_settings.filtered_pathways,
                    orient="index",
                )

                # Format FDR q-values
                pathways_df["FDR q-val"] = pathways_df[
                    "FDR q-val"
                ].apply(
                    lambda x: (
                        decimal.Decimal(x) if x is not None else x
                    )
                )

                pathways_df.sort_values(by="FDR q-val", inplace=True)
                pathways_df["Rank"] = range(1, len(pathways_df) + 1)

                # Save results
                pathways_df.to_excel(
                    condition_settings.output_path, index=False
                )
    except Exception as e:
        logger.error(
            f"An error occurred during enrichment analysis: {e}"
        )
        raise

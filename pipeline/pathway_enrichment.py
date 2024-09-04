import pandas as pd
from .settings import Settings
from scipy.stats import rankdata
import gseapy as gp
from .utils import load_pathways_and_propagation_scores
from statsmodels.stats.multitest import multipletests
from .statistical_methods import compute_mw_python, global_gene_ranking, kolmogorov_smirnov_test_with_ranking
import numpy as np
import logging
import decimal


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a logger
logger = logging.getLogger(__name__)
def perform_statist_ks(settings :Settings, genes_by_pathway: dict, scores: dict):
    """
    Perform statistical enrichment analysis on pathways.

    Parameters:
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    # Populate a set with all genes from the filtered pathways
    for genes in genes_by_pathway.values():
        settings.filtered_genes.update(genes)

    # Create a global ranking of genes
    global_ranking = global_gene_ranking(scores)
    significant_pathways = list(genes_by_pathway.keys())

    # Perform the Kolmogorov-Smirnov test using global ranking
    ks_p_values = []
    for pathway in significant_pathways:
        genes = genes_by_pathway[pathway]
        ks_p_values.append(kolmogorov_smirnov_test_with_ranking(genes, global_ranking))

    if not ks_p_values:
        logger.info("No significant pathways found after hypergeometric test. Skipping KS test.")
        return
    # Apply Benjamini-Hochberg correction to the KS P-values
    adjusted_p_values = multipletests(ks_p_values, method='fdr_bh')[1]

    # Filter significant pathways based on adjusted KS P-values
    settings.ks_significant_pathways_with_genes = {
        pathway: (genes_by_pathway[pathway], adjusted_p_values[i])
        for i, pathway in enumerate(significant_pathways)
        if adjusted_p_values[i] < 0.05
    }
    if not settings.ks_significant_pathways_with_genes:
        logger.info("No significant pathways found after KS test.")


def perform_statist_mann_whitney(settings: Settings, scores: dict):
    """
    Perform Mann-Whitney U test on pathways that passed the KS test and filter significant pathways.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - args (GeneralArgs): General arguments and settings.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    mw_p_values = []  # List to store Mann-Whitney p-values

    # Use filtered_genes for ranking and background scores
    filtered_scores = [scores[gene_id][0] for gene_id in settings.filtered_genes]

    # Rank the scores only for the filtered genes and reverse the ranks
    ranks = rankdata(filtered_scores)
    scores_rank = {gene_id: rank for gene_id, rank in zip(settings.filtered_genes, ranks)}

    # Iterate over pathways that passed the KS test to perform the Mann-Whitney U test
    for pathway, genes_info in settings.ks_significant_pathways_with_genes.items():
        # Extract the set of genes and other information from the tuple
        pathway_genes, p_value = genes_info  # Unpack the tuple

        # Ensure pathway_genes is a set of integers
        pathway_genes = set(pathway_genes)
        pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes]
        background_genes = settings.filtered_genes - pathway_genes

        pathway_ranks = [scores_rank[gene_id] for gene_id in pathway_genes]
        background_ranks = [scores_rank[gene_id] for gene_id in background_genes]

        # Compute the Mann-Whitney U test p-value using scores
        _, rmw_pval = compute_mw_python(pathway_ranks, background_ranks)
        mw_p_values.append(rmw_pval)

        # Collect the significant pathways with their corresponding details
        settings.filtered_pathways[pathway] = {
            'Rank': len(settings.filtered_pathways) + 1,
            'Term': pathway,
            'ES': np.mean(pathway_scores),
            'NES': None,
            'NOM p-val': rmw_pval,
            'FDR q-val': None,
            'FWER p-val': None,
            'Tag %': len(pathway_genes) / len(settings.filtered_genes) * 100,
            'Gene %': len(pathway_genes) / len(scores) * 100,
            'Lead_genes': ','.join(map(str, pathway_genes))
        }

    # Apply Benjamini-Hochberg correction to adjust the p-values
    adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]

    # Update the FDR q-values in the filtered pathways dictionary
    for i, (pathway, _) in enumerate(settings.filtered_pathways.items()):
        # Store the FDR q-val with full precision
        settings.filtered_pathways[pathway]['FDR q-val'] = decimal.Decimal(adjusted_mw_p_values[i])

def perform_enrichment(settings: Settings, scores_df):
    """
    Perform pathway enrichment analysis for a given test.

    Parameters:
    - settings (Settings): Settings object containing the experiment settings.
    - scores_df (pd.DataFrame): DataFrame containing the gene scores.

    Returns:
    - None
    """
    genes_by_pathway, scores = load_pathways_and_propagation_scores(settings, scores_df)
    if settings.run_gsea:
        # Prepare data for GSEA
        gene_ids = list(scores.keys())
        logfc_scores = [score[0] for score in scores.values()]
        gene_expression_data = pd.DataFrame({'gene': gene_ids, 'logFC': logfc_scores})
        gene_expression_data['gene'] = gene_expression_data['gene'].astype(str)
        gene_expression_data = gene_expression_data.sort_values(by='logFC', ascending=False)
        gsea_results = gp.prerank(rnk=gene_expression_data, gene_sets=genes_by_pathway, verbose=True, no_plot=True)
        gsea_results.res2d.to_excel(settings.output_path)

    else:
        perform_statist_ks(settings, genes_by_pathway, scores)
        if settings.ks_significant_pathways_with_genes:
            perform_statist_mann_whitney(settings, scores)
        else:
            logger.info("Skipping Mann-Whitney test.")

        if settings.filtered_pathways:
            pathways_df = pd.DataFrame.from_dict(settings.filtered_pathways, orient='index')

            def format_fdr_q_val(x):
                try:
                    return decimal.Decimal(x)
                except (ValueError, TypeError):
                    return x

            pathways_df['FDR q-val'] = pathways_df['FDR q-val'].apply(format_fdr_q_val)
            pathways_df.sort_values(by='FDR q-val', inplace=True)
            pathways_df['Rank'] = range(1, len(pathways_df) + 1)
            pathways_df.to_excel(settings.output_path,
                                 index=False)



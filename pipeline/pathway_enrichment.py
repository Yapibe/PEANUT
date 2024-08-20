import pandas as pd
from os import path, listdir
from args import EnrichTask, GeneralArgs
from scipy.stats import rankdata
import gseapy as gp
from utils import load_pathways_and_propagation_scores
from statsmodels.stats.multitest import multipletests
from visualization_tools import print_enriched_pathways_to_file
from statistical_methods import jaccard_index, kolmogorov_smirnov_test, compute_mw_python, run_hyper, \
    global_gene_ranking, kolmogorov_smirnov_test_with_ranking
import numpy as np
import logging


# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Create a logger
logger = logging.getLogger(__name__)
def perform_statist(task: EnrichTask, general_args, genes_by_pathway: dict, scores: dict):
    """
    Perform statistical enrichment analysis on pathways.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - general_args (GeneralArgs): General arguments and settings.
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    # Populate a set with all genes from the filtered pathways
    for genes in genes_by_pathway.values():
        task.filtered_genes.update(genes)

    # Create a global ranking of genes
    global_ranking = global_gene_ranking(scores)

    significant_pathways_hyper = list(genes_by_pathway.keys())

    # # Perform the Kolmogorov-Smirnov test to compare distributions of scores between pathway genes and background
    # ks_p_values = []
    # for pathway in significant_pathways_hyper:
    #     genes = genes_by_pathway[pathway]
    #     pathway_scores = [scores[gene_id][0] for gene_id in genes if gene_id in scores]
    #     background_genes = scores_keys - genes
    #     background_scores = [scores[gene_id][0] for gene_id in background_genes]
    #     ks_p_values.append(kolmogorov_smirnov_test(pathway_scores, background_scores))

    # Perform the Kolmogorov-Smirnov test using global ranking
    ks_p_values = []
    for pathway in significant_pathways_hyper:
        genes = genes_by_pathway[pathway]
        ks_p_values.append(kolmogorov_smirnov_test_with_ranking(genes, global_ranking))

    if not ks_p_values:
        logger.info("No significant pathways found after hypergeometric test. Skipping KS test.")
        return
    # Apply Benjamini-Hochberg correction to the KS P-values
    adjusted_p_values = multipletests(ks_p_values, method='fdr_bh')[1]

    # Filter significant pathways based on adjusted KS P-values
    task.ks_significant_pathways_with_genes = {
        pathway: (genes_by_pathway[pathway], adjusted_p_values[i])
        for i, pathway in enumerate(significant_pathways_hyper)
        if adjusted_p_values[i] < 0.05
    }
    if not task.ks_significant_pathways_with_genes:
        logger.info("No significant pathways found after KS test.")


def perform_statist_mann_whitney(task: EnrichTask, args, scores: dict):
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
    filtered_scores = [scores[gene_id][0] for gene_id in task.filtered_genes]

    # Rank the scores only for the filtered genes and reverse the ranks
    ranks = rankdata(filtered_scores)
    scores_rank = {
        gene_id: rank for gene_id, rank in zip(task.filtered_genes, ranks)
    }

    # Iterate over pathways that passed the KS test to perform the Mann-Whitney U test
    for pathway, genes_info in task.ks_significant_pathways_with_genes.items():
        pathway_genes = set(genes_info)
        pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes]
        background_genes = task.filtered_genes - pathway_genes

        pathway_ranks = [scores_rank[gene_id] for gene_id in pathway_genes]
        background_ranks = [scores_rank[gene_id] for gene_id in background_genes]

        # Compute the Mann-Whitney U test p-value using scores
        _, rmw_pval = compute_mw_python(pathway_ranks, background_ranks)
        mw_p_values.append(rmw_pval)

        # Collect the significant pathways with their corresponding details
        task.filtered_pathways[pathway] = {
            'Rank': len(task.filtered_pathways) + 1,
            'Name': pathway,
            'Term': pathway,
            'ES': np.mean(pathway_scores),
            'NES': None,
            'NOM p-val': None,
            'FDR q-val': None,
            'FWER p-val': None,
            'Tag %': len(pathway_genes) / len(task.filtered_genes) * 100,
            'Gene %': len(pathway_genes) / len(scores) * 100,
            'Lead_genes': ','.join(map(str, pathway_genes))
        }

    # Apply Benjamini-Hochberg correction to adjust the p-values
    adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]

    # Update the FDR q-values in the filtered pathways dictionary
    for i, (pathway, _) in enumerate(task.filtered_pathways.items()):
        task.filtered_pathways[pathway]['FDR q-val'] = adjusted_mw_p_values[i]


def perform_enrichment(test_name: str, general_args: GeneralArgs, output_path: str = None):
    """
    Perform pathway enrichment analysis for a given test.

    Parameters:
    - test_name (str): Name of the test for which enrichment analysis is performed.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """

    test_name = test_name.split('.')[0]
    # run enrichment
    propagation_folder = path.join(general_args.propagation_folder, test_name)
    if general_args.run_propagation:
        propagation_file = path.join(f'{propagation_folder}', '{}_{}'.format(general_args.alpha, general_args.date))
        enrich_task = EnrichTask(name=test_name, create_scores=True, target_field='gene_prop_scores',
                                statistic_test=kolmogorov_smirnov_test, propagation_file=propagation_file)
    else:
        files_in_dir = listdir(propagation_folder)

        if len(files_in_dir) == 1:
            propagation_file = path.join(propagation_folder, files_in_dir[0])
        else:
            raise ValueError(
                "Expected exactly one file in the propagation directory, but found {} files.".format(len(files_in_dir)))

        enrich_task = EnrichTask(name=test_name, create_scores=True, target_field='gene_prop_scores',
                                statistic_test=kolmogorov_smirnov_test, propagation_file=propagation_file)

    genes_by_pathway, scores = load_pathways_and_propagation_scores(general_args, enrich_task.propagation_file)

    if general_args.run_gsea:
        # Prepare data for GSEA
        # Unpack the scores dictionary into separate lists for GeneID and Score
        gene_ids = list(scores.keys())
        logfc_scores = [score[0] for score in scores.values()]
        # Create DataFrame for GSEA with string gene identifiers
        gene_expression_data = pd.DataFrame({'gene': gene_ids, 'logFC': logfc_scores})
        gene_expression_data['gene'] = gene_expression_data['gene'].astype(str)
        # Rank the data by logFC in descending order
        gene_expression_data = gene_expression_data.sort_values(by='logFC', ascending=False)
        # Run GSEA
        gsea_results = gp.prerank(rnk=gene_expression_data, gene_sets=genes_by_pathway, outdir=general_args.gsea_out,
                                  verbose=True, permutation_num=1000, no_plot=True)
        # save xlsx
        gsea_results.res2d.to_excel(output_path)

    else:
        # # Stage 1 - calculate nominal p-values and directions
        # perform_statist(enrich_task, general_args, genes_by_pathway, scores)
        # Check if there are significant pathways after the KS test
        enrich_task.ks_significant_pathways_with_genes = genes_by_pathway
        # Populate a set with all genes from the filtered pathways
        for genes in genes_by_pathway.values():
            enrich_task.filtered_genes.update(genes)
        if enrich_task.ks_significant_pathways_with_genes:
            # Further statistical test using Mann-Whitney U test
            perform_statist_mann_whitney(enrich_task, general_args, scores)
        else:
            logger.info("Skipping Mann-Whitney test.")

        # If there are filtered pathways, save them to Excel
        if enrich_task.filtered_pathways:
            # Convert the filtered pathways dictionary to a DataFrame
            pathways_df = pd.DataFrame.from_dict(enrich_task.filtered_pathways, orient='index')
            pathways_df = pathways_df[['Rank', 'Name', 'Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val',
                                       'FWER p-val', 'Tag %', 'Gene %', 'Lead_genes']]
            pathways_df.sort_values(by='FDR q-val', inplace=True)
            # Save the DataFrame to an Excel file
            pathways_df.to_excel(output_path, index=False)


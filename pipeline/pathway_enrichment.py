import pandas as pd
import numpy as np
import logging
from os import path, listdir
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from utils import load_pathways_and_propagation_scores
from args import EnrichTask, GeneralArgs
from statistical_methods import (
    kolmogorov_smirnov_test_with_ranking,
    compute_mw_python,
    global_gene_ranking,
    jaccard_index
)
import gseapy as gp

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def perform_ks_test(task: EnrichTask, genes_by_pathway: dict, scores: dict):
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
    task.filtered_genes = set()
    for genes in genes_by_pathway.values():
        task.filtered_genes.update(genes)

    # Create a global ranking of genes
    global_ranking = global_gene_ranking(scores)

    # Perform the Kolmogorov-Smirnov test for each pathway
    ks_p_values = []
    pathways = list(genes_by_pathway.keys())
    for pathway in pathways:
        genes = genes_by_pathway[pathway]
        p_value = kolmogorov_smirnov_test_with_ranking(genes, global_ranking)
        ks_p_values.append(p_value)

    if not ks_p_values:
        logger.info("No pathways to test. Skipping KS test.")
        return
    
    # Apply Benjamini-Hochberg correction to the KS P-values
    adjusted_p_values = multipletests(ks_p_values, method='fdr_bh')[1]


    # Filter significant pathways based on adjusted KS P-values
    task.ks_significant_pathways_with_genes = {}
    for i, pathway in enumerate(pathways):
        if adjusted_p_values[i] < 0.05:
            task.ks_significant_pathways_with_genes[pathway] = {
                'genes': genes_by_pathway[pathway],
                'ks_p_value': adjusted_p_values[i]
            }

    if not task.ks_significant_pathways_with_genes:
        logger.info("No significant pathways found after KS test.")


def remove_similar_pathways(filtered_pathways: dict, jac_threshold: float) -> dict:
    """
    Remove pathways that are too similar based on Jaccard index.

    Parameters:
    - filtered_pathways (dict): Dictionary of pathways with their details.
    - jac_threshold (float): Jaccard index threshold for removing similar pathways.

    Returns:
    - dict: Filtered pathways after removing similar ones.
    """
    keep_list = [
        'KEGG_THYROID_CANCER',
        'KEGG_NON_SMALL_CELL_LUNG_CANCER',
        'KEGG_ACUTE_MYELOID_LEUKEMIA',
        'KEGG_COLORECTAL_CANCER',
        'KEGG_GLIOMA',
        'KEGG_RENAL_CELL_CARCINOMA',
        'KEGG_PANCREATIC_CANCER',
        'KEGG_PROSTATE_CANCER',
        'KEGG_DILATED_CARDIOMYOPATHY',
        'KEGG_PARKINSONS_DISEASE',
        'KEGG_ALZHEIMERS_DISEASE',
        'KEGG_HUNTINGTONS_DISEASE'
    ]
    # Sort pathways by FDR q-value (ascending)
    sorted_pathways = sorted(filtered_pathways.items(), key=lambda x: x[1]['FDR q-val'])
    pathways_genes = {pathway: set(data['Lead_genes'].split(',')) for pathway, data in sorted_pathways}
    removed_pathways = set()

    for i, (pathway_i, data_i) in enumerate(sorted_pathways):
        if pathway_i in removed_pathways:
            continue
        genes_i = pathways_genes[pathway_i]

        # Compare pathway_i with all subsequent pathways in the list
        for j in range(i + 1, len(sorted_pathways)):
            pathway_j, data_j = sorted_pathways[j]
            if pathway_j in removed_pathways or pathway_j in keep_list:
                # Skip removing pathway_j if it's already removed or in the keep list
                continue
            
            genes_j = pathways_genes[pathway_j]
            ji = jaccard_index(genes_i, genes_j)
            
            if ji > jac_threshold:
                # Remove pathway_j (less significant) only if it's not in the keep list
                if pathway_j not in keep_list:
                    removed_pathways.add(pathway_j)

    # Construct the filtered pathways after removing similar ones, always keeping pathways in the keep list
    filtered_pathways_after_jaccard = {
        pathway: data
        for pathway, data in filtered_pathways.items()
        if pathway not in removed_pathways or pathway in keep_list
    }

    return filtered_pathways_after_jaccard


def perform_mann_whitney_test(task: EnrichTask, args: GeneralArgs, scores: dict):
    """
    Perform Mann-Whitney U test on pathways that passed the KS test and filter significant pathways.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - args (GeneralArgs): General arguments and settings.
    - scores (dict): Mapping of gene IDs to their scores and p-values.

    Returns:
    - None
    """
    if not task.ks_significant_pathways_with_genes:
        logger.info("No significant pathways to test with Mann-Whitney U test.")
        return

    # Use filtered_genes for ranking and background scores
    filtered_scores = [scores[gene_id][0] for gene_id in task.filtered_genes]

    # Rank the scores for the filtered genes
    ranks = rankdata(filtered_scores)
    scores_rank = {gene_id: rank for gene_id, rank in zip(task.filtered_genes, ranks)}

    mw_p_values = []  # List to store Mann-Whitney p-values

    # Initialize filtered pathways
    task.filtered_pathways = {}

    # Iterate over pathways that passed the KS test to perform the Mann-Whitney U test
    for pathway, data in task.ks_significant_pathways_with_genes.items():
        pathway_genes = set(data['genes'])
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
            'FDR q-val': None,  # Placeholder, to be filled after correction
            'FWER p-val': None,
            'Tag %': len(pathway_genes) / len(task.filtered_genes) * 100,
            'Gene %': len(pathway_genes) / len(scores) * 100,
            'Lead_genes': ','.join(map(str, pathway_genes))
        }

    if not mw_p_values:
        logger.info("No pathways to test. Skipping FDR correction.")
        return

    # Apply Benjamini-Hochberg correction to adjust the p-values
    adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]

    # Update the FDR q-values in the filtered pathways dictionary
    for i, pathway in enumerate(task.filtered_pathways.keys()):
        task.filtered_pathways[pathway]['FDR q-val'] = float(adjusted_mw_p_values[i])

    # Remove similar pathways based on Jaccard index threshold
    task.filtered_pathways = remove_similar_pathways(task.filtered_pathways, args.JAC_THRESHOLD)


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

    # Set up propagation file path
    propagation_folder = path.join(general_args.propagation_folder, test_name)
    if general_args.run_propagation:
        propagation_file = path.join(propagation_folder, f"{general_args.alpha}_{general_args.date}")
    else:
        files_in_dir = listdir(propagation_folder)
        if len(files_in_dir) == 1:
            propagation_file = path.join(propagation_folder, files_in_dir[0])
        else:
            raise ValueError(f"Expected exactly one file in the propagation directory, but found {len(files_in_dir)} files.")

    # Create EnrichTask
    enrich_task = EnrichTask(
        name=test_name,
        create_scores=True,
        target_field='gene_prop_scores',
        statistic_test=kolmogorov_smirnov_test_with_ranking,
        propagation_file=propagation_file
    )

    genes_by_pathway, scores = load_pathways_and_propagation_scores(general_args, enrich_task.propagation_file)
    logger.info(f"Running enrichment analysis for test: {test_name}")

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
        # Perform Kolmogorov-Smirnov test
        perform_ks_test(enrich_task, genes_by_pathway, scores)
        if enrich_task.ks_significant_pathways_with_genes:
            # Perform Mann-Whitney U test
            perform_mann_whitney_test(enrich_task, general_args, scores)
        else:
            logger.info("No significant pathways after KS test. Skipping Mann-Whitney U test.")

        # If there are filtered pathways, save them to Excel
        if enrich_task.filtered_pathways:
            # Convert the filtered pathways dictionary to a DataFrame
            pathways_df = pd.DataFrame.from_dict(enrich_task.filtered_pathways, orient='index')
            # Ensure that FDR q-val is stored as float
            pathways_df['FDR q-val'] = pathways_df['FDR q-val'].astype(float)
            # Sort by FDR q-val
            pathways_df.sort_values(by='FDR q-val', inplace=True)
            # Save to Excel
            pathways_df.to_excel(output_path, index=False)
        else:
            logger.info("No significant pathways after Mann-Whitney U test.")



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
    compute_mw,
    global_gene_ranking,
    jaccard_index
)
import gseapy as gp
# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def perform_permutation_test(filtered_pathways: dict, genes_by_pathway: dict, scores: dict, n_iterations: int = 1000, pval_threshold: float = 0.05) -> dict:
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
    if not filtered_pathways:
        logger.info("No pathways to test with the permutation test.")
        return {}
    
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
            logger.warning(f"No valid genes found for pathway '{pathway}'. Skipping permutation test.")
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
        logger.info("No valid p-values obtained from permutation test.")
        # Set permutation_significant to False for all pathways
        for data in filtered_pathways.values():
            data['permutation_significant'] = False
    
    return filtered_pathways


def remove_similar_pathways(filtered_pathways: dict, jac_threshold: float, genes_by_pathway: dict, associated_pathway_name: str) -> dict:
    """
    Identify pathways that are too similar based on the Jaccard index, adding a flag to indicate removal.
    
    Parameters:
    - filtered_pathways (dict): Dictionary of pathways with their details.
    - genes_by_pathway (dict): Mapping of pathway names to their gene sets.
    - jac_threshold (float): Jaccard index threshold for considering pathways similar.
    
    Returns:
    - dict: Updated filtered_pathways with Jaccard similarity flags.
    """
    # Initialize the Jaccard flag for all pathways
    for pathway in filtered_pathways:
        filtered_pathways[pathway]['jaccard_removed'] = False

    # Sort pathways by FDR q-value (ascending)
    sorted_pathways = sorted(filtered_pathways.keys(), key=lambda p: filtered_pathways[p]['FDR q-val'])

    # Iterate over pathways based on FDR ranking
    for pathway_i in sorted_pathways:
        if filtered_pathways[pathway_i]['jaccard_removed']:
            continue  # Skip pathways already flagged for removal

        genes_i = genes_by_pathway[pathway_i]

        for pathway_j in sorted_pathways:
            if pathway_j == pathway_i:
                continue  # Skip self-comparison
            if filtered_pathways[pathway_j]['jaccard_removed']:
                continue  # Skip already removed pathways

            genes_j = genes_by_pathway[pathway_j]
            ji = jaccard_index(genes_i, genes_j)

            if ji > jac_threshold:
                # Remove the pathway with higher FDR q-val
                pval_i = filtered_pathways[pathway_i]['FDR q-val']
                pval_j = filtered_pathways[pathway_j]['FDR q-val']

                if pval_i <= pval_j:
                    filtered_pathways[pathway_j]['jaccard_removed'] = True
                else:
                    filtered_pathways[pathway_i]['jaccard_removed'] = True
                    break  # No need to compare pathway_i with others

    return filtered_pathways



def perform_ks_test(task: EnrichTask, genes_by_pathway: dict, scores: dict, test_name: str):
    """
    Perform statistical enrichment analysis on pathways using the Kolmogorov-Smirnov test.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores (dict): Mapping of gene IDs to their scores.
    - prior_set (pd.DataFrame): Prior set of genes (not used in this context but included for consistency).

    Returns:
    - None
    """
    # Extract pathway name from test name
    test_name = '_'.join(test_name.split('_')[1:])
    # Populate a set with all genes from the pathways
    task.filtered_genes = set()
    for genes in genes_by_pathway.values():
        task.filtered_genes.update(genes)

    # Create a global ranking of genes using rankings
    global_ranking = global_gene_ranking(scores)

    # Perform the Kolmogorov-Smirnov test for each pathway
    ks_p_values = []
    pathways = list(genes_by_pathway.keys())
    for pathway in pathways:
        genes = genes_by_pathway[pathway]
        p_value = kolmogorov_smirnov_test_with_ranking(genes, global_ranking, alternative='greater')
        ks_p_values.append(p_value)
    if not ks_p_values:
        logger.info("No pathways to test. Skipping KS test.")
        return

    # Apply Benjamini-Hochberg correction to the KS P-values
    adjusted_p_values = multipletests(ks_p_values, method='fdr_bh')[1]

    # Store all pathways with their adjusted KS p-values and significance status
    task.pathways_with_ks_pvalues = {}
    for i, pathway in enumerate(pathways):
        adjusted_p_value = adjusted_p_values[i]
        is_significant = adjusted_p_value < 0.05
        task.pathways_with_ks_pvalues[pathway] = {
            'genes': genes_by_pathway[pathway],
            'ks_p_value': adjusted_p_value,
            'ks_significant': is_significant
        }

    significant_pathways = sum(data['ks_significant'] for data in task.pathways_with_ks_pvalues.values())
    if significant_pathways == 0:
        logger.info("No significant pathways found after KS test.")
    else:
        logger.info(f"Found {significant_pathways} significant pathways after KS test.")


def perform_mann_whitney_test(task: EnrichTask, args: GeneralArgs, scores: dict, genes_by_pathway: dict):
    """
    Perform Mann-Whitney U test on all pathways, followed by permutation testing and Jaccard similarity flagging.
    
    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - args (GeneralArgs): General arguments and settings.
    - scores (dict): Mapping of gene IDs to their scores.
    
    Returns:
    - None
    """
    if not task.pathways_with_ks_pvalues:
        logger.info("No pathways to test with Mann-Whitney U test.")
        return
    
    # Collect all genes from the pathways
    task.filtered_genes = set()
    for data in task.pathways_with_ks_pvalues.values():
        task.filtered_genes.update(data['genes'])
    
    # Rank the scores for the filtered genes
    filtered_scores = [scores[gene_id][0] for gene_id in task.filtered_genes if gene_id in scores]
    ranks = rankdata(filtered_scores)
    scores_rank = {gene_id: rank for gene_id, rank in zip(task.filtered_genes, ranks)}
    
    # Initialize lists to store MW p-values and pathway names
    mw_p_values = []
    pathways = list(task.pathways_with_ks_pvalues.keys())
    
    # Perform the Mann-Whitney U test
    for pathway in pathways:
        data = task.pathways_with_ks_pvalues[pathway]
        pathway_genes = set(data['genes'])
        pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes if gene_id in scores]
        background_genes = task.filtered_genes - pathway_genes
    
        pathway_ranks = [scores_rank[gene_id] for gene_id in pathway_genes if gene_id in scores_rank]
        background_ranks = [scores_rank[gene_id] for gene_id in background_genes if gene_id in scores_rank]
    
        _, mw_pval = compute_mw(pathway_ranks, background_ranks)
        mw_p_values.append(mw_pval)
    
        sorted_genes = sorted(pathway_genes, key=lambda gene: scores[gene][0], reverse=True)
    
        # Collect pathway details
        data.update({
            'Name': pathway,
            'mw_p_value': mw_pval,
            'ES': np.mean(pathway_scores) if pathway_scores else None,
            'Lead_genes': ','.join(map(str, sorted_genes))
        })
    
    if not mw_p_values:
        logger.info("No pathways to test. Skipping FDR correction.")
        return
    
    # Apply Benjamini-Hochberg correction to adjust the MW p-values
    adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]
    for i, pathway in enumerate(pathways):
        data = task.pathways_with_ks_pvalues[pathway]
        data['FDR q-val'] = adjusted_mw_p_values[i]
        data['mw_significant'] = adjusted_mw_p_values[i] < 0.05
    
    # Update filtered pathways with all pathways (both significant and non-significant)
    task.filtered_pathways = task.pathways_with_ks_pvalues.copy()
    
    # Perform Permutation Test
    task.filtered_pathways = perform_permutation_test(
        task.filtered_pathways,
        genes_by_pathway,
        scores
    )
    
    # Apply Jaccard similarity flagging
    task.filtered_pathways = remove_similar_pathways(
        task.filtered_pathways,
        args.JAC_THRESHOLD,
        genes_by_pathway,
        associated_pathway_name = task.name.split('_', 1)[1]
    )



def perform_enrichment(test_name: str, general_args: GeneralArgs, output_path: str = None):
    """
    Perform pathway enrichment analysis for a given test.
    
    Parameters:
    - test_name (str): Name of the test for which enrichment analysis is performed.
    - general_args (GeneralArgs): General arguments and settings.
    - output_path (str): Path to save the output.
    - related_pathways (list): List of related pathways for Jaccard filtering.
    
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
    
    
    logger.info(f"Running enrichment analysis for test: {test_name} with {len(genes_by_pathway)} pathways")
    
    if general_args.run_gsea:
        # Prepare data for GSEA
        gene_ids = list(scores.keys())
        logfc_scores = [score[0] for score in scores.values()]
        # Count and log the number of genes with negative logFC scores
        negative_genes_count = sum(1 for score in logfc_scores if score < 0)
        gene_expression_data = pd.DataFrame({'gene': gene_ids, 'logFC': logfc_scores})
        gene_expression_data['gene'] = gene_expression_data['gene'].astype(str)
        gene_expression_data = gene_expression_data.sort_values(by='logFC', ascending=False)
        gsea_results = gp.prerank(rnk=gene_expression_data, gene_sets=genes_by_pathway, outdir=general_args.gsea_out,
                                  verbose=False, permutation_num=1000, no_plot=True)
        # Normalize the GSEA output
        res_df = gsea_results.res2d.copy()
        res_df.reset_index(inplace=True)
        res_df.drop('Name', axis=1, errors='ignore', inplace=True)  # Remove Name column if it exists
        res_df.rename(columns={
            "index": "Rank",   # Assuming the rank is in the index; adjust if necessary
            "Term": "Name"     # Rename Term to Name
        }, inplace=True)
        res_df['Rank'] = res_df['Rank'] + 1  # Add 1 to make ranks start at 1
        res_df.to_excel(output_path)
    else:
        # Perform Kolmogorov-Smirnov test
        perform_ks_test(enrich_task, genes_by_pathway, scores, test_name)
        # Perform Mann-Whitney U test on all pathways
        perform_mann_whitney_test(enrich_task, general_args, scores, genes_by_pathway)
    
        # If there are filtered pathways, save them to Excel
        if enrich_task.filtered_pathways:
            # Convert the filtered pathways dictionary to a DataFrame
            pathways_df = pd.DataFrame.from_dict(enrich_task.filtered_pathways, orient='index')
            # Ensure correct data types
            pathways_df['mw_p_value'] = pathways_df['mw_p_value'].astype(float)
            pathways_df['FDR q-val'] = pathways_df['FDR q-val'].astype(float)
            pathways_df['ks_significant'] = pathways_df['ks_significant'].astype(bool)
            pathways_df['mw_significant'] = pathways_df['mw_significant'].astype(bool)
            pathways_df['permutation_significant'] = pathways_df['permutation_significant'].astype(bool)
            pathways_df['jaccard_removed'] = pathways_df['jaccard_removed'].astype(bool)
            pathways_df['Permutation p-val'] = pathways_df.get('Permutation p-val', np.nan).astype(float)
    
            # Define grouping criteria
            def pathway_group(row):
                if row['ks_significant'] and row['mw_significant'] and row['permutation_significant'] and not row['jaccard_removed']:
                    return 1  # Passed all tests
                elif row['ks_significant'] and row['mw_significant'] and row['permutation_significant'] and row['jaccard_removed']:
                    return 2  # Removed by Jaccard
                elif row['ks_significant'] and row['mw_significant'] and not row['permutation_significant']:
                    return 3  # Failed permutation test
                elif row['ks_significant'] and not row['mw_significant']:
                    return 4  # Failed MW test
                else:
                    return 5  # Did not pass KS test
    
            pathways_df['Group'] = pathways_df.apply(pathway_group, axis=1)
    
            # Sort pathways by group and FDR q-val
            pathways_df.sort_values(by=['Group', 'FDR q-val'], ascending=[True, True], inplace=True)
            pathways_df.reset_index(drop=True, inplace=True)
            pathways_df['Rank'] = pathways_df.index + 1
    
            # Reorder columns to put Name and FDR q-val first
            cols = ['Rank', 'Name', 'Group', 'FDR q-val', 'ks_significant', 'mw_significant', 'permutation_significant', 'jaccard_removed', 'Permutation p-val' ] + [col for col in pathways_df.columns if col not in ['Rank', 'Name', 'FDR q-val', 'Group', 'ks_significant', 'mw_significant', 'permutation_significant', 'jaccard_removed', 'Permutation p-val']]
            pathways_df = pathways_df[cols]
    
            # Save to Excel
            pathways_df.to_excel(output_path, index=False)
        else:
            logger.info("No pathways to output.")
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
    jaccard_index,
    kolmogorov_smirnov_test_with_scores
)
import gseapy as gp
# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)



from statsmodels.stats.multitest import multipletests

import os
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

def plot_score_distributions_comparison(
    pathway_prior_scores,
    background_prior_scores,
    pathway_scores,
    background_scores,
    pathway_name,
    output_dir
):
    # Validate input scores
    if not pathway_prior_scores or not background_prior_scores:
        print(f"Skipping plot for {pathway_name}: Missing prior scores.")
        return
    if not pathway_scores or not background_scores:
        print(f"Skipping plot for {pathway_name}: Missing propagated scores.")
        return

    # Plot histograms
    plt.figure(figsize=(12, 6))

    # Before propagation
    plt.subplot(1, 2, 1)
    plt.hist(background_prior_scores, bins=50, alpha=0.5, label='Background Genes', density=True)
    plt.hist(pathway_prior_scores, bins=50, alpha=0.7, label=f'{pathway_name} Genes', density=True)
    plt.xlabel('Gene Scores (Before Propagation)')
    plt.ylabel('Density')
    plt.title(f'Score Distribution Before Propagation: {pathway_name}')
    plt.legend()

    # After propagation
    plt.subplot(1, 2, 2)
    plt.hist(background_scores, bins=50, alpha=0.5, label='Background Genes', density=True)
    plt.hist(pathway_scores, bins=50, alpha=0.7, label=f'{pathway_name} Genes', density=True)
    plt.xlabel('Gene Scores (After Propagation)')
    plt.ylabel('Density')
    plt.title(f'Score Distribution After Propagation: {pathway_name}')
    plt.legend()

    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'score_distribution_comparison_{pathway_name}.png')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_ecdfs_with_scores_comparison(
    pathway_prior_scores,
    background_prior_scores,
    pathway_scores,
    background_scores,
    pathway_name,
    output_dir
):
    import numpy as np
    import matplotlib.pyplot as plt

    # Sort scores
    pathway_prior_scores_sorted = np.sort(pathway_prior_scores)
    background_prior_scores_sorted = np.sort(background_prior_scores)
    pathway_scores_sorted = np.sort(pathway_scores)
    background_scores_sorted = np.sort(background_scores)

    # Calculate ECDFs
    pathway_prior_ecdf = np.arange(1, len(pathway_prior_scores_sorted)+1) / len(pathway_prior_scores_sorted)
    background_prior_ecdf = np.arange(1, len(background_prior_scores_sorted)+1) / len(background_prior_scores_sorted)
    pathway_ecdf = np.arange(1, len(pathway_scores_sorted)+1) / len(pathway_scores_sorted)
    background_ecdf = np.arange(1, len(background_scores_sorted)+1) / len(background_scores_sorted)

    # Plot ECDFs
    plt.figure(figsize=(12, 6))

    # Before propagation
    plt.subplot(1, 2, 1)
    plt.step(background_prior_scores_sorted, background_prior_ecdf, where='post', label='Background Genes')
    plt.step(pathway_prior_scores_sorted, pathway_prior_ecdf, where='post', label=f'{pathway_name} Genes')
    plt.xlabel('Gene Scores (Before Propagation)')
    plt.ylabel('ECDF')
    plt.title(f'ECDF Before Propagation: {pathway_name}')
    plt.legend()
    plt.grid(True)

    # After propagation
    plt.subplot(1, 2, 2)
    plt.step(background_scores_sorted, background_ecdf, where='post', label='Background Genes')
    plt.step(pathway_scores_sorted, pathway_ecdf, where='post', label=f'{pathway_name} Genes')
    plt.xlabel('Gene Scores (After Propagation)')
    plt.ylabel('ECDF')
    plt.title(f'ECDF After Propagation: {pathway_name}')
    plt.legend()
    plt.grid(True)

    # Save the plot
    output_file = os.path.join(output_dir, f'ecdf_comparison_{pathway_name}.png')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


def plot_score_distributions(pathway_scores, background_scores, pathway_name, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Plot histograms
    plt.figure(figsize=(10, 6))
    plt.hist(background_scores, bins=50, alpha=0.5, label='Background Genes', density=True)
    plt.hist(pathway_scores, bins=50, alpha=0.7, label=f'{pathway_name} Genes', density=True)
    plt.xlabel('Gene Scores')
    plt.ylabel('Density')
    plt.title(f'Score Distribution: {pathway_name} vs. Background')
    plt.legend()

    # Save the plot
    output_file = os.path.join(output_dir, f'score_distribution_{pathway_name}.png')
    plt.savefig(output_file)
    plt.close()


def plot_ecdfs_with_scores(pathway_scores, background_scores, pathway_name, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert sets to lists
    pathway_scores = list(pathway_scores)
    background_scores = list(background_scores)
    
    # Sort scores
    pathway_scores_sorted = np.sort(pathway_scores)
    background_scores_sorted = np.sort(background_scores)
    
    # Calculate ECDFs
    pathway_ecdf = np.arange(1, len(pathway_scores_sorted)+1) / len(pathway_scores_sorted)
    background_ecdf = np.arange(1, len(background_scores_sorted)+1) / len(background_scores_sorted)

    # Plot ECDFs
    plt.figure(figsize=(10, 6))
    plt.step(pathway_scores_sorted, pathway_ecdf, where='post', label=f'{pathway_name} Genes')
    plt.step(background_scores_sorted, background_ecdf, where='post', label='Background Genes')
    plt.xlabel('Gene Scores')
    plt.ylabel('ECDF')
    plt.title(f'ECDF Comparison: {pathway_name} vs. Background')
    plt.legend()
    plt.grid(True)

    # Save the plot
    output_file = os.path.join(output_dir, f'ecdf_{pathway_name}.png')
    plt.savefig(output_file)
    plt.close()

def plot_ecdfs(pathway_genes, global_ranking, pathway_name, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Filter pathway genes present in global ranking
    pathway_genes_in_ranking = [gene for gene in pathway_genes if gene in global_ranking.index]
    if not pathway_genes_in_ranking:
        print(f"No genes from {pathway_name} are in the global ranking.")
        return

    # Get ranks
    pathway_ranks = global_ranking[pathway_genes_in_ranking].values
    background_genes = global_ranking.index.difference(pathway_genes_in_ranking)
    background_ranks = global_ranking[background_genes].values

    # Sort ranks
    pathway_ranks_sorted = np.sort(pathway_ranks)
    background_ranks_sorted = np.sort(background_ranks)

    # Calculate ECDFs
    pathway_ecdf = np.arange(1, len(pathway_ranks_sorted)+1) / len(pathway_ranks_sorted)
    background_ecdf = np.arange(1, len(background_ranks_sorted)+1) / len(background_ranks_sorted)

    # Create DataFrames for the ECDF data
    pathway_ecdf_df = pd.DataFrame({
        'Rank': pathway_ranks_sorted,
        'ECDF': pathway_ecdf
    })
    background_ecdf_df = pd.DataFrame({
        'Rank': background_ranks_sorted,
        'ECDF': background_ecdf
    })

    # Save ECDF data to CSV files
    pathway_ecdf_file = os.path.join(output_dir, f'ecdf_data_{pathway_name}_genes.csv')
    background_ecdf_file = os.path.join(output_dir, f'ecdf_data_{pathway_name}_background.csv')
    pathway_ecdf_df.to_csv(pathway_ecdf_file, index=False)
    background_ecdf_df.to_csv(background_ecdf_file, index=False)
    # # Optionally, print the first few rows of the ECDF data
    # print(f"ECDF data for {pathway_name} genes:")
    # print(pathway_ecdf_df.head())
    # print(f"\nECDF data for {pathway_name} background genes:")
    # print(background_ecdf_df.head())

    # Plot ECDFs
    plt.figure(figsize=(10, 6))
    plt.step(pathway_ranks_sorted, pathway_ecdf, where='post', label=f'{pathway_name} Genes')
    plt.step(background_ranks_sorted, background_ecdf, where='post', label='Background Genes')
    plt.xlabel('Gene Rank')
    plt.ylabel('ECDF')
    plt.title(f'ECDF Comparison: {pathway_name} vs. Background')
    plt.legend()
    plt.grid(True)

    # Save the plot
    output_file = os.path.join(output_dir, f'ecdf_{pathway_name}.png')
    plt.savefig(output_file)
    plt.close()

def plot_boxplots(pathway_genes, scores, pathway_name, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Get pathway and background scores
    pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes if gene_id in scores]
    background_genes = set(scores.keys()) - set(pathway_genes)
    background_scores = [scores[gene_id][0] for gene_id in background_genes]

    # Plot boxplots
    data = [pathway_scores, background_scores]
    labels = [f'{pathway_name} Genes', 'Background Genes']
    plt.figure(figsize=(8, 6))
    plt.boxplot(data, labels=labels)
    plt.ylabel('Gene Scores')
    plt.title(f'Boxplot of Scores: {pathway_name} vs. Background')

    # Save the plot
    output_file = os.path.join(output_dir, f'boxplot_{pathway_name}.png')
    plt.savefig(output_file)
    plt.close()


def perform_permutation_test(filtered_pathways: dict, scores: dict, associated_pathway_name: str, n_iterations: int = 1000, pval_threshold: float = 0.05) -> dict:
    """
    Perform a permutation test to assess the significance of the observed pathway scores.
    Updates pathways with adjusted permutation p-values and removes pathways that do not pass the p-value threshold.

    Parameters:
    - filtered_pathways (dict): Dictionary of pathways with their details.
    - scores (dict): Mapping of gene IDs to their scores.
    - n_iterations (int): Number of permutations to perform.
    - pval_threshold (float): Threshold for adjusted p-values to filter pathways.

    Returns:
    - dict: Updated filtered_pathways with permutation p-values set for each pathway,
            and pathways not passing the p-value threshold removed.
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
        try:
            pathway_genes = set(map(int, data['Lead_genes'].split(',')))
        except ValueError:
            logger.warning(f"Invalid gene IDs in pathway '{pathway}'. Skipping permutation test.")
            data['Permutation p-val'] = np.nan
            continue

        # Find the intersection of pathway genes with available scores
        valid_genes = pathway_genes.intersection(scores.keys())
        pathway_size = len(valid_genes)

        if pathway_size == 0:
            logger.warning(f"No valid genes found for pathway '{pathway}'. Skipping permutation test.")
            data['Permutation p-val'] = np.nan
            continue

        # Calculate the observed mean score for the pathway
        observed_mean = np.mean([scores[gene_id][0] for gene_id in valid_genes])

        # Generate null distribution for this pathway size
        permuted_means = np.random.choice(score_values, size=(n_iterations, pathway_size), replace=True).mean(axis=1)
        empirical_p_val = np.mean(permuted_means >= observed_mean)

        # Collect data for multiple testing correction
        observed_means.append(observed_mean)
        raw_p_values.append(empirical_p_val)
        pathways_list.append(pathway)
        valid_pathways.append((pathway, data))

    # Adjust p-values for multiple testing using Benjamini-Hochberg
    if raw_p_values:
        adjusted_pvals = multipletests(raw_p_values, method='fdr_bh')[1]

        # Update pathways with adjusted p-values and remove those that don't pass the threshold
        filtered_pathways_after_permutation = {}
        for i, (pathway, data) in enumerate(valid_pathways):
            adjusted_pval = adjusted_pvals[i]
            data['Permutation p-val'] = adjusted_pval

            # Keep the pathway if it passes the threshold or if it is the associated pathway
            if adjusted_pval <= pval_threshold or pathway == associated_pathway_name:
                filtered_pathways_after_permutation[pathway] = data
            else:
                if pathway == associated_pathway_name:
                    logger.info(f"Pathway '{pathway}' removed due to adjusted permutation p-value {adjusted_pval:.4f} > {pval_threshold}")

        return filtered_pathways_after_permutation
    else:
        logger.info("No valid p-values obtained from permutation test.")
        return {}



def perform_ks_test(task: EnrichTask, genes_by_pathway: dict, scores: dict, output_path: str, test_name: str, prior_set: pd.DataFrame):
    """
    Perform statistical enrichment analysis on pathways.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - genes_by_pathway (dict): Mapping of pathways to their constituent genes.
    - scores (dict): Mapping of gene IDs to their scores.

    Returns:
    - None
    """
    output_dir = path.join(output_path, 'Plots', test_name)
    # Convert prior_set DataFrame to dict format matching scores
    prior_set_dict = {row['GeneID']: (row['Score'], 0) for _, row in prior_set.iterrows()}
    filtered_prior_set_dict = {gene_id: score for gene_id, score in prior_set_dict.items() if gene_id in scores}
    # Extract pathway name from test name by removing the first part (e.g., GSE21354_)
    test_name = '_'.join(test_name.split('_')[1:])
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
        # Get pathway and background scores
        pathway_scores = [scores[gene_id][0] for gene_id in genes if gene_id in scores]
        background_genes = set(scores.keys()) - set(genes)
        background_scores = [scores[gene_id][0] for gene_id in background_genes]

        pathway_prior_scores = [filtered_prior_set_dict[gene_id][0] for gene_id in genes if gene_id in filtered_prior_set_dict]
        background_prior_genes = set(filtered_prior_set_dict.keys()) - set(genes)
        background_prior_scores = [filtered_prior_set_dict[gene_id][0] for gene_id in background_prior_genes]


        if pathway == test_name:
            # Plot and save distributions
            plot_score_distributions_comparison(
                pathway_prior_scores,
                background_prior_scores,
                pathway_scores,
                background_scores,
                pathway,
                output_dir
            )

            # Plot ECDFs before and after propagation
            plot_ecdfs_with_scores_comparison(
                pathway_prior_scores,
                background_prior_scores,
                pathway_scores,
                background_scores,
                pathway,
                output_dir
            )
        # p_value = kolmogorov_smirnov_test_with_ranking(genes, global_ranking, alternative='greater')
        p_value = kolmogorov_smirnov_test_with_scores(genes, prior_set_dict)
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




def remove_similar_pathways(filtered_pathways: dict, jac_threshold: float, associated_pathway_name: str, related_pathways: list) -> dict:
    """
    Remove pathways that are too similar based on the Jaccard index, retaining only the associated pathway.

    Parameters:
    - filtered_pathways (dict): Dictionary of pathways with their details.
    - jac_threshold (float): Jaccard index threshold for removing similar pathways.
    - associated_pathway_name (str): The specific pathway to retain regardless of similarity.

    Returns:
    - dict: Filtered pathways after removing similar ones, keeping the associated pathway.
    """
    # Sort pathways by FDR q-value in ascending order to prioritize significance
    sorted_pathways = sorted(
        filtered_pathways.items(),
        key=lambda x: x[1]['FDR q-val']
    )

    # Map each pathway to its set of associated genes for efficient Jaccard computation
    pathways_genes = {
        pathway: set(data['Lead_genes'].split(','))
        for pathway, data in sorted_pathways
    }

    removed_pathways = set()

    # Iterate through each pathway to compare with others
    for i, (pathway_i, data_i) in enumerate(sorted_pathways):
        if pathway_i in removed_pathways:
            continue
        genes_i = pathways_genes[pathway_i]

        # Compare the current pathway with all subsequent pathways
        for j in range(i + 1, len(sorted_pathways)):
            pathway_j, data_j = sorted_pathways[j]

            # Skip if pathway_j is already removed or is the associated pathway to retain
            if pathway_j in removed_pathways or pathway_j == associated_pathway_name:
                continue

            genes_j = pathways_genes[pathway_j]
            ji = jaccard_index(genes_i, genes_j)

            if ji > jac_threshold:
                # Remove pathway_j if it's not the associated pathway
                removed_pathways.add(pathway_j)

    # Always retain the associated pathway
    filtered_pathways_after_jaccard = {
        pathway: data
        for pathway, data in filtered_pathways.items()
        if pathway not in removed_pathways or pathway == associated_pathway_name
    }

    return filtered_pathways_after_jaccard


def perform_mann_whitney_test(task: EnrichTask, args: GeneralArgs, scores: dict, related_pathways: list):
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
        _, rmw_pval = compute_mw(pathway_ranks, background_ranks)
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

    associated_pathway_name = task.name.split('_', 1)[1]

    # Remove similar pathways based on Jaccard index threshold
    task.filtered_pathways = remove_similar_pathways(task.filtered_pathways, args.JAC_THRESHOLD, associated_pathway_name, related_pathways)

    # Perform Permutation Test
    task.filtered_pathways = perform_permutation_test(task.filtered_pathways, scores, associated_pathway_name)




def perform_enrichment(test_name: str, general_args: GeneralArgs, output_path: str = None, related_pathways: list = None):
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

    genes_by_pathway, scores, prior_set = load_pathways_and_propagation_scores(general_args, enrich_task.propagation_file)

    logger.info(f"Running enrichment analysis for test: {test_name}")

    if general_args.run_gsea:
        # Prepare data for GSEA
        gene_ids = list(scores.keys())
        logfc_scores = [score[0] for score in scores.values()]
        gene_expression_data = pd.DataFrame({'gene': gene_ids, 'logFC': logfc_scores})
        gene_expression_data['gene'] = gene_expression_data['gene'].astype(str)
        gene_expression_data = gene_expression_data.sort_values(by='logFC', ascending=False)
        gsea_results = gp.prerank(rnk=gene_expression_data, gene_sets=genes_by_pathway, outdir=general_args.gsea_out,
                                  verbose=False, permutation_num=1000, no_plot=True)
        gsea_results.res2d.to_excel(output_path)
    else:
        # Perform Kolmogorov-Smirnov test
        perform_ks_test(enrich_task, genes_by_pathway, scores, general_args.output_dir, test_name, prior_set)
        if enrich_task.ks_significant_pathways_with_genes:
            # Perform Mann-Whitney U test
            perform_mann_whitney_test(enrich_task, general_args, scores, related_pathways)
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



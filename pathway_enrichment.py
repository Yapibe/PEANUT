import pandas as pd
from os import path
from args import EnrichTask
from scipy.stats import rankdata
from utils import load_network_and_pathways
from statsmodels.stats.multitest import multipletests
from visualization_tools import print_enriched_pathways_to_file
from statistic_methods import hypergeometric_sf, jaccard_index , kolmogorov_smirnov_test, compute_mw_python


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
    # Identify genes with P-values below the significance threshold
    significant_p_vals = {gene_id: p_value for gene_id, (score, p_value) in scores.items()
                          if p_value < general_args.FDR_threshold}

    # Retrieve keys (gene IDs) with scores
    scores_keys = set(scores.keys())

    # Filter pathways by those having gene counts within the specified range and that intersect with scored genes
    pathways_with_many_genes = { pathway: set(genes).intersection(scores_keys)
                                 for pathway, genes in genes_by_pathway.items()
                                 if general_args.minimum_gene_per_pathway <=
                                 len(set(genes).intersection(scores_keys)) <= general_args.maximum_gene_per_pathway}

    # Populate a set with all genes from the filtered pathways
    for genes in pathways_with_many_genes.values():
        task.filtered_genes.update(genes)

    # Total number of scored genes
    M = len(scores_keys)
    # Number of genes with significant P-values
    n = len(significant_p_vals)

    # Prepare lists to hold the hypergeometric P-values and corresponding pathway names
    hypergeom_p_values = []
    pathway_names = []

    # Calculate hypergeometric P-values for each pathway with enough genes
    for pathway_name, pathway_genes in pathways_with_many_genes.items():
        N = len(pathway_genes)  # Number of genes in the current pathway
        x = len(set(pathway_genes).intersection(significant_p_vals.keys()))  # Enriched genes in the pathway
        # Apply hypergeometric test; if fewer than 5 enriched genes, assign a P-value of 1 (non-significant)
        pval = hypergeometric_sf(x, M, N, n) if x >= 5 else 1
        hypergeom_p_values.append(pval)
        pathway_names.append(pathway_name)

    # Identify pathways with significant hypergeometric P-values
    significant_pathways = [
        pathway for i, pathway in enumerate(pathway_names) if hypergeom_p_values[i] < 0.05
    ]

    # Perform the Kolmogorov-Smirnov test to compare distributions of scores between pathway genes and background
    ks_p_values = []
    for pathway in significant_pathways:
        genes = pathways_with_many_genes[pathway]
        pathway_scores = [scores[gene_id][0] for gene_id in genes if gene_id in scores]
        background_genes = scores_keys - genes
        background_scores = [scores[gene_id][0] for gene_id in background_genes]
        ks_p_values.append(kolmogorov_smirnov_test(pathway_scores, background_scores))

    # Apply Benjamini-Hochberg correction to the KS P-values
    adjusted_p_values = multipletests(ks_p_values, method='fdr_bh')[1]

    # Filter significant pathways based on adjusted KS P-values
    task.ks_significant_pathways_with_genes = {
        pathway: (pathways_with_many_genes[pathway], adjusted_p_values[i])
        for i, pathway in enumerate(significant_pathways)
        if adjusted_p_values[i] < 0.05
    }


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
        pathway_genes = set(genes_info[0])
        pathway_scores = [scores[gene_id][0] for gene_id in pathway_genes]
        background_genes = task.filtered_genes - pathway_genes
        background_scores = [scores[gene_id][0] for gene_id in background_genes]

        pathway_ranks = [scores_rank[gene_id] for gene_id in pathway_genes]
        background_ranks = [scores_rank[gene_id] for gene_id in background_genes]

        # Compute the Mann-Whitney U test p-value using scores
        # mw_pval = wilcoxon_rank_sums_test(pathway_scores, background_scores)
        # mw_p_values.append(mw_pval)
        _, rmw_pval = compute_mw_python(pathway_ranks, background_ranks)
        mw_p_values.append(rmw_pval)

    # Apply Benjamini-Hochberg correction to adjust the p-values
    adjusted_mw_p_values = multipletests(mw_p_values, method='fdr_bh')[1]


    # Collect significant pathways after adjustment
    filtered_pathways = []
    for i, (pathway, genes) in enumerate(task.ks_significant_pathways_with_genes.items()):
        if adjusted_mw_p_values[i] < args.FDR_threshold:
            filtered_pathways.append({
                'Pathway': pathway,
                'Adjusted_p_value': adjusted_mw_p_values[i],
                'Genes': genes[0]
            })

    # Convert the list of filtered pathways to a DataFrame and sort by p-value
    pathways_df = pd.DataFrame(filtered_pathways)
    pathways_df.sort_values(by='Adjusted_p_value', inplace=True)

    # Filter out pathways with high overlap using the Jaccard index
    for i, row in pathways_df.iterrows():
        current_genes = set(row['Genes'])
        if not any(jaccard_index(current_genes, set(filtered_row['Genes'])) > args.JAC_THRESHOLD
                   for filtered_row in task.filtered_pathways.values()):
            task.filtered_pathways[row['Pathway']] = row


def perform_enrichment(test_name: str, general_args):
    """
    Perform pathway enrichment analysis for a given test.

    Parameters:
    - test_name (str): Name of the test for which enrichment analysis is performed.
    - general_args (GeneralArgs): General arguments and settings.

    Returns:
    - None
    """
    # run enrichment
    propagation_folder = path.join(general_args.propagation_folder, test_name)
    if general_args.run_propagation:
        propagation_file = path.join(f'{propagation_folder}', '{}_{}_{}'.format(test_name, general_args.alpha, general_args.date))
        enrich_task = EnrichTask(name=test_name, create_scores=True, target_field='gene_prop_scores',
                                 statistic_test=kolmogorov_smirnov_test, propagation_file=propagation_file)
    else:
        propagation_file = path.join(f'{propagation_folder}', f'{test_name}_0.1_29_05_2024__15_03_46')
        enrich_task = EnrichTask(name=test_name, create_scores=True, target_field='gene_prop_scores',
                                 statistic_test=kolmogorov_smirnov_test, propagation_file=propagation_file)

    genes_by_pathway, scores = load_network_and_pathways(general_args, enrich_task.propagation_file)

    # Stage 1 - calculate nominal p-values and directions
    perform_statist(enrich_task, general_args, genes_by_pathway, scores)
    if enrich_task.ks_significant_pathways_with_genes:
        # Further statistical test using Mann-Whitney U test
        perform_statist_mann_whitney(enrich_task, general_args, scores)
    # Output the enriched pathways to files
    print_enriched_pathways_to_file(enrich_task, general_args.FDR_threshold)


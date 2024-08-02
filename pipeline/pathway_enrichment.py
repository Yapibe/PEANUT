import pandas as pd
from os import path
from args import EnrichTask, GeneralArgs
from scipy.stats import rankdata
import gseapy as gp
from utils import load_pathways_and_propagation_scores
from statsmodels.stats.multitest import multipletests
from visualization_tools import print_enriched_pathways_to_file
from statistical_methods import jaccard_index , kolmogorov_smirnov_test, compute_mw_python, run_hyper, global_gene_ranking, kolmogorov_smirnov_test_with_ranking


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

    # Populate a set with all genes from the filtered pathways
    for genes in genes_by_pathway.values():
        task.filtered_genes.update(genes)

    # Create a global ranking of genes
    global_ranking = global_gene_ranking(scores)

    if general_args.run_hyper:
        hyper_results_path = path.join(general_args.output_dir, f'hyper_results.txt')
        significant_pathways_hyper = run_hyper(genes_by_pathway, set(scores.keys()), significant_p_vals)
        # Save the significant pathways to a file
        with open(hyper_results_path, 'w') as f:
            for pathway in significant_pathways_hyper:
                f.write(f"{pathway}\n")
    else:
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
        print("No significant pathways found after hypergeometric test. Skipping KS test.")
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
        print("No significant pathways found after KS test.")


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

    if not filtered_pathways:
        print("No significant pathways found after Mann-Whitney U test.")
        return

    # Convert the list of filtered pathways to a DataFrame and sort by p-value
    pathways_df = pd.DataFrame(filtered_pathways)
    pathways_df.sort_values(by='Adjusted_p_value', inplace=True)

    # Filter out pathways with high overlap using the Jaccard index
    for i, row in pathways_df.iterrows():
        current_genes = set(row['Genes'])
        if not any(jaccard_index(current_genes, set(filtered_row['Genes'])) > args.JAC_THRESHOLD
                   for filtered_row in task.filtered_pathways.values()):
            task.filtered_pathways[row['Pathway']] = row


def perform_enrichment(test_name: str, general_args: GeneralArgs, gsea_output_path: str = None):
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

    propagation_file = path.join(f'{propagation_folder}', '{}_{}_{}'.format(test_name, general_args.alpha, general_args.date))
    enrich_task = EnrichTask(name=test_name, create_scores=True, target_field='gene_prop_scores',
                                 statistic_test=kolmogorov_smirnov_test, propagation_file=propagation_file)

    genes_by_pathway, scores = load_pathways_and_propagation_scores(general_args, enrich_task.propagation_file)


    # # Prepare data for GSEA
    # # Unpack the scores dictionary into separate lists for GeneID and Score
    # gene_ids = list(scores.keys())
    # logfc_scores = [score[0] for score in scores.values()]
    # # Create DataFrame for GSEA with string gene identifiers
    # gene_expression_data = pd.DataFrame({'gene': gene_ids, 'logFC': logfc_scores})
    # gene_expression_data['gene'] = gene_expression_data['gene'].astype(str)
    # # Rank the data by logFC in descending order
    # gene_expression_data = gene_expression_data.sort_values(by='logFC', ascending=False)
    # # Run GSEA
    # gsea_results = gp.prerank(rnk=gene_expression_data, gene_sets=genes_by_pathway, outdir=general_args.gsea_out, verbose=True, permutation_num=1000, no_plot=True)


    # # save csv
    # gsea_results.res2d.to_csv(gsea_output_path)

    return scores

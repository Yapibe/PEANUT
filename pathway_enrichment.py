from propagation_routines import generate_propagation_scores
from args import PathwayResults
from statistic_methods import wilcoxon_test, calculate_empirical_p_values, paired_sample_t_test, sign_test
from utils import shuffle_scores, load_network_and_pathways, save_filtered_pathways_to_tsv
from typing import List, Dict
import pandas as pd
import csv
from os import path


def process_tasks(task, general_args, genes_by_pathway, all_experiment_genes_scores):
    """
    Processes a list of tasks for pathway enrichment analysis.
    Args:
        task (Task): tasks to be processed.
        general_args (GeneralArgs): General configuration settings.
        genes_by_pathway (Dict[str, List[int]]): Dictionary mapping pathways to their gene IDs.
        all_experiment_genes_scores (Dict[int, float]): Dictionary mapping all gene IDs in the experiment to their score
    Returns:
        Set[str]: Set of genes by pathway to be included in the analysis.
    """
    significant_pathways_with_genes = {}

    # Filter genes for each pathway
    genes_by_pathway_filtered = {
        pathway: [gene_id for gene_id in genes if gene_id in all_experiment_genes_scores]
        for pathway, genes in genes_by_pathway.items()
    }

    # Print the gene count for 'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS' before filtering
    if 'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS' in genes_by_pathway:
        print(
            f"Before filtering, '{'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'}' has {len(genes_by_pathway['KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'])} genes.")

    # Print the gene count after filtering but before applying gene count criteria
    if 'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS' in genes_by_pathway_filtered:
        print(
            f"After filtering (but before applying gene count criteria), '{'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'}' has {len(genes_by_pathway_filtered['KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'])} genes.")

    # Filter pathways based on gene count criteria
    pathways_with_many_genes = [
        pathway for pathway, genes in genes_by_pathway_filtered.items()
        if general_args.minimum_gene_per_pathway <= len(genes) <= general_args.maximum_gene_per_pathway
    ]

    # Print the gene count after applying all criteria
    if 'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS' in pathways_with_many_genes:
        print(
            f"After all criteria, '{'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'}' is included with {len(genes_by_pathway_filtered['KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS'])} genes.")

    # Manually add a specific pathway
    pathways_with_many_genes.append('REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION')
    # check if pathway KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS is in pathways_with_many_genes
    if 'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS' in pathways_with_many_genes:
        print("KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS included")

    output_file_path = path.join(general_args.output_path, f'{task.name}_filtered_pathways.tsv')
    save_filtered_pathways_to_tsv(pathways_with_many_genes, genes_by_pathway_filtered, output_file_path)
    # Perform statistical tests
    print('After filtering:', len(pathways_with_many_genes))

    # Perform statistical tests
    for pathway in pathways_with_many_genes:
        pathway_scores = [all_experiment_genes_scores[gene_id] for gene_id in genes_by_pathway_filtered[pathway]]

        # Background scores from all genes in the experiment, excluding those in the current pathway
        background_genes = set(all_experiment_genes_scores.keys()) - set(genes_by_pathway_filtered[pathway])
        background_scores = [all_experiment_genes_scores[gene_id] for gene_id in background_genes]
        result = task.statistic_test(pathway_scores, background_scores)
        if result.p_value < general_args.FDR_threshold:
            task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)
            significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[pathway]

    return significant_pathways_with_genes


def run(task1, general_args, num_shuffles=1000):
    """
    Main function to run the pathway enrichment analysis.

    Args:
        task1 (Task): task to be processed.
        general_args (GeneralArgs): General configuration settings.
    Side Effects:
        Executes the entire pathway enrichment analysis workflow, including data loading, processing,
        analysis, and plotting results.

    Returns:
        None
    """
    print("uploding data")
    network_graph, genes_by_pathway, scores, prior_data, input_dict = load_network_and_pathways(general_args, task1)

    # Stage 1 - calculate nominal p-values and directions
    significant_pathways_with_genes = (process_tasks(task1, general_args, genes_by_pathway, scores))

    if task1.alpha == 1:
        # Stage 2 - shuffle scores
        print("shuffling scores")
        updated_prior_data = shuffle_scores(prior_data, 'Score', num_shuffles=num_shuffles)
        # Skip propagation step and directly calculate empirical p-values
        print("Calculating empirical p-values")
        updated_prior_data_float = updated_prior_data.select_dtypes(include=['number']).astype(float)
        empirical_p_values = calculate_empirical_p_values(updated_prior_data_float, significant_pathways_with_genes,
                                                          sign_test, num_simulations=num_shuffles)
    else:
        if task1.create_scores:
            # Stage 2 - shuffle scores
            print("shuffling scores")
            # Step 1: Get all unique gene IDs from the network graph
            all_network_genes = set(network_graph.nodes)

            # Step 2: Find intersected genes and calculate mean score
            intersected_genes = set.intersection(set(prior_data.GeneID), all_network_genes)
            mean_score = prior_data[prior_data.GeneID.isin(intersected_genes)].Score.mean()
            # Step 3: Create new DataFrame for propagation
            propagation_data = pd.DataFrame({'GeneID': list(all_network_genes)})
            propagation_data = propagation_data.merge(prior_data[['GeneID', 'Score']], on='GeneID', how='left')
            propagation_data['Score'].fillna(mean_score, inplace=True)
            # Step 4: Shuffle scores for propagation analysis
            updated_prior_data = shuffle_scores(propagation_data, 'Score', num_shuffles=num_shuffles)
            # Perform propagation for other alpha values

            column_scores_df = generate_propagation_scores(task1, updated_prior_data, network_graph, num_shuffles)
        else:
            # load dataframe from csv
            column_scores_df = pd.read_csv(f'Outputs\\propagation_matrix\\alpha:{task1.alpha}_{num_shuffles}.csv')

        # Stage 3 - calculate empirical p-values
        print("calculating empirical p values")
        empirical_p_values = calculate_empirical_p_values(column_scores_df, significant_pathways_with_genes, wilcoxon_test)

    # save empirical p values to csv
    # get path
    main_dir = path.dirname(path.realpath(__file__))
    empirical_p_values_path = path.join(main_dir, 'Outputs', 'empirical_p_values',
                                        f'empirical_p_values_{task1.alpha}_{num_shuffles}.csv')
    try:
        with open(empirical_p_values_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Pathway', 'Empirical P-Value'])
            for pathway, pval in empirical_p_values.items():
                writer.writerow([pathway, pval])
    except Exception as e:
        print(f"Error occurred while saving CSV: {e}")

    specific_pathways = [
        "WP_PARKINSONS_DISEASE_PATHWAY",
        "KEGG_PARKINSONS_DISEASE",
        "REACTOME_ELASTIC_FIBRE_FORMATION",
        "REACTOME_NUCLEAR_EVENTS_KINASE_AND_TRANSCRIPTION_FACTOR_ACTIVATION",
        "REACTOME_BASIGIN_INTERACTIONS",
        "REACTOME_NGF_STIMULATED_TRANSCRIPTION",
        "WP_PHOTODYNAMIC_THERAPYINDUCED_HIF1_SURVIVAL_SIGNALING",
        "WP_HAIR_FOLLICLE_DEVELOPMENT_CYTODIFFERENTIATION_PART_3_OF_3",
        "WP_OLIGODENDROCYTE_SPECIFICATION_AND_DIFFERENTIATION_LEADING_TO_MYELIN_COMPONENTS_FOR_CNS",
        "WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV",
        "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
        "WP_PHOTODYNAMIC_THERAPYINDUCED_UNFOLDED_PROTEIN_RESPONSE",
        "PID_AP1_PATHWAY",
        "PID_HIF1_TFPATHWAY"
    ]
    for pathway, pval in empirical_p_values.items():
        if pathway in specific_pathways:
            print(f'{pathway}: {pval}')

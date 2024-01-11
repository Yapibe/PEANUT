from propagation_routines import generate_propagation_scores
from args import PathwayResults
from statistic_methods import calculate_empirical_p_values, sign_test, bh_correction, wilcoxon_rank_sums_test, jaccard_index
from utils import shuffle_scores, load_network_and_pathways
import numpy as np
import csv
from os import path


def perform_statist(task, general_args, genes_by_pathway, all_experiment_genes_scores):
    significant_pathways_with_genes = {}
    ks_p_values = []
    pathway_names = []
    mw_p_values = []
    # Filter genes for each pathway, only includes genes that are in the experiment and in the pathway file
    genes_by_pathway_filtered = {
        pathway: [gene_id for gene_id in genes if gene_id in all_experiment_genes_scores]
        for pathway, genes in genes_by_pathway.items()
    }

    # Filter pathways based on gene count criteria
    pathways_with_many_genes = [
        pathway for pathway, genes in genes_by_pathway_filtered.items()
        if general_args.minimum_gene_per_pathway <= len(genes) <= general_args.maximum_gene_per_pathway
    ]

    # Perform statistical tests
    print('After filtering:', len(pathways_with_many_genes))
    for pathway in pathways_with_many_genes:
        pathway_scores = [all_experiment_genes_scores[gene_id] for gene_id in genes_by_pathway_filtered[pathway]]
        background_genes = set(all_experiment_genes_scores.keys()) - set(genes_by_pathway_filtered[pathway])
        background_scores = [all_experiment_genes_scores[gene_id] for gene_id in background_genes]
        result = task.statistic_test(pathway_scores, background_scores)

        task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)
        ks_p_values.append(result.p_value)
        pathway_names.append(pathway)

    # Apply BH correction
    adjusted_p_values = bh_correction(np.array(ks_p_values))

    # Filter significant pathways based on adjusted p-values
    for i, pathway in enumerate(pathway_names):
        if adjusted_p_values[i] < 0.05:  # Using a significance threshold of 0.05
            significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[pathway]

    specific_pathways = [
        "WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV",
        "WP_HIPPOCAMPAL_SYNAPTOGENESIS_AND_NEUROGENESIS",
        "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
        "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION"
    ]
    # print the p-values of specific pathways
    for pathway in specific_pathways:
        if pathway in task.results:
            print(f'{pathway}: {task.results[pathway].p_value}')

    # Mann-Whitney U test and FDR
    for pathway in pathways_with_many_genes:
        pathway_scores = [all_experiment_genes_scores[gene_id] for gene_id in genes_by_pathway_filtered[pathway]]
        background_genes = set(all_experiment_genes_scores.keys()) - set(genes_by_pathway_filtered[pathway])
        background_scores = [all_experiment_genes_scores[gene_id] for gene_id in background_genes]

        # Perform Mann-Whitney U Test
        u_stat, mw_pval =  wilcoxon_rank_sums_test(pathway_scores, background_scores, alternative='two-sided')
        mw_p_values.append(mw_pval)
        pathway_names.append(pathway)

    # Apply BH correction to Mann-Whitney p-values
    adjusted_mw_p_values = bh_correction(np.array(mw_p_values))

    # Filter significant pathways based on adjusted Mann-Whitney p-values
    for i, pathway in enumerate(pathway_names):
        if adjusted_mw_p_values[i] < 0.05:  # Using a significance threshold of 0.05
            significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[pathway]

    # Filter pathways based on adjusted p-values and Jaccard index
    filtered_pathways = {}
    JAC_THRESHOLD = 0.05  # Set your Jaccard threshold
    for i, pathway_i in enumerate(pathway_names):
        if adjusted_mw_p_values[i] > 0.05:
            continue

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
    genes_by_pathway, scores = load_network_and_pathways(task1)

    # Stage 1 - calculate nominal p-values and directions
    significant_pathways_with_genes = perform_statist(task1, general_args, genes_by_pathway, scores)



    print("shuffling scores")
    updated_prior_data = shuffle_scores(scores, 'Score', num_shuffles=num_shuffles)
    print("Calculating empirical p-values")
    updated_prior_data_float = updated_prior_data.select_dtypes(include=['number']).astype(float)
    empirical_p_values = calculate_empirical_p_values(updated_prior_data_float, significant_pathways_with_genes,
                                                          sign_test, num_simulations=num_shuffles)

    # save empirical p values to csv
    # get path
    main_dir = path.dirname(path.realpath(__file__))
    empirical_p_values_path = path.join(main_dir, 'Outputs', 'empirical_p_values',
                                        f'empirical_p_values_ks.csv')
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

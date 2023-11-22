from args import PropagationTask, EnrichTask, GeneralArgs, PathwayResults
from utils import read_network, load_pathways_genes, load_propagation_scores, read_prior_set
from pathway_enrichment import load_network_and_pathways, get_scores
from propagation_routines import propagate_network, generate_similarity_matrix, read_sparse_matrix_txt
import utils
from statistic_methods import wilcoxon_rank_sums_test, students_t_test
import numpy as np
import time
from scipy.stats import rankdata
import pandas as pd
from scipy.stats import ttest_rel


def bh_correction(empirical_values):
    # Extract p-values and their corresponding pathway names
    p_values = np.array(list(empirical_values.values()))
    pathways = list(empirical_values.keys())

    # Apply Benjamini-Hochberg correction
    p_vals_rank = rankdata(p_values, 'max') - 1
    p_vals_rank_ord = rankdata(p_values, 'ordinal') - 1

    p_values_sorted = np.zeros_like(p_vals_rank)
    p_values_sorted[p_vals_rank_ord] = np.arange(len(p_vals_rank_ord))

    p_vals = p_values * (len(p_values) / (p_vals_rank + 1))
    adj_p_vals_by_rank = p_vals[p_values_sorted]

    p_vals_ordered = np.minimum(adj_p_vals_by_rank, np.minimum.accumulate(adj_p_vals_by_rank[::-1])[::-1])
    adj_p_vals = p_vals_ordered[p_vals_rank]

    # Map the adjusted p-values back to the corresponding pathways
    adjusted_values = {pathway: adj_p_val for pathway, adj_p_val in zip(pathways, adj_p_vals)}
    return adjusted_values


# def get_fold_change_ranks(prior_data, score_column):
#     fold_values = {gene: score for gene, score in zip(prior_data['GeneID'], prior_data[score_column])}
#     return {gene: rank for gene, rank in zip(fold_values.keys(), rankdata(list(fold_values.values())))}
#
#
# def compare_pathway_to_randomized(real_scores, shuffled_scores):
#     # Wilcoxon signed-rank test logic
#     differences = real_scores - shuffled_scores
#     signed_ranks = rankdata(np.abs(differences)) * np.sign(differences)
#     return np.sum(signed_ranks)
#
#
# def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
#     empirical_p_values = {}
#
#     # Convert to float once
#     real_data_float = real_data.select_dtypes(include=['number']).astype(float)
#
#     for pathway_name, pathway_genes in pathways.items():
#         count_negative_T = 0
#         count_positive_T = 0
#
#         # Filter once outside the loop
#         pathway_mask = real_data['GeneID'].isin(pathway_genes)
#         pathway_real_scores = real_data_float[pathway_mask]['Score']
#         start = time.time()
#         for i in range(1, num_simulations + 1):
#             shuffled_score_column = f'Shuffled_Score_{i}'
#             # Filter the shuffled scores for pathway genes
#             pathway_shuffled_scores = real_data_float[pathway_mask][shuffled_score_column]
#
#             T = compare_pathway_to_randomized(pathway_real_scores, pathway_shuffled_scores)
#             if T < 0:
#                 count_negative_T += 1
#             elif T > 0:
#                 count_positive_T += 1
#         end = time.time()
#         print(f'pathway {pathway_name} took {end - start} seconds')
#         # Choose the smaller count to calculate the empirical p-value and double it for a two-tailed test
#         min_count = min(count_negative_T, count_positive_T)
#         empirical_p_values[pathway_name] = 2 * (min_count + 1) / (num_simulations + 1)
#
#     return empirical_p_values


# def calculate_pathway_score(ranks, pathway_genes):
#     return sum(ranks.get(gene, 0) for gene in pathway_genes)
#
#
# def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
#     empirical_p_values = {}
#     real_gene_rank = get_fold_change_ranks(real_data, 'Score')
#     shuffled_gene_ranks = {}
#
#     # Precompute shuffled gene ranks for each simulation
#     for i in range(1, num_simulations + 1):
#         shuffled_score_column = f'Shuffled_Score_{i}'
#         shuffled_gene_ranks[i] = get_fold_change_ranks(real_data, shuffled_score_column)
#
#     for pathway_name, pathway_genes in pathways.items():
#         score = 0
#         observed_score = calculate_pathway_score(real_gene_rank, pathway_genes)
#         for i in range(1, num_simulations + 1):
#             shuffled_scores = calculate_pathway_score(shuffled_gene_ranks[i], pathway_genes)
#             if observed_score < shuffled_scores:
#                 score += 1
#         empirical_p_values[pathway_name] = (score + 1) / (num_simulations + 1)
#     return empirical_p_values


def paired_sample_t_test(real_scores, shuffled_scores):
    # Ensure that real_scores and shuffled_scores are NumPy arrays of type float
    real_scores = np.asarray(real_scores, dtype=float)
    shuffled_scores = np.asarray(shuffled_scores, dtype=float)

    # Perform the paired sample t-test
    t_statistic, p_value = ttest_rel(real_scores, shuffled_scores)
    return t_statistic


def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
    empirical_p_values = {}

    for pathway_name, pathway_genes in pathways.items():
        count_negative_t = 0
        count_positive_t = 0
        # Filter the real scores for pathway genes
        pathway_real_scores = real_data[real_data['GeneID'].isin(pathway_genes)]['Score'].astype(float)

        for i in range(1, num_simulations + 1):
            shuffled_score_column = f'Shuffled_Score_{i}'
            # Filter the shuffled scores for pathway genes
            pathway_shuffled_scores = real_data[real_data['GeneID'].isin(pathway_genes)][shuffled_score_column].astype(
                float)

            T = paired_sample_t_test(pathway_real_scores, pathway_shuffled_scores)
            if T < 0:
                count_negative_t += 1
            elif T > 0:
                count_positive_t += 1

        # Choose the smaller count to calculate the empirical p-value and double it for a two-tailed test
        min_count = min(count_negative_t, count_positive_t)
        empirical_p_values[pathway_name] = 2 * (min_count + 1) / (num_simulations + 1)

    return empirical_p_values


def load_and_prepare_data(task):
    data = read_prior_set(task.experiment_file_path)
    data = data[data['GeneID'].apply(lambda x: str(x).isdigit())]
    data['GeneID'] = data['GeneID'].astype(int)
    # remove p_values column
    data = data.drop(columns=['P-value'])
    data = data.reset_index(drop=True)

    return data


def shuffle_scores(dataframe, shuffle_column, num_shuffles=1000):
    if shuffle_column not in dataframe.columns:
        raise ValueError(f"Column '{shuffle_column}' not found in the DataFrame")

    # Dictionary to store all shuffled columns
    shuffled_columns = {}

    # Generate shuffled columns
    for i in range(1, num_shuffles + 1):
        shuffled_column = dataframe[shuffle_column].sample(frac=1).reset_index(drop=True)
        shuffled_columns[f'Shuffled_Score_{i}'] = shuffled_column

    # Convert the dictionary to a DataFrame
    shuffled_df = pd.DataFrame(shuffled_columns)

    # Concatenate the original dataframe with the new shuffled scores DataFrame
    result_df = pd.concat([dataframe, shuffled_df], axis=1)

    return result_df


def process_tasks(task_list, network_graph, general_args, genes_by_pathway):
    pathways_to_display = set()
    significant_pathways_with_genes = {}

    # Create a set to hold all genes that are in some pathway and are in the network.
    all_genes_in_filtered_pathways_and_network = set()

    for task in task_list:
        scores = get_scores(task)
        # Filter genes for each pathway, contains only genes that are in the experiment and in the pathway file
        genes_by_pathway_filtered = {pathway: [id for id in genes if id in scores]
                                     for pathway, genes in genes_by_pathway.items()}

        # keep only pathway with certain amount of genes
        pathways_with_many_genes = [pathway_name for pathway_name in genes_by_pathway_filtered.keys() if
                                    (len(genes_by_pathway_filtered[
                                             pathway_name]) >= general_args.minimum_gene_per_pathway and len(
                                        genes_by_pathway_filtered[
                                            pathway_name]) <= general_args.maximum_gene_per_pathway)]
        # manually add REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION'
        pathways_with_many_genes.append('REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION')
        # pathways_with_many_genes.append('WP_PARKINSONS_DISEASE_PATHWAY')
        # pathways_with_many_genes.append('KEGG_PARKINSONS_DISEASE')

        # Update all_genes_in_filtered_pathways_and_network
        for pathway in pathways_with_many_genes:
            all_genes_in_filtered_pathways_and_network.update(genes_by_pathway_filtered[pathway])
            pathways_to_display.add(pathway)
        # Intersect with network nodes to refine the set
        all_genes_in_filtered_pathways_and_network &= set(network_graph.nodes)

        # Perform statistical tests
        print('after filtering', len(pathways_to_display))
        for pathway in pathways_to_display:
            pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
            background_scores = [scores[id] for id in all_genes_in_filtered_pathways_and_network if
                                 id not in genes_by_pathway_filtered[pathway]]
            result = task.statistic_test(pathway_scores, background_scores)
            print(f'nominal: {pathway}: {result.p_value}')
            if result.p_value < general_args.FDR_threshold:
                task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)
                significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[pathway]

    return np.sort(list(pathways_to_display)), significant_pathways_with_genes


def run(task_list, general_args):
    """
    takes a list of tasks, a general args object and an optional dataset type and runs the pathway enrichment analysis
    :param task_list: list of tasks
    :param general_args: general args object
    :return: None
    """
    network_graph, genes_by_pathway = load_network_and_pathways(general_args)

    pathways_to_display, genes_by_pathway_filtered = (process_tasks(task_list, network_graph, general_args,
                                                                    genes_by_pathway))

    # task2 = EnrichTask(name='TvN', propagation_file='TvN_Score_0.1_20_11_2023__16_43_14',
    #                    propagation_folder=f'Outputs\\propagation_scores\\TvN',
    #                    statistic_test=wilcoxon_rank_sums_test, target_field='gene_prop_scores')
    #
    # prop_scores = get_scores(task2)
    prior_data = load_and_prepare_data(prop_task)
    # # Create a copy of prior_data to avoid SettingWithCopyWarning
    # prior_data = prior_data.copy()
    # all_genes_ids = set.intersection(set(prior_data.GeneID), set(network_graph.nodes))
    # prior_data = prior_data[prior_data['GeneID'].isin(all_genes_ids)]
    #
    # # Update 'Score' column with propagation scores
    # prior_data.loc[:, 'Score'] = prior_data['GeneID'].apply(lambda x: prop_scores[x])

    # Shuffle the scores and add to the DataFrame
    updated_prior_data = shuffle_scores(prior_data, 'Score')
    empirical_p_values = calculate_empirical_p_values(updated_prior_data, genes_by_pathway_filtered)
    for pathway, pval in empirical_p_values.items():
        if pval <= 0.05:
            print(f'{pathway}: {pval}')
    # send to bh correction
    # corrected_p_values = bh_correction(empirical_p_values)
    # filter by FDR 0.05
    # corrected_p_values = {pathway: p_val for pathway, p_val in corrected_p_values.items() if
    #                       p_val < general_args.FDR_threshold}
    # print(f'number of pathways: {len(corrected_p_values)}')
    # print('significant pathways:')
    # for pathway in corrected_p_values:
    #     print(f'{pathway}: {empirical_p_values[pathway]}')


def perform_enrichment():
    # run enrichment
    print("running enrichment")
    tasks = []
    # propagation_scores_file = '{}_{}_{}_{}'.format(task.experiment_name, task.propagation_input_type,
    #                                                task.alpha, task.date)
    task1 = EnrichTask(name='TvN', propagation_file='TvN_Score_1_20_11_2023__16_20_09',
                       propagation_folder=f'Outputs\\propagation_scores\\TvN',
                       statistic_test=wilcoxon_rank_sums_test, target_field='gene_prop_scores')

    FDR_threshold = 0.05

    figure_name = prop_task.experiment_name + '-alpha' + str(
        prop_task.alpha) + '-Threshold' + str(FDR_threshold) + '.pdf'

    general_args = GeneralArgs(prop_task.network_file_path, genes_names_path=prop_task.genes_names_file_path,
                               pathway_members_path=prop_task.pathway_file_dir, FDR_threshold=FDR_threshold,
                               figure_name=figure_name)

    tasks += [task1]
    print('running')
    run(tasks, general_args)


if __name__ == '__main__':
    start = time.time()
    prop_task = PropagationTask(experiment_name='pipeline_ttest', create_similarity_matrix=False)
    perform_enrichment()
    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))

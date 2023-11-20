from args import PropagationTask, EnrichTask, GeneralArgs, PathwayResults
from utils import read_network, load_pathways_genes, load_propagation_scores, read_prior_set
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


def get_fold_change_ranks(prior_data, score_column):
    fold_values = {gene: score for gene, score in zip(prior_data['GeneID'], prior_data[score_column])}
    return {gene: rank for gene, rank in zip(fold_values.keys(), rankdata(list(fold_values.values())))}


# def compare_pathway_to_randomized(real_scores, shuffled_scores):
#     # Ensure that real_scores and shuffled_scores are NumPy arrays of type float
#     real_scores = np.asarray(real_scores, dtype=float)
#     shuffled_scores = np.asarray(shuffled_scores, dtype=float)
#
#     # Wilcoxon signed-rank test logic
#     differences = real_scores - shuffled_scores
#     abs_diff = np.abs(differences)
#     ranks = rankdata(abs_diff)
#     signed_ranks = ranks * np.sign(differences)
#     T = np.sum(signed_ranks)
#     return T
#
#
# def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
#     empirical_p_values = {}
#
#     for pathway_name, pathway_genes in pathways.items():
#         count_negative_T = 0
#         # Filter the real scores for pathway genes
#         pathway_real_scores = real_data[real_data['GeneID'].isin(pathway_genes)]['Score'].astype(float)
#
#         for i in range(1, num_simulations + 1):
#             shuffled_score_column = f'Shuffled_Score_{i}'
#             # Filter the shuffled scores for pathway genes
#             pathway_shuffled_scores = real_data[real_data['GeneID'].isin(pathway_genes)][shuffled_score_column].astype(
#                 float)
#
#             T = compare_pathway_to_randomized(pathway_real_scores, pathway_shuffled_scores)
#             if T < 0:
#                 count_negative_T += 1
#
#         empirical_p_values[pathway_name] = (count_negative_T + 1) / (num_simulations + 1)
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
        count_negative_T = 0
        # Filter the real scores for pathway genes
        pathway_real_scores = real_data[real_data['GeneID'].isin(pathway_genes)]['Score'].astype(float)

        for i in range(1, num_simulations + 1):
            shuffled_score_column = f'Shuffled_Score_{i}'
            # Filter the shuffled scores for pathway genes
            pathway_shuffled_scores = real_data[real_data['GeneID'].isin(pathway_genes)][shuffled_score_column].astype(
                float)

            T = paired_sample_t_test(pathway_real_scores, pathway_shuffled_scores)
            if T < 0:
                count_negative_T += 1

        empirical_p_values[pathway_name] = (count_negative_T + 1) / (num_simulations + 1)

    return empirical_p_values


def load_and_prepare_data(task):
    data = read_prior_set(task.experiment_file_path)
    data['Score'] = data['Score'].apply(lambda x: abs(x))
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


def get_scores(task):
    """
    takes a task, a network graph and a general args object and returns a dictionary of scores
    :param task: EnrichTask or RawScoreTask
    :return: dictionary of scores
    """

    if isinstance(task, EnrichTask):
        result_dict = load_propagation_scores(task, propagation_file_name=task.propagation_file)
        gene_id_to_idx = {xx: x for x, xx in result_dict['gene_idx_to_id'].items()}
        scores = result_dict[task.target_field]
        if task.constrain_to_experiment_genes:
            gene_id_to_idx = {id: idx for id, idx in gene_id_to_idx.items() if id in result_dict['propagation_input']}
        scores = {id: scores[idx][0] for id, idx in gene_id_to_idx.items()}
    else:
        raise ValueError('Invalid task')
    return scores


def process_tasks(task_list, network_graph, general_args, interesting_pathways, genes_by_pathway):
    pathways_to_display = set()
    significant_pathways_with_genes = {}

    # Create a set to hold all genes that are in some pathway and are in the network.
    all_genes_in_filtered_pathways_and_network = set()

    for task in task_list:
        scores = get_scores(task)
        # Filter genes for each pathway, contains only genes that are in the experiment and in the pathway file
        genes_by_pathway_filtered = {pathway: [id for id in genes_by_pathway[pathway] if id in scores]
                                     for pathway in interesting_pathways}

        # keep only pathway with certain amount of genes
        pathways_with_many_genes = [pathway_name for pathway_name in genes_by_pathway_filtered.keys() if
                                    (len(genes_by_pathway_filtered[
                                             pathway_name]) >= general_args.minimum_gene_per_pathway and len(
                                        genes_by_pathway_filtered[
                                            pathway_name]) <= general_args.maximum_gene_per_pathway)]
        # manually add REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION
        pathways_with_many_genes.append('REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION')
        pathways_with_many_genes.append('WP_PARKINSONS_DISEASE_PATHWAY')
        pathways_with_many_genes.append('KEGG_PARKINSONS_DISEASE')

        # Update all_genes_in_filtered_pathways_and_network
        for pathway in pathways_with_many_genes:
            all_genes_in_filtered_pathways_and_network.update(genes_by_pathway_filtered[pathway])

        # Intersect with network nodes to refine the set
        all_genes_in_filtered_pathways_and_network &= set(network_graph.nodes)

        for pathway in pathways_with_many_genes:
            pathways_to_display.add(pathway)

        # Perform statistical tests
        print('after filtering', len(pathways_to_display))
        for pathway in pathways_to_display:
            pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
            background_scores = [scores[id] for id in all_genes_in_filtered_pathways_and_network if
                                 id not in genes_by_pathway_filtered[pathway]]
            result = task.statistic_test(pathway_scores, background_scores)
            print(f'nominal:{pathway}: {result.p_value}')
            if result.p_value < general_args.FDR_threshold:
                task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)
                significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[
                    pathway]  # Store the genes for significant pathways

    pathways_to_display = np.sort(list(pathways_to_display))

    return pathways_to_display, significant_pathways_with_genes


def load_network_and_pathways(general_args):
    network_graph = read_network(general_args.network_file_path)
    genes_by_pathway = load_pathways_genes(general_args.pathway_members_path)
    interesting_pathways = list(genes_by_pathway.keys())

    return network_graph, interesting_pathways, genes_by_pathway


def run(task_list, general_args, prior_data):
    """
    takes a list of tasks, a general args object and an optional dataset type and runs the pathway enrichment analysis
    :param task_list: list of tasks
    :param general_args: general args object
    :return: None
    """
    network_graph, interesting_pathways, genes_by_pathway = load_network_and_pathways(general_args)

    pathways_to_display, genes_by_pathway_filtered = (process_tasks(task_list, network_graph, general_args,
                                                                    interesting_pathways,
                                                                    genes_by_pathway))

    # Shuffle the scores and add to the DataFrame
    updated_prior_data = shuffle_scores(prior_data, 'Score')
    empirical_p_values = calculate_empirical_p_values(updated_prior_data, genes_by_pathway_filtered)
    # send to bh correction
    corrected_p_values = bh_correction(empirical_p_values)
    # filter by FDR 0.05
    # corrected_p_values = {pathway: p_val for pathway, p_val in corrected_p_values.items() if
    #                       p_val < general_args.FDR_threshold}
    # print
    print(f'number of pathways: {len(corrected_p_values)}')
    print('significant pathways:')
    for pathway in corrected_p_values:
        print(f'{pathway}: {empirical_p_values[pathway]}')


def perform_enrichment(task):
    # run enrichment
    print("running enrichment")
    tasks = []
    task1 = EnrichTask(name='TvN', propagation_file='TvN_abs_Score_1_18_11_2023__20_05_16',
                       propagation_folder=f'Outputs\\propagation_scores\\TvN',
                       statistic_test=wilcoxon_rank_sums_test,
                       target_field='gene_prop_scores', constrain_to_experiment_genes=True)

    FDR_threshold = 0.05

    figure_name = task.experiment_name + '-alpha' + str(
        task.alpha) + '-Threshold' + str(FDR_threshold) + '.pdf'

    general_args = GeneralArgs(task.network_file_path, genes_names_path=task.genes_names_file_path,
                               pathway_members_path=task.pathway_file_dir, FDR_threshold=FDR_threshold,
                               figure_name=figure_name)

    tasks += [task1]
    print('running')
    prior_data = load_and_prepare_data(task)
    run(tasks, general_args, prior_data)


if __name__ == '__main__':
    start = time.time()
    task = PropagationTask(experiment_name='TvN', create_similarity_matrix=False)
    perform_enrichment(task)
    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))

from args import PropagationTask, EnrichTask, GeneralArgs, PathwayResults
from utils import read_prior_set, get_propagation_input
from statistic_methods import wilcoxon_rank_sums_test, students_t_test
import numpy as np
import time
from scipy.stats import rankdata
import pandas as pd
from scipy.stats import ttest_rel
from pathway_enrichment import load_network_and_pathways, get_scores
from propagation_routines import read_sparse_matrix_txt, propagate_network


def column_propagation(column_data, all_genes_ids, matrix, genes):
    # Extract the score column name
    score_column = column_data.columns[1]  # Assuming it's the second column

    # Create a dictionary mapping from GeneID to Score
    gene_id_to_score = column_data.set_index('GeneID')[score_column].to_dict()

    # Filter the dictionary to only include relevant gene IDs
    propagation_input = {int(gene_id): gene_id_to_score.get(gene_id, 0) for gene_id in all_genes_ids}
    gene_index = dict([(gene, index) for (index, gene) in enumerate(genes)])
    score_gene_scores_inverse, gene_score_dict = propagate_network(propagation_input, matrix, gene_index)
    return score_gene_scores_inverse.flatten()


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
def compare_pathway_to_randomized(real_scores, shuffled_scores):
    # Wilcoxon signed-rank test logic
    differences = real_scores - shuffled_scores
    signed_ranks = rankdata(np.abs(differences)) * np.sign(differences)
    return np.sum(signed_ranks)


def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
    empirical_p_values = {}

    # Convert to float once
    real_data_float = real_data.select_dtypes(include=['number']).astype(float)

    for pathway_name, pathway_genes in pathways.items():
        count_negative_T = 0
        count_positive_T = 0

        # Filter once outside the loop
        pathway_mask = real_data['GeneID'].isin(pathway_genes)
        pathway_real_scores = real_data_float[pathway_mask]['Score']
        for i in range(1, num_simulations + 1):
            shuffled_score_column = f'Shuffled_Score_{i}'
            # Filter the shuffled scores for pathway genes
            pathway_shuffled_scores = real_data_float[pathway_mask][shuffled_score_column]

            T = compare_pathway_to_randomized(pathway_real_scores, pathway_shuffled_scores)
            if T < 0:
                count_negative_T += 1
            elif T > 0:
                count_positive_T += 1
        # Choose the smaller count to calculate the empirical p-value and double it for a two-tailed test
        min_count = min(count_negative_T, count_positive_T)
        empirical_p_values[pathway_name] = 2 * (min_count + 1) / (num_simulations + 1)

    return empirical_p_values


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


# def paired_sample_t_test(real_scores, shuffled_scores):
#     # Ensure that real_scores and shuffled_scores are NumPy arrays of type float
#     real_scores = np.asarray(real_scores, dtype=float)
#     shuffled_scores = np.asarray(shuffled_scores, dtype=float)
#
#     # Perform the paired sample t-test
#     t_statistic, p_value = ttest_rel(real_scores, shuffled_scores)
#     return t_statistic
#
#
# def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
#     empirical_p_values = {}
#     all_pathway_genes = set.union(*[set(genes) for genes in pathways.values()])
#     filtered_real_data = real_data[real_data['GeneID'].isin(all_pathway_genes)]
#
#     for pathway_name, pathway_genes in pathways.items():
#         # Filter the real scores for pathway genes
#         pathway_real_scores = filtered_real_data[filtered_real_data['GeneID'].isin(pathway_genes)]['Score'].to_numpy(
#             dtype=float)
#
#         count_negative_t = 0
#         count_positive_t = 0
#         for i in range(1, num_simulations + 1):
#             shuffled_score_column = f'Shuffled_Score_{i}'
#             # Convert to NumPy array for faster processing
#             pathway_shuffled_scores = filtered_real_data[filtered_real_data['GeneID'].isin(pathway_genes)][
#                 shuffled_score_column].to_numpy(dtype=float)
#
#             T = paired_sample_t_test(pathway_real_scores, pathway_shuffled_scores)
#             if T < 0:
#                 count_negative_t += 1
#             elif T > 0:
#                 count_positive_t += 1
#
#         # Choose the smaller count to calculate the empirical p-value and double it for a two-tailed test
#         min_count = min(count_negative_t, count_positive_t)
#         empirical_p_values[pathway_name] = 2 * (min_count + 1) / (num_simulations + 1)
#
#     return empirical_p_values


def load_and_prepare_data(task):
    data = read_prior_set(task.experiment_file_path)
    data = data[data['GeneID'].apply(lambda x: str(x).isdigit())]
    data['GeneID'] = data['GeneID'].astype(int)
    # remove p_values column
    data = data.drop(columns=['P-value'])
    data = data.reset_index(drop=True)

    return data


def shuffle_scores(dataframe, shuffle_column, num_shuffles=10000):
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

def get_gene_p_values(gene_id):
    """
        Retrieves the p-value for a specific gene from a DataFrame.
        Returns:
        - float: a dict of gene id: p value
        """
    # Load the DataFrame
    df = pd.read_excel('Inputs/experiments_data/scores_T_v_N.xlsx')
    # Check if the gene_id is in the DataFrame and return the corresponding p-value
    p_value_row = df[df['GeneID'] == gene_id]
    if not p_value_row.empty:
        return p_value_row.iloc[0]['P-value']
    else:
        return None




def process_tasks(task_list, network_graph, general_args, genes_by_pathway):
    pathways_to_display = set()
    significant_pathways_with_genes = {}

    # Create a set to hold all genes that are in some pathway and are in the network.
    all_genes_in_filtered_pathways_and_network = set()

    for task in task_list:
        scores = get_scores(task)

        # Filter genes for each pathway
        genes_by_pathway_filtered = {
            pathway: [gene_id for gene_id in genes if gene_id in scores]
            for pathway, genes in genes_by_pathway.items()
        }

        # Define the pathways of interest
        pathways_of_interest = [
            'WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER',
            'WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV'
        ]

        # Iterate through the pathways of interest
        for pathway in pathways_of_interest:
            print(f"Pathway: {pathway}")
            filtered_genes = genes_by_pathway_filtered[pathway]
            print(f"Number of genes in pathway: {len(filtered_genes)}")
            print("Gene IDs, Scores, and p-values:")

            # Print each gene ID, its score, and p-value if p-value is less than 0.05
            for gene_id in filtered_genes:
                gene_score = scores.get(gene_id, "Not available")
                gene_p_value = get_gene_p_values(gene_id)

                if gene_p_value and gene_p_value < 0.05:
                    print(f"Gene ID: {gene_id}, Score: {gene_score}, p-value: {gene_p_value}")

            print()  # Add an empty line for better readability

        # keep only pathway with certain amount of genes
        pathways_with_many_genes = [pathway_name for pathway_name in genes_by_pathway_filtered.keys() if
                                    (len(genes_by_pathway_filtered[
                                             pathway_name]) >= general_args.minimum_gene_per_pathway and len(
                                        genes_by_pathway_filtered[
                                            pathway_name]) <= general_args.maximum_gene_per_pathway)]
        pathways_with_many_genes.append('REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION')

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

            # Check if the pathway is significant
            if result.p_value < general_args.FDR_threshold:
                task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)
                significant_pathways_with_genes[pathway] = genes_by_pathway_filtered[pathway]

                # Print additional information for specified pathways
                if pathway in pathways_of_interest:
                    print(f"\nPathway: {pathway}")
                    print(f"Number of genes in pathway: {len(genes_by_pathway_filtered[pathway])}")
                    print(f"Mean score of pathway: {np.mean(pathway_scores):.2f}")
                    print(f"P-value: {result.p_value:.2e}")
                    print(f"Directionality: {result.directionality}")

    return np.sort(list(pathways_to_display)), significant_pathways_with_genes


def run(task_list, general_args):
    """
    takes a list of tasks, a general args object and an optional dataset type and runs the pathway enrichment analysis
    :param task_list: list of tasks
    :param general_args: general args object
    :return: None
    """
    # Stage 1 - Load the network and pathways and get nominal p values
    network_graph, genes_by_pathway = load_network_and_pathways(general_args)
    # print genes of the pathway WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER
    print(len(genes_by_pathway['WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER']))
    print(genes_by_pathway['WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER'])
    print(len(genes_by_pathway['WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV']))
    print(genes_by_pathway['WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV'])

    pathways_to_display, genes_by_pathway_filtered = (process_tasks(task_list, network_graph, general_args,
                                                                    genes_by_pathway))
    # Stage 2 - Calculate empirical p values
    print("uploding data")
    # prior_data = load_and_prepare_data(task)
    # # Shuffle the scores and add to the DataFrame
    # print("shuffling scores")
    # updated_prior_data = shuffle_scores(prior_data, 'Score')
    # all_genes_ids = set.intersection(set(updated_prior_data.GeneID), set(network_graph.nodes))
    # # Read or create similarity matrix
    # print("reading similarity matrix")
    # matrix, genes = read_sparse_matrix_txt(network_graph, task.similarity_matrix_path)
    # column_scores = {}
    # # Iterate over shuffled columns
    # print("propagating")
    # for column in updated_prior_data.columns:
    #     if column != 'GeneID' and column != 'Human_Name':
    #         # Create sub dataframe of GeneID and current column
    #         column_data = updated_prior_data[['GeneID', column]]
    #         # Perform propagation
    #         column_scores[column] = column_propagation(column_data, all_genes_ids, matrix, genes)
    # del matrix
    # # turn column_scores to dataframe and add gene name
    # column_scores_df = pd.DataFrame.from_dict(column_scores)
    # column_scores_df['Human_Name'] = updated_prior_data['Human_Name']
    # column_scores_df['GeneID'] = updated_prior_data['GeneID']
    # # change order of columns
    # column_scores_df = column_scores_df[
    #     ['GeneID', 'Human_Name'] + [col for col in column_scores_df.columns if col not in ['GeneID', 'Human_Name']]]
    # # save to csv
    # column_scores_df.to_csv(f'Outputs\\propagation_file\\{task.experiment_name}_10,000.csv', index=False)
    #load dataframe from csv
    column_scores_df = pd.read_csv(f'Outputs\\propagation_file\\matrix_with_prop_10,000.csv')
    print("calculating empirical p values")
    empirical_p_values = calculate_empirical_p_values(column_scores_df, genes_by_pathway_filtered)
    # save empirical p values to csv
    pd.DataFrame.from_dict(empirical_p_values, orient='index').to_csv(
        f'Outputs\\propagation_scores\\scores_with_prop_10,000.csv')
    # for pathway, pval in empirical_p_values.items():
    #     if pval <= 0.05:
    #         print(f'{pathway}: {pval}')
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
        "WP_PHOTODYNAMIC_THERAPYINDUCED_UNFOLDED_PROTEIN_RESPONSE",
        "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
        "PID_AP1_PATHWAY",
        "PID_HIF1_TFPATHWAY"
    ]
    for pathway, pval in empirical_p_values.items():
        if pathway in specific_pathways:
            print(f'{pathway}: {pval}')


def perform_enrichment(task):
    # run enrichment
    print("running enrichment")
    tasks = []
    task1 = EnrichTask(name='TvN', propagation_file='TvN_Score_1_23_11_2023__13_37_07',
                       propagation_folder=f'Outputs\\propagation_scores\\TvN',
                       statistic_test=wilcoxon_rank_sums_test, target_field='gene_prop_scores')

    FDR_threshold = 0.05

    figure_name = task.experiment_name + '-alpha' + str(
        task.alpha) + '-Threshold' + str(FDR_threshold) + '.pdf'

    general_args = GeneralArgs(task.network_file_path, genes_names_path=task.genes_names_file_path,
                               pathway_members_path=task.pathway_file_dir, FDR_threshold=FDR_threshold,
                               figure_name=figure_name)

    tasks += [task1]
    print('running')
    run(tasks, general_args)


if __name__ == '__main__':
    start = time.time()
    task = PropagationTask(experiment_name='TvN', create_similarity_matrix=False)
    perform_enrichment(task)
    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))

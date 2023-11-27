from args import PropagationTask, EnrichTask, GeneralArgs
from utils import save_propagation_score, read_prior_set, get_propagation_input, filter_network_by_prior_data
from propagation_routines import propagate_network, generate_similarity_matrix, read_sparse_matrix_txt
from statistic_methods import wilcoxon_rank_sums_test, students_t_test
from pathway_enrichment import run
import time
import numpy as np


def perform_propagation(prop_task):
    """
    Executes the gene score propagation process using a provided network and prior data.

    This function reads the network graph, loads the prior data, and performs gene score propagation.
    It involves normalization of propagation scores and saving these scores for further analysis.

    Parameters:
    - task (PropagationTask): An object containing all necessary parameters and file paths for propagation.

    Returns:
    - None: This function does not return a value but saves the propagation scores to a specified directory.
    """

    # Load and prepare prior set
    prior_data = read_prior_set(prop_task.experiment_file_path)

    if prop_task.alpha == 1:
        print("Skipping propagation, saving sorted scores directly")

        # Sort prior_data by GeneID
        sorted_prior_data = prior_data.sort_values(by='GeneID').reset_index(drop=True)

        # Create experiment_gene_index based on sorted GeneID
        genes = sorted_prior_data['GeneID']
        experiment_gene_index = {gene_id: idx for idx, gene_id in enumerate(genes)}

        # Create propagation_input as a dictionary
        propagation_input = {gene_id: score for gene_id, score in
                             zip(sorted_prior_data['GeneID'], sorted_prior_data['Score'])}

        # Create gene_scores as a ndarray of scores
        gene_scores = sorted_prior_data['Score'].values
        # make gene_scores a 2 dimensional array
        gene_scores = gene_scores.reshape((len(gene_scores), 1))

        save_propagation_score(propagation_scores=gene_scores, prior_set=sorted_prior_data,
                               propagation_input=propagation_input, genes_id_to_idx=experiment_gene_index,
                               task=prop_task, save_dir=prop_task.output_folder)

        return

    # reads the network graph from a file
    filtered_network = filter_network_by_prior_data(prop_task.network_file_path, prior_data)
    # Intersection with network nodes
    all_genes_ids = filtered_network.nodes()

    # create or upload similarity matrix
    if prop_task.create_similarity_matrix:
        print("generating similarity matrix")
        matrix, network_gene_index = generate_similarity_matrix(filtered_network, prop_task.similarity_matrix_path,
                                                                prop_task.alpha)
    else:
        print("reading similarity matrix")
        matrix, network_gene_index = read_sparse_matrix_txt(filtered_network, prop_task.similarity_matrix_path)
        print("uploaded similarity matrix")

    # Propagate network
    print("propagating network")
    propagation_input = get_propagation_input(all_genes_ids, prior_data, prop_task.propagation_input_type)
    propagation_score, gene_score_dict = propagate_network(propagation_input, matrix, network_gene_index)

    # print("getting ones input")
    ones_input = get_propagation_input(all_genes_ids, prior_data, 'ones')
    ones_gene_scores_inverse, ones_gene_score_dict = propagate_network(ones_input, matrix, network_gene_index)

    del matrix

    # Identify genes with zero normalization score but non-zero propagation score
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(propagation_score == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))

    # Set the normalization score of these genes to 1
    ones_gene_scores_inverse[genes_to_delete] = 1

    # Perform the normalization
    non_zero_indices = np.nonzero(propagation_score != 0)[0]
    propagation_score[non_zero_indices] = propagation_score[non_zero_indices] / np.abs(ones_gene_scores_inverse[non_zero_indices])

    # save propagation score
    print("saving propagation score")
    save_propagation_score(propagation_scores=propagation_score, prior_set=prior_data, propagation_input=propagation_input,
                           genes_id_to_idx=network_gene_index, task=prop_task, save_dir=prop_task.output_folder)


def perform_enrichment(prop_task):
    """
    Executes the enrichment analysis on propagated gene scores.
    This function sets up tasks for enrichment analysis, including defining parameters and file paths.
    It then runs the enrichment analysis and processes the results.
    Parameters:
    - task (EnrichTask): An object containing parameters and file paths for running the enrichment analysis.
    Returns:
    - None: This function does not return a value but may generate output files like plots or data summaries.
    """
    # run enrichment
    print("running enrichment")
    FDR_threshold = 0.05
    num_shuffles = 1000
    alpha = 0.2
    if run_propagation_flag:

        propagation_scores_file = '{}_{}_{}_{}'.format(prop_task.experiment_name, prop_task.propagation_input_type,
                                                       prop_task.alpha, prop_task.date)
        task1 = EnrichTask(name=prop_task.experiment_name, propagation_file=propagation_scores_file,
                           propagation_folder=f'Outputs\\propagation_scores\\{prop_task.experiment_name}',
                           statistic_test=wilcoxon_rank_sums_test, target_field='gene_prop_scores', alpha=alpha,
                           create_propagation_matrix=False, create_scores=True)
    else:
        task1 = EnrichTask(name='TvN', propagation_file='TvN_Score_1_26_11_2023__20_07_31',
                           propagation_folder=f'Outputs\\propagation_scores\\TvN',
                           statistic_test=wilcoxon_rank_sums_test, target_field='gene_prop_scores', alpha=alpha,
                           create_propagation_matrix=False, create_scores=True)

    general_args = GeneralArgs(prop_task.network_file_path, genes_names_path=prop_task.genes_names_file_path,
                               pathway_members_path=prop_task.pathway_file_dir, FDR_threshold=FDR_threshold)
    print('running')
    run(task1, general_args, num_shuffles=num_shuffles)


def main(run_propagation=True, run_enrichment=True):
    """
    Main function to execute propagation and enrichment analysis based on specified flags.
    This function initializes tasks for propagation and enrichment and executes them based on the
    provided flags. It serves as the entry point for running the gene score propagation and enrichment analysis pipeline.
    Parameters:
    - run_propagation (bool): Flag to determine whether to run propagation (default: True).
    - run_enrichment (bool): Flag to determine whether to run enrichment analysis (default: True).
    Returns:
    - None: This function orchestrates the execution of other functions but does not return a value.
    """

    # create a propagation task
    prop_task = PropagationTask(experiment_name='TvN')

    if run_propagation:
        perform_propagation(prop_task)

    if run_enrichment:
        perform_enrichment(prop_task)


if __name__ == '__main__':
    start = time.time()

    # Set these flags to control the tasks to run
    run_propagation_flag = True
    run_enrichment_flag = True

    main(run_propagation=run_propagation_flag, run_enrichment=run_enrichment_flag)

    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))

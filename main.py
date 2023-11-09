from args import PropagationTask, EnrichTask, GeneralArgs
import utils
from propagation_routines import propagate_network, generate_similarity_matrix, read_sparse_matrix_txt
from statistic_methods import wilcoxon_rank_sums_test
from pathway_enrichment import run
import numpy as np
import time


def perform_propagation(task):
    # reads the network graph from a file
    network = utils.read_network(task.network_file_path)

    # load prior set
    prior_data_df = utils.read_prior_set(task.experiment_file_path)
    prior_data = set.intersection(set(prior_data_df.GeneID), set(network.nodes))

    all_genes_ids = set.intersection(set(list(prior_data)), set(network.nodes))

    print("getting propagation input")
    propagation_input = utils.get_propagation_input(all_genes_ids, prior_data_df, task.propagation_input_type,
                                                    network=network)
    # print("getting ones input")
    ones_input = utils.get_propagation_input(all_genes_ids, prior_data_df, 'ones', network=network)

    # create or upload similarity matrix
    if task.create_similarity_matrix:
        print("generating similarity matrix")
        matrix, genes = generate_similarity_matrix(network, task.similarity_matrix_path,
                                                   task.alpha)
    else:
        print("reading similarity matrix")
        matrix, genes = read_sparse_matrix_txt(network, task.similarity_matrix_path)
        print("uploaded similarity matrix")

    print("propagating")
    # Propagate network
    score_genes_id_to_idx, score_gene_scores_inverse, gene_score_dict = propagate_network(propagation_input, matrix, genes)
    ones_genes_id_to_idx, ones_gene_scores_inverse, ones_gene_score_dict = propagate_network(ones_input, matrix, genes)
    score_genes_idx_to_id = {xx: x for x, xx in score_genes_id_to_idx.items()}

    # Identify genes with zero normalization score but non-zero propagation score
    zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    zero_propagation_genes = np.nonzero(score_gene_scores_inverse == 0)[0]
    genes_to_delete = list(set(zero_normalization_genes).difference(zero_propagation_genes))

    # Set the normalization score of these genes to 1
    ones_gene_scores_inverse[genes_to_delete] = 1

    # Perform the normalization
    non_zero_indices = np.nonzero(score_gene_scores_inverse != 0)[0]
    score_gene_scores_inverse[non_zero_indices] = score_gene_scores_inverse[non_zero_indices] / np.abs(
        ones_gene_scores_inverse[non_zero_indices])

    # # Identify genes with zero normalization score but non-zero propagation score
    # zero_normalization_genes = np.nonzero(ones_gene_scores_inverse == 0)[0]
    # non_zero_normalization_genes = np.nonzero(ones_gene_scores_inverse != 0)[0]
    #
    # # Calculate the mean of non-zero normalization scores
    # mean_non_zero_normalization = np.mean(ones_gene_scores_inverse[non_zero_normalization_genes])
    #
    # # Set the normalization score of zero normalization genes to the mean of non-zero genes
    # ones_gene_scores_inverse[zero_normalization_genes] = mean_non_zero_normalization
    #
    # # Perform the normalization
    # non_zero_indices = np.nonzero(score_gene_scores_inverse != 0)[0]
    # normalized_scores = score_gene_scores_inverse.copy()  # Make a copy to keep the original scores
    # normalized_scores[non_zero_indices] /= ones_gene_scores_inverse[non_zero_indices]

    # save propagation score
    print("saving propagation score")
    utils.save_propagation_score(propagation_scores=score_gene_scores_inverse, prior_set=prior_data,
                                 propagation_input=propagation_input, genes_idx_to_id=score_genes_idx_to_id,
                                 task=task, save_dir=task.output_folder, date=task.date)


def perform_enrichment(task):
    # run enrichment
    print("running enrichment")
    tasks = []
    propagation_scores_file = '{}_{}_{}_{}'.format(task.experiment_name, task.propagation_input_type,
                                                   task.alpha, task.date)
    task1 = EnrichTask(name=task.experiment_name, propagation_file=propagation_scores_file,
                       propagation_folder=f'Outputs\\propagation_scores\\{task.experiment_name}',
                       statistic_test=wilcoxon_rank_sums_test,
                       target_field='gene_prop_scores', constrain_to_experiment_genes=True)

    # task1 = EnrichTask(name='500nm_v_T', propagation_file='500nm_v_T_Score_0.1_08_11_2023__12_22_59',
    #                    propagation_folder=f'Outputs\\propagation_scores\\500nm_v_T',
    #                    statistic_test=wilcoxon_rank_sums_test,
    #                    target_field='gene_prop_scores', constrain_to_experiment_genes=True)
    #
    # task2 = EnrichTask(name='TvN', propagation_file='TvN_Score_0.1_07_11_2023__16_17_29',
    #                    propagation_folder=f'Outputs\\propagation_scores\\TvN',
    #                    statistic_test=wilcoxon_rank_sums_test,
    #                    target_field='gene_prop_scores', constrain_to_experiment_genes=True)

    FDR_threshold = 0.01

    figure_name = task.experiment_name + '-alpha' + str(
        task.alpha) + '-Threshold' + str(FDR_threshold) + '.pdf'

    general_args = GeneralArgs(task.network_file_path, genes_names_path=task.genes_names_file_path,
                               pathway_members_path=task.pathway_file_dir, FDR_threshold=FDR_threshold,
                               figure_name=figure_name)

    tasks += [task1]
    print('running')
    run(tasks, general_args)


def main(run_propagation=True, run_enrichment=True):
    # create a propagation task
    task = PropagationTask(experiment_name='TvN', create_similarity_matrix=False)

    if run_propagation:
        perform_propagation(task)

    if run_enrichment:
        perform_enrichment(task)


if __name__ == '__main__':
    start = time.time()

    # Set these flags to control the tasks to run
    run_propagation_flag = True
    run_enrichment_flag = True

    main(run_propagation=run_propagation_flag, run_enrichment=run_enrichment_flag)

    end = time.time()
    print("Time elapsed: {} seconds".format(end - start))

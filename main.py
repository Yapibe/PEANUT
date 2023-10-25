from args import PropagationTask, EnrichTask, GeneralArgs
import utils
from propagation_routines import propagate_network, generate_similarity_matrix, read_sparse_matrix_txt
from statistic_methods import wilcoxon_rank_sums_test


def main(full_calculation=True):
    # create a propagation task
    task = PropagationTask(create_similarity_matrix=False)

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
    score_genes_id_to_idx, score_gene_scores = propagate_network(network, propagation_input, task, matrix, genes)
    ones_genes_id_to_idx, ones_gene_scores = propagate_network(network, ones_input, task, matrix, genes)
    del matrix
    # genes_idx_to_id = {xx: x for x, xx in genes_id_to_idx.items()}
    # Replace zeros in ones_gene_scores with a small constant
    ones_gene_scores[ones_gene_scores == 0] = 1e-10

    # Normalize the actual propagated scores
    normalized_gene_scores = score_gene_scores / ones_gene_scores

    # save propagation score
    print("saving propagation score")
    utils.save_propagation_score(propagation_scores=normalized_gene_scores, prior_set=prior_data,
                                 propagation_input=propagation_input, genes_idx_to_id=score_genes_id_to_idx,
                                 task=task, save_dir=task.output_folder, date=task.date)

    if not full_calculation:
        return

    # run enrichment
    print("running enrichment")
    file_name = f"{task.experiment_name}_{task.propagation_input_type}_{task.alpha}_{task.date}"

    propagation_scores_file = '{}_{}_{}_{}_IPN'.format(task.experiment_name, task.propagation_input_type,
                                                       task.alpha, task.date)
    task1 = EnrichTask(name=task.experiment_name, propagation_file=propagation_scores_file,
                       propagation_folder='propagation_scores', statistic_test=wilcoxon_rank_sums_test,
                       target_field='randomization_ranks', constrain_to_experiment_genes=True)

    figure_name = task.experiment_file.split('.')[
                      0] + '_sortedEnrichment_Prop-' + 'alpha' + task.alpha + '-' + '.pdf'
    general_args = GeneralArgs(task.network_file_path, genes_names_path=task.genes_names_file_path,
                               pathway_members_path=task.pathway_file_dir, figure_name=figure_name)

    print('running')


if __name__ == '__main__':
    main()

from os import path
import numpy as np
from args import EnrichTask, PathwayResults
from statistic_methods import bh_correction
from utils import load_propagation_scores
from utils import read_network, load_pathways_genes
from visualization_tools import plot_enrichment_table


def get_scores(task):
    """
    takes a task, a network graph and a general args object and returns a dictionary of scores
    :param task:
    :return:
    """

    if isinstance(task, EnrichTask):
        result_dict = load_propagation_scores(task, propagation_file_name=task.propagation_file)
        gene_id_to_idx = {xx: x for x, xx in result_dict['gene_idx_to_id'].items()}
        scores = result_dict[task.target_field]
        if task.constrain_to_experiment_genes:
            gene_id_to_idx = {id: idx for id, idx in gene_id_to_idx.items() if id in result_dict['propagation_input']}
        scores = {id: scores[idx] for id, idx in gene_id_to_idx.items()}
    else:
        raise ValueError('Invalid task')
    return scores


def load_network_and_pathways(general_args):
    network_graph = read_network(general_args.network_file_path)
    genes_by_pathway = load_pathways_genes(general_args.pathway_members_path)
    interesting_pathways = list(genes_by_pathway.keys())

    return network_graph, interesting_pathways, genes_by_pathway


def process_tasks(task_list, network_graph, general_args, interesting_pathways, genes_by_pathway):
    pathways_to_display = set()
    all_genes_by_pathway_filtered = {}

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

        # Update all_genes_in_filtered_pathways_and_network
        for pathway in pathways_with_many_genes:
            all_genes_in_filtered_pathways_and_network.update(genes_by_pathway_filtered[pathway])

        # Intersect with network nodes to refine the set
        all_genes_in_filtered_pathways_and_network &= set(network_graph.nodes)


        # Perform statistical tests
        print('after filtering', len(pathways_with_many_genes))
        for pathway in pathways_with_many_genes:
            pathways_to_display.add(pathway)
            pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
            # background_scores1 = [scores[id] for id in scores.keys() if id not in genes_by_pathway_filtered[pathway]]
            background_scores = [scores[id] for id in all_genes_in_filtered_pathways_and_network if
                                   id not in genes_by_pathway_filtered[pathway]]
            result = task.statistic_test(pathway_scores, background_scores)
            task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality,
                                                   score=np.mean(pathway_scores))

        # Store the filtered genes by pathway for this task
        all_genes_by_pathway_filtered[task.name] = genes_by_pathway_filtered

    pathways_to_display = np.sort(list(pathways_to_display))

    return pathways_to_display, all_genes_by_pathway_filtered, pathways_with_many_genes


def create_matrices(task_list, pathways_to_display):
    pathways_to_display = np.sort(list(pathways_to_display))
    p_vals_mat = np.ones((len(pathways_to_display), len(task_list)))
    adj_p_vals_mat = np.ones_like(p_vals_mat)
    directions_mat = np.zeros_like(p_vals_mat)
    coll_names_in_heatmap = []

    for t, task in enumerate(task_list):
        indexes = []
        for p, pathway in enumerate(pathways_to_display):
            if pathway in task.results:
                indexes.append(p)
                p_vals_mat[p, t] = task.results[pathway].p_value
                directions_mat[p, t] = task.results[pathway].direction
        adj_p_vals_mat[indexes, t] = bh_correction(p_vals_mat[indexes, t])
        coll_names_in_heatmap.append(task.name)

    return p_vals_mat, adj_p_vals_mat, directions_mat, coll_names_in_heatmap


def filter_and_adjust_matrices(p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display, general_args):
    if general_args.display_only_significant_pathways:
        keep_rows = np.nonzero(np.any(adj_p_vals_mat <= general_args.significant_pathway_threshold, axis=1))[0]
        pathways_to_display = list(pathways_to_display)  # Convert set to list
        pathways_to_display = [pathways_to_display[x] for x in keep_rows]  # Now this should work

        p_vals_mat = p_vals_mat[keep_rows, :]
        adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
        directions_mat = directions_mat[keep_rows, :]

    return p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display


def prepare_row_names(pathways_to_display, genes_by_pathway_filtered):
    row_names = ['({}) {}'.format(len(genes_by_pathway_filtered[pathway]), pathway) for pathway in pathways_to_display]
    row_names = [x.replace("-", "-\n") for x in row_names]
    return row_names


def filter_by_minimum_p_values(p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display, general_args):
    candidates = np.min(p_vals_mat, axis=1)
    ind = np.sort(
        np.argpartition(candidates, general_args.maximum_number_of_pathways)[:general_args.maximum_number_of_pathways])
    p_vals_mat = p_vals_mat[ind, :]
    adj_p_vals_mat = adj_p_vals_mat[ind, :]
    directions_mat = directions_mat[ind, :]
    pathways_to_display = [pathways_to_display[x] for x in ind]
    return p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display


def plot_results(p_vals_mat, adj_p_vals_mat, directions_mat, row_names, general_args, coll_names_in_heatmap,
                 dataset_type, n_pathways_before):
    res = -np.log10(p_vals_mat)
    fig_out_dir = path.join(general_args.output_path, general_args.figure_name)
    row_names = [x.replace("_", " ") for x in row_names]
    plot_enrichment_table(res, adj_p_vals_mat, directions_mat, row_names, fig_out_dir,
                          experiment_names=coll_names_in_heatmap,
                          title=general_args.figure_title + ' {} {}/{}'.format(dataset_type,
                                                                               len(row_names), n_pathways_before),
                          res_type='-log10(p_val)', adj_p_value_threshold=general_args.significant_pathway_threshold)


def run(task_list, general_args, dataset_type=''):
    """
    takes a list of tasks, a general args object and an optional dataset type and runs the pathway enrichment analysis
    :param task_list: list of tasks
    :param general_args: general args object
    :param dataset_type: optional dataset type
    :return: None
    """
    network_graph, interesting_pathways, genes_by_pathway = load_network_and_pathways(general_args)

    pathways_to_display, genes_by_pathway_filtered, pathways_with_many_genes = (process_tasks(task_list, network_graph,
                                                                                             general_args,
                                                                                             interesting_pathways,
                                                                                             genes_by_pathway))

    # Create matrices for p-values, adjusted p-values, and directions
    p_vals_mat, adj_p_vals_mat, directions_mat, coll_names_in_heatmap = create_matrices(task_list, pathways_to_display)

    # Filter and adjust matrices and pathway lists based on various conditions
    p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display = filter_and_adjust_matrices(p_vals_mat,
                                                                                                 adj_p_vals_mat,
                                                                                                 directions_mat,
                                                                                                 pathways_to_display,
                                                                                                 general_args)

    row_names = ['{}'.format(pathway) for pathway in pathways_to_display]
    # Plot the results
    plot_results(p_vals_mat, adj_p_vals_mat, directions_mat, row_names, general_args, coll_names_in_heatmap,
                 dataset_type, len(pathways_with_many_genes))

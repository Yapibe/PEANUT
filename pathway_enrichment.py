from os import path
import numpy as np
from args import EnrichTask, PathwayResults
from statistic_methods import bh_correction
from utils import load_propagation_scores
from utils import read_network, load_pathways_genes
from visualization_tools import plot_enrichment_table
import gseapy as gp
from typing import List, Dict
import pandas as pd
import csv


def run_gsea(task, ranked_genes: Dict[str, float], gene_sets: Dict[str, List[str]]):
    """
    Executes Gene Set Enrichment Analysis (GSEA) using the gseapy library.

    Args:
        task (Task): The task object containing relevant settings and results.
        ranked_genes (Dict[str, float]): Dictionary mapping gene names to their ranking scores.
        gene_sets (Dict[str, List[str]]): Dictionary where keys are pathway names and values are lists of gene names.

    Modifies:
        task.results (Dict): Updates the task's results dictionary with GSEA results for each pathway.

    Side Effects:
        Writes GSEA results to a CSV file.

    Returns:
        None
    """

    # Convert dictionary to Series, ensuring it's of float data type
    ranked_gene_series = pd.Series(ranked_genes, dtype=float)

    ranked_gene_series.index = ranked_gene_series.index.astype(str)
    gene_sets = {str(k): [str(g) for g in v] for k, v in gene_sets.items()}

    # Run GSEA
    enrichment_results = gp.prerank(rnk=ranked_gene_series, gene_sets=gene_sets, min_size=1, max_size=5000, outdir=None,
                                    no_plot=True)
    results_df = pd.DataFrame(enrichment_results.res2d)
    # save to csv
    results_df.to_csv('gsea_results.csv')

    # Extract and return the p-values for each pathway
    for index, row in results_df.iterrows():
        p_value = row['NOM p-val']
        adj_p_value = row['FDR q-val']
        direction = row['ES'] > 0
        # add to task results by pathway
        task.results[row['Term']] = PathwayResults(p_value=p_value, direction=direction, adj_p_value=adj_p_value)


def get_scores(task):
    """
    Retrieves gene scores based on the given task configuration.

    Args:
        task (EnrichTask or RawScoreTask): Task object specifying the type of analysis and associated parameters.

    Raises:
        ValueError: If the task type is not supported.

    Returns:
        Dict[int, float]: Dictionary mapping gene IDs to their corresponding scores.
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


def load_network_and_pathways(general_args):
    """
    Loads the network graph and pathways based on the provided configuration.
    Args:
        general_args (GeneralArgs): Object containing general configuration settings.
    Returns:
        tuple: A tuple containing the network graph, a list of interesting pathways, and a dictionary mapping
               pathways to their genes.
    """
    network_graph = read_network(general_args.network_file_path)
    genes_by_pathway = load_pathways_genes(general_args.pathway_members_path)
    interesting_pathways = list(genes_by_pathway.keys())

    return network_graph, interesting_pathways, genes_by_pathway


def process_tasks(task_list, network_graph, general_args, interesting_pathways, genes_by_pathway):
    """
    Processes a list of tasks for pathway enrichment analysis.
    Args:
        task_list (List[Task]): List of tasks to be processed.
        network_graph (Graph): Graph representing the network.
        general_args (GeneralArgs): General configuration settings.
        interesting_pathways (List[str]): List of pathways of interest.
        genes_by_pathway (Dict[str, List[int]]): Dictionary mapping pathways to their gene IDs.
    Returns:
        Set[str]: Set of pathways to be displayed in the analysis.
    """
    pathways_to_display = set()
    pathways_with_many_genes = []

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
        # manually add REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION'
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

        if general_args.use_gsea:  # Assume you add this flag to general_args or similar
            # Replace with how you get your ranked genes
            run_gsea(task, scores, genes_by_pathway)
        else:
            # Perform statistical tests
            print('after filtering', len(pathways_with_many_genes))
            for pathway in pathways_with_many_genes:
                pathway_scores = [scores[id] for id in genes_by_pathway_filtered[pathway]]
                # background_scores1 = [scores[id] for id in scores.keys() if id not in genes_by_pathway_filtered[pathway]]
                background_scores = [scores[id] for id in all_genes_in_filtered_pathways_and_network if
                                     id not in genes_by_pathway_filtered[pathway]]
                result = task.statistic_test(pathway_scores, background_scores)
                task.results[pathway] = PathwayResults(p_value=result.p_value, direction=result.directionality)

    pathways_to_display = np.sort(list(pathways_to_display))

    return pathways_to_display


def process_matrices(task_list, pathways_to_display, general_args):
    """
    Processes the matrices for p-values, adjusted p-values, and directions based on task results.
    Args:
        task_list (List[Task]): List of tasks with their respective results.
        pathways_to_display (Set[str]): Set of pathways to be included in the analysis.
        general_args (GeneralArgs): General configuration settings.
    Returns:
        tuple: A tuple containing matrices for adjusted p-values, directions, pathway names, and column names
               for the heatmap.
    """
    p_vals_mat = np.ones((len(pathways_to_display), len(task_list)))
    adj_p_vals_mat = np.ones_like(p_vals_mat)
    directions_mat = np.zeros_like(p_vals_mat)
    coll_names_in_heatmap = []

    # Create matrices
    for t, task in enumerate(task_list):
        indexes = []
        for p, pathway in enumerate(pathways_to_display):
            if pathway in task.results:
                if general_args.use_gsea:
                    # Directly use the adjusted p-value from GSEA results
                    adj_p_vals_mat[p, t] = task.results[pathway].adj_p_value
                    p_vals_mat[p, t] = task.results[pathway].p_value
                else:
                    # Use p-value for non-GSEA results and collect indexes for BH correction
                    p_vals_mat[p, t] = task.results[pathway].p_value
                    indexes.append(p)
                    # Direction should be collected in both cases
                directions_mat[p, t] = task.results[pathway].direction
        if not general_args.use_gsea:
            # Apply BH correction outside the inner loop
            adj_p_vals_mat[indexes, t] = bh_correction(p_vals_mat[indexes, t])
            # Append task name to heatmap column names
        coll_names_in_heatmap.append(task.name)
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
    ]  # List of specific pathways
    # Create a list for CSV rows
    csv_rows = []
    for pathway in specific_pathways:
        if pathway in pathways_to_display:
            idx = np.where(pathways_to_display == pathway)[0][0]  # Find the index of the pathway
            for t in range(len(task_list)):
                # Append information to the rows list
                csv_rows.append([
                    pathway,
                    task_list[t].name,
                    p_vals_mat[idx, t],
                    adj_p_vals_mat[idx, t]
                ])

    # Writing to a CSV file
    with open('pathway_analysis_results_TvN_mann.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(['Pathway', 'Task Name', 'P-Value', 'Adjusted P-Value'])
        # Write data
        writer.writerows(csv_rows)

    # Count pathways with adjusted p-value < 0.05
    count_significant_pathways = np.sum(adj_p_vals_mat < 0.05)

    # Filter and adjust matrices
    keep_rows = np.nonzero(np.any(adj_p_vals_mat <= general_args.significant_pathway_threshold, axis=1))[0]
    pathways_to_display = list(pathways_to_display)  # Convert set to list
    pathways_to_display = [pathways_to_display[x] for x in keep_rows]
    #

    adj_p_vals_mat = adj_p_vals_mat[keep_rows, :]
    directions_mat = directions_mat[keep_rows, :]

    return adj_p_vals_mat, directions_mat, pathways_to_display, coll_names_in_heatmap


def filter_by_minimum_p_values(p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display, general_args):
    """
    Filters pathways based on minimum p-values and limits the number of pathways to display.
    Args:
        p_vals_mat (numpy.ndarray): Matrix of p-values.
        adj_p_vals_mat (numpy.ndarray): Matrix of adjusted p-values.
        directions_mat (numpy.ndarray): Matrix indicating the direction of changes.
        pathways_to_display (List[str]): List of pathway names to be considered.
        general_args (args): General configuration settings.
    Returns:
        tuple: Filtered matrices for p-values, adjusted p-values, directions, and a list of pathway names.
    """
    candidates = np.min(p_vals_mat, axis=1)
    ind = np.sort(
        np.argpartition(candidates, general_args.maximum_number_of_pathways)[:general_args.maximum_number_of_pathways])
    p_vals_mat = p_vals_mat[ind, :]
    adj_p_vals_mat = adj_p_vals_mat[ind, :]
    directions_mat = directions_mat[ind, :]
    pathways_to_display = [pathways_to_display[x] for x in ind]
    return p_vals_mat, adj_p_vals_mat, directions_mat, pathways_to_display


def plot_results(adj_p_vals_mat, directions_mat, row_names, general_args, coll_names_in_heatmap, dataset_type,
                 n_pathways_before):
    """
    Plots the results of pathway enrichment analysis.
    Args:
        adj_p_vals_mat (numpy.ndarray): Matrix of adjusted p-values.
        directions_mat (numpy.ndarray): Matrix indicating the direction of changes.
        row_names (List[str]): Names of the pathways (rows).
        general_args (GeneralArgs): General configuration settings.
        coll_names_in_heatmap (List[str]): Names of the columns in the heatmap.
        dataset_type (str): Type of the dataset.
        n_pathways_before (int): Number of pathways before filtering.
    Side Effects:
        Generates and saves a heatmap plot.
    Returns:
        None
    """
    # Set a small value to represent the minimum p-value
    min_p_val = np.finfo(adj_p_vals_mat.dtype).tiny

    # Replace 0 with the minimum p-value
    adj_p_vals_mat_no_zeros = np.where(adj_p_vals_mat == 0, min_p_val, adj_p_vals_mat)

    # Calculate -log10, avoiding log10(0) which gives inf
    res = -np.log10(adj_p_vals_mat_no_zeros)
    fig_out_dir = path.join(general_args.output_path, general_args.figure_name)
    plot_enrichment_table(res, directions_mat, row_names, fig_out_dir, experiment_names=coll_names_in_heatmap,
                          title=general_args.figure_title + ' {} {}/{}'.format(dataset_type,
                                                                               len(row_names), n_pathways_before),
                          res_type='-log10(p_val)', adj_p_value_threshold=general_args.significant_pathway_threshold)


def run(task_list, general_args, dataset_type=''):
    """
    Main function to run the pathway enrichment analysis.

    Args:
        task_list (List[Task]): List of tasks to be processed.
        general_args (args): General configuration settings.
        dataset_type (str, optional): Type of the dataset. Defaults to an empty string.

    Side Effects:
        Executes the entire pathway enrichment analysis workflow, including data loading, processing,
        analysis, and plotting results.

    Returns:
        None
    """
    network_graph, interesting_pathways, genes_by_pathway = load_network_and_pathways(general_args)

    pathways_to_display = (
        process_tasks(task_list, network_graph, general_args, interesting_pathways, genes_by_pathway))

    # Create matrices for p-values, adjusted p-values, and directions
    adj_p_vals_mat, directions_mat, pathways_to_display, coll_names_in_heatmap = process_matrices(task_list,
                                                                                                  pathways_to_display,
                                                                                                  general_args)

    row_names = ['{}'.format(pathway.replace("_", " ")) for pathway in pathways_to_display]

    # Plot the results
    plot_results(adj_p_vals_mat, directions_mat, row_names, general_args, coll_names_in_heatmap, dataset_type,
                 len(pathways_to_display))

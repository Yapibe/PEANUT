import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle, Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D
from itertools import combinations
import networkx as nx
from matplotlib import colors
from os import path


def sort_y_ticks_by_dir_and_pvalue_of_one_column(enrichment_table, y_ticks, annotation_map, direction, adj_p_mat,
                                                 dir_column_to_sort_by=0):
    # 1. sort 0s from 1s in direction table
    sorted_indexes1 = np.argsort(direction[:, dir_column_to_sort_by])
    enrichment_table = enrichment_table[sorted_indexes1, :]
    y_ticks = list(np.array(y_ticks)[sorted_indexes1])
    annotation_map = annotation_map[sorted_indexes1, :]
    direction = direction[sorted_indexes1, :]
    adj_p_mat = adj_p_mat[sorted_indexes1, :]

    # 2. sort p_value within indexes with sorted direction (negative first, positives after)
    negatives = len(direction[direction[:, dir_column_to_sort_by] == 0, dir_column_to_sort_by])

    negative_enrichment_table = enrichment_table[:negatives, ]
    positive_enrichment_table = enrichment_table[negatives:, ]
    negative_annotation_map = annotation_map[:negatives, ]
    positive_annotation_map = annotation_map[negatives:, ]
    negative_direction = direction[:negatives, ]
    positive_direction = direction[negatives:, ]
    negative_adj_p_mat = adj_p_mat[:negatives, ]
    positive_adj_p_mat = adj_p_mat[negatives:, ]
    y_ticks_n = list(np.array(y_ticks)[:negatives])
    y_ticks_p = list(np.array(y_ticks)[negatives:])

    sorted_positive_indexes = np.argsort(positive_enrichment_table[:, dir_column_to_sort_by])
    sorted_negative_indexes = np.argsort(-negative_enrichment_table[:, dir_column_to_sort_by])  # sort descending

    sorted_positive_enrichment_table = positive_enrichment_table[sorted_positive_indexes, :]
    sorted_negative_enrichment_table = negative_enrichment_table[sorted_negative_indexes, :]
    sorted_positive_annotation_map = positive_annotation_map[sorted_positive_indexes, :]
    sorted_negative_annotation_map = negative_annotation_map[sorted_negative_indexes, :]
    sorted_positive_direction = positive_direction[sorted_positive_indexes, :]
    sorted_negative_direction = negative_direction[sorted_negative_indexes, :]
    sorted_positive_adj_p_mat = positive_adj_p_mat[sorted_positive_indexes, :]
    sorted_negative_adj_p_mat = negative_adj_p_mat[sorted_negative_indexes, :]
    y_ticks_ns = list(np.array(y_ticks_n)[sorted_negative_indexes])
    y_ticks_ps = list(np.array(y_ticks_p)[sorted_positive_indexes])

    # 3. return joined objects, with negatives first.
    return np.concatenate((sorted_negative_enrichment_table, sorted_positive_enrichment_table),
                          axis=0), y_ticks_ns + y_ticks_ps, np.concatenate(
        (sorted_negative_annotation_map, sorted_positive_annotation_map), axis=0), np.concatenate(
        (sorted_negative_direction, sorted_positive_direction), axis=0), np.concatenate(
        (sorted_negative_adj_p_mat, sorted_positive_adj_p_mat), axis=0)


def sort_y_ticks_by_dir(enrichment_table, y_ticks, annotation_map, direction):
    # warning NO ADJ_P_MAT
    simplified_dir = []
    for row in direction:
        if np.sum(row) > 1:
            simplified_dir.append(1)
        else:
            simplified_dir.append(-1)
    sorted_indexes = np.argsort(simplified_dir)

    return enrichment_table[sorted_indexes, :], list(np.array(y_ticks)[sorted_indexes]), annotation_map[sorted_indexes,
                                                                                         :], direction[sorted_indexes,
                                                                                             :]


def sort_y_ticks(enrichment_table, y_ticks, annotation_map):
    # works (tested reversing the orded)
    # warning NO ADJ_P_MAT

    def custom_sort(list_of_names):
        sorted_list = list_of_names
        return sorted_list

    old_ind_y_tick = {tick: y_ticks.index(tick) for tick in y_ticks}
    nonum = [x.split()[1] + ' ' + x.split()[0] for x in y_ticks]
    srtd_y_ticks = [x.split()[1] + ' ' + x.split()[0] for x in custom_sort(nonum)]
    new_ind = [old_ind_y_tick[tick] for tick in srtd_y_ticks]
    return enrichment_table[new_ind, :], srtd_y_ticks, annotation_map[new_ind, :]


def plot_enrichment_table(enrichment_table, adj_p_mat, direction, interesting_pathways, save_dir=None,
                          experiment_names=None,
                          title=None, res_type=None, adj_p_value_threshold=0.03):
    fig, ax = plt.subplots()

    enriched_clusters = np.nonzero(np.sum(enrichment_table, axis=0) != 0)[0]
    found_pathways = np.nonzero(np.sum(enrichment_table, axis=1) != 0)[0]
    enrichment_table = enrichment_table[:, enriched_clusters][found_pathways, :]
    interesting_pathways_filtered = {x: xx for x, xx in enumerate(interesting_pathways) if x in found_pathways}
    annotation_map = (np.round(enrichment_table, 3)).astype(str)
    annotation_map[annotation_map == '0.0'] = ''
    y_ticks = [x[:60] for x in interesting_pathways_filtered.values()]
    enrichment_table, y_ticks, annotation_map, direction, adj_p_mat = sort_y_ticks_by_dir_and_pvalue_of_one_column(
        enrichment_table, y_ticks, annotation_map, direction, adj_p_mat)
    important_indexes = np.where(adj_p_mat < 0.05)

    if direction is not None and res_type != 'z_score':
        print("Shape of enrichment_table:", enrichment_table.shape)
        print("Shape of direction:", direction.shape)

        # set low propagation scores to be negative in order to color them blue
        enrichment_table[np.logical_not(direction)] = -enrichment_table[np.logical_not(direction)]
        # enrichment_table[np.logical_not(direction[:, 0]), 0] = -enrichment_table[np.logical_not(direction[:, 0]), 0]
        # enrichment_table[np.logical_not(direction[:, 2]), 2] = -enrichment_table[np.logical_not(direction[:, 2]), 2]
        # enrichment_table[np.logical_not(direction[:, 3]), 3] = -enrichment_table[np.logical_not(direction[:, 3]), 3]
        # enrichment_table[np.logical_not(direction[:, 5]), 5] = -enrichment_table[np.logical_not(direction[:, 5]), 5]

    # set color bar size and location
    cax = inset_axes(ax,
                     width="100%",  # width: 40% of parent_bbox width
                     height="70%",  # height: 10% of parent_bbox height
                     loc='lower left',
                     bbox_to_anchor=(1.1, 0.2, 0.2, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )
    colorbar_edge = np.maximum(np.abs(np.min(enrichment_table)), np.abs(np.max(enrichment_table)))
    heatmap = sns.heatmap(enrichment_table, fmt=".4s", yticklabels=y_ticks,
                          cbar_ax=cax, annot=annotation_map, cmap="coolwarm",
                          linewidths=.1, linecolor='gray',
                          cbar_kws={'label': res_type}, ax=ax, vmin=-20,
                          vmax=20)  # vmin=np.min([0, -colorbar_edge]), vmax=colorbar_edge)
    ax.set_xticklabels(experiment_names, fontsize=24)

    # circle significant scores (<0.05)
    for i in range(len(important_indexes[0])):
        heatmap.add_patch(
            Rectangle((important_indexes[1][i], important_indexes[0][i]), 1, 1, fill=False, edgecolor='black', lw=3))
    # set font size of colobar ticks
    heatmap.figure.axes[-1].yaxis.label.set_size(20)

    # set colorbar tick values to be positive in both direction
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(np.round(np.abs(cbar.get_ticks()), 1))

    # rotate test names
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')

    # add square color legend
    legend_handles = [Patch(fill=False, edgecolor='black', label=f'adj_p_value<{adj_p_value_threshold}')]
    if direction is not None:
        legend_handles.append(Patch(fill=True, color='red', label='Upwards'))
        legend_handles.append(Patch(fill=True, color='blue', label='Downwards'))

    # locate legend
    ax.legend(handles=legend_handles, bbox_to_anchor=[1, 0, 1, 1], ncol=1, loc='lower left', fontsize=14)
    ax.set_title(title, fontsize=28)
    sns.set(font_scale=0.9)
    fig.set_size_inches(20, 40)
    plt.savefig(save_dir, bbox_inches='tight')
    sns.set(font_scale=1 / 0.9)


def visualise_pathway(network_graph, pathway_genes, reference_scores, propagation_scores, pathway_name,
                      id_to_label_dict=None,
                      mark_second_neighbors=True, significant_genes=None, save_dir=None):
    """

    :param network_graph: networkx Graph of the whole network
    :param pathway_genes: a list of gene ids in the pathway
    :param reference_scores: a dictionary of the format gene_id: score of the reference score
    :param propagation_scores: a dictionary of the format gene_id: score of the propagation score
    :param pathway_name: pathway name, string
    :param id_to_label_dict: a diction of the format gene_id: label
    :param mark_second_neighbors: whether to mark second order neighbors with dashed lines
    :param significant_genes: a list of significant gene ids
    :param save_dir: directory to save figure
    """
    # generate pathway subgraph
    G = nx.Graph()
    labels = dict()
    for gene in pathway_genes:
        G.add_node(gene)
        if id_to_label_dict:
            labels[gene] = id_to_label_dict[gene]
        else:
            labels[gene] = gene

    edges = []
    for edge in combinations(pathway_genes, 2):
        if edge[0] < edge[1]:
            if edge in network_graph.edges():
                G.add_edge(edge[0], edge[1], weight=1)
                edges.append((edge[0], edge[1]))

    if mark_second_neighbors:
        secondary_edges = set()
        for node in G.nodes:
            neighbors = nx.single_source_shortest_path_length(network_graph, node, cutoff=2)
            neighbors = set(neighbors).intersection(G.nodes)
            for neighbor in neighbors:
                if neighbor in G.nodes and (node < neighbor) and ((node, neighbor) not in G.edges):
                    G.add_edge(node, neighbor, weight=1.5)
                    secondary_edges.add((node, neighbor))

    for node in list(G.nodes):
        if len(nx.single_source_shortest_path_length(G, node, cutoff=1)) == 1:
            G.remove_node(node)
            labels.pop(node)

    pos = nx.kamada_kawai_layout(G, )
    pos_labels = {id: np.array([coord[0] + 0.075, coord[1]]) for id, coord in pos.items()}
    # nx.drawing.lay
    gene_propagation_scores = [propagation_scores[id] for id in G.nodes]
    path_reference_scores = {id: reference_scores[id] for id in G.nodes if id in reference_scores}

    if significant_genes:
        significant_nodes = [node for node in G.nodes if id_to_label_dict[node] in significant_genes]
    else:
        significant_genes = []

    # draw graph
    colorbar_edge = np.maximum(np.abs(np.min(gene_propagation_scores)), np.abs(np.max(gene_propagation_scores)))
    vmin = np.min([0, -colorbar_edge])
    vmax = colorbar_edge

    sm = plt.cm.ScalarMappable(cmap='coolwarm', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    border_colors = [sm.to_rgba(path_reference_scores[id]) if id in reference_scores else colors.to_rgba('lightgreen')
                     for id in G.nodes]

    # draw networkx
    fig, ax = plt.subplots()
    nx.draw_networkx(G, pos, node_color=gene_propagation_scores, node_size=600, cmap='coolwarm', vmin=vmin, vmax=vmax,
                     linewidths=6, edgecolors=border_colors, with_labels=False, ax=ax, edgelist=edges)
    nx.draw_networkx_labels(G, pos_labels, labels=labels)

    # mark significant genes
    nx.draw_networkx(G, pos, nodelist=significant_nodes, node_size=600, node_color='none',
                     labels={id: 'X' for id in significant_nodes}, with_labels=True, ax=ax, edgelist=[])

    if mark_second_neighbors:
        nx.draw_networkx(G, pos, edgelist=secondary_edges, node_size=600, node_color='none', edgecolors='none',
                         style='dashed', with_labels=False, alpha=0.6)

    # draw color bar
    sm._A = []
    cax = inset_axes(ax,
                     width="50%",  # width: 40% of parent_bbox width
                     height="70%",  # height: 10% of parent_bbox height
                     loc='lower left',
                     bbox_to_anchor=(1.025, 0.2, 0.2, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0,
                     )
    cb = plt.colorbar(sm, cax=cax, )
    cb.set_label(label='Normalized Score Rank', fontsize=16)

    legend_handles = [Line2D([0], [0], marker='o', color='w', label='Exp absent gene',
                             markerfacecolor=None, markeredgecolor='lightgreen', markersize=15, markeredgewidth=2),
                      Line2D([0], [0], marker='o', color='w', label='log2FC score',
                             markerfacecolor=None, markeredgecolor='blue', markersize=15, markeredgewidth=2),
                      Line2D([0], [0], marker='o', color='w', label='Propagation score',
                             markerfacecolor='blue', markeredgecolor='None', markersize=15),
                      Line2D([0], [0], marker='X', color='w', label='Exp significant gene',
                             markerfacecolor='black', markeredgecolor='None', markersize=15),
                      Line2D([0], [0], linestyle='-', color='black', label='Network edge',
                             markerfacecolor='black', markeredgecolor='None', markersize=15),
                      Line2D([0], [0], linestyle='--', color='black', label='Second order neighbor',
                             markerfacecolor='black', markeredgecolor='None', markersize=15),
                      ]

    ax.legend(handles=legend_handles, bbox_to_anchor=[1, 0, 1, 1], ncol=1, loc='lower left', fontsize=12)
    ax.set_title('{}'.format(pathway_name), fontsize=24)
    fig.set_size_inches(18, 14)
    plt.savefig(save_dir, bbox_inches='tight', format='pdf')

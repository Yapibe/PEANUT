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


def sort_y_ticks_by_dir_and_pvalue_of_one_column(enrichment_table, y_ticks, annotation_map, direction,
                                                 dir_column_to_sort_by=0):
    """
    Sorts the rows of the enrichment table by direction and p-value.
    Args:
        enrichment_table (numpy.ndarray): The enrichment table to sort.
        y_ticks (list): The y-axis labels corresponding to rows in the enrichment table.
        annotation_map (numpy.ndarray): Annotations for each cell in the enrichment table.
        direction (numpy.ndarray): Array indicating the direction of enrichment for each row.
        dir_column_to_sort_by (int): Index of the column in the direction array to sort by.
    Returns:
        tuple: Tuple containing sorted enrichment table, y-ticks, annotation map, and direction array.
    """
    # 1. sort 0s from 1s in direction table
    sorted_indexes1 = np.argsort(direction[:, dir_column_to_sort_by])
    enrichment_table = enrichment_table[sorted_indexes1, :]
    y_ticks = list(np.array(y_ticks)[sorted_indexes1])
    annotation_map = annotation_map[sorted_indexes1, :]
    direction = direction[sorted_indexes1, :]

    # 2. sort p_value within indexes with sorted direction (negative first, positives after)
    negatives = len(direction[direction[:, dir_column_to_sort_by] == 0, dir_column_to_sort_by])

    negative_enrichment_table = enrichment_table[:negatives, ]
    positive_enrichment_table = enrichment_table[negatives:, ]
    negative_annotation_map = annotation_map[:negatives, ]
    positive_annotation_map = annotation_map[negatives:, ]
    negative_direction = direction[:negatives, ]
    positive_direction = direction[negatives:, ]
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
    y_ticks_ns = list(np.array(y_ticks_n)[sorted_negative_indexes])
    y_ticks_ps = list(np.array(y_ticks_p)[sorted_positive_indexes])

    # 3. return joined objects, with negatives first.
    return np.concatenate((sorted_negative_enrichment_table, sorted_positive_enrichment_table),
                          axis=0), y_ticks_ns + y_ticks_ps, np.concatenate(
        (sorted_negative_annotation_map, sorted_positive_annotation_map), axis=0), np.concatenate(
        (sorted_negative_direction, sorted_positive_direction), axis=0)


def sort_y_ticks_by_dir(enrichment_table, y_ticks, annotation_map, direction):
    """
    Sorts the rows of the enrichment table by the overall direction of enrichment.
    Args:
        enrichment_table (numpy.ndarray): The enrichment table to sort.
        y_ticks (list): The y-axis labels corresponding to rows in the enrichment table.
        annotation_map (numpy.ndarray): Annotations for each cell in the enrichment table.
        direction (numpy.ndarray): Array indicating the direction of enrichment for each row.
    Returns:
        tuple: Tuple containing sorted enrichment table, y-ticks, annotation map, and direction array.
    """
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
    """
    Sorts the rows of the enrichment table based on a custom sorting of y-ticks.
    Args:
        enrichment_table (numpy.ndarray): The enrichment table to sort.
        y_ticks (list): The y-axis labels corresponding to rows in the enrichment table.
        annotation_map (numpy.ndarray): Annotations for each cell in the enrichment table.
    Returns:
        tuple: Tuple containing sorted enrichment table, y-ticks, and annotation map.
    """
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


def plot_enrichment_table(enrichment_table, direction, interesting_pathways, save_dir=None, experiment_names=None,
                          title=None, res_type=None, adj_p_value_threshold=0.01):
    """
    Plots a heatmap of enrichment analysis results, highlighting significant pathways.
    Args:
        enrichment_table (numpy.ndarray): 2D array of enrichment scores for each pathway and experiment.
        direction (numpy.ndarray): 2D array indicating the direction of enrichment.
        interesting_pathways (list): List of pathway names corresponding to rows in the enrichment table.
        save_dir (str, optional): Directory to save the generated plot.
        experiment_names (list, optional): Names of experiments corresponding to columns in the enrichment table.
        title (str, optional): Title for the plot.
        res_type (str, optional): Type of results being plotted, e.g., 'z_score'.
        adj_p_value_threshold (float, optional): Threshold for adjusted p-value to highlight significant pathways.
    Returns:
        None: The function saves the plot to the specified directory.
    """

    # Initialize the plot
    fig, ax = plt.subplots()

    # Find non-zero columns and rows in the enrichment table
    enriched_clusters = np.nonzero(np.sum(enrichment_table, axis=0) != 0)[0]
    found_pathways = np.nonzero(np.sum(enrichment_table, axis=1) != 0)[0]

    # Filter the enrichment table to include only non-zero rows and columns
    enrichment_table = enrichment_table[:, enriched_clusters][found_pathways, :]
    interesting_pathways_filtered = {i: path for i, path in enumerate(interesting_pathways) if i in found_pathways}

    # Prepare the annotations for the heatmap
    annotation_map = (np.round(enrichment_table, 3)).astype(str)
    annotation_map[annotation_map == '0.0'] = ''
    y_ticks = [path[:60] for path in interesting_pathways_filtered.values()]

    # Sort pathways for visualization
    enrichment_table, y_ticks, annotation_map, direction = sort_y_ticks_by_dir_and_pvalue_of_one_column(
        enrichment_table, y_ticks, annotation_map, direction)

    # Find indexes of important pathways (below the adj. p-value threshold)
    important_indexes = np.where(enrichment_table)

    # Color negative enrichment values if direction is given and res_type is not 'z_score'
    if direction is not None and res_type != 'z_score':
        enrichment_table[np.logical_not(direction)] = -enrichment_table[np.logical_not(direction)]
        # Create a color axis for the heatmap
    colorbar_edge = np.max(np.abs(enrichment_table))
    cax = inset_axes(ax, width="100%", height="70%", loc='lower left',
                     bbox_to_anchor=(1.1, 0.2, 0.2, 1), bbox_transform=ax.transAxes, borderpad=0)

    # Plot the heatmap
    heatmap = sns.heatmap(enrichment_table, fmt=".4s", yticklabels=y_ticks, cbar_ax=cax,
                          annot=annotation_map, cmap="coolwarm", linewidths=.1, linecolor='gray',
                          cbar_kws={'label': res_type}, ax=ax, vmin=-colorbar_edge, vmax=colorbar_edge)
    ax.set_xticklabels(experiment_names, fontsize=24)

    # Highlight significant pathways with a rectangle
    for i in range(len(important_indexes[0])):
        heatmap.add_patch(Rectangle((important_indexes[1][i], important_indexes[0][i]), 1, 1,
                                    fill=False, edgecolor='black', lw=3))
    # set font size of colobar ticks
    heatmap.figure.axes[-1].yaxis.label.set_size(20)

    # Set colorbar properties
    cbar = heatmap.collections[0].colorbar
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(np.round(np.abs(cbar.get_ticks()), 1))

    # Adjust x-tick labels
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')

    # Add legend for significance and direction
    legend_handles = [Patch(fill=False, edgecolor='black', label=f'adj_p_value<{adj_p_value_threshold}')]
    if direction is not None:
        legend_handles.append(Patch(fill=True, color='red', label='Upwards'))
        legend_handles.append(Patch(fill=True, color='blue', label='Downwards'))
    ax.legend(handles=legend_handles, bbox_to_anchor=[1, 0, 1, 1], ncol=1, loc='lower left', fontsize=14)

    # Set plot title
    ax.set_title(title, fontsize=28)

    # Adjust figure size and save the plot
    fig.set_size_inches(20, 40)
    plt.savefig(save_dir, bbox_inches='tight')
    print(f"plot saved to {save_dir}")

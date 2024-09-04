import numpy as np
import pandas as pd
from .settings import Settings
from os import path, makedirs
import matplotlib.pyplot as plt

def print_aggregated_pathway_information(args: Settings, all_pathways: dict) -> str:
    """
    Print aggregated pathway information including P-values, trends, and significant genes
    for each pathway to a text file based on a given experiment.

    Parameters:
    - args (GeneralArgs): General arguments and settings.
    - all_pathways (dict): Dictionary containing pathway information across different conditions.

    Returns:
    - file_path (str): Path to the output file.
    """
    # Define the path for the output file
    file_path = path.join(args.output_path, 'Text', f'{args.experiment_name}_{args.pathway_file}'
                                                   f'_{args.alpha}_aggregated.txt')

    if not all_pathways:
        with open(file_path, 'w') as file:
            file.write("No significant pathways found.\n")
        print("No significant pathways found. File written with message.")
        return

    # Create a list of (pathway, best_p_value) tuples
    pathways_p_values = []
    for pathway, conditions in all_pathways.items():
        # Find the minimum P-value for each pathway across all conditions
        best_p_value = min(condition_data['P-value'] for condition_data in conditions.values())
        pathways_p_values.append((pathway, best_p_value))

    # Sort pathways by the best (lowest) P-value
    pathways_sorted = sorted(pathways_p_values, key=lambda x: x[1])

    # Write to the output file
    with open(file_path, 'w') as file:
        for pathway, best_p_value in pathways_sorted:
            file.write(f"Pathway: {pathway} P-value: {best_p_value:.5f}\n")  # Print pathway and its best P-value

            # Aggregate and write trends for all conditions for the pathway
            trends = [f"{condition_name}: {all_pathways[pathway][condition_name]['Trend']}"
                      for condition_name in all_pathways[pathway]]
            file.write(f"Trends: {', '.join(trends)}\n")

            # Aggregate and write significant genes across all conditions
            file.write("Significant Genes:\n")
            gene_scores_across_conditions = {}
            for condition_name, condition_data in all_pathways[pathway].items():
                for gene_id, gene_info in condition_data.get('significant_genes', {}).items():
                    if gene_id not in gene_scores_across_conditions:
                        gene_scores_across_conditions[gene_id] = {'Symbol': gene_info['Symbol'], 'Scores': []}
                    gene_scores_across_conditions[gene_id]['Scores'].append(gene_info['Score'])

            # List each gene with its scores across conditions
            for gene_id, gene_data in gene_scores_across_conditions.items():
                scores_str = ', '.join(map(str, gene_data['Scores']))
                file.write(f"    {gene_data['Symbol']}: {scores_str}\n")

            file.write("\n")
    print(f"Aggregated pathway information written to {file_path}")
    return file_path

def print_enriched_pathways_to_file(task: EnrichTask, FDR_threshold: float) -> None:
    """
    Write enriched pathways that pass the FDR threshold to a text file.

    Parameters:
    - task (EnrichTask): Enrichment task containing task-specific settings.
    - FDR_threshold (float): False Discovery Rate threshold for significance.

    Returns:
    - None
    """
    output_file_path = path.join(task.temp_output_folder, f'{task.name}.txt')
    significant_count = 0

    with open(output_file_path, 'w') as file:
        for pathway, details in task.filtered_pathways.items():
            p_value = details.get('Adjusted_p_value')
            if p_value is not None and p_value < FDR_threshold:
                file.write(f"{pathway} {p_value:.5f}\n")
                significant_count += 1

    print(f"Total significant pathways written: {significant_count}")

def plot_pathways_mean_scores(args: Settings, all_pathways: dict) -> str:
    """
    Plot mean scores of pathways across all conditions and save the plot as a PNG file.

    Parameters:
    - general_args (GeneralArgs): General arguments and settings.
    - all_pathways (dict): Dictionary containing pathway information across different conditions.

    Returns:
    - output_file_path (str): Path to the output plot file.
    """
    if not all_pathways:
        print("No pathways to plot. Exiting function.")
        return

    # Initialize dictionaries to store mean scores and p-values for each condition
    mean_scores_data = {}
    p_values_data = {}
    for pathway, conditions in all_pathways.items():
        for condition_name, condition_data in conditions.items():
            mean_scores_data.setdefault(condition_name, {})[pathway] = condition_data.get('Mean', 0)
            p_values_data.setdefault(condition_name, {})[pathway] = condition_data.get('P-value', 1)

    # Create DataFrames from the dictionaries
    data_df = pd.DataFrame(mean_scores_data)
    p_values_df = pd.DataFrame(p_values_data)

    if data_df.empty:
        print("Data for plotting is empty. Exiting function.")
        return

    # Sort pathways alphabetically
    sorted_pathways = sorted(data_df.index)

    # Reorder the DataFrames based on the sorted pathways
    data_df = data_df.loc[sorted_pathways]
    p_values_df = p_values_df.loc[sorted_pathways]

    # Create a large figure to accommodate the potentially large number of pathways
    plt.figure(figsize=(20, 60))  # This may need adjustment based on the actual data
    ax = plt.subplot(111)

    # Prepare data for plotting
    conditions = list(mean_scores_data.keys())
    total_pathways = data_df.index
    num_conditions = len(conditions)
    bar_height = 0.8 / num_conditions  # Calculate bar height based on the number of conditions
    positions = np.arange(len(total_pathways))

    # Generate a color map for the conditions
    colors = plt.colormaps['viridis'](np.linspace(0, 1, num_conditions))

    # Define keywords for bold formatting
    keywords = ['NEURO', 'SYNAP']

    # Plot each condition's mean scores for each pathway
    for i, condition in enumerate(conditions):
        mean_scores = data_df[condition].values
        p_values = p_values_df[condition].values

        # Plot bars with different styles based on p-value significance
        for j, (score, p_value) in enumerate(zip(mean_scores, p_values)):
            bar_style = {"color": "white", "edgecolor": colors[i], "hatch": "//"} if p_value > args.FDR_threshold else {
                "color": colors[i]}
            ax.barh(positions[j] + bar_height * i, score, height=bar_height, **bar_style)

    # Set y-axis labels to be pathway names, replace underscores with spaces for readability
    ax.set_yticks(positions + bar_height * (num_conditions / 2) - bar_height / 2)
    formatted_pathways = [pathway.replace('_', ' ') for pathway in total_pathways]
    ax.set_yticklabels(formatted_pathways, fontsize=12)

    # Bold labels for specific keywords
    for i, label in enumerate(ax.get_yticklabels()):
        if any(keyword in label.get_text().upper() for keyword in keywords):
            label.set_fontweight('bold')

    # Label axes and set title
    ax.set_ylabel('Pathways', fontsize=16)
    ax.set_xlabel('Mean Scores', fontsize=16)
    ax.set_title('Pathway Mean Scores Across Different Conditions', fontsize=20)

    # Create a legend for the conditions
    plt.legend([plt.Rectangle((0, 0), 1, 1, color=colors[i]) for i in range(num_conditions)], conditions,
               prop={'size': 14})

    # Adjust subplot layout to avoid clipping of tick-labels
    plt.subplots_adjust(left=0.4)

    # Save the figure to a PDF file in the specified output directory
    output_file_path = path.join(args.output_path, 'Plots', f"{args.experiment_name}_{args.pathway_file}"
                                                           f"_{args.alpha}_plot.pdf")
    makedirs(path.dirname(output_file_path), exist_ok=True)
    plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
    plt.close()
    return output_file_path


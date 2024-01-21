import pandas as pd
import matplotlib.pyplot as plt
import os


def read_scores(file_name):
    """Reads scores from a given file into a dictionary."""
    return pd.read_csv(file_name, sep=' ', header=None, names=['Pathway', 'Score'], index_col='Pathway')[
        'Score'].to_dict()


def read_pathways(file_name):
    """Reads pathways from a given file into a dictionary."""
    with open(file_name, 'r') as file:
        return {line.split()[0]: [int(gene) for gene in line.split()[2:]] for line in file}


def process_experiment(condition_file, experiment_file, pathways_file):
    """Processes an experiment and returns scores, enriched pathway genes, and mean scores."""
    scores_dict = read_scores(condition_file)
    experiment_data = pd.read_excel(experiment_file)
    experiment_data_filtered = experiment_data[experiment_data['Score'] != 0]

    pathways = read_pathways(pathways_file)
    enriched_pathway_genes = {}
    pathway_mean_scores = {}

    for pathway, genes in pathways.items():
        filtered_genes = experiment_data_filtered[experiment_data_filtered['GeneID'].isin(genes)]
        if not filtered_genes.empty:
            enriched_pathway_genes[pathway] = filtered_genes.set_index('GeneID')[['Human_Name', 'Score']].to_dict(
                orient='index')
            pathway_mean_scores[pathway] = filtered_genes['Score'].mean()

    return scores_dict, enriched_pathway_genes, pathway_mean_scores


def calculate_trend(pathway_mean_scores):
    """Calculates the trend (up or down) of each pathway based on its mean score."""
    trends = {}
    for pathway, mean_score in pathway_mean_scores.items():
        trend = "Up" if mean_score > 0 else "Down"
        trends[pathway] = trend
    return trends


def bold_keywords(text, keywords):
    """Returns text with keywords in bold."""
    for keyword in keywords:
        if keyword in text:
            text = text.replace(text, f"**{text}**")
    return text


def print_pathway_information(condition, scores_dict, pathway_genes_dict, pathway_trends, output_dir, experiment_name):
    """

    :param condition:
    :param scores_dict:
    :param pathway_genes_dict:
    :param pathway_trends:
    :param output_dir:
    :param experiment_name:
    :return:
    """
    file_path = os.path.join(output_dir, 'Text', f'{experiment_name}.txt')
    with open(file_path, 'a') as file:
        file.write(f"Condition: {condition}\n------------------------\n")
        for pathway, p_value in scores_dict.items():
            trend = pathway_trends.get(pathway, "N/A")
            file.write(f"{pathway} (Trend: {trend}), P-Value: {p_value}\n")

            for gene_id, details in pathway_genes_dict.get(pathway, {}).items():
                if details['Score'] > 1 or details['Score'] < -1:
                    file.write(f"{details['Human_Name']} (ID: {gene_id}), Score: {details['Score']}\n")

            file.write("\n")


def plot_pathways_mean_scores(pathway_mean_scores_data, scores_dict, output_dir, experiment_name):
    """Plots mean scores of pathways across experiments in a grouped bar chart."""
    data_df = pd.DataFrame(pathway_mean_scores_data)

    # Increase the figure size significantly
    plt.figure(figsize=(60, 20))  # Adjust the size as needed
    ax = plt.subplot(111)

    conditions = list(pathway_mean_scores_data.keys())
    num_conditions = len(conditions)

    total_pathways = set()
    for condition in conditions:
        total_pathways.update(scores_dict[condition].keys())
    total_pathways = list(total_pathways)
    num_pathways = len(total_pathways)

    bar_width = 0.8 / num_conditions
    positions = list(range(num_pathways))

    for i, condition in enumerate(conditions):
        mean_scores = data_df[condition].reindex(total_pathways)
        ax.bar([p + bar_width * i for p in positions], mean_scores, width=bar_width, label=condition)

    ax.set_xticks([p + bar_width * (num_conditions / 2) - bar_width / 2 for p in positions])

    # Removing underscores and applying bold formatting
    keywords = ['NEURO', 'SYNAP']
    formatted_pathways = [pathway.replace('_', ' ') for pathway in total_pathways]
    ax.set_xticklabels(formatted_pathways, rotation=90, fontsize=14)  # Increase fontsize for readability

    # Apply bold formatting based on keywords
    for label, pathway in zip(ax.get_xticklabels(), total_pathways):
        if any(keyword in pathway.upper() for keyword in keywords):
            label.set_weight('bold')

    ax.set_xlabel('Pathways', fontsize=16)
    ax.set_ylabel('Mean Scores', fontsize=16)
    ax.set_title('Pathway Mean Scores Across Different Conditions', fontsize=20)
    ax.legend(prop={'size': 14})

    plt.subplots_adjust(bottom=0.4)  # Adjust for layout

    # Define the output file path
    output_file_path = os.path.join(output_dir,"Plots", f'{experiment_name}.pdf')

    # Save the plot to the specified directory
    plt.savefig(output_file_path, format='pdf', bbox_inches='tight')
    plt.show()

experiment_name = 'Parkinson'
condition_files = ['scores_of_T_v_N.txt', 'scores_of_500nm.txt']
test_file_paths = ['Inputs/experiments_data/Parkinson/roded_T_v_N.xlsx', 'Inputs/experiments_data/Parkinson/roded_500nm.xlsx']
pathway_file_dir = 'Data/H_sapiens/pathways/pathway_file'
output_plot_dir = 'Outputs'
all_enriched_genes = {}
all_mean_scores = {}

for condition_file, experiment_file in zip(condition_files, test_file_paths):
    scores_dict, pathway_genes_dict, pathway_mean_scores = process_experiment(condition_file, experiment_file, pathway_file_dir)
    pathway_trends = calculate_trend(pathway_mean_scores)

    print_pathway_information(condition_file, scores_dict, pathway_genes_dict, pathway_trends, output_plot_dir, experiment_name)

    all_enriched_genes[condition_file] = scores_dict
    all_mean_scores[condition_file] = pathway_mean_scores

plot_pathways_mean_scores(all_mean_scores, all_enriched_genes, output_plot_dir, experiment_name)




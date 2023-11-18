import numpy as np
import random
from os import path
from collections import defaultdict
from scipy.stats import rankdata, wilcoxon
from statsmodels.stats.multitest import multipletests
from args import PropagationTask
import utils
from pathway_enrichment import load_pathways_genes
import pandas as pd
import time


def shuffle_scores(dataframe, shuffle_column, num_shuffles=1000):
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


def filter_significant_pathways(pathways, p_values, threshold=0.05):
    return {pathway: pval for pathway, pval in zip(pathways, p_values) if pval < threshold}


def save_results(filename, data):
    with open(filename, 'w') as file:
        for pathway, pval in data.items():
            file.write(f'{pathway}: {pval}\n')


# Step 1: Data Preparation
def load_and_prepare_data(task):
    data = utils.read_prior_set(task.experiment_file_path)
    data['Score'] = data['Score'].apply(lambda x: abs(x))
    data = data[data['GeneID'].apply(lambda x: str(x).isdigit())]
    data['GeneID'] = data['GeneID'].astype(int)
    # remove p_values column
    data = data.drop(columns=['P-value'])
    data = data.reset_index(drop=True)

    return data


# Extracting experiment_gene_set
def get_experiment_gene_set(prior_data):
    return set(prior_data['GeneID'])


def get_fold_change_ranks(prior_data):
    fold_values = {gene: score for gene, score in zip(prior_data['GeneID'], prior_data['Score'])}
    return {gene: rank for gene, rank in zip(fold_values.keys(), rankdata(list(fold_values.values())))}


def filter_pathways_genes(genes_by_pathway, gene_set):
    # Filter out genes not in the experiment and pathways not within the size limits
    filtered_pathways = {}
    for pathway, genes in genes_by_pathway.items():
        filtered_genes = [gene for gene in genes if gene in gene_set]
        if 10 <= len(filtered_genes) <= 60:
            filtered_pathways[pathway] = filtered_genes
    return filtered_pathways


# Step 2: Randomization of Scores
def shuffle_gene_scores(data):
    shuffled_data = data.copy()
    shuffled_scores = shuffled_data['Score'].sample(frac=1).reset_index(drop=True)
    shuffled_data['Score'] = shuffled_scores
    return shuffled_data


# Step 4: Empirical P-Value Calculation
def compare_pathway_to_randomized(real_scores, shuffled_scores):
    # Ensure that real_scores and shuffled_scores are NumPy arrays of type float
    real_scores = np.asarray(real_scores, dtype=float)
    shuffled_scores = np.asarray(shuffled_scores, dtype=float)

    # Wilcoxon signed-rank test logic
    differences = real_scores - shuffled_scores
    abs_diff = np.abs(differences)
    ranks = rankdata(abs_diff)
    signed_ranks = ranks * np.sign(differences)
    T = np.sum(signed_ranks)
    return T


# Step 5: Processing Results
def calculate_empirical_p_values(real_data, pathways, num_simulations=1000):
    empirical_p_values = {}

    for pathway_name, pathway_genes in pathways.items():
        start = time.time()
        count_negative_T = 0
        # Filter the real scores for pathway genes
        pathway_real_scores = real_data[real_data['GeneID'].isin(pathway_genes)]['Score'].astype(float)

        for i in range(1, num_simulations + 1):
            shuffled_score_column = f'Shuffled_Score_{i}'
            # Filter the shuffled scores for pathway genes
            pathway_shuffled_scores = real_data[real_data['GeneID'].isin(pathway_genes)][shuffled_score_column].astype(
                float)

            T = compare_pathway_to_randomized(pathway_real_scores, pathway_shuffled_scores)
            if T < 0:
                count_negative_T += 1

        empirical_p_values[pathway_name] = (count_negative_T + 1) / (num_simulations + 1)
        end = time.time()
        print(f'{pathway_name}: {empirical_p_values[pathway_name]} ({end - start} seconds)')

    return empirical_p_values


def benjamini_hochberg_correction(p_values):
    n = len(p_values)
    sorted_p_values = np.sort(p_values)
    sorted_indices = np.argsort(p_values)
    ranks = np.arange(1, n + 1)

    # Calculate corrected p-values
    corrected_p_values = np.minimum(1, sorted_p_values * n / ranks)

    # Ensure monotonicity
    for i in range(n - 1, 0, -1):
        corrected_p_values[i - 1] = min(corrected_p_values[i - 1], corrected_p_values[i])

    # Back to original order
    corrected_p_values_back_to_original = np.empty(n)
    corrected_p_values_back_to_original[sorted_indices] = corrected_p_values

    # _, corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')

    return corrected_p_values_back_to_original


# Main analysis workflow
task = PropagationTask(experiment_name='random', create_similarity_matrix=False)
prior_data = load_and_prepare_data(task)

# Shuffle the scores and add to the DataFrame
updated_prior_data = shuffle_scores(prior_data, 'Score')
experiment_gene_set = get_experiment_gene_set(prior_data)

dir_path = path.dirname(path.realpath(__file__))
pathway_file_path = path.join(dir_path, '../Data', 'H_sapiens', 'pathways', 'pathway_file')
genes_by_pathway = load_pathways_genes(pathway_file_path)

# Filter pathways to include only genes in the experiment and with sizes 10-60
filtered_genes_by_pathway = filter_pathways_genes(genes_by_pathway, experiment_gene_set)

# random_scores = generate_random_scores(gene_ranks, filtered_genes_by_pathway)

empirical_p_values = calculate_empirical_p_values(updated_prior_data, filtered_genes_by_pathway)

p_values = np.array(list(empirical_p_values.values()))
corrected_p_values = benjamini_hochberg_correction(p_values)
significant_pathways = filter_significant_pathways(empirical_p_values.keys(), corrected_p_values)
# index correct p-values by pathway name
corrected_p_values = {pathway: pval for pathway, pval in zip(empirical_p_values.keys(), corrected_p_values)}

specific_pathways = [
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

csv_data = []
for pathway in specific_pathways:
    # print the pathway name, p-value, and corrected p-value
    print(f'{pathway}: {empirical_p_values[pathway]} {corrected_p_values[pathway]}')

save_results('../Mann.txt', significant_pathways)

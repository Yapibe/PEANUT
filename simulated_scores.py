import pandas as pd
import numpy as np
from os import path, remove, listdir, rmdir
from scipy.stats import norm
from pipeline.utils import load_pathways_genes
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc

# Function to simulate scores and p-values
def simulate_scores(pathways, delta=1.0, percentage=0.5, num_decoy_pathways=10, seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Create a list of all unique genes in all pathways
    all_genes = list(set(gene for genes in pathways.values() for gene in genes))

    # Create a dictionary mapping gene IDs to tuples of (Score, P-value)
    scores = np.random.normal(0, 1, len(all_genes))

    # Calculate two-tailed p-values based on the scores
    abs_scores = np.abs(scores)
    tail_probabilities = norm.cdf(-abs_scores)
    p_values = tail_probabilities * 2

    gene_score_pvalue_dict = {gene: (score, p_value) for gene, score, p_value in zip(all_genes, scores, p_values)}

    # Filter pathways to include only those with sizes between 20 and 200 genes for decoys
    valid_pathways = {k: v for k, v in pathways.items() if 20 <= len(v) <= 200}

    # Randomly select decoy pathways from the valid pathways
    selected_pathways = np.random.choice(list(valid_pathways.keys()), num_decoy_pathways, replace=False)

    # Create a set of all genes in the selected decoy pathways
    decoy_genes_set = set(gene for pathway in selected_pathways for gene in valid_pathways[pathway])

    # Calculate the number of genes to change based on the percentage
    num_genes_to_change = int(len(decoy_genes_set) * percentage)

    # Randomly select the genes to change
    genes_to_change = np.random.choice(list(decoy_genes_set), num_genes_to_change, replace=False)

    # Assign scores to the selected decoy pathway genes from a different distribution
    decoy_scores = np.random.normal(delta, 1, len(genes_to_change))
    for gene, score in zip(genes_to_change, decoy_scores):
        decoy_p_value = norm.cdf(-np.abs(score)) * 2
        gene_score_pvalue_dict[gene] = (score, decoy_p_value)

    return gene_score_pvalue_dict, selected_pathways

def run_pipeline(alpha):
    from pipeline.main import main
    main(run_propagation=True, alpha=alpha)

# Function to calculate precision, recall, and AUPR
def calculate_metrics(true_decoys, identified_pathways, all_pathways):
    true_positives = len(set(true_decoys).intersection(set(identified_pathways)))
    false_positives = len(set(identified_pathways) - set(true_decoys))
    false_negatives = len(set(true_decoys) - set(identified_pathways))

    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0

    y_true = [1 if pathway in true_decoys else 0 for pathway in all_pathways]
    y_scores = [1 if pathway in identified_pathways else 0 for pathway in all_pathways]

    precision_curve, recall_curve, _ = precision_recall_curve(y_true, y_scores)
    aupr = auc(recall_curve, precision_curve)

    return precision, recall, aupr

# Main script
root_dir = path.dirname(path.abspath(__file__))
pathways_file = path.join(root_dir, 'pipeline', 'Data', 'H_sapiens', 'pathways', 'bio_pathways')
pathways = load_pathways_genes(pathways_file)

deltas = [0.3, 1, 3, 10]
percentages = [0.1, 0.3, 0.75, 1]  # Percentages of genes to change
alphas = [1]  # Different alpha values to loop through
num_decoy_pathways_list = [10, 30, 100]  # Different numbers of decoy pathways to loop through
n_runs = 10  # Number of runs for each combination

results_dict = {
    'delta': [],
    'percentage': [],
    'alpha': [],
    'num_decoy_pathways': [],
    'precision': [],
    'recall': [],
    'aupr': [],
    'aupr_hyper': []
}

for delta in deltas:
    for percentage in percentages:
        for alpha in alphas:
            for num_decoy_pathways in num_decoy_pathways_list:
                precision_values = []
                recall_values = []
                aupr_values = []
                aupr_hyper_values = []

                for _ in range(n_runs):
                    gene_score_pvalue_dict, true_decoys = simulate_scores(pathways, delta=delta, percentage=percentage, num_decoy_pathways=num_decoy_pathways, seed=None)
                    simulated_scores_path = path.join(root_dir, 'pipeline', 'Inputs', 'Simulated', f'simulated_scores_{delta}_{percentage}_{alpha}_{num_decoy_pathways}.xlsx')
                    result_file_path = path.join(root_dir, 'pipeline', 'Outputs', 'Temp', f'simulated_scores_{delta}_{percentage}_{alpha}_{num_decoy_pathways}.txt')
                    hyper_results_path = path.join(root_dir, 'pipeline', 'Outputs', 'hyper_results.txt')

                    # Save the simulated scores
                    results = pd.DataFrame.from_dict(gene_score_pvalue_dict, orient='index', columns=['Score', 'P-value'])
                    results.reset_index(inplace=True)
                    results.rename(columns={'index': 'GeneID'}, inplace=True)
                    results.to_excel(simulated_scores_path, index=False)

                    # Run the pipeline
                    run_pipeline(alpha=alpha)

                    with open(result_file_path, 'r') as f:
                        identified_pathways = [line.split()[0] for line in f.readlines()]

                    # Calculate metrics for the whole run
                    precision, recall, aupr = calculate_metrics(true_decoys, identified_pathways, list(pathways.keys()))
                    precision_values.append(precision)
                    recall_values.append(recall)
                    aupr_values.append(aupr)

                    # Read hypergeometric test results from the file
                    with open(hyper_results_path, 'r') as f:
                        significant_pathways_hyper = [line.strip() for line in f.readlines()]

                    # Calculate AUPR for hypergeometric test
                    _, _, aupr_hyper = calculate_metrics(true_decoys, significant_pathways_hyper, list(pathways.keys()))
                    aupr_hyper_values.append(aupr_hyper)

                    # Delete temporary files
                    if path.exists(simulated_scores_path):
                        remove(simulated_scores_path)
                    if path.exists(hyper_results_path):
                        remove(hyper_results_path)
                    if path.exists(result_file_path):
                        remove(result_file_path)
                    # Erase each folder and its contents in Outputs/Propagation_Scores that start with simulated_scores
                    for folder in [f for f in listdir(path.join(root_dir, 'pipeline', 'Outputs', 'Propagation_Scores')) if f.startswith('simulated_scores')]:
                        folder_path = path.join(root_dir, 'pipeline', 'Outputs', 'Propagation_Scores', folder)
                        for file in listdir(folder_path):
                            remove(path.join(folder_path, file))
                        rmdir(folder_path)

                # Average the precision, recall, and AUPR values
                results_dict['delta'].append(delta)
                results_dict['percentage'].append(percentage)
                results_dict['alpha'].append(alpha)
                results_dict['num_decoy_pathways'].append(num_decoy_pathways)
                results_dict['precision'].append(np.mean(precision_values))
                results_dict['recall'].append(np.mean(recall_values))
                results_dict['aupr'].append(np.mean(aupr_values))
                results_dict['aupr_hyper'].append(np.mean(aupr_hyper_values))

# Convert results to DataFrame
results_df = pd.DataFrame(results_dict)

# Plot the results for each combination separately
for percentage in percentages:
    for alpha in alphas:
        for num_decoy_pathways in num_decoy_pathways_list:
            subset = results_df[(results_df['percentage'] == percentage) & (results_df['alpha'] == alpha) & (results_df['num_decoy_pathways'] == num_decoy_pathways)]
            plt.figure(figsize=(15, 10))
            plt.plot(subset['delta'].values, subset['precision'].values, marker='o', label='Precision')
            plt.plot(subset['delta'].values, subset['recall'].values, marker='o', label='Recall')
            plt.plot(subset['delta'].values, subset['aupr'].values, marker='o', label='AUPR')
            plt.plot(subset['delta'].values, subset['aupr_hyper'].values, marker='o', label='AUPR Hyper')
            plt.xlabel('Delta')
            plt.ylabel('Value')
            plt.title(f'Precision, Recall, and AUPR (Percentage: {percentage*100}%, Alpha: {alpha}, Decoy Pathways: {num_decoy_pathways})')
            plt.xticks(deltas)
            plt.legend()
            plt.grid(True)
            plt.savefig(path.join(root_dir, f'precision_recall_aupr_chart_{percentage}_{alpha}_{num_decoy_pathways}.png'))
            plt.show()

print(results_df)

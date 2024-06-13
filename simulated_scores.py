import pandas as pd
import numpy as np
from os import path, remove
from scipy.stats import norm
from pipeline.utils import load_pathways_genes
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

# Function to simulate scores and p-values
def simulate_scores(pathways, delta=1.0, num_decoy_pathways=10, seed=None):
    if seed is not None:
        np.random.seed(seed)

    # Create a list of all unique genes in all pathways
    all_genes = list(set(gene for genes in pathways.values() for gene in genes))

    # Create a dictionary mapping gene IDs to tuples of (Score, P-value)
    scores = np.random.normal(0, 1, len(all_genes))

    # Calculate two-tailed p-values based on the scores
    abs_scores = np.abs(scores)  # Step 1: Take the absolute value of the scores
    tail_probabilities = norm.cdf(-abs_scores)  # Step 2: Calculate the tail probability using the CDF
    p_values = tail_probabilities * 2  # Step 3: Multiply by 2 to get the two-tailed p-value

    gene_score_pvalue_dict = {gene: (score, p_value) for gene, score, p_value in zip(all_genes, scores, p_values)}

    # Filter pathways to include only those with sizes between 20 and 200 genes for decoys
    valid_pathways = {k: v for k, v in pathways.items() if 20 <= len(v) <= 200}

    # Randomly select decoy pathways from the valid pathways
    selected_pathways = np.random.choice(list(valid_pathways.keys()), num_decoy_pathways, replace=False)

    # Create a set of all genes in the selected decoy pathways
    decoy_genes_set = set(gene for pathway in selected_pathways for gene in valid_pathways[pathway])

    # Assign scores to decoy pathway genes from a different distribution
    decoy_scores = np.random.normal(delta, 1, len(decoy_genes_set))
    for gene, score in zip(decoy_genes_set, decoy_scores):
        decoy_p_value = norm.cdf(-np.abs(score)) * 2  # Calculate p-value for decoy gene scores
        gene_score_pvalue_dict[gene] = (score, decoy_p_value)

    return gene_score_pvalue_dict, selected_pathways

def run_pipeline():
    from pipeline.main import main
    main(run_propagation=True)

# Function to calculate precision, recall, and AUC
def calculate_metrics(true_decoys, identified_pathways, all_pathways):
    true_positives = len(set(true_decoys).intersection(set(identified_pathways)))
    false_positives = len(set(identified_pathways) - set(true_decoys))
    false_negatives = len(set(true_decoys) - set(identified_pathways))

    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0

    y_true = [1 if pathway in true_decoys else 0 for pathway in all_pathways]
    y_scores = [1 if pathway in identified_pathways else 0 for pathway in all_pathways]

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    return precision, recall, roc_auc

# Main script
root_dir = path.dirname(path.abspath(__file__))
pathways_file = path.join(root_dir, 'Pipeline', 'Data', 'H_sapiens', 'pathways', 'bio_pathways')
pathways = load_pathways_genes(pathways_file)

deltas = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
n_runs = 10  # Number of runs for each delta
precision_list = []
recall_list = []
auc_list = []

for delta in deltas:
    precision_values = []
    recall_values = []
    auc_values = []

    for _ in range(n_runs):
        gene_score_pvalue_dict, true_decoys = simulate_scores(pathways, delta=delta, seed=None)
        simulated_scores_path = path.join(root_dir, 'pipeline', 'Inputs', 'Simulated', f'simulated_scores_{delta}.xlsx')
        result_file_path = path.join(root_dir, 'pipeline', 'Outputs', 'Temp', f'simulated_scores_{delta}.txt')

        # Save the simulated scores
        results = pd.DataFrame.from_dict(gene_score_pvalue_dict, orient='index', columns=['Score', 'P-value'])
        results.reset_index(inplace=True)
        results.rename(columns={'index': 'GeneID'}, inplace=True)
        results.to_excel(simulated_scores_path, index=False)

        # Run the pipeline
        run_pipeline()

        with open(result_file_path, 'r') as f:
            identified_pathways = [line.split()[0] for line in f.readlines()]

        precision, recall, roc_auc = calculate_metrics(true_decoys, identified_pathways, list(pathways.keys()))
        precision_values.append(precision)
        recall_values.append(recall)
        auc_values.append(roc_auc)

        # Delete simulated scores files
        if path.exists(simulated_scores_path):
            remove(simulated_scores_path)

    # Average the precision, recall, and AUC values
    precision_list.append(np.mean(precision_values))
    recall_list.append(np.mean(recall_values))
    auc_list.append(np.mean(auc_values))

plt.figure(figsize=(10, 6))
plt.plot(deltas, precision_list, marker='o', label='Precision')
plt.plot(deltas, recall_list, marker='o', label='Recall')
plt.plot(deltas, auc_list, marker='o', label='AUC')
plt.xlabel('Delta')
plt.ylabel('Value')
plt.title('Precision, Recall, and AUC for Different Delta Values, Alpha 1 SGEA')
plt.xticks(deltas)
plt.legend()
plt.grid(True)
plt.savefig(path.join(root_dir, 'precision_recall_auc_chart_alpha1_SGEA.png'))
plt.show()

print("Precision values:", precision_list)
print("Recall values:", recall_list)
print("AUC values:", auc_list)

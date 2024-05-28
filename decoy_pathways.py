import random
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

def generate_and_save_combined_pathways(real_pathways, all_genes, output_file):
    """
    Generate decoy pathways by randomly selecting nodes while preserving the topological structure,
    and save both real and decoy pathways to a file.

    Parameters:
    - real_pathways (dict): Dictionary containing real pathways with their gene lists.
    - all_genes (list): List of all possible genes.
    - output_file (str): Path to the file where the combined pathways will be saved.

    Returns:
    - None
    """
    # Convert all genes to integers
    all_genes = list(map(int, all_genes))

    # Generate decoy pathways
    decoy_pathways = {}
    for pathway, genes in real_pathways.items():
        decoy_genes = random.sample(all_genes, len(genes))
        decoy_pathways['Decoy_' + pathway] = decoy_genes

    # Combine real and decoy pathways
    combined_pathways = {**real_pathways, **decoy_pathways}

    # Prepare data for saving
    with open(output_file, 'w') as file:
        file.write("Pathway\tLength\tGenes\n")
        for pathway, genes in combined_pathways.items():
            genes = list(map(int, genes))  # Ensure genes are integers
            length = len(genes)
            genes_str = ' '.join(map(str, genes))
            file.write(f"{pathway}\t{length}\t{genes_str}\n")

    print(f"Combined pathways saved to {output_file}")


def evaluate_performance(real_results, decoy_results):
    """
    Evaluate the performance of the pathway analysis using ROC curves.

    Parameters:
    - real_results (dict): Dictionary containing the results for real pathways.
    - decoy_results (dict): Dictionary containing the results for decoy pathways.

    Returns:
    - float: The area under the ROC curve (AUC).
    """
    y_true = [1] * len(real_results) + [0] * len(decoy_results)
    y_scores = [result['Adjusted_p_value'] for result in real_results.values()] + \
               [result['Adjusted_p_value'] for result in decoy_results.values()]

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:0.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC)')
    plt.legend(loc="lower right")
    plt.show()

    return roc_auc
import random
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from os import path


def generate_decoy_pathways(pathway_file, output_file):
    """
    Generate decoy pathways by randomly shuffling genes in each pathway
    from the given file and save them to a separate file.

    Parameters:
      - pathway_file (str): Path to the two-tab text file containing pathways.
      - output_file (str): Path to the file where decoy pathways will be saved.

    Returns:
      - None
    """

    # Read pathways and genes
    pathways = {}
    genes_pool = set()  # Use set for faster membership checks
    with open(pathway_file, "r") as file:
        for line in file:
            pathway, *genes = line.strip().split("\t")
            genes = [int(gene) for gene in genes]
            pathways[pathway] = genes
            genes_pool.update(genes)  # Add genes to the pool

    # Generate decoy pathways
    decoy_pathways = {}
    for pathway, genes in pathways.items():
        decoy_genes = random.sample(
            genes_pool, len(genes)
        )  # Shuffle from the gene pool
        decoy_pathways[f"Decoy_{pathway}"] = decoy_genes

    # Save decoy pathways
    with open(output_file, "w") as file:
        for pathway, genes in decoy_pathways.items():
            genes_str = " ".join(map(str, genes))
            file.write(f"{pathway}\t{genes_str}\n")

    print(f"Decoy pathways saved to {output_file}")


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
    y_scores = [result["Adjusted_p_value"] for result in real_results.values()] + [
        result["Adjusted_p_value"] for result in decoy_results.values()
    ]

    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(
        fpr, tpr, color="darkorange", lw=2, label=f"ROC curve (area = {roc_auc:0.2f})"
    )
    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Receiver Operating Characteristic (ROC)")
    plt.legend(loc="lower right")
    plt.show()

    return roc_auc


# Usage
root_folder = path.dirname(path.abspath(__file__))
pathway_file = path.join(root_folder, "Data", "Anat", "pathways", "kegg")
output_file = path.join(root_folder, "Data", "Anat", "pathways", "kegg_decoy")
generate_decoy_pathways(pathway_file, output_file)

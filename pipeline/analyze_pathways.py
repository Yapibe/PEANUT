import os
import time
import logging
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml
import matplotlib.pyplot as plt
import seaborn as sns

# Project and Directory Configuration
def get_project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PROJECT_ROOT = get_project_root()

config_path = os.path.join(PROJECT_ROOT, 'pipeline', 'config.yaml')
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

logging.basicConfig(level=config['logging']['level'], format=config['logging']['format'])
logger = logging.getLogger(__name__)

DIRECTORIES = {key: os.path.join(PROJECT_ROOT, value) for key, value in config['directories'].items()}
os.makedirs(DIRECTORIES['summary_base'], exist_ok=True)
NETWORKS = config['analysis_settings']['networks']
PATHWAY_FILES = config['analysis_settings']['pathway_files']
METHODS = ["PEANUT", "GSEA"]
MAX_WORKERS = config['analysis_settings']['max_workers']
PRE_CALCULATED_DATA = config['pre_calculated_data']

def get_file_list(input_dir):
    return [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]

from scipy.stats import wilcoxon

def perform_wilcoxon_test(analysis_dfs):
    """
    Perform a Wilcoxon signed-rank test to compare PEANUT and GSEA ranks.

    Args:
        analysis_dfs (dict): A dictionary with method names as keys and DataFrames with 'Highest Related Rank' as values.

    Returns:
        None: Prints the test results to the log.
    """
    if "PEANUT" not in analysis_dfs or "GSEA" not in analysis_dfs:
        logger.error("Both PEANUT and GSEA results are required for the Wilcoxon test.")
        return

    # Prepare paired ranks
    peanut_ranks = analysis_dfs["PEANUT"].sort_values("Dataset")["Highest Related Rank"].values
    gsea_ranks = analysis_dfs["GSEA"].sort_values("Dataset")["Highest Related Rank"].values

    # Ensure the datasets match
    peanut_datasets = analysis_dfs["PEANUT"].sort_values("Dataset")["Dataset"].values
    gsea_datasets = analysis_dfs["GSEA"].sort_values("Dataset")["Dataset"].values

    if not all(peanut_datasets == gsea_datasets):
        logger.error("Datasets do not match between PEANUT and GSEA. Ensure proper alignment.")
        return

    # Perform Wilcoxon signed-rank test
    stat, p_value = wilcoxon(peanut_ranks, gsea_ranks, alternative="less")

    # Log the results
    logger.info(f"Wilcoxon test statistic: {stat}")
    logger.info(f"P-value: {p_value}")
    if p_value < 0.05:
        logger.info("PEANUT significantly outperforms GSEA (lower ranks for related pathways).")
    else:
        logger.info("No significant difference between PEANUT and GSEA.")

def calculate_highest_related_rank(file_results, disease_to_pathways, associated_pathway):
    """
    Calculate the rank of the highest related pathway, including the associated pathway.

    Parameters:
    - file_results (pd.DataFrame): DataFrame containing pathway ranking results with columns 'Name' and 'Rank'.
    - disease_to_pathways (dict): Mapping of associated pathways to their related pathways.
    - associated_pathway (str): The pathway associated with the current dataset.

    Returns:
    - int: The rank of the highest related pathway, guaranteed to include the associated pathway.
    """
    # Ensure related_pathways is a set and includes the associated pathway
    related_pathways = set(disease_to_pathways.get(associated_pathway, []))
    related_pathways.add(associated_pathway)

    # Get the ranks of related pathways
    related_pathways_ranks = file_results[file_results['Name'].isin(related_pathways)]['Rank'].tolist()

    # Return the highest related rank (minimum rank)
    highest_related_rank = min(related_pathways_ranks)
    return highest_related_rank


def process_single_file(pathway_file, network_name, method, file_name, disease_to_pathways):
    """
    Process a single pathway file to calculate the rank of the highest related pathway.

    Parameters:
    - pathway_file (str): Pathway file name.
    - network_name (str): Name of the network.
    - method (str): Method used.
    - file_name (str): Name of the file.
    - disease_to_pathways (dict): Mapping of associated pathways to related pathways.

    Returns:
    - dict: Dictionary containing dataset, original pathway, and the highest related rank.
    """
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    output_dir = os.path.join(DIRECTORIES['output_base'], method, network_name, pathway_file, "alpha_0.2", "filtered")
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    try:
        results_df = pd.read_excel(output_file_path)
    except Exception as e:
        logger.error(f"Error reading {output_file_path}: {e}")
        return None

    # Get the rank of the highest related pathway
    highest_related_rank = calculate_highest_related_rank(results_df, disease_to_pathways, pathway_name)

    # Return a simplified dictionary with only the required data
    return {
        'Dataset': dataset_name,
        'Original Pathway': pathway_name,
        'Highest Related Rank': highest_related_rank
    }


def plot_comparative_false_positives(analysis_dfs, output_dir, methods, all_datasets):
    """
    Generate and save a comparative plot for 'Highest Related Rank' across methods.

    Args:
        analysis_dfs (dict): A dictionary where keys are methods and values are DataFrames containing analysis results.
        output_dir (str): Directory to save the plot.
        methods (list): List of methods to include in the plot.
        all_datasets (set): Set of all dataset names expected.
    """
    combined_data = []
    for method, analysis_df in analysis_dfs.items():

        # Prepare data for plotting
        filtered_df = analysis_df[analysis_df['Highest Related Rank'] != 'N/A'].copy()
        filtered_df['Method'] = method
        combined_data.append(filtered_df)

    combined_df = pd.concat(combined_data, ignore_index=True)
    # Create the plot
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        data=combined_df,
        x='Dataset',
        y='Highest Related Rank',
        hue='Method',
        style='Method',
        markers=True,
        s=100  # Adjust marker size as needed
    )
    plt.title('Highest Related Rank Across Methods')
    plt.ylabel('Highest Related Pathway Rank')
    plt.xlabel('Dataset')
    plt.xticks(rotation=90)
    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Save the plot
    plot_path = os.path.join(output_dir, 'highest_related_rank_comparison.png')
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()
    logger.info(f"Comparative Highest Related Rank plot saved to {plot_path}")



def main():
    start_time = time.time()
    disease_to_pathways = PRE_CALCULATED_DATA

    higher_ranked_records = {}
    analysis_dfs = {}
    all_datasets = set()

    FILE_LIST = get_file_list(DIRECTORIES["input"])
    logger.info(f"Found {len(FILE_LIST)} files.")

    for method in METHODS:
        logger.info(f"Processing method: {method}")
        higher_ranked_records[method] = []

        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = []
            for file_name in FILE_LIST:
                dataset_name, pathway_name = file_name.replace(".xlsx", "").split("_", 1)
                all_datasets.add(dataset_name)  # Collect dataset names
                futures.append(
                    executor.submit(
                        process_single_file,
                        "kegg", NETWORKS[0], method, file_name, disease_to_pathways
                    )
                )

            for future in as_completed(futures):
                result = future.result()
                if result:
                    higher_ranked_records[method].append(result)

        analysis_dfs[method] = pd.DataFrame(higher_ranked_records[method])

    # Plot the comparative results
    plot_comparative_false_positives(analysis_dfs, DIRECTORIES["summary_base"], METHODS, all_datasets)

    # Perform the Wilcoxon signed-rank test
    perform_wilcoxon_test(analysis_dfs)

    logger.info(f"Total time taken: {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
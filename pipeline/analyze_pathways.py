import os
import time
import logging
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import yaml
from utils import (
    load_pathways_genes,
    read_network
)
import matplotlib.pyplot as plt
import seaborn as sns
from statistical_methods import jaccard_index
logging.getLogger('matplotlib').setLevel(logging.INFO)

# Project and Directory Configuration
def get_project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PROJECT_ROOT = get_project_root()

config_path = os.path.join(PROJECT_ROOT, 'pipeline', 'config.yaml')
with open(config_path, 'r') as config_file:
    config = yaml.safe_load(config_file)

logging.basicConfig(level=config['logging']['level'], format=config['logging']['format'])
logger = logging.getLogger(__name__)

def get_file_list(input_dir, limit=None):
    files = [f for f in os.listdir(input_dir) if f.endswith('.xlsx')]
    if limit:
        files = files[:limit]
    return files

def load_data(networks, pathway_files, pathways_dir):
    networks_data = {
        name: read_network(os.path.join(PROJECT_ROOT, 'pipeline', 'Data', 'H_sapiens', 'network', name))
        for name in networks
    }
    pathways_data = {
        file: load_pathways_genes(os.path.join(pathways_dir, file))
        for file in pathway_files
    }
    return networks_data, pathways_data

DIRECTORIES = {key: os.path.join(PROJECT_ROOT, value) for key, value in config['directories'].items()}
os.makedirs(DIRECTORIES['summary_base'], exist_ok=True)
NETWORKS = config['analysis_settings']['networks']
PATHWAY_FILES = config['analysis_settings']['pathway_files']
PROP_METHODS = config['analysis_settings']['prop_methods']
ALPHAS = config['analysis_settings']['alphas']
MAX_WORKERS = config['analysis_settings']['max_workers']

# Pre-Calculated Data
PRE_CALCULATED_DATA = config['pre_calculated_data']

def get_pre_calculated_data(disease_name):
    """Retrieve pre-calculated related pathways for a given disease."""
    return PRE_CALCULATED_DATA.get(disease_name, [])


def get_pathway_rank(output_path, pathway_name):
    try:
        results_df = pd.read_excel(output_path)
        if 'FDR q-val' not in results_df.columns:
            logger.error(f"'FDR q-val' column not found in {output_path}.")
            return None, None, None

        results_df = results_df.sort_values(by='FDR q-val', ascending=True).reset_index(drop=True)
        results_df['Rank'] = results_df.index + 1

        pathway_row = results_df[results_df['Term'] == pathway_name]
        if not pathway_row.empty:
            rank = pathway_row['Rank'].values[0]
            fdr_q_val = pathway_row['FDR q-val'].values[0]
            return rank, fdr_q_val, results_df
        else:
            logger.error(f"Pathway '{pathway_name}' not found in {output_path}.")
    except Exception as e:
        logger.error(f"Error reading {output_path}: {e}")
    return None, None, None

def plot_analysis_results(analysis_df, classification_df, output_dir):
    """
    Generate and save analysis plots, including the new 'Overlapping' category.

    Args:
        analysis_df (pd.DataFrame): DataFrame containing analysis results.
        classification_df (pd.DataFrame): DataFrame containing classification counts.
        output_dir (str): Directory to save plots.
    """
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    axes = axes.flatten()

    # Plot 1: True Positives vs. False Positives vs. Overlapping
    tp_fp_counts = classification_df[['Associated', 'Related', 'Overlapping', 'Unrelated']].sum().dropna()
    tp_fp_counts = tp_fp_counts[tp_fp_counts > 0]  # Exclude zero counts
    sns.barplot(x=tp_fp_counts.index, y=tp_fp_counts.values, ax=axes[0], palette='viridis')
    if axes[0].legend_:
        axes[0].legend_.remove()
    axes[0].set_title('Pathway Classifications')
    axes[0].set_ylabel('Count')
    axes[0].set_xlabel('Classification')
    for i, v in enumerate(tp_fp_counts.values):
        axes[0].text(i, v + 0.5, str(v), ha='center')

    # Plot 2: Distribution of Ranks for Each Classification
    sns.histplot(data=analysis_df, x='Rank', hue='Classification', element='step', kde=True, ax=axes[1], palette='viridis', bins=20)
    axes[1].set_title('Distribution of Ranks by Classification')
    axes[1].set_ylabel('Frequency')
    axes[1].set_xlabel('Rank')

    # Plot 3: False Positives Above the Highest Ranked Related Pathway
    filtered_df = analysis_df[analysis_df['Highest Related Rank'] != 'N/A']
    sorted_df = filtered_df.sort_values(by='False Positives Above Top Related')
    sns.lineplot(data=sorted_df, x='Dataset', y='False Positives Above Top Related', marker='o', ax=axes[2])
    axes[2].set_title('False Positives Above Top Related Pathway')
    axes[2].set_ylabel('False Positives Count')
    axes[2].set_xlabel('Dataset')
    axes[2].tick_params(axis='x', rotation=90)

    # Plot 4: Counts of Each Classification per Dataset
    # Compute total counts per dataset
    classification_df['Total'] = classification_df[['Associated', 'Related', 'Overlapping', 'Unrelated']].sum(axis=1)

    # Sort datasets based on total counts
    dataset_order = classification_df.sort_values('Total', ascending=False)['Dataset']

    # Melt the classification_df for plotting
    classification_melted = classification_df.melt(
        id_vars='Dataset',
        value_vars=['Associated', 'Related', 'Overlapping', 'Unrelated'],
        var_name='Classification',
        value_name='Count'
    )

    # Plot using the specified dataset order
    sns.barplot(
        data=classification_melted,
        x='Dataset',
        y='Count',
        hue='Classification',
        ax=axes[3],
        palette='viridis',
        order=dataset_order
    )
    axes[3].set_title('Pathway Classifications per Dataset')
    axes[3].set_ylabel('Count')
    axes[3].set_xlabel('Dataset')
    axes[3].tick_params(axis='x', rotation=90)
    axes[3].legend(title='Classification', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plot_path = os.path.join(output_dir, 'analysis_results.png')
    plt.savefig(plot_path)
    plt.close(fig)
    logger.info(f"Analysis plots saved to {plot_path}")


def calculate_fp_above_related_per_file(file_results, disease_to_pathways, associated_pathway):
    related_pathways = disease_to_pathways.get(associated_pathway, set())
    if not related_pathways:
        related_pathways = {associated_pathway}
    
    file_results = file_results.sort_values(by='Rank').reset_index(drop=True)
    
    # Get the rank of the associated pathway
    associated_pathway_rank = file_results[file_results['Term'] == associated_pathway]['Rank'].min()
    
    # Get the ranks of the related pathways
    related_pathways_ranks = file_results[file_results['Term'].isin(related_pathways)]['Rank']
    
    # Combine the ranks, ensuring associated_pathway_rank is included
    combined_ranks = related_pathways_ranks.tolist()
    if pd.notna(associated_pathway_rank):
        combined_ranks.append(associated_pathway_rank)
    
    if combined_ranks:
        highest_related_rank = min(combined_ranks)
    else:
        # If no ranks are found, assign a high rank (e.g., max rank + 1)
        highest_related_rank = file_results['Rank'].max() + 1
    
    # Count the false positives above the highest related rank
    false_positives_above_related = file_results[
        (file_results['Rank'] < highest_related_rank) &
        (~file_results['Term'].isin(related_pathways))
    ].shape[0]
    
    return false_positives_above_related, highest_related_rank


def process_single_file(pathway_file, network_name, alpha, method, file_name, disease_to_pathways, jac_threshold=0.3):
    """
    Process a single pathway file and classify pathways, including the new 'Overlapping' category
    based on Jaccard similarity.

    Args:
        pathway_file (str): Pathway file name.
        network_name (str): Name of the network.
        alpha (float): Alpha value.
        method (str): Method used.
        file_name (str): Name of the file.
        disease_to_pathways (dict): Mapping of diseases to pathways.
        jac_threshold (float): Jaccard index threshold for overlap.

    Returns:
        dict or None: Processed results or None if errors occur.
    """
    dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
    output_dir = os.path.join(DIRECTORIES['output_base'], method, network_name, pathway_file, f"alpha_{alpha}")
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    rank, fdr_q_val, results_df = get_pathway_rank(output_file_path, pathway_name)
    if results_df is None:
        return None

    # Calculate False Positives Above Related Pathways per File
    fp_above_related, highest_related_rank = calculate_fp_above_related_per_file(
        results_df, disease_to_pathways, pathway_name
    )

    # Collect higher-ranked pathways for analysis
    higher_ranked_paths = results_df[results_df['Rank'] < rank]
    higher_ranked_details = []
    for _, row in higher_ranked_paths.iterrows():
        higher_ranked_details.append({
            'Dataset': dataset_name,
            'Original Pathway': pathway_name,
            'Higher-Ranked Pathway': row['Term'],
            'Rank': row['Rank'],
            'FDR q-val': row['FDR q-val'],
            'Highest Related Rank': highest_related_rank,
            'FP Above Related': fp_above_related
        })

    # Classify all pathways in the results, including Overlapping
    associated_pathway = pathway_name
    related_pathways = disease_to_pathways.get(associated_pathway, set())
    if not isinstance(related_pathways, set):
        related_pathways = set(related_pathways)
    related_pathways.discard(associated_pathway)  # Exclude the associated pathway if present

    # Retrieve gene sets for related pathways
    related_gene_sets = []
    for rel_pathway in related_pathways:
        # Attempt to get genes from results_df
        rel_genes_series = results_df.loc[results_df['Term'] == rel_pathway, 'Lead_genes']
        if not rel_genes_series.empty:
            genes = set(rel_genes_series.iloc[0].split(','))
        else:
            # Fetch genes from PATHWAYS_DATA if not present in results_df
            genes = get_lead_genes_from_pathways_data(rel_pathway, pathway_file)
            if not genes:
                logger.warning(f"Genes for related pathway '{rel_pathway}' not found in PATHWAYS_DATA.")
        related_gene_sets.append(genes)

    pathway_classifications = []
    for _, row in results_df.iterrows():
        pathway_term = row['Term']
        pathway_genes = set(row['Lead_genes'].split(',')) if pd.notna(row['Lead_genes']) else set()

        if pathway_term == associated_pathway:
            classification = 'Associated'
        elif pathway_term in related_pathways:
            classification = 'Related'
        else:
            # Calculate Jaccard similarity with all related pathways
            jaccard_indices = [jaccard_index(pathway_genes, rel_genes) for rel_genes in related_gene_sets]
            max_jaccard = max(jaccard_indices) if jaccard_indices else 0.0
            if max_jaccard >= jac_threshold:
                classification = 'Overlapping'
            else:
                classification = 'Unrelated'

        pathway_classifications.append({
            'Dataset': dataset_name,
            'Pathway': pathway_term,
            'Classification': classification
        })

    # Count the number of pathways in each classification
    classification_counts = pd.DataFrame(pathway_classifications)['Classification'].value_counts().to_dict()

    return {
        'Dataset': dataset_name,
        'Original Pathway': pathway_name,
        'Rank': rank,
        'FDR q-val': fdr_q_val,
        'Significant': int(fdr_q_val <= 0.05) if fdr_q_val is not None else 0,
        'Higher-Ranked Details': higher_ranked_details,
        'Classification Counts': classification_counts  # Include counts for plotting
    }

def get_lead_genes_from_pathways_data(rel_pathway: str, pathway_file: str) -> set:
    """
    Retrieve the Lead_genes for a given related pathway from PATHWAYS_DATA.

    Args:
        rel_pathway (str): Related pathway name.
        pathway_file (str): Pathway file name ('kegg' in this context).

    Returns:
        set: Set of lead genes or empty set if not found.
    """
    try:
        genes_list = PATHWAYS_DATA['kegg'].get(rel_pathway, [])
        return set(genes_list)
    except KeyError:
        logger.error(f"Pathway file '{pathway_file}' not found in PATHWAYS_DATA.")
        return set()


def analyze_disease_pathway_rankings(disease_to_pathways, higher_ranked_df, classification_df, summary_dir, method):
    """
    Analyze pathway rankings and classifications, including the new 'Overlapping' category.

    Args:
        disease_to_pathways (dict): Mapping of diseases to pathways.
        higher_ranked_df (pd.DataFrame): DataFrame of higher-ranked pathways.
        classification_df (pd.DataFrame): DataFrame of classification counts.
        summary_dir (str): Directory to save summaries.
        method (str): Method used for analysis.
    """
    analysis_results = []

    # Process each higher-ranked pathway
    higher_ranked_df = higher_ranked_df.sort_values(by=['Dataset', 'Original Pathway', 'Rank']).reset_index(drop=True)

    for _, row in higher_ranked_df.iterrows():
        dataset = row['Dataset']
        original_pathway = row['Original Pathway']
        higher_ranked_pathway = row['Higher-Ranked Pathway']
        rank = row['Rank']
        highest_related_rank = row['Highest Related Rank']
        fp_above_related = row['FP Above Related']

        related_pathways = disease_to_pathways.get(original_pathway, set())
        is_related = any(term.lower() in higher_ranked_pathway.lower() for term in related_pathways)
        classification = 'Related' if is_related else 'False Positive'

        # Append results for this row
        analysis_results.append({
            'Dataset': dataset,
            'Original Pathway': original_pathway,
            'Higher-Ranked Pathway': higher_ranked_pathway,
            'Rank': rank,
            'Classification': classification,
            'Highest Related Rank': highest_related_rank,
            'False Positives Above Top Related': fp_above_related
        })

    # Convert the analysis results into a DataFrame
    analysis_df = pd.DataFrame(analysis_results)

    # Generate plots including the new 'Overlapping' category
    plot_analysis_results(analysis_df, classification_df, summary_dir)

    # Save analysis results to a CSV file
    analysis_output_path = os.path.join(summary_dir, f'higher_ranked_pathway_analysis_{method}.csv')
    analysis_df.to_csv(analysis_output_path, index=False)
    logger.info(f"Higher-ranked pathway analysis saved to {analysis_output_path}")



FILE_LIST = get_file_list(DIRECTORIES['input'])
logger.info(f"Found {len(FILE_LIST)} .xlsx files in {DIRECTORIES['input']}")

NETWORKS_DATA, PATHWAYS_DATA = load_data(NETWORKS, PATHWAY_FILES, DIRECTORIES['pathways'])
logger.info("Networks loaded and pathway densities calculated.")

# Additional Diagnostic Function
def inspect_pathways_data():
    """
    Inspect the structure of PATHWAYS_DATA for each pathway_file.

    Helps in debugging the structure and contents of PATHWAYS_DATA.
    """
    for pathway_file, pathways in PATHWAYS_DATA.items():
        logger.info(f"Inspecting PATHWAYS_DATA for '{pathway_file}'")
        logger.info(f"Type: {type(pathways)}")
        if isinstance(pathways, list):
            if pathways:
                first_item = pathways[0]
                logger.info(f"First item type: {type(first_item)}")
                if isinstance(first_item, dict):
                    logger.info(f"First item keys: {first_item.keys()}")
                else:
                    logger.info(f"First item: {first_item}")
            else:
                logger.info("No pathways found in the list.")
        else:
            logger.error(f"Unexpected type for PATHWAYS_DATA[{pathway_file}]: {type(pathways)}")

def main():
    """
    Main function to execute the pathway analysis pipeline.
    """
    start_time = time.time()
    disease_to_pathways = PRE_CALCULATED_DATA
    logger.info("Loaded disease to pathway mappings.")
    
    # Inspect PATHWAYS_DATA structure
    inspect_pathways_data()
    
    futures = []
    higher_ranked_records = []  # List to collect all higher-ranked pathways across files
    classification_records = []  # List to collect classification counts for plotting

    network_name = NETWORKS[0]
    
    # Retrieve Jaccard threshold from config or set a default
    jac_threshold = config['analysis_settings'].get('jaccard_threshold', 0.2)

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        for file_name in FILE_LIST:
            dataset_name, pathway_name = file_name.replace('.xlsx', '').split('_', 1)
            if pathway_name in PATHWAYS_DATA['kegg']:
                futures.append(
                    executor.submit(
                        process_single_file,
                        'kegg', network_name, 0.2, 'PEANUT', file_name, disease_to_pathways, jac_threshold
                    )
                )

        for future in as_completed(futures):
            result = future.result()
            if result is None:
                continue
            higher_ranked_details = result.pop('Higher-Ranked Details', [])
            higher_ranked_records.extend(higher_ranked_details)

            # Collect classification counts
            classification_counts = result.pop('Classification Counts', {})
            dataset_name = result['Dataset']
            classification_counts['Dataset'] = dataset_name
            classification_records.append(classification_counts)

    # Convert higher-ranked pathways to DataFrame
    higher_ranked_df = pd.DataFrame(higher_ranked_records)

    # Convert classification counts to DataFrame
    classification_df = pd.DataFrame(classification_records).fillna(0)
    # Ensure all expected classifications are present
    for col in ['Associated', 'Related', 'Overlapping', 'Unrelated']:
        if col not in classification_df.columns:
            classification_df[col] = 0

    # Now we can analyze the higher-ranked pathways
    analyze_disease_pathway_rankings(disease_to_pathways, higher_ranked_df, classification_df, DIRECTORIES['summary_base'], 'PEANUT')

    elapsed_time = time.time() - start_time
    logger.info(f"Total time taken: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main()
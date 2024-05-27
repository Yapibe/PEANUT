import zlib
import pickle
import numpy as np
import pandas as pd
import networkx as nx
from os import path, makedirs


def filter_network_by_prior_data(network_filename, prior_data):
    """
    Filters a network to only include nodes present in the prior_data DataFrame.

    Args:
        network_filename (str): Path to the network file.
        prior_data (pd.DataFrame): DataFrame containing gene information.

    Returns:
        networkx.Graph: A filtered graph object.
    """
    # Read the network
    network = read_network(network_filename)

    # Get the set of gene IDs from prior_data
    gene_ids_in_prior_data = set(prior_data['GeneID'])

    # Filter the network
    filtered_network = network.copy()
    for node in network.nodes:
        if node not in gene_ids_in_prior_data:
            filtered_network.remove_node(node)

    return filtered_network


# noinspection PyTypeChecker
def read_network(network_filename):
    """
    Reads a network from a file and returns a NetworkX graph.

    Args:
        network_filename (str): Path to the file containing the network data.

    Returns:
        networkx.Graph: A graph object representing the network.
    """
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


def read_prior_set(excel_dir):
    """
    Reads prior data set from an Excel file and applies preprocessing.
    Args:
        excel_dir (str): Path to the Excel file containing the prior data.
    Returns:
        pandas.DataFrame: DataFrame containing the preprocessed prior data.
    """
    prior_data = pd.read_excel(excel_dir, engine='openpyxl')

    # Drop duplicate GeneID values
    # print all duplicates
    # print how many genes pre editing
    print(f"Number of genes before editing: {len(prior_data)}")
    prior_data = prior_data.drop_duplicates(subset='GeneID')
    # remove any row with no value in Score column
    prior_data = prior_data[prior_data['Score'].notna()]
    # remove any row with "?" in Score column
    prior_data = prior_data[prior_data['Score'] != '?']
    # Filter out GeneIDs that are not purely numeric (to exclude concatenated IDs)
    prior_data = prior_data[prior_data['GeneID'].apply(lambda x: str(x).isdigit())]
    # Convert GeneID to integer
    prior_data['GeneID'] = prior_data['GeneID'].astype(int)

    # Reset the DataFrame index
    prior_data = prior_data.reset_index(drop=True)
    print(f"Number of genes after editing: {len(prior_data)}")
    return prior_data


def get_propagation_input(prior_gene_ids, prior_data, input_type='Score'):
    """
    Generates propagation inputs based on specified type and network.

    Args:
        prior_gene_ids (set): List of gene IDs for the prior set.
        prior_data (pandas.DataFrame): DataFrame containing all experimental data.
        input_type (str): Type of input to generate (e.g., 'ones', 'abs_Score', etc.).

    Returns:
        dict: A dictionary mapping gene IDs to their corresponding input values.
    """
    inputs = dict()
    if input_type == 'ones' or input_type is None:
        inputs = {int(x): 1 for x in prior_gene_ids}
    elif input_type == 'abs_Score':
        for id in prior_gene_ids:
            try:
                # Your original code for processing each id
                inputs[int(id)] = np.abs(float(prior_data[prior_data.GeneID == id]['Score'].values[0]))
            except TypeError:
                print(f"Error processing ID: {id}")
                print(f"Data for this ID: {prior_data[prior_data.GeneID == id]}")
                # Optionally, you can break or continue based on your needs
                break

    elif input_type == 'Score':
        inputs = {int(id): float(prior_data[prior_data.GeneID == id]['Score'].values[0]) for id in prior_gene_ids}
    else:
        assert 0, '{} is not a valid input type'.format(input_type)

    inputs = {id: np.round(input_score, 3) for id, input_score in inputs.items()}
    return inputs


def save_file(obj, save_dir=None, compress=True):
    """
    Saves an object to a file, with optional compression.
    Args:
        obj (object): The object to be saved.
        save_dir (str, optional): The directory where the file will be saved.
        compress (bool, optional): Whether to compress the file.
    Returns:
        None
    """
    obj = pickle.dumps(obj)
    if compress:
        obj = zlib.compress(obj)
    with open(save_dir, 'wb') as f:
        pickle.dump(obj, f)
    print('File was saved in {}'.format(save_dir))


def save_propagation_score(propagation_scores, prior_set, propagation_input, genes_id_to_idx, task, general_args, save_dir=None):
    """
    Saves the propagation scores to a file.

    Args:
        propagation_scores (ndarray): The propagation scores.
        prior_set (dataframe): The set of prior genes.
        propagation_input (dict): The input used for propagation.
        genes_id_to_idx (dict): Mapping from gene indices to gene IDs.
        task (PropagationTask): The propagation task object containing program arguments.
        save_dir (str, optional): Directory to save the results.

    Returns:
        dict: A dictionary containing the saved data.
    """
    file_name = f"{task.test_name}_{general_args.alpha}_{general_args.date}"
    save_dir = save_dir
    makedirs(save_dir, exist_ok=True)
    propagation_results_path = path.join(save_dir, file_name)

    save_dict = {
        'args': task, 'prior_set': prior_set, 'propagation_input': propagation_input,
        'gene_id_to_idx': genes_id_to_idx, 'gene_prop_scores': propagation_scores,
    }

    save_file(save_dict, propagation_results_path)
    return save_dict


def load_pathways_genes(pathways_dir):
    """
    Loads the pathways and their associated genes from a file.

    Args:
        pathways_dir (str): Path to the file containing the pathway data.

    Returns:
        dict: A dictionary mapping pathway names to lists of genes in each pathway.
    """
    pathways = {}
    # Open the file containing pathway data
    try:
        with open(pathways_dir, 'r') as file:
            for line in file:
                # Process each line, normalize case, and split by tab
                parts = line.strip().upper().split('\t')
                # Skip lines that don't have at least 3 parts or where the second part isn't a digit
                if len(parts) < 3 or not parts[1].isdigit():
                    continue

                # Parse pathway name and expected size
                pathway_name, pathway_size = parts[0], int(parts[1])

                # Collect gene IDs ensuring they are numeric and don't exceed the pathway size
                genes = [int(gene) for gene in parts[2].split()[:pathway_size] if gene.isdigit()]
                pathways[pathway_name] = genes

    except FileNotFoundError:
        print(f"File not found: {pathways_dir}")
    except Exception as e:
        print(f"An error occurred while loading pathways: {e}")

    return pathways


def load_file(file_path, decompress=True):
    """
    Loads the propagation score data from a file.

    Args:
        file_path (str): The path to the file to be loaded.

    Returns:
        dict: The loaded data containing propagation scores and other information.
    """
    with open(file_path, 'rb') as file:
        decompressed_data = pickle.load(file)
    if decompress:
        try:
            decompressed_data = pickle.loads(zlib.decompress(decompressed_data))
        except zlib.error:
            print('Entered an uncompressed file but asked to decompress it')
    return decompressed_data


def load_propagation_scores(propagation_file_path):
    propagation_result_dict: dict = load_file(propagation_file_path, decompress=True)
    return propagation_result_dict


def load_network_and_pathways(general_args, propagation_file_path):
    """
    Loads the network graph and pathways based on the provided configuration.
    Returns:
        tuple: A tuple containing the network graph, a list of interesting pathways, and a dictionary mapping
               pathways to their genes.
    """
    # network_graph = read_network(general_args.network_file_path)
    genes_by_pathway = load_pathways_genes(general_args.pathway_file_dir)
    scores = get_scores(propagation_file_path)

    return genes_by_pathway, scores


def get_scores(score_path):
    """
        Loads gene scores and P-values from a file and returns a dictionary mapping gene IDs to their scores and P-values.

        Args:
            score_path (str): The path to the file containing the propagation scores.

        Returns:
            dict: A dictionary with gene IDs as keys and tuples of (Score, P-value) as values.
        """
    try:
        # Load propagation results from the specified file
        propagation_results = load_propagation_scores(score_path)

        # Sort the propagation scores by GeneID
        sorted_scores = propagation_results['gene_prop_scores'].sort_values(by='GeneID').reset_index(drop=True)

        # Create a dictionary mapping gene IDs to tuples of (Score, P-value)
        gene_scores = {gene_id: (score, pvalue) for gene_id, score, pvalue
                       in zip(sorted_scores['GeneID'], sorted_scores['Score'], sorted_scores['P-value'])}

        return gene_scores

    except FileNotFoundError:
        print(f"File not found: {score_path}")
        return {}
    except Exception as error:
        print(f"An error occurred: {error}")
        return {}


def read_pathways(file_name):
    """
    Reads pathways from a file into a dictionary mapping each pathway to a list of gene IDs.

    Parameters:
    - file_name (str): Path to the file containing the pathways.

    Returns:
    dict: A dictionary mapping pathway names to lists of gene IDs.
    """
    try:
        with open(file_name, 'r') as file:
            return {line.split()[0]: [int(gene) for gene in line.split()[2:]] for line in file if line.strip()}
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return {}
    except Exception as e:
        print(f"Error reading pathways from {file_name}: {e}")
        return {}


def read_temp_scores(file_name):
    """
    Read scores from a file into a dictionary.

    Parameters:
    - file_name (str): Path to the file containing the scores.

    Returns:
    dict: A dictionary mapping pathways to their scores.
    """
    try:
        scores = pd.read_csv(file_name, sep=' ', header=None, names=['Pathway', 'Score'], index_col='Pathway')
        return scores['Score'].to_dict()
    except FileNotFoundError:
        print(f"File not found: {file_name}")
        return {}
    except Exception as e:
        print(f"Error reading scores from {file_name}: {e}")
        return {}


def process_condition(condition_file, experiment_file, pathways_file, all_pathways, P_VALUE_THRESHOLD=0.05):
    """
    Processes experimental data to determine the enrichment and trends of pathways based on specified conditions.

    Parameters:
    - condition_file (str): Path to the file containing score conditions.
    - experiment_file (str): Path to the file containing experimental data.
    - pathways_file (str): Path to the file containing pathways data.

    Returns:
    dict: A dictionary containing scores, enriched pathway genes, and mean scores for pathways under the given condition.
    """
    # Read scores for the condition, mapping pathways to their scores
    enriched_pathway_dict = read_temp_scores(condition_file)

    # Load experiment data and filter out entries where the score is zero. experiment file is xlsx
    condition_data_df = pd.read_excel(experiment_file)
    experiment_data_filtered_df = condition_data_df[condition_data_df['Score'] != 0]

    # Extract the condition name from the file name
    condition_name = path.basename(condition_file).split('.')[-2]

    # Load pathway data mapping pathway names to lists of gene IDs
    homo_sapien_pathway_dict = read_pathways(pathways_file)

    # Dictionary to store enriched pathway genes
    enriched_pathway_genes = {}

    # Loop through each pathway
    for pathway in all_pathways:
        # Initialize a dictionary for the pathway under the current condition
        all_pathways[pathway][condition_name] = {}

        # List of genes associated with the current pathway
        pathway_genes = homo_sapien_pathway_dict[pathway]

        # Filter the experiment data to only include genes that are part of the current pathway
        pathway_filtered_genes = experiment_data_filtered_df[experiment_data_filtered_df['GeneID'].isin(pathway_genes)]

        # Store details of filtered genes in a dictionary
        enriched_pathway_genes[pathway] = pathway_filtered_genes.set_index('GeneID')[['Symbol', 'Score', 'P-value']].to_dict(
            orient='index')

        # Filter to find significant genes based on the P-value threshold
        significant_genes = {gene_id: gene_details for gene_id, gene_details in
                             enriched_pathway_genes[pathway].items() if
                             gene_details['P-value'] <= P_VALUE_THRESHOLD}

        # Calculate the mean score of significant genes or set to 0 if none are significant
        mean_score = np.mean(
            [gene_details['Score'] for gene_details in significant_genes.values()]) if significant_genes else 0

        # Store the mean score and significant genes for the pathway under the condition
        all_pathways[pathway][condition_name]['Mean'] = mean_score
        all_pathways[pathway][condition_name]['significant_genes'] = significant_genes

        # Check if the pathway is in the enriched pathway dictionary to assign a P-value and trend
        if pathway in enriched_pathway_dict:
            all_pathways[pathway][condition_name]['P-value'] = enriched_pathway_dict[pathway]
            all_pathways[pathway][condition_name]['Trend'] = "Up*" if mean_score > 0 else "Down*"
        else:
            all_pathways[pathway][condition_name]['P-value'] = 1  # Default P-value if not in enriched dictionary
            all_pathways[pathway][condition_name]['Trend'] = "Up" if mean_score > 0 else "Down"
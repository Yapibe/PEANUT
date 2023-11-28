from args import PropagationTask, PathwayResults
import utils
from visualization_tools import plot_enrichment_table
from os import path
from scipy.stats import hypergeom
from decimal import Decimal
from scipy.stats import rankdata
import numpy as np
from statsmodels.stats.multitest import multipletests


def get_experiment_gene_set(prior_data):
    return set(prior_data['GeneID'])


def filter_pathways_genes(genes_by_pathway, gene_set):
    # Filter out genes not in the experiment and pathways not within the size limits
    filtered_pathways = {}
    for pathway, genes in genes_by_pathway.items():
        filtered_genes = [gene for gene in genes if gene in gene_set]
        if pathway == 'REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION':
            filtered_pathways[pathway] = filtered_genes
        elif 10 <= len(filtered_genes) <= 60:
            filtered_pathways[pathway] = filtered_genes
    return filtered_pathways


def create_matrices(filtered_results):
    num_pathways = len(filtered_results)
    adj_p_vals_mat = np.zeros(num_pathways)
    directions_mat = np.zeros(num_pathways, dtype=bool)

    for i, (pathway_name, result) in enumerate(filtered_results.items()):
        adj_p_vals_mat[i] = result.adj_p_value
        directions_mat[i] = result.direction

    # Reshape the arrays to 2D column vectors
    adj_p_vals_mat = adj_p_vals_mat.reshape(-1, 1)
    directions_mat = directions_mat.reshape(-1, 1)

    return adj_p_vals_mat, directions_mat


def iterative_factorial(n):
    """Compute factorial iteratively to avoid recursion depth issues."""
    result = 1
    for i in range(2, n + 1):
        result *= i
    return result


def binomial_coefficient(n, k):
    """Calculate the binomial coefficient 'n choose k'."""
    return iterative_factorial(n) // (iterative_factorial(k) * iterative_factorial(n - k))


def hypergeometric_pmf(K, k, N, n):
    """Calculate the probability mass function for the hypergeometric distribution."""
    # # Number of ways to choose k successes from K possible successes
    # success_ways = binomial_coefficient(K, k)
    # # Number of ways to choose n-k failures from N-K possible failures
    # failure_ways = binomial_coefficient(N - K, n - k)
    # # Total number of ways to choose n draws out of N possible draws
    # total_ways = binomial_coefficient(N, n)
    #
    # # Probability calculation
    # probability = success_ways * failure_ways / total_ways

    # use scipy.stats.hypergeom
    probability = hypergeom.pmf(k, N, K, n)
    return probability


def hypergeom_bh_correction(p_values_dict):
    """
    Perform the Benjamini-Hochberg correction on a dictionary of p-values.

    :param p_values_dict: Dictionary of p-values with pathway keys
    :return: Dictionary of adjusted p-values with pathway keys
    """
    # Extract p-values and keys
    pathways = list(p_values_dict.keys())
    p_values = np.array(list(p_values_dict.values()))

    # Perform BH correction using statsmodels
    _, adj_p_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Re-associate the adjusted p-values with their corresponding pathways
    adj_p_values_dict = dict(zip(pathways, adj_p_values))

    return adj_p_values_dict


def get_p_values(prior_data, threshold=Decimal('0.01')):
    """
    Get p-values for genes and filter based on the threshold.

    :param prior_data: DataFrame with GeneID and P-value columns
    :param threshold: Decimal p-value threshold for filtering
    :return: Dictionary of gene IDs with their corresponding p-values
             if they are below the threshold, as Decimal objects.
    """
    # Convert the P-value column to Decimals to avoid precision loss
    prior_data['P-value'] = prior_data['P-value'].apply(Decimal)

    # Filter the DataFrame based on the threshold
    filtered_data = prior_data[prior_data['P-value'] <= threshold]

    # Create a dictionary with gene IDs and their corresponding p-values
    inputs = {
        int(gene_id): p_value
        for gene_id, p_value in zip(filtered_data['GeneID'], filtered_data['P-value'])
    }
    return inputs


# run a hypergeometric test on the results of the propagation
# load prior set
# create a propagation task
task = PropagationTask(experiment_name='hyper', create_similarity_matrix=False)
prior_data_df = utils.read_prior_set(task.experiment_file_path)
# get the p-values
enriched_p_vals = get_p_values(prior_data_df)

# get pathway gene sets
dir_path = path.dirname(path.realpath(__file__))
genes_by_pathway = utils.load_pathways_genes(path.join(dir_path, '../Data', 'H_sapiens', 'pathways', 'pathway_file'))
experiment_gene_set = get_experiment_gene_set(prior_data_df)
filtered_genes_by_pathway = filter_pathways_genes(genes_by_pathway, experiment_gene_set)


# Total number of genes in the background
M = len(prior_data_df['GeneID'])
# Number of enriched genes
n = len(enriched_p_vals)

# Perform the hypergeometric test for each pathway
for pathway_name, pathway_genes in filtered_genes_by_pathway.items():
    # Number of genes in the pathway
    N = len(pathway_genes)

    # Number of enriched genes in the pathway, i.e. the intersection of the pathway and the enriched genes
    x = len(set(pathway_genes).intersection(set(enriched_p_vals.keys())))

    # Calculate the hypergeometric p-value
    pval = hypergeometric_pmf(N, x, M, n)
    # calculate the direction by the mean of the gene scores of the pathway
    direction = np.mean(prior_data_df[prior_data_df['GeneID'].isin(pathway_genes)]['Score']) > 0
    # pathway_pvals[pathway_name] = pval
    task.results[pathway_name] = PathwayResults(p_value=pval, direction=direction)

specific_pathways = [
    "WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV",
    "WP_HIPPOCAMPAL_SYNAPTOGENESIS_AND_NEUROGENESIS",
    "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
    "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION"

]

for pathway in specific_pathways:
    if pathway in task.results:
        print(f'{pathway}: {task.results[pathway].p_value}')

# Perform the Benjamini-Hochberg correction
adj_pvals = hypergeom_bh_correction({k: v.p_value for k, v in task.results.items()})

# keep only the pathways with adjusted p-values below 0.01
adj_pvals = {k: v for k, v in adj_pvals.items() if v <= 0.05}

# add the adjusted p-values to the results
for pathway_name, adj_pval in adj_pvals.items():
    task.results[pathway_name].adj_p_value = adj_pval

# Filter PathwayResults based on adjusted p-values
filtered_results = {k: task.results[k] for k in adj_pvals.keys()}

# Create matrices
adj_p_vals_mat, directions_mat = create_matrices(filtered_results)

for pathway in specific_pathways:
    if pathway in adj_pvals:
        print(f'{pathway}: {adj_pvals[pathway]}')

# Get dict of genes and their names from the priordata[Human_Name]
genes_names = {int(gene_id): gene_name for gene_id, gene_name in zip(prior_data_df['GeneID'], prior_data_df['Human_Name'])}

# Save the results to a file
with open('11_28_23_HG.txt', 'w') as f:
    for pathway, pval in adj_pvals.items():
        # Retrieve gene names for the pathway and ensure they are strings
        gene_names_in_pathway = [str(genes_names[gene_id]) for gene_id in filtered_genes_by_pathway[pathway] if gene_id in genes_names]

        # Filter out any non-string values (like NaNs) and join the names
        gene_names_str = ', '.join([name for name in gene_names_in_pathway if isinstance(name, str)])

        # Write pathway, p-value, and gene names to the file
        f.write(f'{pathway}: {pval}\nPath Genes: {gene_names_str}\n\n')

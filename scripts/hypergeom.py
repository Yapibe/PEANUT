from args import PropagationTask
import utils
from pathway_enrichment import load_pathways_genes
from os import path
from scipy.stats import hypergeom
from decimal import Decimal
from scipy.stats import rankdata
import numpy as np
from statsmodels.stats.multitest import multipletests


def iterative_factorial(n):
    """Compute factorial iteratively to avoid recursion depth issues."""
    result = 1
    for i in range(2, n+1):
        result *= i
    return result


def binomial_coefficient(n, k):
    """Calculate the binomial coefficient 'n choose k'."""
    return iterative_factorial(n) // (iterative_factorial(k) * iterative_factorial(n - k))


def hypergeometric_pmf(K, k, N, n):
    """Calculate the probability mass function for the hypergeometric distribution."""
    # Number of ways to choose k successes from K possible successes
    success_ways = binomial_coefficient(K, k)
    # Number of ways to choose n-k failures from N-K possible failures
    failure_ways = binomial_coefficient(N - K, n - k)
    # Total number of ways to choose n draws out of N possible draws
    total_ways = binomial_coefficient(N, n)

    # Probability calculation
    probability = success_ways * failure_ways / total_ways
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
genes_by_pathway = load_pathways_genes(path.join(dir_path, '../Data', 'H_sapiens', 'pathways', 'pathway_file'))

# Total number of genes in the background
M = len(prior_data_df['GeneID'])
# Number of enriched genes
n = len(enriched_p_vals)

# Dictionary to hold pathway p-values
pathway_pvals = {}

# Perform the hypergeometric test for each pathway
for pathway_name, pathway_genes in genes_by_pathway.items():
    # Number of genes in the pathway
    N = len(pathway_genes)

    # Number of enriched genes in the pathway, i.e. the intersection of the pathway and the enriched genes
    x = len(set(pathway_genes).intersection(set(enriched_p_vals.keys())))

    # Calculate the hypergeometric p-value
    pval = hypergeometric_pmf(N, x, M, n)

    pathway_pvals[pathway_name] = pval

# Perform the Benjamini-Hochberg correction
adj_pvals = hypergeom_bh_correction(pathway_pvals)

# keep only the pathways with adjusted p-values below 0.01
adj_pvals = {k: v for k, v in adj_pvals.items() if v <= 0.01}

specific_pathways = [
        "REACTOME_PRESYNAPTIC_DEPOLARIZATION_AND_CALCIUM_CHANNEL_OPENING",
        "REACTOME_NEUROTRANSMITTER_RECEPTORS_AND_POSTSYNAPTIC_SIGNAL_TRANSMISSION",
        "REACTOME_DOPAMINE_CLEARANCE_FROM_THE_SYNAPTIC_CLEFT",
        "REACTOME_ACTIVATION_OF_NMDA_RECEPTORS_AND_POSTSYNAPTIC_EVENTS",
        "REACTOME_PRESYNAPTIC_FUNCTION_OF_KAINATE_RECEPTORS",
        "REACTOME_HIGHLY_SODIUM_PERMEABLE_POSTSYNAPTIC_ACETYLCHOLINE_NICOTINIC_RECEPTORS",
        "REACTOME_HIGHLY_CALCIUM_PERMEABLE_POSTSYNAPTIC_NICOTINIC_ACETYLCHOLINE_RECEPTORS",
        "REACTOME_SYNAPTIC_ADHESION_LIKE_MOLECULES",
        "WP_SPLICING_FACTOR_NOVA_REGULATED_SYNAPTIC_PROTEINS",
        "WP_DISRUPTION_OF_POSTSYNAPTIC_SIGNALING_BY_CNV",
        "WP_SYNAPTIC_VESICLE_PATHWAY",
        "WP_SYNAPTIC_SIGNALING_PATHWAYS_ASSOCIATED_WITH_AUTISM_SPECTRUM_DISORDER",
        "KEGG_LYSOSOME",
        "REACTOME_LYSOSPHINGOLIPID_AND_LPA_RECEPTORS",
        "REACTOME_LYSOSOME_VESICLE_BIOGENESIS",
        "REACTOME_PREVENTION_OF_PHAGOSOMAL_LYSOSOMAL_FUSION",
        "PID_LYSOPHOSPHOLIPID_PATHWAY",
        "BIOCARTA_CARM_ER_PATHWAY",
        "REACTOME_SYNTHESIS_OF_PIPS_AT_THE_ER_MEMBRANE",
        "REACTOME_ER_TO_GOLGI_ANTEROGRADE_TRANSPORT",
        "REACTOME_N_GLYCAN_TRIMMING_IN_THE_ER_AND_CALNEXIN_CALRETICULIN_CYCLE",
        "REACTOME_COPI_DEPENDENT_GOLGI_TO_ER_RETROGRADE_TRAFFIC",
        "REACTOME_COPI_INDEPENDENT_GOLGI_TO_ER_RETROGRADE_TRAFFIC",
        "REACTOME_INTRA_GOLGI_AND_RETROGRADE_GOLGI_TO_ER_TRAFFIC",
        "REACTOME_GOLGI_TO_ER_RETROGRADE_TRANSPORT",
        "REACTOME_ER_QUALITY_CONTROL_COMPARTMENT_ERQC",
        "WP_METABOLISM_OF_SPHINGOLIPIDS_IN_ER_AND_GOLGI_APPARATUS",
        "PID_ER_NONGENOMIC_PATHWAY",
        "WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING",
        "WP_RELATIONSHIP_BETWEEN_INFLAMMATION_COX2_AND_EGFR",
        "WP_RESISTIN_AS_A_REGULATOR_OF_INFLAMMATION",
        "WP_APOE_AND_MIR146_IN_INFLAMMATION_AND_ATHEROSCLEROSIS",
        "WP_SUPRESSION_OF_HMGB1_MEDIATED_INFLAMMATION_BY_THBD",
        "WP_RESOLVIN_E1_AND_RESOLVIN_D1_SIGNALING_PATHWAYS_PROMOTING_INFLAMMATION_RESOLUTION",
        "WP_NEUROINFLAMMATION"
    ]  # List of specific pathways

# check if any of the specific pathways are in the results
for pathway in specific_pathways:
    if pathway in adj_pvals:
        print(f'{pathway}: {adj_pvals[pathway]}')

# save the results to a file
with open('../old.txt', 'w') as f:
    for pathway, pval in adj_pvals.items():
        f.write(f'{pathway}: {pval}\n')






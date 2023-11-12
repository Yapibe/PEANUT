import numpy as np
import random
from os import path
from collections import defaultdict
from decimal import Decimal
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from args import PropagationTask
import utils
from pathway_enrichment import load_pathways_genes
import pandas as pd


def load_and_prepare_data(task):
    prior_data = utils.read_prior_set(task.experiment_file_path)
    prior_data['Score'] = prior_data['Score'].apply(lambda x: abs(Decimal(x)))
    prior_data = prior_data[prior_data['GeneID'].apply(lambda x: str(x).isdigit())].copy()
    prior_data['GeneID'] = prior_data['GeneID'].astype(int)
    return prior_data


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


def shuffle_gene_ranks(ranks):
    # Shuffle the ranks and return a new dictionary with the same keys but shuffled values
    shuffled_values = list(ranks.values())
    random.shuffle(shuffled_values)
    return dict(zip(ranks.keys(), shuffled_values))


def generate_random_scores(ranks, pathways, num_sim=10000):
    random_scores = defaultdict(list)
    for _ in range(num_sim):
        shuffled_ranks = shuffle_gene_ranks(ranks)
        for pathway_name, pathway_genes in pathways.items():
            score = 0
            for gene in pathway_genes:
                if gene not in shuffled_ranks:
                    print(f"Pathway: {pathway_name}, Missing Gene: {gene}")
                score += shuffled_ranks.get(gene, 0)
            random_scores[pathway_name].append(score)
    return random_scores


def calculate_pathway_score(ranks, pathway_genes):
    return sum(ranks.get(gene, 0) for gene in pathway_genes)


def compare_with_random_scores(observed_score, random_scores):
    score = sum(score < observed_score for score in random_scores) + 1
    return score


def calculate_empirical_p_value(observed_score, random_scores):
    return compare_with_random_scores(observed_score, random_scores) / (len(random_scores) + 1)


# Replace the original calculate_empirical_p_values function
def calculate_empirical_p_values(ranks, pathways, random_scores):
    empirical_p_values = {}
    for pathway_name, pathway_genes in pathways.items():
        observed_score = calculate_pathway_score(ranks, pathway_genes)
        random_scores_list = random_scores[pathway_name]
        empirical_p_values[pathway_name] = calculate_empirical_p_value(observed_score, random_scores_list)
        if empirical_p_values[pathway_name] < 0.005:
            print(f'Pathway: {pathway_name}, Empirical P-Value: {empirical_p_values[pathway_name]}')
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


def filter_significant_pathways(pathways, p_values, threshold=0.05):
    return {pathway: pval for pathway, pval in zip(pathways, p_values) if pval < threshold}


def save_results(filename, data):
    with open(filename, 'w') as file:
        for pathway, pval in data.items():
            file.write(f'{pathway}: {pval}\n')


# Main analysis workflow
task = PropagationTask(experiment_name='random', create_similarity_matrix=False)
prior_data = load_and_prepare_data(task)
gene_ranks = get_fold_change_ranks(prior_data)
experiment_gene_set = set(gene_ranks.keys())

dir_path = path.dirname(path.realpath(__file__))
pathway_file_path = path.join(dir_path, 'Data', 'H_sapiens', 'pathways', 'pathway_file')
genes_by_pathway = load_pathways_genes(pathway_file_path)

# Filter pathways to include only genes in the experiment and with sizes 10-60
filtered_genes_by_pathway = filter_pathways_genes(genes_by_pathway, experiment_gene_set)

random_scores = generate_random_scores(gene_ranks, filtered_genes_by_pathway)
empirical_p_values = calculate_empirical_p_values(gene_ranks, filtered_genes_by_pathway, random_scores)

p_values = np.array(list(empirical_p_values.values()))
corrected_p_values = benjamini_hochberg_correction(p_values)
significant_pathways = filter_significant_pathways(empirical_p_values.keys(), corrected_p_values)
# index correct p-values by pathway name
corrected_p_values = {pathway: pval for pathway, pval in zip(empirical_p_values.keys(), corrected_p_values)}

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

csv_data = []
for pathway in specific_pathways:
    if pathway in significant_pathways:
        print(f'{pathway}: {significant_pathways[pathway]}')

for pathway in specific_pathways:
    if pathway in empirical_p_values:
        csv_data.append({
            'Pathway': pathway,
            'Empirical P-Value': empirical_p_values[pathway],
            'Corrected P-Value': corrected_p_values[pathway],
            'Gene Count': len(genes_by_pathway[pathway])
        })

# Create DataFrame
df = pd.DataFrame(csv_data)
df.to_csv('Mann.csv', index=False)

save_results('Mann.txt', significant_pathways)

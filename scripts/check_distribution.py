from os import path
from collections import Counter, defaultdict

# Open pathway file and check how many genes in each pathway
pathfile = path.join(path.dirname(path.dirname(path.realpath(__file__))), 'Data', 'H_sapiens', 'pathways',
                     'pathway_file')

with open(pathfile, 'r') as f:
    pathway = f.readlines()

# Dictionary to hold the count of genes in each pathway
pathway_gene_count = {}

# Iterate through each line in the pathway file content
for line in pathway:
    parts = line.strip().split('\t')  # Split the line by tab
    pathway_name = parts[0]  # Get the pathway name
    genes = parts[1:]  # Get the genes (all elements after the first)
    num_genes = len(genes)  # Count the number of genes
    pathway_gene_count[pathway_name] = num_genes  # Store the count in the dictionary

# Extract the number of genes for each pathway
plot_data = list(pathway_gene_count.values())

# Use Counter to get the frequency of each gene count
gene_count_frequency = Counter(plot_data)

# Sort the dictionary by gene count for easier interpretation
sorted_gene_count_frequency = {k: gene_count_frequency[k] for k in sorted(gene_count_frequency)}

print("Sorted gene count frequency:", sorted_gene_count_frequency)

# Create a dictionary to store the frequency of gene counts in bins of 10
binned_gene_count_frequency = defaultdict(int)

# Iterate through the sorted_gene_count_frequency to populate the binned dictionary
for gene_count, freq in sorted_gene_count_frequency.items():
    bin_index = gene_count // 10  # Integer division to get the bin index
    bin_key = f"{bin_index * 10}-{(bin_index + 1) * 10 - 1}"  # Create the bin key, e.g., "0-9", "10-19", etc.
    binned_gene_count_frequency[bin_key] += freq  # Add the frequency to the bin

# Convert the defaultdict to a regular dict for easier reading
binned_gene_count_frequency = dict(binned_gene_count_frequency)

print("Binned gene count frequency:", binned_gene_count_frequency)

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
]

# Create a dictionary to hold the count of genes in each pathway
specific_pathway_gene_count = {}

# Iterate through each line in the pathway file content
for specific_pathway in specific_pathways:
    # Get the genes for the pathway from the genes_by_pathway dictionary
    num_genes = pathway_gene_count[specific_pathway]
    specific_pathway_gene_count[specific_pathway] = num_genes  # Store the count in the dictionary

# sort the dictionary based on smallest to largest gene count
# sort the dictionary based on the gene count (smallest to largest)
sorted_specific_pathway_gene_count = {k: v for k, v in sorted(specific_pathway_gene_count.items(), key=lambda item: item[1])}
print("Sorted specific pathway gene count:", sorted_specific_pathway_gene_count)

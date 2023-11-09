import pandas as pd
import numpy as np
from os import path
import utils

main_dir = path.dirname(path.dirname(path.realpath(__file__)))

network = utils.read_network(path.join(main_dir, 'Data', 'H_sapiens', 'network', 'H_sapiens.net'))
# Create a list of all GeneIDs
all_gene_ids = list(network.nodes)
print(len(all_gene_ids))

# Randomly select 20000 GeneIDs
selected_gene_ids = np.random.choice(all_gene_ids, 15000, replace=False)

# Simulate scores from a standard normal distribution (N(0, 1)) for all selected GeneIDs
scores_standard = np.random.normal(10, 1, len(selected_gene_ids))

# Simulate a new set of scores from a standard normal distribution (N(0, 1)) for 'Score_Unique'
scores_unique_full_set = np.random.normal(10, 1, len(selected_gene_ids))

# Create a DataFrame for the selected GeneIDs and their scores
df = pd.DataFrame({
    'GeneID': selected_gene_ids,
    'Score_Standard': scores_standard,
    'Score_Unique': scores_unique_full_set
})

# Randomly select 50 genes to belong to the decoy pathway
decoy_genes = np.random.choice(df['GeneID'], 50, replace=False)

# Generate scores for the decoy genes from a normal distribution with mean=1 and std=1
decoy_scores = np.random.normal(4, 1, len(decoy_genes))

# Update the scores for the selected decoy genes in the DataFrame
df.loc[df['GeneID'].isin(decoy_genes), 'Score_Unique'] = decoy_scores

# Calculate log2 fold-change
# Adding a small constant to avoid division by zero and log of zero
constant = 1e-9
df['Log2FoldChange'] = np.log2((df['Score_Unique'] + constant) / (df['Score_Standard'] + constant))


# Prepare the decoy genes string
decoy_genes_str = 'AAAA_PATHWAY\t' + '\t'.join(map(str, decoy_genes))

#reorder columns
df = df[['GeneID', 'Log2FoldChange', 'Score_Unique', 'Score_Standard']]

# change the column names
df.columns = ['GeneID', 'Score', 'Score_Unique', 'Score_Standard']
# Save the DataFrame to a file in the Data folder
df.to_excel(path.join(main_dir, 'Inputs', 'experiments_data', 'simulated_scores.xlsx'), index=False)

print(decoy_genes_str)
df.head()

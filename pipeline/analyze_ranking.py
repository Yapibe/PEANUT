import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the pathways file with correct headers
pathways_file_path = 'Data/H_sapiens/pathways/kegg.gmt'
pathways_df = pd.read_csv(pathways_file_path, header=None)  # Read without header

# Set the correct column names
pathways_df.columns = ['Pathway', 'Genes']

# Ensure 'Genes' is split into individual gene entries if needed
pathways_df['Genes'] = pathways_df['Genes'].apply(lambda x: x.split('\\t'))

# Load the rankings summary
rankings_path = 'Outputs/NGSEA/rankings_summary.xlsx'
rankings_df = pd.read_excel(rankings_path)

# Calculate performance differences
rankings_df['Diff_vs_GSEA'] = rankings_df['PROP'] - rankings_df['GSEA']
rankings_df['Diff_vs_NGSEA'] = rankings_df['PROP'] - rankings_df['NGSEA']
rankings_df['Diff_vs_ABS_PROP'] = rankings_df['PROP'] - rankings_df['ABS-PROP']

# Classify pathways
rankings_df['Better_than_GSEA'] = rankings_df['Diff_vs_GSEA'] < 0
rankings_df['Better_than_NGSEA'] = rankings_df['Diff_vs_NGSEA'] < 0
rankings_df['Better_than_ABS_PROP'] = rankings_df['Diff_vs_ABS_PROP'] < 0
rankings_df['Worse_than_GSEA'] = rankings_df['Diff_vs_GSEA'] > 0
rankings_df['Worse_than_NGSEA'] = rankings_df['Diff_vs_NGSEA'] > 0
rankings_df['Worse_than_ABS_PROP'] = rankings_df['Diff_vs_ABS_PROP'] > 0

# Add pathway sizes
pathway_sizes = pathways_df.explode('Genes').groupby('Pathway').size().reset_index(name='Size')
rankings_df = rankings_df.merge(pathway_sizes, on='Pathway', how='left')

# Summary statistics
better_than_gsea = rankings_df[rankings_df['Better_than_GSEA']]
worse_than_gsea = rankings_df[rankings_df['Worse_than_GSEA']]
better_than_ngsea = rankings_df[rankings_df['Better_than_NGSEA']]
worse_than_ngsea = rankings_df[rankings_df['Worse_than_NGSEA']]
better_than_abs_prop = rankings_df[rankings_df['Better_than_ABS_PROP']]
worse_than_abs_prop = rankings_df[rankings_df['Worse_than_ABS_PROP']]

# Print summary statistics
print("Summary Statistics:")
print(f"Number of pathways better than GSEA: {len(better_than_gsea)}")
print(f"Number of pathways worse than GSEA: {len(worse_than_gsea)}")
print(f"Number of pathways better than NGSEA: {len(better_than_ngsea)}")
print(f"Number of pathways worse than NGSEA: {len(worse_than_ngsea)}")
print(f"Number of pathways better than ABS-PROP: {len(better_than_abs_prop)}")
print(f"Number of pathways worse than ABS-PROP: {len(worse_than_abs_prop)}")

# Perform statistical analysis (e.g., mean difference, t-test)
mean_diff_gsea = np.mean(rankings_df['Diff_vs_GSEA'])
mean_diff_ngsea = np.mean(rankings_df['Diff_vs_NGSEA'])
mean_diff_abs_prop = np.mean(rankings_df['Diff_vs_ABS_PROP'])
std_diff_gsea = np.std(rankings_df['Diff_vs_GSEA'])
std_diff_ngsea = np.std(rankings_df['Diff_vs_NGSEA'])
std_diff_abs_prop = np.std(rankings_df['Diff_vs_ABS_PROP'])

print(f"Mean difference vs GSEA: {mean_diff_gsea} (std: {std_diff_gsea})")
print(f"Mean difference vs NGSEA: {mean_diff_ngsea} (std: {std_diff_ngsea})")
print(f"Mean difference vs ABS-PROP: {mean_diff_abs_prop} (std: {std_diff_abs_prop})")

# Visualize results
plt.figure(figsize=(18, 6))

# Histogram of differences
plt.subplot(1, 3, 1)
plt.hist(rankings_df['Diff_vs_GSEA'], bins=30, alpha=0.7, label='Diff vs GSEA')
plt.hist(rankings_df['Diff_vs_NGSEA'], bins=30, alpha=0.7, label='Diff vs NGSEA')
plt.hist(rankings_df['Diff_vs_ABS_PROP'], bins=30, alpha=0.7, label='Diff vs ABS-PROP')
plt.xlabel('Difference in Rank')
plt.ylabel('Frequency')
plt.legend()
plt.title('Histogram of Rank Differences')

# Boxplot of differences
plt.subplot(1, 3, 2)
plt.boxplot([rankings_df['Diff_vs_GSEA'], rankings_df['Diff_vs_NGSEA'], rankings_df['Diff_vs_ABS_PROP']],
            labels=['Diff vs GSEA', 'Diff vs NGSEA', 'Diff vs ABS-PROP'])
plt.ylabel('Difference in Rank')
plt.title('Boxplot of Rank Differences')

# Size of pathways
better_than_gsea_sizes = better_than_gsea['Size']
worse_than_gsea_sizes = worse_than_gsea['Size']
better_than_ngsea_sizes = better_than_ngsea['Size']
worse_than_ngsea_sizes = worse_than_ngsea['Size']
better_than_abs_prop_sizes = better_than_abs_prop['Size']
worse_than_abs_prop_sizes = worse_than_abs_prop['Size']

# Boxplot of pathway sizes
plt.subplot(1, 3, 3)
plt.boxplot([better_than_gsea_sizes, worse_than_gsea_sizes,
             better_than_ngsea_sizes, worse_than_ngsea_sizes,
             better_than_abs_prop_sizes, worse_than_abs_prop_sizes],
            labels=['Better GSEA', 'Worse GSEA', 'Better NGSEA', 'Worse NGSEA', 'Better ABS-PROP', 'Worse ABS-PROP'])
plt.ylabel('Pathway Size')
plt.title('Boxplot of Pathway Sizes')

plt.tight_layout()
plt.savefig('Outputs/NGSEA/rank_differences_analysis.png')
plt.close()

print("Rank differences analysis saved to Outputs/NGSEA/rank_differences_analysis.png")

# Function to extract genes from pathway
def extract_genes_from_pathway(pathway_name):
    # Extract genes based on your data
    genes = pathways_df[pathways_df['Pathway'] == pathway_name]['Genes'].values[0]
    return genes

# Additional analyses for important genes
def analyze_important_genes(rankings_df, condition):
    pathways = rankings_df[condition]['Pathway'].unique()
    important_genes = []
    for pathway in pathways:
        genes = extract_genes_from_pathway(pathway)
        important_genes.extend(genes)
    return pd.Series(important_genes).value_counts().head(10)

print("Top 10 important genes where pipeline is better than GSEA:")
print(analyze_important_genes(rankings_df, rankings_df['Better_than_GSEA']))

print("Top 10 important genes where pipeline is worse than GSEA:")
print(analyze_important_genes(rankings_df, rankings_df['Worse_than_GSEA']))

print("Top 10 important genes where pipeline is better than NGSEA:")
print(analyze_important_genes(rankings_df, rankings_df['Better_than_NGSEA']))

print("Top 10 important genes where pipeline is worse than NGSEA:")
print(analyze_important_genes(rankings_df, rankings_df['Worse_than_NGSEA']))

print("Top 10 important genes where pipeline is better than ABS-PROP:")
print(analyze_important_genes(rankings_df, rankings_df['Better_than_ABS_PROP']))

print("Top 10 important genes where pipeline is worse than ABS-PROP:")
print(analyze_important_genes(rankings_df, rankings_df['Worse_than_ABS_PROP']))

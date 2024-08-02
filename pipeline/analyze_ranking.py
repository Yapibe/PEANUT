import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon, pearsonr

# Load the rankings summary
file_path = 'Outputs/NGSEA/Summary/H_sapiens/rankings_summary_H_sapiens.xlsx'
rankings_df = pd.read_excel(file_path)

# Function to calculate the changes in rank
def calculate_rank_changes(df, method1, method2):
    df[f'{method1}_vs_{method2}'] = df[method2] - df[method1]
    return df

# Calculate rank changes between methods
rankings_df = calculate_rank_changes(rankings_df, 'ABS-PROP', 'GSEA')
rankings_df = calculate_rank_changes(rankings_df, 'ABS-PROP', 'GSE')
rankings_df = calculate_rank_changes(rankings_df, 'ABS-PROP', 'PROP')

# Save the processed data to CSV
rankings_df.to_csv('Outputs/GSE/processed_rankings.csv', index=False)

# Perform Wilcoxon signed-rank test
w_stat_gsea, p_val_gsea = wilcoxon(rankings_df['ABS-PROP'], rankings_df['GSEA'])
w_stat_ngsea, p_val_ngsea = wilcoxon(rankings_df['ABS-PROP'], rankings_df['GSE'])
w_stat_prop, p_val_prop = wilcoxon(rankings_df['ABS-PROP'], rankings_df['PROP'])

# Generate rank distribution plot with significance
plt.figure(figsize=(10, 6))
sns.boxplot(data=rankings_df[['ABS-PROP', 'GSEA', 'GSE', 'PROP']], palette="Set2")
plt.title('Rank Distribution of Matched KEGG Pathway Terms')
plt.ylabel('Rank')
plt.xlabel('Method')

# Add significance asterisks
def add_significance_annotation(ax, p_value, x1, x2, y, h, col):
    if p_value < 0.05:
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        ax.text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color=col)

ax = plt.gca()
y_max = rankings_df[['ABS-PROP', 'GSEA', 'GSE', 'PROP']].values.max()
h = 2000  # height of the annotation
col = 'k'
add_significance_annotation(ax, p_val_gsea, 0, 1, y_max, h, col)
add_significance_annotation(ax, p_val_ngsea, 0, 2, y_max + h, h, col)
add_significance_annotation(ax, p_val_prop, 0, 3, y_max + 2*h, h, col)

plt.savefig('Outputs/GSE/rank_distribution.png')
plt.show()

# Function to compute PCC
def compute_pcc(df, method1, method2):
    pcc_values = []
    for pathway in df['Pathway'].unique():
        subset = df[df['Pathway'] == pathway]
        if len(subset) > 1:  # Ensure there are at least 2 points to calculate PCC
            pcc, _ = pearsonr(subset[method1], subset[method2])
            pcc_values.append(pcc)
    return pcc_values

# Compute PCC for same diseases and different diseases
pcc_abs_prop_vs_gsea = compute_pcc(rankings_df, 'ABS-PROP', 'GSEA')
pcc_abs_prop_vs_ngsea = compute_pcc(rankings_df, 'ABS-PROP', 'GSE')

# Generate PCC distribution plot
plt.figure(figsize=(10, 6))
sns.boxplot(data=[pcc_abs_prop_vs_gsea, pcc_abs_prop_vs_ngsea], palette="Set2")
plt.title('PCC Distribution of NES')
plt.ylabel('PCC')
plt.xticks([0, 1], ['ABS-PROP vs GSEA', 'ABS-PROP vs GSE'])
plt.savefig('Outputs/GSE/pcc_distribution.png')
plt.show()

# Remove "KEGG_" prefix from Pathway names and underscores
rankings_df['Pathway'] = rankings_df['Pathway'].str.replace('KEGG_', '').str.replace('_', ' ')

# Create combined labels for y-axis
rankings_df['Label'] = rankings_df['Pathway'] + ' (' + rankings_df['Dataset'] + ')'

# Sort the data by Pathway names to group similar pathways together
rankings_df = rankings_df.sort_values(by=['Pathway', 'Dataset'])

# Plotting the comparison
plt.figure(figsize=(12, 8))

# Plot ranks for ABS-PROP and GSEA
for i, row in rankings_df.iterrows():
    color = 'orange' if row['ABS-PROP'] < row['GSEA'] else 'black'
    plt.scatter(row['ABS-PROP'], i, color='black', label='ABS-PROP' if i == 0 else "", marker='o', edgecolor='black')
    plt.scatter(row['GSEA'], i, color='white', edgecolor='black', label='GSEA' if i == 0 else "", marker='o')

# Custom y-ticks with colored labels
y_ticks_labels = rankings_df['Label'].tolist()
y_ticks_colors = ['orange' if row['ABS-PROP'] < row['GSEA'] else 'black' for _, row in rankings_df.iterrows()]

plt.yticks(range(len(rankings_df)), labels=y_ticks_labels, fontsize=8, color='black')
ax = plt.gca()
for tick_label, tick_color in zip(ax.get_yticklabels(), y_ticks_colors):
    tick_label.set_color(tick_color)

plt.xlabel('Rank')
plt.ylabel('Gene expression data set (GSE ID)')
plt.title('Rank Comparison of Matched Pathways between ABS-PROP and GSEA')
plt.legend(loc='upper right')
ax.grid(False)  # Disable grid lines
plt.tight_layout()

# Save the plot
plt.savefig('Outputs/GSE/rank_comparison_ABS_PROP_GSEA.png')
plt.show()

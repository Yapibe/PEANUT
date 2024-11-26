import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon, pearsonr, ranksums
import numpy as np

# Define configurations
data_dir = 'pipeline/Outputs/NGSEA/Summary'
network = 'H_sapiens'
alphas = [0.2]
methods = ['PEANUT', 'ABS SCORE', 'GSEA', 'NGSEA']

# Define specific method pairs to compare and plot significance bars
methods_to_compare = [
    ('GSEA', 'NGSEA'),
    ('GSEA', 'PEANUT'),
    ('NGSEA', 'PEANUT'),
    ('PEANUT', 'ABS SCORE')
]
run_type = "filtered"
# Output directory for saving plots
output_plot_dir = "pipeline/Outputs/Plots"
os.makedirs(output_plot_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Map datasets to diseases (you need to provide this mapping)
disease_mapping = {
    'GSE1297': 'Alzheimer\'s',
    'GSE14762': 'Renal Cancer',
    'GSE15471': 'Pancreatic Cancer',
    'GSE16515': 'Pancreatic Cancer',
    'GSE18842': 'Non-Small Cell Lung Cancer',
    'GSE19188': 'Non-Small Cell Lung Cancer',
    'GSE19728': 'Glioma',
    'GSE20153': 'Parkinson\'s',
    'GSE20291': 'Parkinson\'s',
    'GSE21354': 'Glioma',
    'GSE3467': 'Thyroid Cancer',
    'GSE3585': 'Dilated Cardiomyopathy',
    'GSE3678': 'Thyroid Cancer',
    'GSE4107': 'Colorectal Cancer',
    'GSE5281EC': 'Alzheimer\'s',
    'GSE5281HIP': 'Alzheimer\'s',
    'GSE5281VCX': 'Alzheimer\'s',
    'GSE6956AA': 'Prostate Cancer',
    'GSE6956C': 'Prostate Cancer',
    'GSE781': 'Renal Cancer',
    'GSE8671': 'Colorectal Cancer',
    'GSE8762': 'Huntington\'s Disease',
    'GSE9348': 'Colorectal Cancer',
    'GSE9476': 'Acute Myeloid Leukemia'
}

def compute_pcc_for_diseases(df):
    """
    Compute PCC values for both same disease and different disease comparisons.

    Args:
        df (pd.DataFrame): DataFrame containing dataset, pathway, and rank information.

    Returns:
        pd.DataFrame: DataFrame with PCC values and comparison type ("Same Disease" or "Different Disease").
    """
    # Add disease column
    df['Disease'] = df['Dataset'].map(disease_mapping)

    # Prepare results
    pcc_results = []

    # Group by pathway
    grouped = df.groupby('Pathway')

    for pathway, group in grouped:
        datasets = group['Dataset'].unique()
        ranks = group[['Dataset', 'Rank PEANUT']].set_index('Dataset')['Rank PEANUT']

        # Same Disease Comparisons
        for disease, disease_group in group.groupby('Disease'):
            disease_datasets = disease_group['Dataset'].unique()
            if len(disease_datasets) > 1:  # Only compare within diseases with multiple datasets
                for i, dataset1 in enumerate(disease_datasets):
                    for j, dataset2 in enumerate(disease_datasets):
                        if i >= j:
                            continue
                        try:
                            # Ensure both datasets exist in the index
                            if dataset1 not in ranks.index or dataset2 not in ranks.index:
                                print(f"Skipping comparison: {dataset1} or {dataset2} not in ranks index.")
                                continue

                            # Extract ranks
                            ranks1 = ranks.loc[dataset1]
                            ranks2 = ranks.loc[dataset2]

                            # Ensure sufficient data for correlation
                            if len(ranks1) < 2 or len(ranks2) < 2:
                                print(f"Skipping pair {dataset1} and {dataset2}: insufficient data for correlation.")
                                continue

                            # Compute Pearson Correlation Coefficient (PCC)
                            pcc, _ = pearsonr(ranks1, ranks2)
                            pcc_results.append({
                                'Pathway': pathway,
                                'Dataset1': dataset1,
                                'Dataset2': dataset2,
                                'Disease1': disease,
                                'Disease2': disease,
                                'PCC': pcc,
                                'Comparison': 'Same Disease'
                            })
                        except Exception as e:
                            print(f"Error comparing {dataset1} and {dataset2} in pathway {pathway}: {e}")

        # Different Disease Comparisons
        diseases = group['Disease'].unique()
        for i, disease1 in enumerate(diseases):
            for j, disease2 in enumerate(diseases):
                if i >= j or disease1 == disease2:
                    continue
                disease1_datasets = group[group['Disease'] == disease1]['Dataset'].unique()
                disease2_datasets = group[group['Disease'] == disease2]['Dataset'].unique()
                for dataset1 in disease1_datasets:
                    for dataset2 in disease2_datasets:
                        try:
                            # Ensure both datasets exist in the index
                            if dataset1 not in ranks.index or dataset2 not in ranks.index:
                                print(f"Skipping comparison: {dataset1} or {dataset2} not in ranks index.")
                                continue

                            # Extract ranks
                            ranks1 = ranks.loc[dataset1]
                            ranks2 = ranks.loc[dataset2]

                            # Ensure sufficient data for correlation
                            if len(ranks1) < 2 or len(ranks2) < 2:
                                print(f"Skipping pair {dataset1} and {dataset2}: insufficient data for correlation.")
                                continue

                            # Compute Pearson Correlation Coefficient (PCC)
                            pcc, _ = pearsonr(ranks1, ranks2)
                            pcc_results.append({
                                'Pathway': pathway,
                                'Dataset1': dataset1,
                                'Dataset2': dataset2,
                                'Disease1': disease1,
                                'Disease2': disease2,
                                'PCC': pcc,
                                'Comparison': 'Different Disease'
                            })
                        except Exception as e:
                            print(f"Error comparing {dataset1} and {dataset2} in pathway {pathway}: {e}")

    return pd.DataFrame(pcc_results)



def plot_pcc_results(pcc_df, output_plot_dir):
    """
    Plot the PCC results for same and different disease comparisons.

    Args:
        pcc_df (pd.DataFrame): DataFrame containing PCC results.
        output_plot_dir (str): Directory to save the plot.
    """
    # Boxplot for PCC values
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='Comparison', y='PCC', data=pcc_df, palette='Set2')
    plt.title('PCC Distribution: Same Disease vs Different Disease')
    plt.xlabel('')
    plt.ylabel('Pearson Correlation Coefficient (PCC)')

    # Perform statistical test (Wilcoxon rank-sum test)
    same_disease_pcc = pcc_df[pcc_df['Comparison'] == 'Same Disease']['PCC']
    different_disease_pcc = pcc_df[pcc_df['Comparison'] == 'Different Disease']['PCC']
    stat, p_value = ranksums(same_disease_pcc, different_disease_pcc)

    # Add significance marker
    significance_marker = f'p = {p_value:.2e}'
    plt.text(0.5, 0.9, significance_marker, horizontalalignment='center',
             verticalalignment='center', transform=plt.gca().transAxes, fontsize=12)

    # Save plot
    plot_filepath = os.path.join(output_plot_dir, 'PCC_Comparison.png')
    plt.tight_layout()
    plt.savefig(plot_filepath, dpi=300)
    print(f"PCC plot saved to: {plot_filepath}")
    plt.show()



def plot_comparative_ranks(df, alpha, output_plot_dir):
    """
    Plot comparative ranks for PEANUT vs GSEA, grouped by similar pathways with consistent group shading.
    
    Args:
        df (pd.DataFrame): DataFrame containing rank information.
        alpha (float): Alpha value for the plot title.
        output_plot_dir (str): Directory to save the plot.
    """
    # Remove "KEGG" prefix and underscores from pathway names
    df['Pathway'] = df['Pathway'].str.replace('KEGG', '').str.replace('_', ' ').str.strip()

    # Combine pathway and dataset for y-axis labels
    df['Y Label'] = df['Pathway'] + ', ' + df['Dataset']

    # Sort pathways by similarity (custom logic to group similar pathways)
    df.sort_values(by=['Pathway', 'Dataset'], inplace=True)

    # Assign groups based on similarity
    unique_groups = df['Pathway'].str.extract(r'(^[^\s]+)')[0]  # Extract group identifier (e.g., first word)
    df['Group'] = unique_groups
    df['Group ID'] = df['Group'].factorize()[0]  # Assign numeric group IDs for alternating shades

    # Prepare data for scatter plot
    comparison_df = df[['Y Label', 'Rank PEANUT', 'Rank GSEA']]
    plot_df = comparison_df.melt(id_vars=['Y Label'], var_name='Method', value_name='Rank')

    # Create the plot
    plt.figure(figsize=(8, len(df) * 0.4))  # Adjusted for compact spacing
    ax = plt.gca()

    # Add consistent background shades for each group (use only valid data rows)
    group_shades = ['#f2f2f2', '#ffffff']  # Two alternating colors
    for idx, (group_id, color) in enumerate(zip(df['Group ID'], group_shades * len(df))):
        ax.axhspan(idx - 0.5, idx + 0.5, facecolor=color, alpha=0.7)  # Limit shading to valid rows only

    # Scatter plot for the ranks
    sns.scatterplot(data=plot_df, x='Rank', y='Y Label', hue='Method', style='Method', s=100, palette=['blue', 'orange'])

    # Customize plot
    ax.set_title(f'', fontsize=14)  # Slightly smaller title
    ax.set_xlabel('Rank', fontsize=12)  # Smaller font for the x-axis
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=10)  # Smaller font for y-axis labels
    ax.tick_params(axis='x', labelsize=10)
    plt.legend(title='Method', fontsize=10)

    # Adjust layout for compactness
    plt.tight_layout()

    # Save the plot
    plot_filename = f'Comparative_Ranks_Alpha_{alpha}_{run_type}.png'
    plot_filepath = os.path.join(output_plot_dir, plot_filename)
    plt.savefig(plot_filepath, format='png', dpi=300)
    print(f"Plot saved to: {plot_filepath}")

def load_data(alpha, data_dir):
    file_name = f'rankings_summary_{run_type}.xlsx'
    file_path = os.path.join(data_dir, file_name)
    if os.path.exists(file_path):
        df = pd.read_excel(file_path)
        
        # Remove the second row
        df = df.drop(df.index[0])
        
        # Ensure the new expected columns based on updated file structure
        expected_columns = ['Dataset'] + [f'Rank {method}' for method in methods] + [f'Significant {method}' for method in methods]
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            print(f"Missing expected columns: {missing_columns}")
            return None
        
        # Convert rank columns to numeric
        rank_columns = [f'Rank {method}' for method in methods]
        for col in rank_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Drop rows with NaN in rank columns
        df.dropna(subset=rank_columns, inplace=True)
        
        df['Alpha'] = alpha
        return df
    else:
        print(f"File not found: {file_path}")
        return None

def significance_marker(p_value):
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'  # Not significant

def add_significance_bars(ax, method_pairs, p_values, y_max):
    """Draw significance bars between specified pairs of methods."""
    bar_height = y_max * 0.02  # Adjust the bar height relative to y_max
    y_increment = y_max * 0.05  # Extra space for each bar to make the markers visible
    for i, (pair_indices, p_value) in enumerate(zip(method_pairs, p_values)):
        if p_value < 0.05:  # Draw only if p-value is significant
            # Get x positions for the pair of methods
            x1, x2 = pair_indices
            x1_pos, x2_pos = x1, x2
            y = y_max + y_increment * (i + 1)  # Increment y position for each new bar
            # Draw the horizontal bar and add the significance marker
            ax.plot([x1_pos, x1_pos, x2_pos, x2_pos], [y, y + bar_height, y + bar_height, y], color='black')
            ax.text((x1_pos + x2_pos) * 0.5, y + 1.5 * bar_height, significance_marker(p_value),
                    ha='center', va='bottom', fontsize=14, color='red')

def analyze_data(df, alpha, output_plot_dir):
    # Extract ranks for all methods (focusing on GSEA, NGSEA, and PEANUT)
    ranks = {method: df[f'Rank {method}'] for method in methods}
    
    # Perform Wilcoxon signed-rank tests between specified method pairs
    p_values = []
    test_results = []
    mean_rank_diffs = []
    
    for method1, method2 in methods_to_compare:
        rank1 = ranks[method1]
        rank2 = ranks[method2]
        
        # Compute the rank differences
        rank_diff = rank1 - rank2
        
        # Compute mean rank difference to determine which method ranks better
        mean_diff = rank_diff.mean()
        mean_rank_diffs.append(mean_diff)  # Store the mean difference
        
        # Perform the Wilcoxon signed-rank test
        w_result = wilcoxon(rank1, rank2)
        
        # Determine which method ranks better based on sign of the differences
        better_method = method1 if mean_diff < 0 else method2
        
        p_values.append(w_result.pvalue)
        test_results.append((method1, method2, w_result.statistic, w_result.pvalue, better_method))
    
    # Print test results with the signed-rank statistic, better method, and mean rank difference
    print(f"\nAlpha {alpha} - Wilcoxon signed-rank test results:")
    for i, (method1, method2, stat, p_value, better_method) in enumerate(test_results):
        print(f"{method1} vs {method2}: Statistic={stat}, p-value={p_value}, Better Method={better_method}, Mean Rank Difference={mean_rank_diffs[i]}")
    
    # Prepare data for plotting
    plot_df = df.melt(id_vars=['Dataset'], value_vars=[f'Rank {method}' for method in methods],
                      var_name='Method', value_name='Rank')
    plot_df['Method'] = plot_df['Method'].str.replace('Rank ', '')
    
    # Plot boxplots
    plt.figure(figsize=(14, 8))
    ax = sns.boxplot(x='Method', y='Rank', data=plot_df, palette='Set2')
    ax.set_title(f'', fontsize=16)
    ax.set_ylabel('Rank')
    ax.set_xlabel('')
    ax.set_xticklabels(methods)
    plt.gca().invert_yaxis()  # Lower ranks at the top

    # Get y-axis maximum to position bars and markers
    y_max = plot_df['Rank'].max()
    
    # Adjust method indices for x positions
    method_indices = {method: idx for idx, method in enumerate(methods)}
    # Adjust method pairs to indices for plotting
    method_pairs_indices = [(method_indices[method1], method_indices[method2]) for method1, method2 in methods_to_compare]

    # Adjust p_values to match the selected method pairs
    p_values_to_plot = [p_values[methods_to_compare.index(pair)] for pair in methods_to_compare]

    # Add significance bars and asterisks
    add_significance_bars(ax, method_pairs_indices, p_values_to_plot, y_max)

    # Save the plot as a PNG file
    plot_filename = f'Rank_Comparison_Alpha_{alpha}_{run_type}.png'
    plot_filepath = os.path.join(output_plot_dir, plot_filename)
    plt.savefig(plot_filepath, format='png', dpi=300)
    print(f"Plot saved to: {plot_filepath}")

def main():
    for alpha in alphas:
        # Load data for the current alpha
        df = load_data(alpha, data_dir)
        
        if df is not None and not df.empty:
            # Analyze and plot existing data
            analyze_data(df, alpha, output_plot_dir)
            plot_comparative_ranks(df, alpha, output_plot_dir)
            
            # # Compute PCC between datasets grouped by pathway
            # pcc_df = compute_pcc_for_diseases(df)
            
            # # Analyze and plot PCC results
            # plot_pcc_results(pcc_df, output_plot_dir)
        else:
            print(f"No data available for alpha {alpha}")

if __name__ == "__main__":
    main()


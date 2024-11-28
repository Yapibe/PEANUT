import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon, pearsonr, ranksums
import numpy as np
from matplotlib.ticker import MultipleLocator  # Import MultipleLocator

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
    plt.figure(figsize=(6, len(df) * 0.2))  # Adjusted for smaller size
    ax = plt.gca()

    # Add consistent background shades for each group
    group_shades = ['#f2f2f2', '#ffffff']  # Two alternating colors
    group_ids = df['Group ID'].unique()
    for idx, group_id in enumerate(group_ids):
        color = group_shades[idx % len(group_shades)]
        group_df = df[df['Group ID'] == group_id]
        y_positions = group_df['Y Label']
        y_start = group_df.index[0] - 0.5
        y_end = group_df.index[-1] + 0.5
        ax.axhspan(y_start, y_end, facecolor=color, alpha=0.7)

    # Scatter plot for the ranks
    sns.scatterplot(data=plot_df, x='Rank', y='Y Label', hue='Method', style='Method', s=100, palette=['blue', 'orange'])

    # Customize plot
    ax.set_title(f'', fontsize=12)  # Reduced font size
    ax.set_xlabel('Rank', fontsize=10)  # Reduced font size
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=8)  # Reduced font size
    ax.tick_params(axis='x', labelsize=8)
    plt.legend(title='Method', fontsize=8, title_fontsize=10)

    # Set x-axis ticks every 50 ranks
    min_rank = plot_df['Rank'].min()
    max_rank = plot_df['Rank'].max()
    x_ticks = np.arange(0, max_rank + 50, 50)
    ax.set_xticks(x_ticks)

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
    plt.figure(figsize=(8, 6))  # Adjusted figure size as needed
    ax = sns.boxplot(x='Method', y='Rank', data=plot_df, palette='Set2')
    ax.set_title(f'', fontsize=16)
    ax.set_ylabel('Rank', fontsize=14)  # Font size for y-axis label
    ax.set_xlabel('')
    ax.set_xticklabels(methods, fontsize=12)  # Font size for x-tick labels
    ax.tick_params(axis='y', labelsize=12)  # Font size for y-tick labels

    # Invert y-axis to have lower ranks at the top
    ax.invert_yaxis()

    # Add more y-axis ticks every 10 units
    ax.yaxis.set_major_locator(MultipleLocator(10))

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
    plt.savefig(plot_filepath, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to: {plot_filepath}")

def main():
    for alpha in alphas:
        # Load data for the current alpha
        df = load_data(alpha, data_dir)
        
        if df is not None and not df.empty:
            # Analyze and plot existing data
            analyze_data(df, alpha, output_plot_dir)
            plot_comparative_ranks(df, alpha, output_plot_dir)
        else:
            print(f"No data available for alpha {alpha}")

if __name__ == "__main__":
    main()

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon
from collections import Counter

# Define configurations
data_dir = 'pipeline/Outputs/NGSEA/Summary'
network = 'H_sapiens'
alphas = [0.2]
methods = ['GSEA', 'NGSEA', 'ABS_PROP', 'PROP', 'MW']

# Define specific method pairs to compare and plot significance bars
methods_to_compare = [
    ('GSEA', 'ABS_PROP'),
    ('GSEA', 'MW'),
    ('ABS_PROP', 'MW')
]

# Output directory for saving plots
output_plot_dir = "pipeline/Outputs/Plots"
os.makedirs(output_plot_dir, exist_ok=True)  # Create the directory if it doesn't exist

def load_data(alpha, data_dir, network):
    file_name = f'combined_rankings_summary_all_methods.xlsx'
    file_path = os.path.join(data_dir, file_name)
    if os.path.exists(file_path):
        df = pd.read_excel(file_path)
        
        # Remove any rows with 'Average' or 'percent' in 'Dataset' column
        df = df[~df['Dataset'].str.contains('average|percent', case=False, na=False)]
        
        # Ensure relevant columns are present
        expected_columns = ['Dataset'] + [f'{method} Rank' for method in methods]
        missing_columns = [col for col in expected_columns if col not in df.columns]
        if missing_columns:
            print(f"Missing expected columns: {missing_columns}")
            return None
        
        # Convert rank columns to numeric
        rank_columns = [f'{method} Rank' for method in methods]
        for col in rank_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Drop rows with NaN in rank columns
        df.dropna(subset=rank_columns, inplace=True)
        
        df['Alpha'] = alpha
        return df
    else:
        print(f"File not found: {file_path}")
        return None


# Annotate the significance level with asterisks
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
    # Extract ranks for all methods
    ranks = {method: df[f'{method} Rank'] for method in methods}
    
    # Perform Wilcoxon signed-rank tests between all pairs of methods
    from itertools import combinations
    method_pairs = list(combinations(methods, 2))
    p_values = []
    test_results = []
    
    for method1, method2 in method_pairs:
        rank1 = ranks[method1]
        rank2 = ranks[method2]
        # Check if there are at least one non-zero difference
        if (rank1 != rank2).any():
            w_result = wilcoxon(rank1, rank2)
            p_values.append(w_result.pvalue)
            test_results.append((method1, method2, w_result))
        else:
            # If all ranks are equal, p-value is 1
            p_values.append(1)
            test_results.append((method1, method2, 'Ranks are identical'))
    
    # Print test results
    print(f"\nAlpha {alpha} - Wilcoxon signed-rank test results:")
    for i, (method1, method2, result) in enumerate(test_results):
        if isinstance(result, str):
            print(f"{method1} vs {method2}: {result}")
        else:
            print(f"{method1} vs {method2}: Statistic={result.statistic}, p-value={result.pvalue}")
    
    # Prepare data for plotting
    plot_df = df.melt(id_vars=['Dataset'], value_vars=[f'{method} Rank' for method in methods],
                      var_name='Method', value_name='Rank')
    plot_df['Method'] = plot_df['Method'].str.replace(' Rank', '')
    
    # Plot boxplots
    plt.figure(figsize=(14, 8))
    ax = sns.boxplot(x='Method', y='Rank', data=plot_df, palette='Set2')
    ax.set_title(f'Rank Comparison for Alpha={alpha}', fontsize=16)
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
    p_values_to_plot = [p_values[method_pairs.index(pair)] for pair in methods_to_compare]

    # Add significance bars and asterisks
    add_significance_bars(ax, method_pairs_indices, p_values_to_plot, y_max)

    # Save the plot as a PNG file
    plot_filename = f'Rank_Comparison_Alpha_{alpha}.png'
    plot_filepath = os.path.join(output_plot_dir, plot_filename)
    plt.savefig(plot_filepath, format='png', dpi=300)
    print(f"Plot saved to: {plot_filepath}")

    plt.tight_layout()
    plt.show()


def main():
    for alpha in alphas:
        df = load_data(alpha, data_dir, network)
        if df is not None and not df.empty:
            analyze_data(df, alpha, output_plot_dir)
        else:
            print(f"No data available for alpha {alpha}")

if __name__ == "__main__":
    main()
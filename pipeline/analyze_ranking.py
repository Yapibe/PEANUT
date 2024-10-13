import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon
from collections import Counter

# Define configurations
data_dir = 'Outputs/NGSEA/Summary'
network = 'H_sapiens'
alphas = [0.1, 0.2]
methods = ['GSEA', 'NGSEA', 'ABS_PROP']

# Output directory for saving plots
output_plot_dir = "Outputs/Plots"
os.makedirs(output_plot_dir, exist_ok=True)  # Create the directory if it doesn't exist

# Load data for each alpha value
def load_data(alpha, data_dir, network):
    file_name = f'rankings_summary_{network}_kegg_alpha_{alpha}.xlsx'
    file_path = os.path.join(data_dir, network, 'kegg', f'alpha {alpha}', file_name)
    if os.path.exists(file_path):
        df = pd.read_excel(file_path)
        
        # Remove any rows with 'Average' or 'percent' in 'Dataset' column
        df = df[~df['Dataset'].str.contains('average|percent', case=False, na=False)]
        
        # Ensure relevant columns are present
        expected_columns = ['Dataset', 'ABS_PROP Rank', 'GSEA Rank', 'NGSEA Rank']
        for col in expected_columns:
            if col not in df.columns:
                print(f"Missing expected column: {col}")
                return None
        
        # Convert rank columns to numeric
        rank_columns = ['ABS_PROP Rank', 'GSEA Rank', 'NGSEA Rank']
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

# Function to draw significance bars on the boxplot
def add_significance_bars(ax, method_pairs, p_values, y_max):
    """Draw significance bars between specified pairs of methods."""
    bar_height = 5  # Distance between bars
    y_increment = 20  # Extra space for each bar to make the markers visible
    for i, (pair, p_value) in enumerate(zip(method_pairs, p_values)):
        if p_value < 0.05:  # Draw only if p-value is significant
            # Get x positions for the pair of methods
            x1, x2 = pair
            # Adjust x positions using ax.get_xticks() to find the actual boxplot positions
            x1_pos, x2_pos = ax.get_xticks()[x1], ax.get_xticks()[x2]
            y = y_max + y_increment * (i + 1)  # Increment y position for each new bar
            # Draw the horizontal bar and add the significance marker
            ax.plot([x1_pos, x1_pos, x2_pos, x2_pos], [y, y + bar_height, y + bar_height, y], color='black')
            ax.text((x1_pos + x2_pos) * 0.5, y + 1.5 * bar_height, significance_marker(p_value),
                    ha='center', va='bottom', fontsize=14, color='red')


# Analyze data and plot results with significance markers
# Analyze data and plot results with significance markers
def analyze_data(df, alpha, output_plot_dir):
    # Extract ranks
    abs_prop_ranks = df['ABS_PROP Rank']
    gsea_ranks = df['GSEA Rank']
    ngsea_ranks = df['NGSEA Rank']
    
    # Perform Wilcoxon signed-rank tests
    w_abs_gsea = wilcoxon(abs_prop_ranks, gsea_ranks)
    w_abs_ngsea = wilcoxon(abs_prop_ranks, ngsea_ranks)
    w_gsea_ngsea = wilcoxon(gsea_ranks, ngsea_ranks)
    
    # Corrected method pairs: (0 -> GSEA), (1 -> NGSEA), (2 -> ABS_PROP)
    method_pairs = [(0, 2), (1, 2), (0, 1)]  # Pairs of method indices (x positions on the plot)
    p_values = [w_abs_gsea.pvalue, w_abs_ngsea.pvalue, w_gsea_ngsea.pvalue]

    # Print test results
    print(f"\nAlpha {alpha} - Wilcoxon signed-rank test results:")
    print(f"ABS_PROP vs GSEA: Statistic={w_abs_gsea.statistic}, p-value={w_abs_gsea.pvalue}")
    print(f"ABS_PROP vs NGSEA: Statistic={w_abs_ngsea.statistic}, p-value={w_abs_ngsea.pvalue}")
    print(f"GSEA vs NGSEA: Statistic={w_gsea_ngsea.statistic}, p-value={w_gsea_ngsea.pvalue}")
    
    # Prepare data for plotting
    plot_df = df.melt(id_vars=['Dataset'], value_vars=['GSEA Rank', 'NGSEA Rank', 'ABS_PROP Rank'],
                      var_name='Method', value_name='Rank')
    plot_df['Method'] = plot_df['Method'].str.replace(' Rank', '')

    # Plot boxplots
    plt.figure(figsize=(10, 7))
    ax = sns.boxplot(x='Method', y='Rank', data=plot_df, palette='Set2')
    ax.set_title(f'Rank Comparison for Alpha={alpha}', fontsize=16)
    ax.set_ylabel('Rank')
    ax.set_xlabel('')
    ax.set_xticklabels(['GSEA', 'NGSEA', 'PROP'])
    plt.gca().invert_yaxis()  # Lower ranks at the top

    # Get y-axis maximum to position bars and markers
    y_max = plot_df['Rank'].max() + 50  # Offset to ensure markers are visible
    
    # Add significance bars and asterisks
    add_significance_bars(ax, method_pairs, p_values, y_max)

    # Save the plot as a PNG file
    plot_filename = f'Rank_Comparison_Alpha_{alpha}.png'
    plot_filepath = os.path.join(output_plot_dir, plot_filename)
    plt.savefig(plot_filepath, format='png', dpi=300)
    print(f"Plot saved to: {plot_filepath}")

    plt.tight_layout()
    plt.show()


# Main function to execute the analysis
def main():
    for alpha in alphas:
        df = load_data(alpha, data_dir, network)
        if df is not None and not df.empty:
            analyze_data(df, alpha, output_plot_dir)
        else:
            print(f"No data available for alpha {alpha}")

if __name__ == "__main__":
    main()

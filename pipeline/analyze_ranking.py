import os
import pandas as pd
from scipy.stats import pearsonr, ttest_rel, wilcoxon
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import numpy as np

# Load the files
h_sapiens_file_path = 'Outputs/NGSEA/Summary/H_sapiens/kegg/alpha 0.2/rankings_summary_H_sapiens_kegg_alpha_0.2.xlsx'
string_file_path = 'Outputs/NGSEA/Summary/String/kegg/alpha 0.2/rankings_summary_String_kegg_alpha_0.2.xlsx'
Hnet_file_path = 'Outputs/NGSEA/Summary/HumanNet/kegg/alpha 0.2/rankings_summary_HumanNet_kegg_alpha_0.2.xlsx'

df_h_sapiens = pd.read_excel(h_sapiens_file_path)
df_string = pd.read_excel(string_file_path)
df_human_net = pd.read_excel(Hnet_file_path)

# Extract the ranks for ABS_PROP and GSEA
abs_prop_ranks_hs = df_h_sapiens['ABS_PROP Rank']
gsea_ranks_hs = df_h_sapiens['GSEA Rank']

abs_prop_ranks_string = df_string['ABS_PROP Rank']
gsea_ranks_string = df_string['GSEA Rank']

abs_prop_ranks_human_net = df_human_net['ABS_PROP Rank']
gsea_ranks_human_net = df_human_net['GSEA Rank']

# Perform the Wilcoxon signed-rank test
result_hs = wilcoxon(abs_prop_ranks_hs, gsea_ranks_hs)
result_string = wilcoxon(abs_prop_ranks_string, gsea_ranks_string)
result_human_net = wilcoxon(abs_prop_ranks_human_net, gsea_ranks_human_net)

# Output the results
print(f"H_sapiens file - Wilcoxon statistic: {result_hs.statistic}, p-value: {result_hs.pvalue}")
print(f"String file - Wilcoxon statistic: {result_string.statistic}, p-value: {result_string.pvalue}")
print(f"HumanNet file - Wilcoxon statistic: {result_human_net.statistic}, p-value: {result_human_net.pvalue}")



# Specify the directory containing the data
data_dir = 'Outputs/NGSEA/Summary'

# List of networks and alpha values
networks = ['String_', 'H_sapiens', 'HumanNet']
alphas = [0.1, 0.2]
methods = ['GSEA', 'NGSEA', 'PROP', 'ABS_PROP']

# Construct the file paths dynamically
file_paths = []
for network in networks:
    for alpha in alphas:
        file_name = f'rankings_summary_{network}_kegg_alpha_{alpha}.xlsx'
        file_path = os.path.join(data_dir, network, 'kegg', f'alpha {alpha}', file_name)
        file_paths.append(file_path)

# Load the datasets
dataframes = []
for file_path in file_paths:
    df = pd.read_excel(file_path)
    df = df[~df['Dataset'].str.contains('average|percent', case=False, na=False)]

    # Extract network and alpha from file name
    filename = os.path.basename(file_path)
    network = filename.split('_')[2]
    alpha = float(filename.split('_')[-1].replace('.xlsx', ''))
    df['Network'] = network
    df['Alpha'] = alpha

    # Correct network name if it's listed as 'H'
    df['Network'].replace({'H': 'H_sapiens'}, inplace=True)
    df['Network'].replace({'String': 'String_'}, inplace=True)

    dataframes.append(df)

# Combine all data into one DataFrame
combined_df = pd.concat(dataframes, ignore_index=True)

# Aggregate data by Pathway
pathway_df = combined_df.groupby('Pathway').agg({
    'Density': 'mean',
    'Num Genes': 'mean',
    'Avg Diameter': 'mean',
    'PROP Rank': 'mean',
    'ABS_PROP Rank': 'mean'
}).reset_index()

# Network metrics
network_metrics = {
    'Network': ['HumanNet', 'H_sapiens', 'String_'],
    'Num Nodes': [18459, 19950, 17991],
    'Num Edges': [977495, 899737, 712446],
    'Network Density': [0.005737883534057266, 0.004521489698480498, 0.004402460633689089],
    'Avg Degree': [105.90985427162902, 90.19919799498747, 79.2002668000667],
    'Connected Components': [3, 2, 4],
    'Largest Component Size': [18454, 19948, 17984],
    'Avg Clustering Coefficient': [0.22184569149353095, 0.13695278983250672, 0.18843538587424444]
}

# Convert to DataFrame
network_metrics_df = pd.DataFrame(network_metrics)

# Merge the network metrics with the combined dataset
combined_df = pd.merge(combined_df, network_metrics_df, on='Network', how='left')

# Group by 'Dataset' to get unique GSEA and NGSEA results (no need to group by network)
unique_gsea_ngsea_df = combined_df.groupby('Dataset').first().reset_index()

# Calculate Pearson correlations using the unique sets for GSEA and NGSEA, and the full data for PROP and ABS_PROP
correlation_results = {}

# Calculate correlation for GSEA and NGSEA across the different networks
for network in networks:
    # Filter the data for the current network
    network_df = combined_df[combined_df['Network'] == network]

    # Calculate correlations
    correlation_results[network] = {
        'GSEA': unique_gsea_ngsea_df[['GSEA Rank', 'GSEA Significant']].corr().iloc[0, 1],
        'NGSEA': unique_gsea_ngsea_df[['NGSEA Rank', 'NGSEA Significant']].corr().iloc[0, 1],
        'PROP': network_df[['PROP Rank', 'PROP Significant']].corr().iloc[0, 1],
        'ABS_PROP': network_df[['ABS_PROP Rank', 'ABS_PROP Significant']].corr().iloc[0, 1],
    }

print("Pearson Correlation Results by Network:", correlation_results)


# Paired t-tests for PROP and ABS_PROP ranks, comparing alpha=0.1 vs. alpha=0.2 within each network
def perform_paired_ttest(network, method):
    alpha_0_1 = combined_df[(combined_df['Network'] == network) & (combined_df['Alpha'] == 0.1)][f'{method} Rank']
    alpha_0_2 = combined_df[(combined_df['Network'] == network) & (combined_df['Alpha'] == 0.2)][f'{method} Rank']

    t_stat, p_value = ttest_rel(alpha_0_1, alpha_0_2)
    return t_stat, p_value


ttest_results = {}
for network in networks:
    ttest_results[network] = {}
    for method in ['PROP', 'ABS_PROP']:
        t_stat, p_value = perform_paired_ttest(network, method)
        ttest_results[network][method] = {'t_stat': t_stat, 'p_value': p_value}

print("Paired T-Test Results:", ttest_results)

# Updated factors to include network metrics
factors = ('Network + Alpha + Density + Q("Num Genes") + Q("Avg Diameter") + '
           'Q("Num Nodes") + Q("Num Edges") + Q("Network Density") + Q("Avg Degree") + '
           'Q("Connected Components") + Q("Largest Component Size") + Q("Avg Clustering Coefficient")')

anova_results = {}

for method in methods:
    formula = f'Q("{method} Rank") ~ {factors}'
    model = smf.ols(formula=formula, data=combined_df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_results[method] = anova_table

for method, result in anova_results.items():
    print(f"ANOVA Results for {method} Rank:\n", result)

# Multivariate regression including interaction effects between network and alpha
interaction_formula_prop = (
    'Q("PROP Rank") ~ Network * Alpha + Density + Q("Num Genes") + Q("Avg Diameter") + '
    'Q("Num Nodes") + Q("Num Edges") + Q("Network Density") + Q("Avg Degree") + '
    'Q("Connected Components") + Q("Largest Component Size") + Q("Avg Clustering Coefficient")'
)

interaction_formula_abs_prop = (
    'Q("ABS_PROP Rank") ~ Network * Alpha + Density + Q("Num Genes") + Q("Avg Diameter") + '
    'Q("Num Nodes") + Q("Num Edges") + Q("Network Density") + Q("Avg Degree") + '
    'Q("Connected Components") + Q("Largest Component Size") + Q("Avg Clustering Coefficient")'
)

prop_model = smf.ols(formula=interaction_formula_prop, data=combined_df).fit()
abs_prop_model = smf.ols(formula=interaction_formula_abs_prop, data=combined_df).fit()

print("PROP Multivariate Regression Summary:\n", prop_model.summary())
print("ABS_PROP Multivariate Regression Summary:\n", abs_prop_model.summary())


# Function to determine the best number of clusters using the elbow method and silhouette score
# Determine the best number of clusters using the elbow method and silhouette score
def determine_best_clusters(data):
    sse = []
    silhouette_scores = []
    k_values = range(2, 10)  # You can adjust the range of k values

    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)  # Explicitly set n_init to 10
        kmeans.fit(data)
        sse.append(kmeans.inertia_)  # Sum of squared distances to the nearest cluster center
        silhouette_scores.append(silhouette_score(data, kmeans.labels_))

    # Choose the best number of clusters based on silhouette score
    best_k = k_values[np.argmax(silhouette_scores)]
    print(f"Best number of clusters according to silhouette score: {best_k}")

    return best_k


# Perform PCA and K-means clustering for each individual dataset
for file_path in file_paths:
    df = pd.read_excel(file_path)
    df = df[~df['Dataset'].str.contains('average|percent', case=False, na=False)]

    # Standardize the data, excluding the ranks
    clustering_data = pathway_df[['Density', 'Num Genes', 'Avg Diameter']].dropna()
    scaler = StandardScaler()
    clustering_data_scaled = scaler.fit_transform(clustering_data)

    # Determine the best number of clusters for this dataset
    best_k = determine_best_clusters(clustering_data_scaled)

    # Perform K-means clustering using the best number of clusters
    kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
    pathway_df['Cluster'] = kmeans.fit_predict(clustering_data_scaled)

    # Print the pathways in each cluster
    print(f"\nClusters for pathways:")
    for cluster in range(best_k):
        pathways_in_cluster = pathway_df[pathway_df['Cluster'] == cluster]['Pathway']
        print(f"Cluster {cluster + 1}:")
        print(pathways_in_cluster.to_list())

    # Perform PCA for dimensionality reduction
    pca = PCA(n_components=2)
    pca_components = pca.fit_transform(clustering_data_scaled)
    pca_df = pd.DataFrame(data=pca_components, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = pathway_df['Cluster']
    pca_df['Pathway'] = pathway_df['Pathway']  # Add the pathway names to the PCA DataFrame
    pca_df['PROP Rank'] = pathway_df['PROP Rank']  # Keep track of the PROP Rank
    pca_df['ABS_PROP Rank'] = pathway_df['ABS_PROP Rank']  # Keep track of the ABS_PROP Rank

    # Visualize the PCA results with clusters and pathway names
    plt.figure(figsize=(10, 6))
    for cluster in pca_df['Cluster'].unique():
        subset = pca_df[pca_df['Cluster'] == cluster]
        plt.scatter(subset['PC1'], subset['PC2'], label=f'Cluster {cluster}')
        for i in range(subset.shape[0]):
            plt.text(subset['PC1'].iloc[i], subset['PC2'].iloc[i], subset['Pathway'].iloc[i], fontsize=8)
    plt.title(f'PCA of Pathways')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Display the principal components and the explained variance
    pca_explained_variance = pca.explained_variance_ratio_
    pca_components_df = pd.DataFrame(pca.components_, columns=clustering_data.columns,
                                     index=[f'PC{i + 1}' for i in range(pca.n_components_)])

    print(f"Explained Variance by Principal Components: {pca_explained_variance}")
    print("\nPrincipal Components:\n", pca_components_df)

    # Optionally: Correlate the PCA components with PROP and ABS_PROP Rank
    corr_prop = pca_df[['PC1', 'PC2']].corrwith(pca_df['PROP Rank'])
    corr_abs_prop = pca_df[['PC1', 'PC2']].corrwith(pca_df['ABS_PROP Rank'])

    print(f"Correlation of PCA Components with PROP Rank:\n", corr_prop)
    print(f"Correlation of PCA Components with ABS_PROP Rank:\n", corr_abs_prop)

# Random Forest model to predict PROP and ABS_PROP ranks
features = ['Density', 'Num Genes', 'Avg Diameter', 'Network', 'Alpha',
            'Num Nodes', 'Num Edges', 'Network Density', 'Avg Degree',
            'Connected Components', 'Largest Component Size', 'Avg Clustering Coefficient']

rf_results = {}

for method in ['PROP Rank', 'ABS_PROP Rank']:
    X = combined_df[features]
    X = pd.get_dummies(X, columns=['Network', 'Alpha'], drop_first=True)
    y = combined_df[method].dropna()

    rf = RandomForestRegressor(random_state=42)
    rf.fit(X, y)

    importance = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
    rf_results[method] = importance

for method, importance in rf_results.items():
    print(f"Random Forest Feature Importance for {method}:\n", importance)

# Visualizing interaction effects between network and pathway characteristics
# Plotting Density vs. PROP Rank by Network
plt.figure(figsize=(10, 6))
for network in networks:
    subset = combined_df[combined_df['Network'] == network]
    plt.scatter(subset['Density'], subset['PROP Rank'], label=network)
    plt.plot(np.unique(subset['Density']),
             np.poly1d(np.polyfit(subset['Density'], subset['PROP Rank'], 1))(np.unique(subset['Density'])))
plt.title('Density vs. PROP Rank by Network')
plt.xlabel('Density')
plt.ylabel('PROP Rank')
plt.legend()
plt.grid(True)
plt.show()

# Plotting Density vs. ABS_PROP Rank by Network
plt.figure(figsize=(10, 6))
for network in networks:
    subset = combined_df[combined_df['Network'] == network]
    plt.scatter(subset['Density'], subset['ABS_PROP Rank'], label=network)
    plt.plot(np.unique(subset['Density']),
             np.poly1d(np.polyfit(subset['Density'], subset['ABS_PROP Rank'], 1))(np.unique(subset['Density'])))
plt.title('Density vs. ABS_PROP Rank by Network')
plt.xlabel('Density')
plt.ylabel('ABS_PROP Rank')
plt.legend()
plt.grid(True)
plt.show()

# Pathway attribute analysis for both PROP and ABS_PROP ranks
# Define the pathway attributes to analyze
pathway_attributes = ['Density', 'Num Genes', 'Avg Diameter']
attribute_results = {}

for method in ['PROP Rank', 'ABS_PROP Rank']:
    formula = f'Q("{method}") ~ ' + ' + '.join([f'Q("{attr}")' for attr in pathway_attributes])
    model = smf.ols(formula=formula, data=combined_df).fit()

    # Store the results
    attribute_results[method] = model.summary()

    # Print the summary
    print(f"\nOLS Regression Results for {method}:\n")
    print(attribute_results[method])

# Additionally, visualize the coefficients for each attribute
for method in ['PROP Rank', 'ABS_PROP Rank']:
    model = smf.ols(f'Q("{method}") ~ ' + ' + '.join([f'Q("{attr}")' for attr in pathway_attributes]),
                    data=combined_df).fit()
    coefficients = model.params.drop('Intercept')

    # Increase the figure size and adjust layout to prevent cutoff
    plt.figure(figsize=(12, 8))
    coefficients.plot(kind='barh')
    plt.title(f'Impact of Pathway Attributes on {method}')
    plt.xlabel('Coefficient Value')
    plt.ylabel('Pathway Attributes')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

import pandas as pd
import numpy as np
import json
#
# # Load the gene name dictionary from a JSON file
# gene_dict_path = 'Data/H_sapiens/gene_names/gene_info.json'
# with open(gene_dict_path, 'r') as f:
#     gene_name_dict = json.load(f)
#
# # Load the new Excel file with raw data results, skipping the first row
# file_path = 'Inputs/experiments_data/DIO3.xlsx'
# df_raw = pd.read_excel(file_path, header=1)
#
# # Extract the first gene name if multiple names are present
# def extract_first_gene_name(gene_name):
#     return gene_name.split(';')[0].strip() if pd.notna(gene_name) else gene_name
#
# # Apply the function to extract the first gene name
# df_raw['First Gene Name'] = df_raw['Gene names'].apply(extract_first_gene_name)
#
# # Count and print rows with duplicate gene names at the start
# initial_duplicates = df_raw[df_raw.duplicated(subset=['First Gene Name'], keep=False)]
# print("Number of rows with duplicate gene names at the start:", initial_duplicates.shape[0])
# print("Duplicate gene names at the start:\n", initial_duplicates['First Gene Name'].drop_duplicates().values)
#
# # Map Gene names to Gene IDs
# df_raw['GeneID'] = df_raw['First Gene Name'].map(gene_name_dict)
#
# # Print list of gene names without an ID
# gene_names_without_id = df_raw[df_raw['GeneID'].isna()]['First Gene Name'].unique()
# print("Gene names without an ID:", gene_names_without_id)
#
# # Drop rows with missing GeneID values
# df_raw = df_raw.dropna(subset=['GeneID'])
#
# # Define columns for control and experiment
# control_cols = ['LFQ intensity Control_H8_scrambled-1_1', 'LFQ intensity Control_H8_scrambled-2_1']
# experiment_cols = ['LFQ intensity exp_C11-D3-KD-1_1', 'LFQ intensity exp_C11-D3-KD-2_1']
#
# # Aggregate duplicates by taking the mean of LFQ intensity columns
# df_raw = df_raw.groupby(['First Gene Name', 'GeneID']).agg({
#     'LFQ intensity Control_H8_scrambled-1_1': 'mean',
#     'LFQ intensity Control_H8_scrambled-2_1': 'mean',
#     'LFQ intensity exp_C11-D3-KD-1_1': 'mean',
#     'LFQ intensity exp_C11-D3-KD-2_1': 'mean',
#     'Gene names': 'first'
# }).reset_index()
#
# # Filter genes with at least one non-zero control and one non-zero experiment value
# df_filtered = df_raw[(df_raw[control_cols].replace(0, np.nan).count(axis=1) > 0) &
#                      (df_raw[experiment_cols].replace(0, np.nan).count(axis=1) > 0)].copy()
#
# # Calculate average LFQ intensities ignoring zeros
# df_filtered['Control_Avg'] = df_filtered[control_cols].replace(0, np.nan).mean(axis=1)
# df_filtered['Experiment_Avg'] = df_filtered[experiment_cols].replace(0, np.nan).mean(axis=1)
#
# # Calculate fold change and log2 fold change
# df_filtered['Fold_Change'] = df_filtered['Experiment_Avg'] / df_filtered['Control_Avg']
# df_filtered['Score'] = np.log2(df_filtered['Fold_Change'])
#
# # Create a P-value column with empty values
# df_filtered['P-value'] = ''
#
# # Rename 'First Gene Name' to 'Symbol'
# df_filtered = df_filtered.rename(columns={'First Gene Name': 'Symbol'})
#
# # Count and print rows with duplicate gene names at the end
# final_duplicates = df_filtered[df_filtered.duplicated(subset=['Symbol'], keep=False)]
# print("Number of rows with duplicate gene names at the end:", final_duplicates.shape[0])
# print("Duplicate gene names at the end:\n", final_duplicates['Symbol'].drop_duplicates().values)
#
# # Check for duplicate Gene IDs
# duplicate_gene_ids = df_filtered[df_filtered.duplicated(subset=['GeneID'], keep=False)]
# print("Number of rows with duplicate Gene IDs at the end:", duplicate_gene_ids.shape[0])
# print("Duplicate Gene IDs at the end:\n", duplicate_gene_ids['GeneID'].drop_duplicates().values)
#
# # Ensure 'GeneID' column exists before attempting to drop duplicates
# if 'GeneID' in df_filtered.columns:
#     df_filtered = df_filtered.drop_duplicates(subset='GeneID')
# else:
#     raise KeyError("The 'GeneID' column is not found in the DataFrame.")
#
# # Select and reorder the required columns
# df_final = df_filtered[['GeneID', 'Symbol', 'Score', 'P-value']]
#
# # Save the final DataFrame to a new Excel file
# output_path = 'Inputs/experiments_data/Metabolic/DIO3.xlsx'
# df_final.to_excel(output_path, index=False)
#
# print("The data has been formatted and saved to", output_path)

# Load the gene name dictionary from a JSON file
gene_dict_path = '../Data/H_sapiens/gene_names/gene_info.json'
with open(gene_dict_path, 'r') as f:
    gene_name_dict = json.load(f)

# Load the Excel file with raw data results
file_path = 'Inputs/experiments_data/Metabolic/Article Tables.xlsx'
df_raw = pd.read_excel(file_path, sheet_name='Table 1')

# Extract the first gene name if multiple names are present
def extract_first_gene_name(gene_name):
    return gene_name.split(';')[0].strip() if pd.notna(gene_name) else gene_name

# Apply the function to extract the first gene name
df_raw['First Gene Name'] = df_raw['Gene name'].apply(extract_first_gene_name)

# Count and print rows with duplicate gene names at the start
initial_duplicates = df_raw[df_raw.duplicated(subset=['First Gene Name'], keep=False)]
print("Number of rows with duplicate gene names at the start:", initial_duplicates.shape[0])
print("Duplicate gene names at the start:\n", initial_duplicates['First Gene Name'].drop_duplicates().values)

# Map Gene names to Gene IDs
df_raw['GeneID'] = df_raw['First Gene Name'].map(gene_name_dict)

# Print list of gene names without an ID
gene_names_without_id = df_raw[df_raw['GeneID'].isna()]['First Gene Name'].unique()
print("Gene names without an ID:", gene_names_without_id)

# Drop rows with missing GeneID values
df_raw = df_raw.dropna(subset=['GeneID'])

# No control or experiment columns in the provided data, skipping this part

# Create a P-value column from the existing p-value data
df_raw['P-value'] = df_raw['p-value']

# Rename 'First Gene Name' to 'Symbol'
df_raw = df_raw.rename(columns={'First Gene Name': 'Symbol'})

# Count and print rows with duplicate gene names at the end
final_duplicates = df_raw[df_raw.duplicated(subset=['Symbol'], keep=False)]
print("Number of rows with duplicate gene names at the end:", final_duplicates.shape[0])
print("Duplicate gene names at the end:\n", final_duplicates['Symbol'].drop_duplicates().values)

# Check for duplicate Gene IDs
duplicate_gene_ids = df_raw[df_raw.duplicated(subset=['GeneID'], keep=False)]
print("Number of rows with duplicate Gene IDs at the end:", duplicate_gene_ids.shape[0])
print("Duplicate Gene IDs at the end:\n", duplicate_gene_ids['GeneID'].drop_duplicates().values)

# Ensure 'GeneID' column exists before attempting to drop duplicates
if 'GeneID' in df_raw.columns:
    df_raw = df_raw.drop_duplicates(subset='GeneID')
else:
    raise KeyError("The 'GeneID' column is not found in the DataFrame.")

# Select and reorder the required columns
df_final = df_raw[['GeneID', 'Symbol', 'P-value']]

# Save the final DataFrame to a new Excel file
output_path = '../Inputs/experiments_data/after_filter.xlsx'
df_final.to_excel(output_path, index=False)

print("The data has been formatted and saved to", output_path)

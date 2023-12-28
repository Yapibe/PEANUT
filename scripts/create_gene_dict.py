import pandas as pd


def load_and_filter_cortical(file_path):
    """
    Load the CorticalOrganoid file, filter out 'lncRNA', and remove the 'gene_biotype' column.
    """
    df = pd.read_excel(file_path)
    return df


def load_idmap(file_path):
    """
    Load the copy_of_idmap file.
    """
    id_file = pd.read_excel(file_path)
    # print how many duplicates in the query column
    print('duplicates in the query column:', id_file['query'].duplicated().sum())
    # drop duplicates from the query column
    id_file.drop_duplicates(subset=['query'], inplace=True)
    # print how duplicated in the query column after dropping duplicates
    print('duplicates in the query column after dropping duplicates:', id_file['query'].duplicated().sum())
    return id_file


def join_dataframes(df_cortical, df_idmap):
    """
    Join CorticalOrganoid and idmap DataFrames on gene and query columns.
    """
    # print how many rows in each DataFrame
    print('df_cortical:', df_cortical.shape[0])
    print('df_idmap:', df_idmap.shape[0])
    merged_df = pd.merge(df_cortical, df_idmap[['query', 'entrezgene']], left_on='gene', right_on='query', how='left')
    merged_df.drop(columns=['query'], inplace=True)
    # print how many rows in the merged DataFrame
    print('merged_df:', merged_df.shape[0])
    # print how many nan values in the entrezgene column
    print('nan values in entrezgene column:', merged_df['entrezgene'].isna().sum())
    # remove all rows with nan values in the entrezgene column
    merged_df.dropna(subset=['entrezgene'], inplace=True)
    # move column entrezgene to the second column
    merged_df = merged_df[['gene', 'entrezgene'] + [col for col in merged_df.columns if col not in ['gene', 'entrezgene']]]
    # print how many rows in the merged DataFrame after removing nan values
    print('merged_df after removing nan values:', merged_df.shape[0])
    return merged_df


def handle_duplicates(df):
    """
    Handle duplicates by grouping on 'gene' and calculating the mean for other columns.
    """
    # print how many duplicates in the gene column
    print('duplicates in the gene column:', df['gene'].duplicated().sum())
    df_grouped = df.groupby('gene', as_index=False).mean()
    # print how many rows in the grouped DataFrame
    # print how many duplicates in the gene column
    print('duplicates in the gene column after grouping:', df_grouped['gene'].duplicated().sum())
    return df_grouped


# Usage
file_path_cortical = '../Inputs/experiments_data/iPSC_filtered.xlsx'
file_path_idmap = '../Inputs/experiments_data/Copy_of_idmap.xlsx'

df_cortical = load_and_filter_cortical(file_path_cortical)
df_idmap = load_idmap(file_path_idmap)
df_joined = join_dataframes(df_cortical, df_idmap)
df_final = handle_duplicates(df_joined)

# Save the final DataFrame
df_final.to_excel('final_output_iPSC.xlsx', index=False)

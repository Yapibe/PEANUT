import pandas as pd


def load_and_filter_data(file_path):
    df = pd.read_excel(file_path)
    # Load the datasheet file, in column 'gene_biotype' only keep rows with value: 'protein_coding'
    df = df[df['gene_biotype'] == 'protein_coding']
    # remove column 'gene_biotype'
    df.drop(columns=['gene_biotype'], inplace=True)

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
data_file_path = '../pipeline/Inputs/experiments_data/Parkinson/Parkinson_t_v_n_500nm_v_t.xlsx'
file_path_idmap = '../pipeline/Data/H_sapiens/gene_names/ID_to_Name_Map.xlsx'

data_file = load_and_filter_data(data_file_path)

df_idmap = load_idmap(file_path_idmap)
df_joined = join_dataframes(data_file, df_idmap)
df_final = handle_duplicates(df_joined)

# Save the final DataFrame
df_final.to_excel('final_output_iPSC.xlsx', index=False)

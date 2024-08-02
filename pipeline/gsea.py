import gseapy as gp
import pandas as pd
from os import path

def add_description_to_pathway_file(input_path, output_path):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            # Split the line by tab
            parts = line.strip().split('\t')
            # Get the gene set name (first column)
            gene_set_name = parts[0]
            # Duplicate the gene set name as the description
            new_line = '\t'.join([gene_set_name, gene_set_name] + parts[1:])
            outfile.write(new_line + '\n')

def run_gsea(data, gene_sets, out_dir):
    """
    Run GSEA analysis.

    Parameters:
    - data: DataFrame with 'gene' and 'logFC' columns
    - gene_sets: Path to gene sets file in GMT format
    - out_dir: Output directory for GSEA results
    """
    pre_res = gp.prerank(rnk=data, gene_sets=gene_sets, outdir=out_dir, verbose=True, min_size=20, max_size=200, permutation_num=1000, no_plot=True)
    return pre_res


def load_gene_expression_data(file_path):
    """
    Load gene expression data.

    Parameters:
    - file_path: Path to the gene expression data file

    Returns:
    DataFrame with 'gene' and 'logFC' columns
    """
    data = pd.read_excel(file_path)
    # return sorted by descending data columns GeneID and Score
    return data.sort_values(by=['GeneID', 'Score'], ascending=False)


root_dir = path.dirname(path.abspath(__file__))
pathway_path = path.join(root_dir, 'Data', 'Human', 'pathways')
# add_description_to_pathway_file(path.join(pathway_path, 'decoy_pathways'), path.join(pathway_path, 'decoy_pathways_gmt.gmt'))
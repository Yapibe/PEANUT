import pandas as pd
import numpy as np
from os import path

def load_pathways(pathway_file):
    with open(pathway_file, 'r') as f:
        lines = [str.upper(x.strip()).split('\t') for x in f]
    pathways = {x[0]: [int(y) for y in x[2:]] for x in lines}
    return pathways


def create_decoy_pathway(pathway, all_genes):
    all_genes = list(all_genes)
    decoy = pathway.copy()
    for index, gene in enumerate(pathway):
        # Replace each gene with a random gene from the set of all possible genes
        decoy[index] = np.random.choice(all_genes)
    return decoy

def main():
    main_dir = path.dirname(path.realpath(__file__))
    pathway_file = path.join(main_dir,'Data', 'H_sapiens', 'pathways', 'pathway_file')
    pathways = load_pathways(pathway_file)

    all_genes = set()
    decoy_pathways = {}
    for pathway in pathways:
        all_genes.update(pathways[pathway])

    for pathway in pathways:
        decoy = create_decoy_pathway(pathways[pathway], all_genes)
        # create name for decoy pathway
        decoy_name = 'decoy_' + pathway
        decoy_pathways.update({decoy_name: decoy})

    # Combine real and decoy pathways
    combined_pathways = {}
    for pathway in pathways:
        combined_pathways.update({pathway: pathways[pathway]})
    for pathway in decoy_pathways:
        combined_pathways.update({pathway: decoy_pathways[pathway]})

    #write combined pathways to tsv file where first column is name of pathway and then the rest of the columns are genes
    with open(path.join(main_dir,'Data', 'H_sapiens', 'pathways', 'combined_pathways.txt'), 'w') as f:
        for pathway in combined_pathways:
            f.write(pathway + '\t' + '\t'.join([str(x) for x in combined_pathways[pathway]]) + '\n')

if __name__ == "__main__":
    main()
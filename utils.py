import pandas as pd
import networkx as nx
import numpy as np
import os
from os import path

def read_network(network_filename):
    """
    read a network from a file
    :param network_filename: the file containing the network
    :return: a networkx graph
    """
    network = pd.read_table(network_filename, header=None, usecols=[0, 1, 2])
    return nx.from_pandas_edgelist(network, 0, 1, 2)


def read_prior_set(excel_dir):
    prior_data = pd.read_excel(excel_dir, engine='openpyxl')
    return prior_data


def get_propagation_input(prior_gene_ids, prior_data, input_type, network):
    """
    :param network:
    :param prior_gene_ids: list of gene_ids
    :param prior_data: all excel file
    :param input_type:
    :return:
    """

    if input_type == 'ones':
        inputs = {int(x): 1 for x in prior_gene_ids if x in network.nodes}
    elif input_type is None:
        inputs = {int(x): 1 for x in prior_gene_ids if x in network.nodes}
    elif input_type == 'abs_Score':
        inputs = {int(id): np.abs(float(prior_data[prior_data.GeneID == id]['Score'])) for id in prior_gene_ids}
    elif input_type == 'Score':
        inputs = {int(id): float(prior_data[prior_data.GeneID == id]['Score'].values[0]) for id in prior_gene_ids}
    elif input_type == 'Score_all':
        inputs = {int(id): float(prior_data[prior_data.GeneID == id]['Score'].values[0])
                  for name, id in prior_gene_ids}
        mean_input = np.mean([x for x in inputs.values()])
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = mean_input

    elif input_type == 'abs_Score_all':
        inputs = {int(id): np.abs(float(prior_data[prior_data.GeneID == id]['Score'].values[0]))
                  for id in prior_gene_ids}
        mean_input = np.mean([x for x in inputs.values()])
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = mean_input
    elif input_type == 'ones_all':
        inputs = dict()
        for id in network.nodes:
            if id not in inputs:
                inputs[id] = 1
    else:
        assert 0, '{} is not a valid input type'.format(input_type)

    inputs = {id: np.round(input_score, 3) for id, input_score in inputs.items() if id in network.nodes}
    return inputs


def save_file(obj, save_dir=None, compress=True):
    import pickle, zlib
    obj = pickle.dumps(obj)
    if compress:
        obj = zlib.compress(obj)
    with open(save_dir, 'wb') as f:
        pickle.dump(obj, f)
    print('File was saved in {}'.format(save_dir))


def save_propagation_score(propagation_scores, prior_set, propagation_input, genes_idx_to_id, task,
                           save_dir=None, date=None):
    """
    this function saves the propagation scores to a file
    :param propagation_scores: the propagation scores
    :param prior_set: the set of prior genes
    :param propagation_input: the propagation input
    :param genes_idx_to_id: a dictionary from gene index to gene id
    :param task: the arguments of the program
    :param date: the date of the experiment
    :param file_name: the name of the file to save the results to
    :param save_dir: the directory to save the results to
    :return: None
    """
    file_name = f"{task.experiment_name}_{task.propagation_input_type}_{task.alpha}_{task.date}"
    save_dir = save_dir or task.propagation_scores_path
    os.makedirs(save_dir, exist_ok=True)
    propagation_results_path = path.join(save_dir, file_name)

    save_dict = {
        'args': task, 'prior_set': prior_set, 'propagation_input': propagation_input,
        'gene_idx_to_id': genes_idx_to_id, 'gene_prop_scores': propagation_scores,
    }

    save_file(save_dict, propagation_results_path)
    return save_dict


def get_root_path():
    return path.dirname(path.realpath(__file__))
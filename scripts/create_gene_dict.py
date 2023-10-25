import json
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Function to add unique keys of dict1 to dict2
def add_unique_keys_to_dict(dict1, dict2):
    unique_keys_dict1 = set(dict1.keys()) - set(dict2.keys())
    for key in unique_keys_dict1:
        dict2[key] = dict1[key]
    return dict2


# Function to find and print differences between two dictionaries based on their keys
def print_dict_differences(dict1, dict2, labels=('Dict1', 'Dict2')):
    set1 = set(dict1.keys())
    set2 = set(dict2.keys())

    # Find differences
    only_in_dict1 = set1 - set2
    only_in_dict2 = set2 - set1
    common_keys = set1 & set2

    print(f"Keys only in {labels[0]}: {len(only_in_dict1)}")
    print(only_in_dict1)
    print(f"\nKeys only in {labels[1]}: {len(only_in_dict2)}")
    print(only_in_dict2)


def load_dict_from_text_file(input_file_path):
    gene_dict = {}
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:-1]:  # Skip the first and last lines (which are "{" and "}")
            key, value = line.strip().strip(",").split(": ")
            key = key.strip('"')
            value = value.strip('"')
            gene_dict[key] = value
    return gene_dict


# Function to create a Venn diagram to compare keys of two dictionaries
def plot_venn(dict1, dict2, labels=('Dict1', 'Dict2')):
    set1 = set(dict1.keys())
    set2 = set(dict2.keys())

    venn = venn2([set1, set2], set_labels=labels)

    plt.title("Venn Diagram Comparing Dictionary Keys")
    plt.show()


# Define a function to create a gene dictionary from a tab-separated file
def create_gene_dict(file_path):
    gene_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 3:  # Check to ensure the line has enough columns
                gene_id, symbol = columns[1:3]
                gene_dict[symbol] = gene_id
    return gene_dict


# Function to save the dictionary to a text file in the specified format
def save_dict_to_text_file(gene_dict, output_file_path):
    with open(output_file_path, 'w') as file:
        file.write("{\n")
        for key, value in gene_dict.items():
            file.write(f'    "{key}": "{value}",\n')
        file.write("}\n")


# save_dict_to_text_file(gene_dict, 'output.txt')


# Create two gene dictionaries using your function (or any other method)
dict1 = load_dict_from_text_file("../Data/H_sapiens/gene_names/H_sapiens.gene_info")
dict2 = load_dict_from_text_file("../Data/H_sapiens/gene_names/H_sapiens.gene_info")

# Create the Venn diagram
plot_venn(dict1, dict2, labels=('output', 'H_sapiens'))

# Print the differences between the dictionaries
print_dict_differences(dict1, dict2, labels=('File 1', 'File 2'))

updated_dict2 = add_unique_keys_to_dict(dict1, dict2)

# Create the Venn diagram
plot_venn(dict1, dict2, labels=('output', 'H_sapiens'))

# save updated_dict2 to a text file
save_dict_to_text_file(updated_dict2, '../Data/H_sapiens/gene_names/H_sapiens.gene_info')

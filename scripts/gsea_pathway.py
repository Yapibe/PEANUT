from os import path


# Function to add a description column to the pathway file by duplicating the gene set name
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


main_dir = path.dirname(path.dirname(path.realpath(__file__)))
# Example usage
input_pathway_file = path.join(main_dir, 'Data', 'H_sapiens', 'pathways', 'pathway_file')
output_pathway_file = path.join(main_dir, 'Data', 'H_sapiens', 'pathways', 'pathway_file_with_description')

add_description_to_pathway_file(input_pathway_file, output_pathway_file)

import json

# Load the input JSON file
input_file_path = "../Data/H_sapiens/gene_names/hgnc_complete_set.json"
output_file_path = "../Data/H_sapiens/gene_names/hgnc.json"

with open(input_file_path, "r") as input_file:
    data = json.load(input_file)

# Initialize the result dictionary
gene_symbol_dict = {}

# Extract symbols and hgnc_ids
for doc in data["response"]["docs"]:
    hgnc_id = doc.get("entrez_id")
    symbol = doc.get("symbol")

    if hgnc_id:
        # Remove the "HGNC:" prefix and convert to integer
        hgnc_id = int(hgnc_id.replace("HGNC:", ""))

        if symbol:
            gene_symbol_dict[symbol] = hgnc_id

# Extract alias symbols (prev_symbol) and add them if they do not already exist
for doc in data["response"]["docs"]:
    hgnc_id = doc.get("entrez_id")
    alias_symbols = doc.get("prev_symbol", [])

    if hgnc_id:
        hgnc_id = int(hgnc_id.replace("HGNC:", ""))

        for alias in alias_symbols:
            if alias not in gene_symbol_dict:
                gene_symbol_dict[alias] = hgnc_id

# Sort the dictionary by keys
gene_symbol_dict = dict(sorted(gene_symbol_dict.items()))

# Save the resulting dictionary to a new JSON file
with open(output_file_path, "w") as output_file:
    json.dump(gene_symbol_dict, output_file, indent=4)

print(f"Gene symbol dictionary saved to {output_file_path}")

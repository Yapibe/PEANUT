import pandas as pd
import json

# Load the first file (edges with weights) using space as the delimiter, and skip the header row
edges_file = "../Data/H_sapiens/network/protein.links"
edges_df = pd.read_csv(
    edges_file,
    sep="\s+",
    names=[
        "protein1",
        "protein2",
        "experimental",
        "database",
        "textmining",
        "combined_score",
    ],
    dtype={
        "protein1": str,
        "protein2": str,
        "experimental": float,
        "database": float,
        "textmining": float,
        "combined_score": float,
    },
    skiprows=1,
)

# Ensure the experimental score is numeric and divide by 1000
edges_df["experimental"] = edges_df["experimental"] / 1000

# Remove rows where the experimental column is 0
edges_df = edges_df[edges_df["experimental"] != 0]

# Remove unnecessary columns
edges_df = edges_df.drop(columns=["database", "textmining", "combined_score"])

# Load the second file (node translation dictionary)
dictionary_file = "../Data/H_sapiens/network/protein.info"
dictionary_df = pd.read_csv(
    dictionary_file,
    sep="\t",
    usecols=["#string_protein_id", "preferred_name", "annotation"],
)

# Remove rows where 'annotation' contains "novel" or "Uncharacterized"
dictionary_df = dictionary_df[
    ~dictionary_df["annotation"].str.contains(
        "novel|Uncharacterized", case=False, na=False
    )
]

# Create a dictionary for translation, ensure the keys and values are stripped of any extra whitespace
translation_dict = pd.Series(
    dictionary_df["preferred_name"].values,
    index=dictionary_df["#string_protein_id"].str.strip(),
).to_dict()

# Adding an attribute column, here we'll set it to 0 as a placeholder
edges_df["attribute"] = 0

# Translate node IDs using the dictionary, stripping whitespace just in case
edges_df["protein1"] = (
    edges_df["protein1"].str.strip().map(translation_dict).fillna(edges_df["protein1"])
)
edges_df["protein2"] = (
    edges_df["protein2"].str.strip().map(translation_dict).fillna(edges_df["protein2"])
)

# Load the gene_info.json file
with open("../Data/H_sapiens/gene_names/gene_info.json", "r") as f:
    gene_info = json.load(f)

# Map gene symbols to their gene IDs using the gene_info dictionary
edges_df["protein1"] = edges_df["protein1"].map(gene_info).fillna(edges_df["protein1"])
edges_df["protein2"] = edges_df["protein2"].map(gene_info).fillna(edges_df["protein2"])

# Convert the gene IDs to integers, and handle cases where the mapping was not successful
edges_df["protein1"] = pd.to_numeric(edges_df["protein1"], errors="coerce").fillna(
    edges_df["protein1"]
)
edges_df["protein2"] = pd.to_numeric(edges_df["protein2"], errors="coerce").fillna(
    edges_df["protein2"]
)

# Identify rows where gene ID has letters
invalid_rows = edges_df[
    (
        edges_df["protein1"].apply(
            lambda x: isinstance(x, str) and any(c.isalpha() for c in x)
        )
    )
    | (
        edges_df["protein2"].apply(
            lambda x: isinstance(x, str) and any(c.isalpha() for c in x)
        )
    )
]

# Print the invalid gene IDs
print("Invalid gene IDs:")
print(invalid_rows[["protein1", "protein2"]])

# Count and print the number of rows with invalid gene IDs
num_invalid_rows = len(invalid_rows)
print(f"Number of rows with gene IDs containing letters: {num_invalid_rows}")

# Remove rows with invalid gene IDs
edges_df = edges_df.drop(invalid_rows.index)

# Convert the remaining gene IDs to integers
edges_df["protein1"] = edges_df["protein1"].astype(int)
edges_df["protein2"] = edges_df["protein2"].astype(int)

# Save the final output to a new file
output_file = "../Data/H_sapiens/network/String_"
edges_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"Network file saved to {output_file}")

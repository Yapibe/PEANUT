import json
import os

import pandas as pd

# Load the gene_info.json file
gene_dict_path = "../Data/H_sapiens/gene_names/hgnc.json"
with open(gene_dict_path, "r") as f:
    gene_name_dict = json.load(f)


def read_gct(file_path):
    # Read the actual data starting from the third line (skipping metadata lines)
    data = pd.read_csv(file_path, sep="\t", skiprows=2, header=None)
    data.columns = ["Symbol", "Description", "Score"]
    # Convert Score to numeric, coerce errors to NaN
    data["Score"] = pd.to_numeric(data["Score"], errors="coerce")
    return data


# Directories for GCT, RNK, and XLSX files
gct_dir = "../Inputs/experiments_data/NGSEA/GCT"
rnk_dir = "../Inputs/experiments_data/NGSEA/RNK"
xlsx_dir = "../Inputs/experiments_data/NGSEA/XLSXN"

# Ensure RNK and XLSX directories exist
os.makedirs(rnk_dir, exist_ok=True)
os.makedirs(xlsx_dir, exist_ok=True)

# Loop through all GCT files in the directory
for file_name in os.listdir(gct_dir):
    if file_name.endswith(".gct"):
        gct_file_path = os.path.join(gct_dir, file_name)

        # Read the GCT file
        gct_data = read_gct(gct_file_path)

        # Create RNK file
        rnk_data = gct_data[["Symbol", "Score"]]
        rnk_file_name = file_name.replace(".gct", ".rnk")
        rnk_file_path = os.path.join(rnk_dir, rnk_file_name)
        rnk_data.to_csv(
            rnk_file_path, sep="\t", header=False, index=False
        )

        # Read back the RNK file and remove the first line
        with open(rnk_file_path, "r") as rnk_file:
            lines = rnk_file.readlines()
        with open(rnk_file_path, "w") as rnk_file:
            rnk_file.writelines(lines[1:])

        # Map Symbols to GeneIDs
        gct_data["GeneID"] = gct_data["Symbol"].map(gene_name_dict)

        # Print rows without a corresponding GeneID with values of more than abs(1)
        missing_gene_ids = gct_data[gct_data["GeneID"].isnull()]
        missing_gene_ids = missing_gene_ids[
            missing_gene_ids["Score"].abs() > 1
        ]
        if not missing_gene_ids.empty:
            print(
                f"Missing GeneIDs for the following symbols in {file_name}:"
            )
            print(missing_gene_ids[["Symbol", "Description"]])

        # Remove rows without a GeneID
        gct_data = gct_data.dropna(subset=["GeneID"])

        # Ensure GeneID is of integer type
        gct_data["GeneID"] = gct_data["GeneID"].astype(int)

        # Add a P-value column with all values set to zero
        gct_data["P-value"] = 0

        # Select only the required columns
        final_data = gct_data[
            ["GeneID", "Symbol", "Score", "P-value"]
        ]

        # Save the updated DataFrame to a new Excel file
        output_file_name = file_name.replace(".gct", ".xlsx")
        output_path = os.path.join(xlsx_dir, output_file_name)
        final_data.to_excel(output_path, index=False)
        print(f"The data has been updated and saved to {output_path}")

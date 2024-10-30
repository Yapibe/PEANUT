import pandas as pd
import os

# Get the project root directory (assuming your script is in the pipeline folder)
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
new_file_path = os.path.join(project_root, "pipeline", "Data", "H_sapiens", "pathways", "c2.cp.v2024.1.Hs.entrez.gmt")
output_folder = os.path.join(project_root, "pipeline", "Data", "H_sapiens", "pathways")
os.makedirs(output_folder, exist_ok=True)  # Create the folder if it doesn't exist
output_file_path = os.path.join(output_folder, "c2cp")

new_data = []

with open(new_file_path, 'r') as file:
    for line in file:
        parts = line.strip().split('\t')
        pathway_name = parts[0]
        genes = parts[2:]  # Skip the URL part
        new_data.append([pathway_name] + genes)

# Convert to DataFrame
new_pathways_df = pd.DataFrame(new_data)

# Save the DataFrame to a TSV file in the specified folder
new_pathways_df.to_csv(output_file_path, index=False, header=False, sep='\t')

# Display the first few rows to verify the structure
new_pathways_df.head()

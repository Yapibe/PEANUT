import os
import pandas as pd
import numpy as np
import networkx as nx
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from settings import Settings
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import load_pathways_genes, read_network, read_prior_set
import time

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Directories and settings
input_dir = os.path.join("../Inputs", "experiments_data", "GSE", "XLSX")
output_base_dir = os.path.join("../Outputs", "NGSEA")
summary_base_dir = os.path.join(output_base_dir, "Summary")
pathways_dir = os.path.join("../Data", "Anat", "pathways")
networks = ["Anat", "String", "HumanNet", "String_"]
pathway_files = ["kegg"]
prop_methods = ["MW", "PROP", "ABS_PROP"]
alphas = [0.1, 0.2]
max_workers = 60

# Ensure output directories exist
os.makedirs(summary_base_dir, exist_ok=True)

# List all files to be processed
file_list = [f for f in os.listdir(input_dir) if f.endswith(".xlsx")]


# Load network and pathways data
def load_data():
    loaded_networks = {
        name: read_network(os.path.join("../Data", "Anat", "network", name))
        for name in networks
    }
    loaded_pathways = {
        file: load_pathways_genes(os.path.join(pathways_dir, file))
        for file in pathway_files
    }
    return loaded_networks, loaded_pathways


def run_analysis(
    test_name,
    prior_data,
    network,
    network_name,
    alpha,
    method,
    output_path,
    pathway_file,
):
    general_args = Settings(
        network=network_name,
        pathway_file=pathway_file,
        method=method,
        alpha=alpha if method in ["PROP", "ABS_PROP"] else 1,
        run_propagation=True,
    )

    if method == "ABS_PROP":
        general_args.input_type = "abs_Score"

    if general_args.run_propagation:
        perform_propagation(test_name, general_args, network, prior_data)

    perform_enrichment(test_name, general_args, output_path)


# Updated get_pathway_rank function
def get_pathway_rank(output_path, pathway_name):
    try:
        results_df = pd.read_excel(output_path)
        pathway_row = results_df[results_df["Term"] == pathway_name]
        if not pathway_row.empty:
            fdr_p_val = np.float64(pathway_row["FDR q-val"].values[0])
            matching_rows = results_df[results_df["FDR q-val"] == fdr_p_val]
            return matching_rows.index[0], fdr_p_val
    except Exception as e:
        logger.error(f"Error reading {output_path}: {e}")
    return None, None


# Process a single file
# Load the pre-calculated data into a dictionary or DataFrame
pre_calculated_data = {
    "KEGG_THYROID_CANCER": (1.954022989, 29, 4),
    "KEGG_NON_SMALL_CELL_LUNG_CANCER": (2.062678063, 54, 4),
    "KEGG_ACUTE_MYELOID_LEUKEMIA": (2.039721946, 57, 5),
    "KEGG_COLORECTAL_CANCER": (2.107879429, 62, 4),
    "KEGG_GLIOMA": (2.039072039, 65, 4),
    "KEGG_RENAL_CELL_CARCINOMA": (2.11032967, 70, 5),
    "KEGG_PANCREATIC_CANCER": (2.140724947, 70, 5),
    "KEGG_PROSTATE_CANCER": (2.106543291, 89, 5),
    "KEGG_DILATED_CARDIOMYOPATHY": (2.454394693, 90, 7),
    "KEGG_PARKINSONS_DISEASE": (1.949191984, 130, 5),
    "KEGG_ALZHEIMERS_DISEASE": (2.150090558, 166, 5),
    "KEGG_HUNTINGTONS_DISEASE": (1.758788428, 182, 4),
}


def get_pre_calculated_data(pathway_name):
    return pre_calculated_data.get(pathway_name, (None, None, None))


def process_single_file(network, pathway_file, network_name, alpha, method, file_name):
    dataset_name, pathway_name = file_name.replace(".xlsx", "").split("_", 1)
    prior_data = read_prior_set(os.path.join(input_dir, file_name))
    output_dir = os.path.join(
        output_base_dir, method, network_name, pathway_file, f"alpha_{alpha}"
    )
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    # Get the pre-calculated values
    pathway_density, num_genes, avg_diameter = get_pre_calculated_data(pathway_name)

    run_analysis(
        file_name,
        prior_data,
        network,
        network_name,
        alpha,
        method,
        output_file_path,
        pathway_file,
    )
    rank, fdr_q_val = get_pathway_rank(output_file_path, pathway_name)

    return {
        "Dataset": dataset_name,
        "Pathway": pathway_name,
        "Network": network_name,
        "Pathway file": pathway_file,
        "Alpha": alpha,
        "Method": method,
        "Rank": rank,
        "FDR q-val": fdr_q_val,
        "Significant": int(float(fdr_q_val) < 0.05) if fdr_q_val else 0,
        "Density": pathway_density,
        "Num Genes": num_genes,
        "Avg Diameter": avg_diameter,
    }


def generate_summary_df(filtered_df):
    # Pivot the DataFrame to organize by dataset, pathway, and method
    pivot_df = filtered_df.pivot_table(
        index=["Dataset", "Pathway", "Density", "Num Genes", "Avg Diameter"],
        columns="Method",
        values=["Rank", "Significant"],
        aggfunc="first",
    ).reset_index()

    # Flatten the multi-level column index created by pivot_table
    pivot_df.columns = [
        "Dataset",
        "Pathway",
        "Density",
        "Num Genes",
        "Avg Diameter",
    ] + [f"{col[1]} {col[0]}" for col in pivot_df.columns[5:]]

    # Sort columns so that each method has its rank and significance together
    ordered_columns = ["Dataset", "Pathway", "Density", "Num Genes", "Avg Diameter"]
    for method in prop_methods:
        ordered_columns += [f"{method} Rank", f"{method} Significant"]

    pivot_df = pivot_df[ordered_columns]

    # Add missing columns for each method if not present
    for method in prop_methods:
        pivot_df[f"{method} Rank"] = pivot_df.get(f"{method} Rank", np.nan)
        pivot_df[f"{method} Significant"] = pivot_df.get(
            f"{method} Significant", np.nan
        )

    # Calculate the average for all columns (both rank and significant)
    avg_row = pivot_df.mean(numeric_only=True).to_frame().T
    avg_row[["Dataset", "Pathway", "Density", "Num Genes", "Avg Diameter"]] = [
        "Average",
        "",
        "",
        "",
        "",
    ]

    # Concatenate the average row to the summary DataFrame
    summary_df = pd.concat([pivot_df, avg_row], ignore_index=True)

    # Sort the summary DataFrame by Dataset A to Z
    summary_df.sort_values(by="Dataset", inplace=True)

    return summary_df


# Main process
start_time = time.time()
loaded_networks, loaded_pathways = load_data()
logger.info("Networks loaded and pathway densities calculated.")

futures = []
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    for network_name in networks:
        for pathway_file in pathway_files:
            for alpha in alphas:
                for file_name in file_list:
                    dataset_name, pathway_name = file_name.replace(".xlsx", "").split(
                        "_", 1
                    )
                    if pathway_name in loaded_pathways[pathway_file]:
                        for prop_method in prop_methods:
                            futures.append(
                                executor.submit(
                                    process_single_file,
                                    loaded_networks[network_name],
                                    pathway_file,
                                    network_name,
                                    alpha,
                                    prop_method,
                                    file_name,
                                )
                            )

results = [future.result() for future in as_completed(futures)]
results_df = pd.DataFrame(results)

# Save summary DataFrames
for network_name in networks:
    for pathway_file in pathway_files:
        for alpha in alphas:
            filtered_df = results_df[
                (results_df["Network"] == network_name)
                & (results_df["Pathway file"] == pathway_file)
                & (results_df["Alpha"] == alpha)
            ]
            summary_df = generate_summary_df(filtered_df)

            summary_output_dir = os.path.join(
                summary_base_dir, network_name, pathway_file, f"alpha {alpha}"
            )
            os.makedirs(summary_output_dir, exist_ok=True)
            rankings_output_path = os.path.join(
                summary_output_dir,
                f"MW_rankings_summary_{network_name}_{pathway_file}_alpha_{alpha}.xlsx",
            )
            summary_df.to_excel(rankings_output_path, index=False)
            logger.info(f"MW Rankings summary saved to {rankings_output_path}")

elapsed_time = time.time() - start_time
logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

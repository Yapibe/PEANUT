from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pandas as pd
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend
import logging
from args import GeneralArgs
from pathway_enrichment import perform_enrichment
from propagation_routines import perform_propagation
from utils import read_network, read_prior_set
import time
from tqdm import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Create a logger
logger = logging.getLogger(__name__)

# Define directories
pipeline_dir = os.path.dirname(os.path.abspath(__file__))
input_dir = os.path.join(pipeline_dir, "Inputs", "experiments_data", "GSE", "XLSX")
output_base_dir = os.path.join(pipeline_dir, "Outputs", "NGSEA")
summary_base_dir = os.path.join(output_base_dir, "Summary")

# Ensure summary output directories exist
os.makedirs(summary_base_dir, exist_ok=True)


def run_propagation_and_enrichment(
    test_name, prior_data, network, network_name, alpha, method, output_path
):
    if method in ["PROP", "ABS_PROP"]:
        # Set alpha before initializing GeneralArgs for PROP and ABS_PROP
        general_args = GeneralArgs(
            network=network_name, pathway_file="kegg_decoy", method=method, alpha=alpha
        )
        if method == "ABS_PROP":
            general_args.input_type = "abs_Score"
    elif method in ["GSEA", "NGSEA"]:
        # Initialize GeneralArgs for GSEA without modifying alpha
        general_args = GeneralArgs(
            network=network_name, pathway_file="kegg_decoy", method=method
        )
        if method == "NGSEA":
            general_args.run_NGSEA = True

    perform_propagation(test_name, general_args, network, prior_data)
    perform_enrichment(test_name, general_args, output_path)


def calculate_significant_pathway_percentage(gsea_output_path):
    try:
        results_df = pd.read_excel(gsea_output_path)
        if "FDR q-val" in results_df.columns:
            significant_pathways = results_df[results_df["FDR q-val"] < 0.05]
            if not results_df.empty:
                percentage_significant = (
                    len(significant_pathways) / len(results_df)
                ) * 100
                return percentage_significant
    except Exception as e:
        logger.error(f"Error reading {gsea_output_path}: {e}")
    return 0.0


def process_file(network, network_name, alpha, prop_method, file_name):
    dataset_name = file_name.replace(".xlsx", "")
    prior_data = read_prior_set(os.path.join(input_dir, file_name))
    output_dir = os.path.join(
        output_base_dir, prop_method, network_name, "Decoy", f"alpha_{alpha}"
    )
    os.makedirs(output_dir, exist_ok=True)
    output_file_path = os.path.join(output_dir, file_name)

    run_propagation_and_enrichment(
        file_name,
        prior_data,
        network,
        network_name,
        alpha,
        prop_method,
        output_file_path,
    )
    significant_percentage = calculate_significant_pathway_percentage(output_file_path)

    return {
        "Dataset": dataset_name,
        "Network": network_name,
        "Alpha": alpha,
        "Method": prop_method,
        "Significant Percentage": significant_percentage,
    }


# Start timing the entire process
start_time = time.time()

network_names = ["HumanNet", "String_", "Anat", "String"]
prop_methods = ["PROP", "GSEA", "NGSEA", "ABS_PROP"]
alpha_values = [0.2, 0.1]

file_list = [f for f in os.listdir(input_dir) if f.endswith(".xlsx")]

logger.info("Loading network and processing files...")

# Process files in parallel using ProcessPoolExecutor
futures = []
with ProcessPoolExecutor(max_workers=60) as executor:
    for network_name in network_names:
        network_file = os.path.join(
            pipeline_dir, "Data", "Anat", "network", network_name
        )
        network = read_network(network_file)

        for alpha in alpha_values:
            for file_name in file_list:
                for prop_method in prop_methods:
                    futures.append(
                        executor.submit(
                            process_file,
                            network,
                            network_name,
                            alpha,
                            prop_method,
                            file_name,
                        )
                    )

results = []
for future in tqdm(as_completed(futures), total=len(futures), desc="Processing Files"):
    results.append(future.result())

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Calculate the average percentage of significant pathways per method
average_sig_percent = (
    results_df.groupby(["Network", "Alpha", "Method"])["Significant Percentage"]
    .mean()
    .reset_index()
)
average_sig_percent.columns = [
    "Network",
    "Alpha",
    "Method",
    "Average Significant Percentage",
]

# Save the summary DataFrame
for network_name in network_names:
    for alpha in alpha_values:
        summary_output_dir = os.path.join(
            summary_base_dir, network_name, "Decoy", f"alpha_{alpha}"
        )
        os.makedirs(summary_output_dir, exist_ok=True)
        summary_output_path = os.path.join(
            summary_output_dir, f"Decoy_Summary_{network_name}_alpha_{alpha}.xlsx"
        )
        average_sig_percent[
            (average_sig_percent["Network"] == network_name)
            & (average_sig_percent["Alpha"] == alpha)
        ].to_excel(summary_output_path, index=False)
        logger.info(f"Significance summary saved to {summary_output_path}")

end_time = time.time()
elapsed_time = end_time - start_time
logger.info(f"Total time taken: {elapsed_time:.2f} seconds")

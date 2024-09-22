from pipeline.settings import Settings
from pipeline.propagation_routines import perform_propagation
from pipeline.utils import read_network
from app.models import SettingsInput
from typing import List
import pandas as pd


async def pipeline_main(prior_data_list: List[pd.DataFrame], settings: SettingsInput) -> tuple:
    """
    Execute propagation and enrichment analysis based on specified flags.

    Initializes tasks for propagation and enrichment, executes them based on the provided flags,
    and processes results for visualization.

    Parameters:
    - prior_data_list (List[pd.DataFrame]): List of DataFrames containing the input data.
    - settings (SettingsInput): SettingsInput object containing the user-defined settings.

    Returns:
    - tuple: A tuple containing the list of output paths for each condition and the path to the visualization file.
    """
    run_settings = Settings(
        experiment_name=settings.test_name,
        species=settings.species,
        alpha=settings.alpha,
        network=settings.network,
        pathway_file=settings.pathway_file,
        minimum_gene_per_pathway=settings.min_genes_per_pathway,
        maximum_gene_per_pathway=settings.max_genes_per_pathway,
        FDR_THRESHOLD=settings.fdr_threshold,
        JAC_THRESHOLD=settings.JAC_threshold
    )
    network = read_network(run_settings.network_file_path)

    for prior_data in prior_data_list:
        perform_propagation(run_settings, network=network, prior_data=prior_data)

    return run_settings.condition_output_paths, run_settings.visualization_file_path


if __name__ == '__main__':
    settings = SettingsInput(
        test_name="test",
        species="H_sapiens",
        alpha=1,
        network="H_sapiens",
        pathway_file="kegg",
        min_genes_per_pathway=15,
        max_genes_per_pathway=500,
        fdr_threshold=0.05,
        JAC_threshold=0.2
    )

    # Example prior data list
    prior_data_list = [
        pd.read_excel("path/to/experiment1.xlsx"),
        pd.read_excel("path/to/experiment2.xlsx"),
        pd.read_excel("path/to/experiment3.xlsx")
    ]

    condition_output_paths, visualization_file_path = pipeline_main(prior_data_list, settings)
    print("Condition output paths:", condition_output_paths)
    print("Visualization file path:", visualization_file_path)
from os import path
from pipeline.settings import Settings
from pipeline.propagation_routines import perform_propagation
from pipeline.utils import read_network, calculate_time, read_prior_set
import logging
from app.models import SettingsInput
import pandas as pd

async def pipeline_main(prior_data: pd.DataFrame, settings: SettingsInput):
    """
    Execute propagation and enrichment analysis based on specified flags.

    Initializes tasks for propagation and enrichment, executes them based on the provided flags,
    and processes results for visualization.

    Parameters:
    - settings: SettingsInput object containing the user-defined settings
    - prior_data: DataFrame containing the input data

    Returns:
    - None
    """

    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(message)s:        ')
    logger = logging.getLogger(__name__)
    run_settings = Settings(
        experiment_name= settings.test_name,
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
    perform_propagation(run_settings, network=network, prior_data=prior_data)

    return run_settings.output_path


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
    prior_data = read_prior_set("pipeline/Inputs/experiments_data/GSE/XLSX/GSE5281VCX_KEGG_ALZHEIMERS_DISEASE.xlsx")
    output_path = pipeline_main(settings, prior_data)
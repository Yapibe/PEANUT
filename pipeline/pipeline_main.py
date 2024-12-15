import asyncio
import logging
import zipfile
from pathlib import Path
from typing import List, Tuple, Dict, Any
import pandas as pd
import yaml
from dataclasses import replace
import sys

# Set up logger
logger = logging.getLogger(__name__)

# Update Python path to include parent directory
current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent
sys.path.append(str(parent_dir))

# Import necessary modules
from app.models import SettingsInput
from pipeline.propagation_routines import perform_propagation
from pipeline.settings import ConditionSettings, Settings
from pipeline.utils import process_condition, read_network
from pipeline.visualization_tools import plot_pathways_mean_scores


async def pipeline_main(
    conditions_data: List[Tuple[str, pd.DataFrame]],
    settings: SettingsInput,
) -> Path:
    """
    Main pipeline execution function for propagation and GSEA analyses.

    Parameters:
    - conditions_data (List[Tuple[str, pd.DataFrame]]): List of condition names and their data.
    - settings (SettingsInput): User-defined pipeline settings.

    Returns:
    - Path: Path to the ZIP archive containing all results and plots.
    """
    try:
        # Initialize settings and load the network
        run_settings = Settings(
            experiment_name=settings.test_name,
            species=settings.species,
            network=settings.network,
            create_similarity_matrix=settings.create_similarity_matrix,
            alpha=settings.alpha,
            pathway_file=settings.pathway_file,
            minimum_gene_per_pathway=settings.min_genes_per_pathway,
            maximum_gene_per_pathway=settings.max_genes_per_pathway,
            FDR_THRESHOLD=settings.fdr_threshold,
        )
        network = read_network(run_settings.network_file_path)

        # Lists to hold settings for each condition
        all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []
        gsea_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []

        # # Process conditions for both propagation and GSEA
        # for condition_name, condition_data_df in conditions_data:
        #     # Propagation run
        #     propagation_settings = ConditionSettings(
        #         condition_name=f"{condition_name}_Propagation",
        #         experiment_name=run_settings.experiment_name,
        #     )
        #     await perform_propagation(
        #         run_settings, propagation_settings, network, condition_data_df
        #     )
        #     all_condition_settings.append((propagation_settings, condition_data_df))

        #     # GSEA run (if enabled)
        #     if settings.run_gsea:
        #         gsea_settings = ConditionSettings(
        #             condition_name=f"{condition_name}_GSEA",
        #             experiment_name=run_settings.experiment_name,
        #         )
        #         gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
        #         await perform_propagation(
        #             gsea_settings_run, gsea_settings, network, condition_data_df
        #         )
        #         gsea_condition_settings.append((gsea_settings, condition_data_df))

        for condition_path in conditions_data:
            condition_data_df = pd.read_excel(condition_path)
            condition_settings = ConditionSettings(
                condition_name=condition_path.stem,
                experiment_name=run_settings.experiment_name,
            )
            await perform_propagation(
                run_settings, condition_settings, network, condition_data_df
            )
            all_condition_settings.append((condition_settings, condition_data_df))

            if settings.run_gsea:
                gsea_settings = ConditionSettings(
                    condition_name=f"{condition_path.stem}_GSEA",
                    experiment_name=run_settings.experiment_name,
                )
                gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
                await perform_propagation(
                    gsea_settings_run, gsea_settings, network, condition_data_df
                )
                gsea_condition_settings.append((gsea_settings, condition_data_df))

        # Aggregate and process results for propagation
        all_pathways: Dict[str, Dict[str, Any]] = {}
        for condition_settings, condition_data_df in all_condition_settings:
            process_condition(
                condition_settings=condition_settings,
                settings=run_settings,
                condition_data_df=condition_data_df,
                all_pathways=all_pathways,
            )

        # Aggregate and process results for GSEA
        gsea_pathways = None
        if settings.run_gsea:
            gsea_pathways: Dict[str, Dict[str, Any]] = {}
            for condition_settings, condition_data_df in gsea_condition_settings:
                process_condition(
                    condition_settings=gsea_settings,
                    settings=gsea_settings_run,
                    condition_data_df=condition_data_df,
                    all_pathways=gsea_pathways,
                )

        # Generate plots
        propagation_plot_path = plot_pathways_mean_scores(run_settings, all_pathways, gsea_pathways)
        # Collect all output files
        output_files = [settings.output_path for settings, _ in all_condition_settings + gsea_condition_settings]
        if propagation_plot_path:
            output_files.append(propagation_plot_path)

        # Create ZIP archive for all results
        output_zip_dir = run_settings.root_folder / "Outputs" / "Archives"
        output_zip_dir.mkdir(parents=True, exist_ok=True)
        output_zip_path = output_zip_dir / f"{run_settings.experiment_name}_results.zip"

        with zipfile.ZipFile(output_zip_path, 'w') as zipf:
            for file_path in output_files:
                zipf.write(file_path, arcname=file_path.name)

        logger.info(f"Results archived at: {output_zip_path}")
        return output_zip_path

    except Exception as e:
        logger.error(f"An error occurred in pipeline_main: {e}", exc_info=True)
        raise


async def main():
    """Entry point for the pipeline."""
    try:
        output_zip_path = await pipeline_main(conditions_data=condition_data_files, settings=settings)
        logger.info(f"Results successfully archived at: {output_zip_path}")
    except Exception as e:
        logger.error(f"An error occurred in main: {e}")
        raise


if __name__ == "__main__":
    try:
        # Load configuration from YAML file
        config_path = Path(__file__).resolve().parent / "config.yaml"
        with config_path.open("r") as f:
            config = yaml.safe_load(f)

        settings = SettingsInput(**config)
        condition_data_files = [Path(file) for file in config["condition_data_files"]]

        asyncio.run(main())
    except Exception as e:
        logger.error(f"An error occurred in __main__: {e}")
        raise

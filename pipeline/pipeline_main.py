import asyncio
import logging
import zipfile
from pathlib import Path
from typing import List, Tuple, Dict, Any
import pandas as pd
import yaml

import sys

logger = logging.getLogger(__name__)

# Add the parent directory to the Python path
current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent
sys.path.append(str(parent_dir))

from app.models import SettingsInput
from pipeline.propagation_routines import perform_propagation
from pipeline.settings import ConditionSettings, Settings
from pipeline.utils import (
    process_condition,
    read_network,
)
from pipeline.visualization_tools import plot_pathways_mean_scores


async def pipeline_main(
    conditions_data: List[Tuple[str, pd.DataFrame]],
    settings: SettingsInput,
) -> Path:
    """
    Execute propagation and enrichment analysis for multiple conditions.

    Parameters:
    - conditions_data: List of tuples (condition_name, pd.DataFrame)
    - settings: SettingsInput object containing user-defined settings

    Returns:
    - output_zip_path: Path to the output ZIP file containing results
    """
    try:
        run_settings = Settings(
            experiment_name=settings.test_name,
            species=settings.species,
            alpha=settings.alpha,
            network=settings.network,
            pathway_file=settings.pathway_file,
            minimum_gene_per_pathway=settings.min_genes_per_pathway,
            maximum_gene_per_pathway=settings.max_genes_per_pathway,
            FDR_THRESHOLD=settings.fdr_threshold
        )
        network = read_network(run_settings.network_file_path)

        all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []

        for condition_name, condition_data_df in conditions_data:
            condition_settings = ConditionSettings(
                condition_name=condition_name,
                experiment_name=run_settings.experiment_name,
            )
        # for condition_path in conditions_data:
        #     condition_data_df = pd.read_excel(condition_path)
        #     condition_settings = ConditionSettings(
        #         condition_name=condition_path.stem,
        #         experiment_name=run_settings.experiment_name,
        #     )

            # Perform propagation and enrichment for this condition
            await perform_propagation(
                run_settings,
                condition_settings,
                network=network,
                condition_data=condition_data_df,
            )

            # Collect condition_settings for later use
            all_condition_settings.append(
                (condition_settings, condition_data_df)
            )

        # Aggregate pathways across conditions for visualization
        all_pathways: Dict[str, Dict[str, Any]] = {}
        for condition_settings, condition_data_df in all_condition_settings:
            # Process condition for plotting
            process_condition(
                condition_settings=condition_settings,
                settings=run_settings,
                condition_data_df= condition_data_df,
                all_pathways=all_pathways,
            )

        # Plot aggregated pathway mean scores
        plot_output_path = plot_pathways_mean_scores(
            run_settings, all_pathways
        )

        # Collect output file paths
        output_files = []
        for condition_settings, _ in all_condition_settings:
            output_files.append(condition_settings.output_path)

        # Add the plot file
        output_files.append(plot_output_path)

        # Create a ZIP archive of the output files
        output_zip_dir = run_settings.root_folder / "Outputs" / "Archives"
        output_zip_dir.mkdir(parents=True, exist_ok=True)
        output_zip_path = output_zip_dir / f"{run_settings.experiment_name}_results.zip"

        with zipfile.ZipFile(output_zip_path, 'w') as zipf:
            for file_path in output_files:
                zipf.write(
                    file_path,
                    arcname=file_path.name
                )
        logger.info(f"Results archived at: {output_zip_path}")

        return output_zip_path

    except Exception as e:
        logger.error(f"An error occurred in pipeline_main: {e}")
        raise



async def main():
    """Initialize the main function for the pipeline."""
    try:
        output_zip_path = await pipeline_main(
            conditions_data=condition_data_files,
            settings=settings
        )
        logger.info(f"Results archived at: {output_zip_path}")
    except Exception as e:
        logger.error(f"An error occurred in main: {e}")
        raise


if __name__ == "__main__":
    try:
        # Changed the config_path to use a relative path
        config_path = Path(__file__).resolve().parent / "config.yaml"
        with config_path.open("r") as f:
            config = yaml.safe_load(f)
        
        settings = SettingsInput(**config)
        condition_data_files = [
            Path(file) for file in config["condition_data_files"]
        ]
        
        asyncio.run(main())
    except Exception as e:
        logger.error(f"An error occurred in __main__: {e}")
        raise

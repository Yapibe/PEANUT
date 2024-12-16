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
    debug: bool = False
) -> Path:
    """
    Main pipeline execution function for propagation and GSEA analyses.

    Parameters:
    - conditions_data (List[Tuple[str, pd.DataFrame]]): List of condition names and their data.
    - settings (SettingsInput): User-defined pipeline settings.
    - debug (bool): Whether to run in debug mode with different data handling.

    Returns:
    - Path: Path to the ZIP archive containing all results and plots.
    """
    try:
        # Initialize settings and load the network
        if debug:
            # Debug mode: Use settings directly
            run_settings = Settings(
                experiment_name=settings.experiment_name,
                species=settings.species,
                network=settings.network,
                create_similarity_matrix=settings.create_similarity_matrix,
                alpha=settings.alpha,
                pathway_file=settings.pathway_file,
                minimum_gene_per_pathway=settings.min_genes_per_pathway,
                maximum_gene_per_pathway=settings.max_genes_per_pathway,
                FDR_THRESHOLD=settings.fdr_threshold,
                figure_title=settings.figure_title
            )
        else:
            # Production mode: Process network path
            network_path = settings.network
            if settings.network == 'custom' and settings.network_file:
                # Save uploaded network file and use its path
                network_path = f"custom/{settings.network_file}"
                # Ensure matrix creation for custom networks
                create_matrix = True
            else:
                create_matrix = settings.create_similarity_matrix

            run_settings = Settings(
                experiment_name=settings.experiment_name,
                species=settings.species,
                network=network_path,
                create_similarity_matrix=create_matrix,
                alpha=settings.alpha,
                pathway_file=settings.pathway_file,
                minimum_gene_per_pathway=settings.min_genes_per_pathway,
                maximum_gene_per_pathway=settings.max_genes_per_pathway,
                FDR_THRESHOLD=settings.fdr_threshold,
                figure_title=settings.figure_title
            )
        network = read_network(run_settings.network_file_path)

        # Lists to hold settings for each condition
        all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []
        gsea_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []

        if debug:
            # Debug mode: Process conditions directly from data
            for condition_name, condition_data_df in conditions_data:
                # Propagation run
                condition_data_df.columns = ['GeneID', 'Score']
                propagation_settings = ConditionSettings(
                    condition_name=f"{condition_name}_Propagation",
                    experiment_name=run_settings.experiment_name,
                )
                await perform_propagation(
                    run_settings, propagation_settings, network, condition_data_df
                )
                all_condition_settings.append((propagation_settings, condition_data_df))

                # GSEA run (if enabled)
                if settings.run_gsea:
                    gsea_settings = ConditionSettings(
                        condition_name=f"{condition_name}_GSEA",
                        experiment_name=run_settings.experiment_name,
                    )
                    gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
                    await perform_propagation(
                        gsea_settings_run, gsea_settings, network, condition_data_df
                    )
                    gsea_condition_settings.append((gsea_settings, condition_data_df))
        else:
            # Production mode: Load and process files
            for condition_name, condition_path in conditions_data:
                logger.debug(f"Processing condition: {condition_name}, file: {condition_path}")
                if isinstance(condition_path, pd.DataFrame):
                    condition_data_df = condition_path
                else:
                    condition_data_df = pd.read_excel(condition_path)
                condition_settings = ConditionSettings(
                    condition_name=condition_name,
                    experiment_name=run_settings.experiment_name,
                )
                await perform_propagation(
                    run_settings, condition_settings, network, condition_data_df
                )
                all_condition_settings.append((condition_settings, condition_data_df))

                if settings.run_gsea:
                    gsea_settings = ConditionSettings(
                        condition_name=f"{condition_name}_GSEA",
                        experiment_name=run_settings.experiment_name,
                    )
                    gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
                    await perform_propagation(
                        gsea_settings_run, gsea_settings, network, condition_data_df
                    )
                    gsea_condition_settings.append((gsea_settings, condition_data_df))

        # Ensure output directory exists
        for condition_settings, _ in all_condition_settings + gsea_condition_settings:
            condition_settings.output_path.parent.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Ensured output directory exists: {condition_settings.output_path.parent}")

        # Aggregate and process results for propagation
        all_pathways: Dict[str, Dict[str, Any]] = {}
        for condition_settings, condition_data_df in all_condition_settings:
            process_condition(
                condition_settings=condition_settings,
                settings=run_settings,
                condition_data_df=condition_data_df,
                all_pathways=all_pathways,
            )
            logger.debug(f"Processed condition: {condition_settings.condition_name}, output path: {condition_settings.output_path}")

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
                logger.debug(f"Processed GSEA condition: {condition_settings.condition_name}, output path: {condition_settings.output_path}")

        # Generate plots
        propagation_plot_path = plot_pathways_mean_scores(run_settings, all_pathways, gsea_pathways)
        logger.debug(f"Generated propagation plot at: {propagation_plot_path}")

        # Collect all output files
        output_files = [settings.output_path for settings, _ in all_condition_settings + gsea_condition_settings]
        if propagation_plot_path:
            output_files.append(propagation_plot_path)

        # Log all output files
        for file_path in output_files:
            logger.debug(f"Output file to be archived: {file_path}")

        # Create ZIP archive for all results
        output_zip_dir = run_settings.root_folder / "Outputs" / "Archives"
        output_zip_dir.mkdir(parents=True, exist_ok=True)
        output_zip_path = output_zip_dir / f"{run_settings.experiment_name}_results.zip"

        with zipfile.ZipFile(output_zip_path, 'w') as zipf:
            for file_path in output_files:
                if file_path.exists():
                    zipf.write(file_path, arcname=file_path.name)
                    logger.debug(f"Archived file: {file_path}")
                else:
                    logger.error(f"File not found, skipping: {file_path}")

        logger.info(f"Results archived at: {output_zip_path}")
        return output_zip_path

    except Exception as e:
        logger.error(f"An error occurred in pipeline_main: {e}", exc_info=True)
        raise


async def main():
    """Entry point for the pipeline."""
    try:
        # Load configuration from YAML file
        config_path = Path(__file__).resolve().parent / "config.yaml"
        with config_path.open("r") as f:
            config = yaml.safe_load(f)

        # Process network configuration
        if config['network_type'] == 'custom':
            network = f"custom/{config['network_name']}"
        else:
            network = config['network_name']

        # Create settings from config
        settings = SettingsInput(
            experiment_name=config['experiment_name'],
            species=config['species'],
            alpha=config['alpha'],
            network=network,  # Use processed network name
            pathway_file=config['pathway_file'],
            min_genes_per_pathway=config['min_genes_per_pathway'],
            max_genes_per_pathway=config['max_genes_per_pathway'],
            fdr_threshold=config['fdr_threshold'],
            run_gsea=config['run_gsea'],
            restrict_to_network=config['restrict_to_network'],
            create_similarity_matrix=config['create_similarity_matrix'],
            figure_title=config['figure_title']
        )

        # Load condition data files
        conditions_data = []
        for file_path in config['condition_data_files']:
            file_path = Path(file_path)
            if not file_path.exists():
                raise FileNotFoundError(f"Input file not found: {file_path}")
            
            df = pd.read_excel(file_path) if file_path.suffix == '.xlsx' else pd.read_csv(file_path, sep='\t')
            conditions_data.append((file_path.stem, df))

        # For debug mode with custom network, ensure the file exists
        if config['network_type'] == 'custom':
            network_path = Path(config['network_path'])
            if not network_path.exists():
                raise FileNotFoundError(f"Network file not found: {network_path}")
            
            # Copy network file to the correct location if needed
            target_path = Path(__file__).resolve().parent / "Data" / config['species'] / "network" / "custom" / f"{config['network_name']}.txt"
            target_path.parent.mkdir(parents=True, exist_ok=True)
            
            if not target_path.exists():
                import shutil
                shutil.copy2(network_path, target_path)
                logger.info(f"Copied network file to: {target_path}")

        # Run pipeline in debug mode
        output_zip_path = await pipeline_main(conditions_data, settings, debug=True)
        logger.info(f"Results successfully archived at: {output_zip_path}")
    except Exception as e:
        logger.error(f"An error occurred in main: {e}")
        raise


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except Exception as e:
        logger.error(f"An error occurred in __main__: {e}")
        raise
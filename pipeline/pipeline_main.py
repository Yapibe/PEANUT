"""
Main pipeline module for the PEANUT (Pathway Enrichment Analysis using Network Propagation) system.

This module contains the main pipeline functions for running propagation analyses,
managing conditions, and creating output files.
"""

import asyncio
import logging
import zipfile
from pathlib import Path
from typing import List, Tuple, Dict, Any, Set, Optional, Union
import pandas as pd
from dataclasses import replace
import sys
import networkx as nx

# Configure logger
logger = logging.getLogger(__name__)

# Update Python path to include parent directory
current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent
sys.path.append(str(parent_dir))

# Import necessary modules
from app.models import SettingsInput
from pipeline.propagation_routines import perform_propagation
from pipeline.settings import ConditionSettings, Settings
from pipeline.utils import process_condition, read_network, load_config
from pipeline.visualization_tools import plot_pathways_mean_scores


async def pipeline_main(
    conditions_data: List[Tuple[str, pd.DataFrame]],
    settings: SettingsInput,
    debug: bool = False
) -> Optional[Path]:
    """
    Main pipeline execution function for propagation and GSEA analyses.

    Parameters:
    - conditions_data (List[Tuple[str, pd.DataFrame]]): List of condition names and their data.
    - settings (SettingsInput): User-defined pipeline settings.
    - debug (bool): Whether to run in debug mode with different data handling.

    Returns:
    - Optional[Path]: Path to the ZIP archive containing all results and plots, or None if an error occurred.
    """
    try:
        logger.info(f"Starting pipeline for experiment: {settings.experiment_name}")
        logger.info(f"Running with settings: alpha={settings.alpha}, network={settings.network}")
        
        # Initialize settings and load the network
        run_settings = _initialize_settings(settings, debug)
        
        # Load network
        try:
            network = read_network(run_settings.network_file_path)
            logger.info(f"Successfully loaded network with {len(network.nodes())} nodes")
        except Exception as e:
            logger.error(f"Failed to load network from {run_settings.network_file_path}: {e}")
            raise
            
        # Check for valid conditions
        if not conditions_data:
            logger.error("No condition data provided. Exiting pipeline.")
            raise ValueError("No condition data provided")
        
        logger.info(f"Processing {len(conditions_data)} conditions")

        # Initialize data containers
        all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []
        gsea_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]] = []
        pathway_sizes: Dict[str, int] = {}
        genes_by_pathway: Dict[str, Set[str]] = {}
        scores: Dict[str, float] = {}

        # Process all conditions
        if debug:
            await _process_debug_conditions(conditions_data, run_settings, network, settings,
                                           all_condition_settings, gsea_condition_settings)
        else:
            await _process_production_conditions(conditions_data, run_settings, network, settings,
                                                all_condition_settings, gsea_condition_settings)
            
        # Create output directories
        for condition_settings, _ in all_condition_settings + gsea_condition_settings:
            condition_settings.output_path.parent.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Ensured output directory exists for {condition_settings.condition_name}")

        # Aggregate and process results for propagation
        all_pathways: Dict[str, Dict[str, Any]] = {}
        logger.info("Processing propagation results")
        for condition_settings, condition_data_df in all_condition_settings:
            process_condition(
                condition_settings=condition_settings,
                settings=run_settings,
                condition_data_df=condition_data_df,
                all_pathways=all_pathways,
                pathway_sizes=pathway_sizes,
                genes_by_pathway=genes_by_pathway,
                scores=scores,
            )
            logger.info(f"Processed propagation for {condition_settings.condition_name}")

        # Aggregate and process results for GSEA
        gsea_pathways: Dict[str, Dict[str, Any]] = {}
        if settings.run_gsea:
            logger.info("Processing GSEA results")
            for condition_settings, condition_data_df in gsea_condition_settings:
                process_condition(
                    condition_settings=condition_settings,
                    settings=run_settings,
                    condition_data_df=condition_data_df,
                    all_pathways=gsea_pathways,
                    pathway_sizes=pathway_sizes,
                    genes_by_pathway=genes_by_pathway,
                    scores=scores,
                )
                logger.info(f"Processed GSEA for {condition_settings.condition_name}")

        # Generate plots
        logger.info("Generating pathway plots")
        propagation_plot_path = plot_pathways_mean_scores(
            settings=run_settings,
            all_pathways=all_pathways,
            gsea_pathways=gsea_pathways,
            pathway_sizes=pathway_sizes,
            genes_by_pathway=genes_by_pathway,
            scores=scores,
        )
        
        if propagation_plot_path:
            logger.info(f"Generated propagation plot at: {propagation_plot_path}")
        else:
            logger.warning("No plot was generated")

        # Collect all output files
        output_files = [settings.output_path for settings, _ in all_condition_settings + gsea_condition_settings]
        if propagation_plot_path and propagation_plot_path.exists():
            output_files.append(propagation_plot_path)

        # Create ZIP archive for all results
        output_zip_path = _create_results_archive(run_settings, output_files)
        
        logger.info(f"Pipeline completed successfully. Results archived at: {output_zip_path}")
        return output_zip_path

    except Exception as e:
        logger.error(f"An error occurred in pipeline_main: {e}", exc_info=True)
        raise  # Re-raise the exception instead of returning None


def _initialize_settings(settings: SettingsInput, debug: bool) -> Settings:
    """
    Initialize Settings object based on provided input and mode.
    
    Parameters:
    - settings (SettingsInput): User-defined settings
    - debug (bool): Whether to run in debug mode
    
    Returns:
    - Settings: Initialized Settings object
    """
    try:
        if debug:
            # Debug mode: Use settings directly
            run_settings = Settings(
                experiment_name=settings.experiment_name,
                species=settings.species,
                network=settings.network,
                create_similarity_matrix=settings.create_similarity_matrix,
                alpha=settings.alpha,
                pathway_file=settings.pathway_file,
                min_genes_per_pathway=settings.min_genes_per_pathway,
                max_genes_per_pathway=settings.max_genes_per_pathway,
                fdr_threshold=settings.fdr_threshold,
                figure_title=settings.figure_title,
                run_gsea=settings.run_gsea,
                restrict_to_network=settings.restrict_to_network
            )
        else:
            # Production mode: Process network path
            network_path = settings.network
            create_matrix = settings.create_similarity_matrix
            
            if settings.network == 'custom' and settings.network_file:
                # Save uploaded network file and use its path
                network_path = f"custom/{settings.network_file.filename}"
                create_matrix = True
                logger.info(f"Using custom network from {network_path}")

            run_settings = Settings(
                experiment_name=settings.experiment_name,
                species=settings.species,
                network=network_path,
                create_similarity_matrix=create_matrix,
                alpha=settings.alpha,
                pathway_file=settings.pathway_file,
                min_genes_per_pathway=settings.min_genes_per_pathway,
                max_genes_per_pathway=settings.max_genes_per_pathway,
                fdr_threshold=settings.fdr_threshold,
                figure_title=settings.figure_title,
                run_gsea=settings.run_gsea,
                restrict_to_network=settings.restrict_to_network
            )
            
        logger.info(f"Settings initialized: experiment={run_settings.experiment_name}, network={run_settings.network}")
        return run_settings
    except Exception as e:
        logger.error(f"Error initializing settings: {e}", exc_info=True)
        raise


async def _process_debug_conditions(
    conditions_data: List[Tuple[str, pd.DataFrame]],
    run_settings: Settings,
    network: nx.Graph,
    settings: SettingsInput,
    all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]],
    gsea_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]]
) -> None:
    """
    Process conditions in debug mode.
    
    Parameters:
    - conditions_data: List of condition names and dataframes
    - run_settings: Settings object
    - network: Network graph
    - settings: User input settings
    - all_condition_settings: List to store condition settings for propagation
    - gsea_condition_settings: List to store condition settings for GSEA
    """
    for condition_name, condition_data_df in conditions_data:
        logger.info(f"Processing debug condition: {condition_name}")
        
        # Ensure correct column names for debug data
        condition_data_df.columns = ['GeneID', 'Score']

        # Process propagation
        propagation_settings = ConditionSettings(
            condition_name=f"{condition_name}_Propagation",
            experiment_name=run_settings.experiment_name,
        )
        await perform_propagation(
            run_settings, propagation_settings, network, condition_data_df
        )
        all_condition_settings.append((propagation_settings, condition_data_df))

        # Process GSEA (if enabled)
        if settings.run_gsea:
            gsea_settings = ConditionSettings(
                condition_name=f"{condition_name}_GSEA",
                experiment_name=run_settings.experiment_name,
            )
            # Create a copy of settings with GSEA enabled
            gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
            await perform_propagation(
                gsea_settings_run, gsea_settings, network, condition_data_df
            )
            gsea_condition_settings.append((gsea_settings, condition_data_df))


async def _process_production_conditions(
    conditions_data: List[Tuple[str, Union[str, pd.DataFrame]]],
    run_settings: Settings,
    network: nx.Graph,
    settings: SettingsInput,
    all_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]],
    gsea_condition_settings: List[Tuple[ConditionSettings, pd.DataFrame]]
) -> None:
    """
    Process conditions in production mode.
    
    Parameters:
    - conditions_data: List of condition names and paths/dataframes
    - run_settings: Settings object
    - network: Network graph
    - settings: User input settings
    - all_condition_settings: List to store condition settings for propagation
    - gsea_condition_settings: List to store condition settings for GSEA
    """
    for condition_name, condition_path in conditions_data:
        logger.info(f"Processing condition: {condition_name}")
        
        # Load condition data if it's a path
        if isinstance(condition_path, pd.DataFrame):
            condition_data_df = condition_path
        else:
            try:
                path = Path(condition_path)
                if path.suffix.lower() == '.xlsx':
                    condition_data_df = pd.read_excel(path)
                elif path.suffix.lower() in ('.csv', '.txt', '.tsv'):
                    condition_data_df = pd.read_csv(path, sep='\t' if path.suffix.lower() == '.tsv' else ',')
                else:
                    logger.warning(f"Unrecognized file extension: {path.suffix}. Attempting to read as Excel.")
                    condition_data_df = pd.read_excel(path)
                    
                logger.info(f"Loaded condition data from {path} with {len(condition_data_df)} entries")
            except Exception as e:
                logger.error(f"Error loading condition data from {condition_path}: {e}", exc_info=True)
                continue

        # Process propagation
        condition_settings = ConditionSettings(
            condition_name=condition_name,
            experiment_name=run_settings.experiment_name,
        )
        await perform_propagation(
            run_settings, condition_settings, network, condition_data_df
        )
        all_condition_settings.append((condition_settings, condition_data_df))

        # Process GSEA (if enabled)
        if settings.run_gsea:
            gsea_settings = ConditionSettings(
                condition_name=f"{condition_name}_GSEA",
                experiment_name=run_settings.experiment_name,
            )
            # Create a copy of settings with GSEA enabled
            gsea_settings_run = replace(run_settings, run_gsea=True, restrict_to_network=settings.restrict_to_network)
            await perform_propagation(
                gsea_settings_run, gsea_settings, network, condition_data_df
            )
            gsea_condition_settings.append((gsea_settings, condition_data_df))


def _create_results_archive(run_settings: Settings, output_files: List[Path]) -> Path:
    """
    Create a ZIP archive containing all result files.
    
    Parameters:
    - run_settings: Settings object
    - output_files: List of files to include in the archive
    
    Returns:
    - Path to the created ZIP archive
    """
    output_zip_dir = run_settings.root_folder / "Outputs" / "Archives"
    output_zip_dir.mkdir(parents=True, exist_ok=True)
    output_zip_path = output_zip_dir / f"{run_settings.experiment_name}_results.zip"

    with zipfile.ZipFile(output_zip_path, 'w') as zipf:
        files_added = 0
        for file_path in output_files:
            if file_path.exists():
                zipf.write(file_path, arcname=file_path.name)
                logger.debug(f"Added file to archive: {file_path}")
                files_added += 1
            else:
                logger.warning(f"File not found, skipping: {file_path}")
                
    logger.info(f"Created archive with {files_added} files at {output_zip_path}")
    return output_zip_path


async def main():
    """Entry point for the pipeline when run as a script."""
    # Configure basic logging 
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    try:
        logger.info("Starting PEANUT pipeline from configuration file")
        
        # Load configuration from YAML file
        config = load_config()
            
        logger.info(f"Loaded configuration for experiment: {config['experiment_name']}")

        # Process network configuration
        if config['network_type'] == 'custom':
            network = f"custom/{config['network_name']}"
            logger.info(f"Using custom network: {network}")
        else:
            network = config['network_name']
            logger.info(f"Using standard network: {network}")

        # Create settings from config
        settings = SettingsInput(
            experiment_name=config['experiment_name'],
            species=config['species'],
            alpha=config['alpha'],
            network=network,
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
        for file_path_str in config['condition_data_files']:
            file_path = Path(file_path_str)
            if not file_path.exists():
                logger.error(f"Input file not found: {file_path}")
                continue
                
            try:
                if file_path.suffix.lower() == '.xlsx':
                    df = pd.read_excel(file_path)
                else:
                    df = pd.read_csv(file_path, sep='\t')
                    
                conditions_data.append((file_path.stem, df))
                logger.info(f"Loaded condition from {file_path}: {len(df)} entries")
            except Exception as e:
                logger.error(f"Error loading condition file {file_path}: {e}")
                continue
                
        if not conditions_data:
            logger.error("No valid condition data files found. Exiting.")
            return

        # For debug mode with custom network, ensure the file exists
        if config['network_type'] == 'custom':
            network_path = Path(config['network_path'])
            if not network_path.exists():
                logger.error(f"Network file not found: {network_path}")
                return
                
            # Copy network file to the correct location if needed
            target_path = Path(__file__).resolve().parent / "Data" / config['species'] / "network" / "custom" / f"{config['network_name']}.txt"
            target_path.parent.mkdir(parents=True, exist_ok=True)
            
            if not target_path.exists():
                import shutil
                shutil.copy2(network_path, target_path)
                logger.info(f"Copied network file to: {target_path}")

        # Run pipeline in debug mode
        logger.info("Running pipeline in debug mode")
        output_zip_path = await pipeline_main(conditions_data, settings, debug=True)
        
        if output_zip_path:
            logger.info(f"Results successfully archived at: {output_zip_path}")
        else:
            logger.error("Pipeline failed to complete successfully")
            
    except Exception as e:
        logger.error(f"An error occurred in main: {e}", exc_info=True)


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except Exception as e:
        logger.error(f"An error occurred in __main__: {e}", exc_info=True)
        sys.exit(1)
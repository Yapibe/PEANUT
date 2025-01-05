import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from .temp_manager import TempFileManager
from .models import SettingsInput
from pipeline.settings import Settings
from pipeline.pipeline_main import pipeline_main as execute_pipeline
import pandas as pd
import traceback

logger = logging.getLogger(__name__)

class PipelineRunner:
    """Handles pipeline execution and temporary file management."""

    def __init__(self, job_id: str):
        self.job_id = job_id
        self.temp_manager = TempFileManager()
        logger.info(f"=== Initialized PipelineRunner for job {job_id} ===")

    async def run(self, settings_input: SettingsInput, file_dfs: List[pd.DataFrame], filenames: List[str]) -> Dict[str, Any]:
        logger.info(f"Starting pipeline for job {self.job_id}")
        try:
            # Detect if temporary files are needed
            use_temp_files = settings_input.is_custom_network() or settings_input.is_custom_pathway()
            logger.info(f"Using temporary files: {use_temp_files}")

            if use_temp_files:
                logger.debug("Creating temporary files")
                with self.temp_manager.create_temp_files(
                    network_filename=settings_input.network_file.filename if settings_input.is_custom_network() else None,
                    pathway_filename=settings_input.pathway_file_upload.filename if settings_input.is_custom_pathway() else None
                ) as temp_mgr:
                    logger.debug("Temporary files created successfully")

                    # Handle custom network file
                    if settings_input.is_custom_network():
                        logger.debug("Handling custom network file")
                        network_content = await settings_input.network_file.read()
                        with open(temp_mgr.network_path, 'wb') as f:
                            f.write(network_content)
                        settings_input.network = str(temp_mgr.network_path)
                        logger.info(f"Custom network saved to {temp_mgr.network_path}")

                    # Handle custom pathway file
                    if settings_input.is_custom_pathway():
                        logger.debug("Handling custom pathway file")
                        pathway_content = await settings_input.pathway_file_upload.read()
                        with open(temp_mgr.pathway_file_path, 'wb') as f:
                            f.write(pathway_content)
                        settings_input.pathway_file = str(temp_mgr.pathway_file_path)
                        logger.info(f"Custom pathway saved to {temp_mgr.pathway_file_path}")

                    # Execute the pipeline
                    logger.debug("Executing pipeline with temporary files")
                    output_file = await execute_pipeline(
                        [(name, df) for name, df in zip(filenames, file_dfs)],
                        settings_input
                    )
                    logger.info(f"Pipeline execution completed with output: {output_file}")

            else:
                # Run pipeline without temporary files
                logger.debug("Executing pipeline with default files")
                output_file = await execute_pipeline(
                    [(name, df) for name, df in zip(filenames, file_dfs)],
                    settings_input
                )
                logger.info(f"Pipeline execution completed with output: {output_file}")

            # Ensure the result is returned as a dictionary
            return {
                "output_file": str(output_file),  # Convert Path object to string
                "status": "success"
            }

        except Exception as e:
            logger.error(f"Pipeline execution failed for job {self.job_id}")
            logger.error(f"Error details: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
        finally:
            logger.debug("Ensuring cleanup after pipeline execution")
            self.temp_manager.cleanup()

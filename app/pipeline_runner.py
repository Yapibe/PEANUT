import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
from .temp_manager import TempFileManager
from .models import SettingsInput
from pipeline.settings import Settings
from pipeline.pipeline_main import pipeline_main as execute_pipeline
import pandas as pd

logger = logging.getLogger(__name__)

class PipelineRunner:
    """Handles pipeline execution and temporary file management."""
    
    def __init__(self, job_id: str):
        self.job_id = job_id
        self.temp_manager = TempFileManager()
        
    async def run(self, settings_input: SettingsInput, file_dfs: List[pd.DataFrame], filenames: List[str]) -> Dict[str, Any]:
        """
        Runs the pipeline with proper temporary file management.
        
        Args:
            settings_input: Input settings for the pipeline
            file_dfs: List of DataFrames containing input data
            filenames: List of input file names
            
        Returns:
            Dictionary containing job results
            
        Raises:
            Exception: If pipeline execution fails
        """
        try:
            logger.info(f"Starting pipeline for job {self.job_id}")
            
            # Handle custom network if provided
            if settings_input.is_custom_network():
                logger.debug("Custom network detected, preparing to handle file")
                
                # Ensure the file is open and seek to the beginning
                if not settings_input.network_file.file.closed:
                    settings_input.network_file.file.seek(0)
                else:
                    logger.error("Network file is closed before reading")
                    raise ValueError("Network file is closed before reading")
                
                # Read file content early to avoid "closed file" issues
                logger.debug("Reading network file content")
                try:
                    network_content = await settings_input.network_file.read()
                    logger.debug(f"Successfully read network file content, size: {len(network_content)} bytes")
                except Exception as e:
                    logger.error(f"Failed to read network file: {str(e)}", exc_info=True)
                    raise ValueError(f"Failed to read network file: {str(e)}")
                
                logger.debug("Creating temporary files")
                with self.temp_manager.create_temp_files(settings_input.network_file.filename) as temp_mgr:
                    logger.debug(f"Temporary files created at: {temp_mgr.network_path}")
                    
                    # Save network content to temporary path
                    logger.debug("Writing network content to temporary file")
                    try:
                        with open(temp_mgr.network_path, 'wb') as f:
                            f.write(network_content)
                        logger.info(f"Saved custom network to {temp_mgr.network_path}")
                    except Exception as e:
                        logger.error(f"Failed to write network file: {str(e)}", exc_info=True)
                        raise ValueError(f"Failed to write network file: {str(e)}")
                
                    # Update network path in settings
                    settings_input.network = str(temp_mgr.network_path)
                    logger.debug(f"Updated network path in settings: {settings_input.network}")
                
                    # Run pipeline while temporary files are valid
                    logger.debug("Starting pipeline execution with custom network")
                    results = await execute_pipeline(
                        [(name, df) for name, df in zip(filenames, file_dfs)],
                        settings_input
                    )
                    logger.debug("Pipeline execution completed")
            else:
                # Run pipeline with default network
                logger.debug("Using default network")
                results = await execute_pipeline(
                    [(name, df) for name, df in zip(filenames, file_dfs)],
                    settings_input
                )
            
            logger.info(f"Pipeline completed successfully for job {self.job_id}")
            return results
            
        except Exception as e:
            logger.error(f"Pipeline failed for job {self.job_id}: {str(e)}", exc_info=True)
            raise
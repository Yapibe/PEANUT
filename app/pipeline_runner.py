from typing import Optional, Dict, Any, List, AsyncGenerator
from contextlib import asynccontextmanager
from .temp_manager import TempFileManager
from .models import SettingsInput
from pipeline.pipeline_main import pipeline_main as execute_pipeline
import pandas as pd
import traceback
from datetime import datetime
from .utils import setup_logging
import asyncio

logger = setup_logging()

class PipelineRunner:
    """Handles pipeline execution and temporary file management for PEANUT.
    
    This class manages the execution of the PEANUT pipeline, including:
    - Handling custom network and pathway files
    - Managing temporary files
    - Executing the pipeline with proper error handling
    - Cleaning up resources
    
    Attributes:
        job_id: Unique identifier for the pipeline run
        temp_manager: Manages temporary files
        start_time: When the pipeline started
        end_time: When the pipeline completed
    """
    
    def __init__(self, job_id: str):
        """Initialize the pipeline runner.
        
        Args:
            job_id: Unique identifier for the pipeline run
        """
        self.job_id = job_id
        self.temp_manager = TempFileManager()
        self.start_time: Optional[datetime] = None
        self.end_time: Optional[datetime] = None
        logger.info(f"Initialized PipelineRunner for job {job_id}")

    @asynccontextmanager
    async def _handle_temp_files(self, settings_input: SettingsInput) -> AsyncGenerator[None, None]:
        """Context manager for handling temporary files.
        
        Args:
            settings_input: Pipeline settings that may require temporary files
            
        Yields:
            None
        """
        use_temp_files = settings_input.is_custom_network() or settings_input.is_custom_pathway()
        logger.info(f"Using temporary files: {use_temp_files}")

        if not use_temp_files:
            yield
            return

        try:
            # For path-based files, we'll use those directly
            if settings_input.network_file_path:
                logger.debug(f"Using provided network file path: {settings_input.network_file_path}")
                settings_input.network = settings_input.network_file_path
                
            if settings_input.pathway_file_path:
                logger.debug(f"Using provided pathway file path: {settings_input.pathway_file_path}")
                settings_input.pathway_file = settings_input.pathway_file_path
                
            # For UploadFile objects, we need to create temporary files
            if settings_input.network_file or settings_input.pathway_file_upload:
                with self.temp_manager.create_temp_files(
                    network_filename=settings_input.network_file.filename if settings_input.network_file else None,
                    pathway_filename=settings_input.pathway_file_upload.filename if settings_input.pathway_file_upload else None
                ) as temp_mgr:
                    # Handle custom network file
                    if settings_input.network_file:
                        logger.debug("Handling custom network file from upload")
                        network_content = await settings_input.network_file.read()
                        with open(temp_mgr.network_path, 'wb') as f:
                            f.write(network_content)
                        settings_input.network = str(temp_mgr.network_path)
                        logger.info(f"Custom network saved to {temp_mgr.network_path}")

                    # Handle custom pathway file
                    if settings_input.pathway_file_upload:
                        logger.debug("Handling custom pathway file from upload")
                        pathway_content = await settings_input.pathway_file_upload.read()
                        with open(temp_mgr.pathway_file_path, 'wb') as f:
                            f.write(pathway_content)
                        settings_input.pathway_file = str(temp_mgr.pathway_file_path)
                        logger.info(f"Custom pathway saved to {temp_mgr.pathway_file_path}")

            yield
        except Exception as e:
            logger.error(f"Error handling temporary files: {str(e)}")
            raise

    async def run(self, settings_input: SettingsInput, file_dfs: List[pd.DataFrame], filenames: List[str]) -> Dict[str, Any]:
        """Run the PEANUT pipeline with the given settings and input files.
        
        Args:
            settings_input: Pipeline settings
            file_dfs: List of input DataFrames
            filenames: List of corresponding filenames
            
        Returns:
            Dictionary containing pipeline results
            
        Raises:
            Exception: If pipeline execution fails
        """
        self.start_time = datetime.now()
        logger.info(f"Starting pipeline for job {self.job_id}")
        
        try:
            # Log matrix generation needs
            if settings_input.needs_matrix_generation():
                if settings_input.is_custom_network():
                    logger.info(f"Custom network detected - matrix generation will take approximately 120 minutes")
                else:
                    logger.info(f"Custom alpha {settings_input.alpha} with default network - matrix generation will take approximately 120 minutes")
            
            # Set create_similarity_matrix flag based on needs
            if settings_input.needs_matrix_generation():
                logger.info("Setting create_similarity_matrix=True for computation")
                settings_input.create_similarity_matrix = True
            
            # Real pipeline execution
            async with self._handle_temp_files(settings_input):
                logger.debug("Executing pipeline")
                output_file = await execute_pipeline(
                    [(name, df) for name, df in zip(filenames, file_dfs)],
                    settings_input
                )
                logger.info(f"Pipeline execution completed with output: {output_file}")

            self.end_time = datetime.now()
            duration = (self.end_time - self.start_time).total_seconds()
            logger.info(f"Pipeline completed in {duration:.2f} seconds")

            return {
                "output_file": str(output_file),
                "status": "success",
                "duration_seconds": duration
            }

        except Exception as e:
            self.end_time = datetime.now()
            duration = (self.end_time - self.start_time).total_seconds()
            logger.error(f"Pipeline execution failed for job {self.job_id} after {duration:.2f} seconds")
            logger.error(f"Error details: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            raise
        finally:
            logger.debug("Ensuring cleanup after pipeline execution")
            self.temp_manager.cleanup()

    def get_execution_time(self) -> Optional[float]:
        """Get the pipeline execution time in seconds.
        
        Returns:
            Execution time in seconds, or None if not completed
        """
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return None

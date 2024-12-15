import logging
from typing import List
import pandas as pd
from pathlib import Path

from .models import SettingsInput
from pipeline.pipeline_main import pipeline_main

logger = logging.getLogger(__name__)


async def run_pipeline(
    file_dfs: List[pd.DataFrame],
    filenames: List[str],
    settings: SettingsInput,
    job_storage: dict,
    job_code: str,
):
    """
    Run the pipeline with the given settings.
    Parameters:
    - file_dfs: List of validated DataFrames (one for each input file).
    - filenames: List of filenames corresponding to the DataFrames.
    - settings: User-provided settings for the pipeline.
    - job_storage: Dictionary to track job statuses.
    - job_code: Unique code for the job.
    """
    try:
        # Prepare conditions_data
        conditions_data = []
        for filename, df in zip(filenames, file_dfs):
            condition_name = Path(filename).stem
            conditions_data.append((condition_name, df))

        # Execute pipeline
        output_file_path = await pipeline_main(conditions_data, settings)

        logger.debug(f"Pipeline main returned output_file_path: {output_file_path}")

        # Handle relative path
        relative_output_path = Path(output_file_path).relative_to(Path.cwd())

        # Update job status
        job_storage[job_code]["status"] = "Finished"
        job_storage[job_code]["output_file"] = str(relative_output_path)

        logger.info(f"Job {job_code} finished successfully.")

    except Exception as e:
        logger.error(f"Error in pipeline execution for job_code {job_code}: {str(e)}", exc_info=True)
        job_storage[job_code]["status"] = "Failed"
        job_storage[job_code]["error"] = str(e)
        raise

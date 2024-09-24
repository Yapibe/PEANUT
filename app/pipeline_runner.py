import logging
from typing import List
from pathlib import Path

from .models import SettingsInput
from .utils import validate_file
from pipeline.pipeline_main import pipeline_main

logger = logging.getLogger(__name__)


async def run_pipeline(
    files: List[bytes],
    filenames: List[str],
    settings: SettingsInput,
    job_storage: dict,
    job_code: str,
):
    """Run the pipeline with the given settings."""
    try:
        conditions_data = []
        for file_content, filename in zip(files, filenames):
            df = await validate_file(file_content, filename)
            condition_name = Path(filename).stem
            conditions_data.append((condition_name, df))

        # Actual pipeline execution
        output_file_path = await pipeline_main(conditions_data, settings)

        # Update job status to Finished
        job_storage[job_code]["status"] = "Finished"
        job_storage[job_code]["output_file"] = str(output_file_path)

        logger.info(f"Job {job_code} finished successfully.")

    except Exception as e:
        logger.error(f"Error in pipeline execution: {str(e)}")
        job_storage[job_code]["status"] = "Failed"
        job_storage[job_code]["error"] = str(e)
        raise

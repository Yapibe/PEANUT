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

        # Actual pipeline execution using the loaded configuration
        output_file_path = await pipeline_main(conditions_data, settings)

        logger.debug(f"Pipeline main returned output_file_path: {output_file_path}")

        # Update job status to Finished with relative path
        try:
            # Ensure the output_file_path is relative to Path.cwd()
            relative_output_path = Path(output_file_path).relative_to(Path.cwd())
            logger.debug(f"Relative output path: {relative_output_path}")
        except ValueError:
            # If output_file_path is not under Path.cwd(), handle accordingly
            # This removes the leading slash if present
            relative_output_path = Path(output_file_path).as_posix().lstrip('/')
            logger.debug(f"Relative output path after lstrip: {relative_output_path}")

        job_storage[job_code]["status"] = "Finished"
        job_storage[job_code]["output_file"] = str(relative_output_path)

        logger.info(f"Job {job_code} finished successfully.")

    except Exception as e:
        logger.error(f"Error in pipeline execution for job_code {job_code}: {str(e)}", exc_info=True)
        job_storage[job_code]["status"] = "Failed"
        job_storage[job_code]["error"] = str(e)
        raise


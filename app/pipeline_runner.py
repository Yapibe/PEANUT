from .models import PipelineOutput, SettingsInput
from .utils import validate_file
from pipeline.pipeline_main import pipeline_main
import pandas as pd
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


async def run_pipeline(
    file_content: bytes, filename: str, settings: SettingsInput, job_storage, job_code
):
    try:
        # Validate and convert the uploaded file to a DataFrame
        df = await validate_file(file_content, filename)

        # Actual pipeline execution
        output_file_path = await pipeline_main(df, settings)

        # Update job status to Finished
        job_storage[job_code]["status"] = "Finished"
        job_storage[job_code]["output_file"] = output_file_path

        return output_file_path
    except Exception as e:
        logger.error(f"Error in pipeline execution: {str(e)}")
        job_storage[job_code]["status"] = "Failed"
        job_storage[job_code]["error"] = str(e)
        raise

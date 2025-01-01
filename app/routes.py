import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
from fastapi import APIRouter, UploadFile, HTTPException, BackgroundTasks, Request, Form, File
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.templating import Jinja2Templates
import pandas as pd
import uuid
from datetime import datetime
import tempfile

from .models import SettingsInput, JobStatus
from .utils import validate_network_file, validate_input_file, convert_gmt_to_pathway_format
from .pipeline_runner import PipelineRunner
from .temp_manager import TempFileManager

logger = logging.getLogger(__name__)
router = APIRouter()

# Set up templates
templates = Jinja2Templates(directory="app/templates")

# In-memory job storage
job_storage: Dict[str, Dict] = {}

@router.post("/run-pipeline")
async def run_pipeline(
    background_tasks: BackgroundTasks,
    experiment_name: str = Form(...),
    species: str = Form(...),
    alpha: float = Form(...),
    network: str = Form(...),
    pathway_file: str = Form(...),
    pathway_file_upload: Optional[UploadFile] = None,
    min_genes_per_pathway: Optional[int] = Form(15),
    max_genes_per_pathway: Optional[int] = Form(500),
    fdr_threshold: float = Form(0.05),
    run_gsea: bool = Form(False),
    restrict_to_network: bool = Form(False),
    plot_title: Optional[str] = Form("Pathway Enrichment"),
    network_file: Optional[UploadFile] = None,
    files: List[UploadFile] = File(...),
):
    """
    Endpoint to run the pipeline with the provided settings and input files.
    """
    try:
        # Log initial inputs
        logger.info("New job request received.")
        logger.info(f"Pathway file: {pathway_file}")
        logger.info(f"Custom pathway uploaded: {pathway_file_upload is not None and pathway_file_upload.filename}")

        # Generate unique job code
        job_code = str(uuid.uuid4())
        logger.info(f"Job code: {job_code}")

        # Initialize job storage
        job_storage[job_code] = {
            "status": "Processing",
            "output_file": None,
            "error": None
        }

        # Read and store content of custom network and pathway files in memory
        if network == 'custom' and network_file:
            logger.debug("Processing custom network file")
            network_content = await network_file.read()
            network_file = UploadFile(
                filename=network_file.filename,
                file=pd.io.common.BytesIO(network_content)
            )

        if pathway_file == 'custom' and pathway_file_upload:
            logger.debug("Processing custom pathway file")
            pathway_content = await pathway_file_upload.read()
            pathway_file_upload = UploadFile(
                filename=pathway_file_upload.filename,
                file=pd.io.common.BytesIO(pathway_content)
            )

        # Validate and process input files
        file_dfs, filenames = [], []
        for file in files:
            file_type = file.filename.split('.')[-1].lower()
            if file_type not in ['rnk', 'xlsx']:
                raise HTTPException(status_code=400, detail=f"Unsupported file type: {file_type}")

            content = await file.read()
            temp_path = Path(f"/tmp/{file.filename}")
            temp_path.write_bytes(content)
            try:
                df = validate_input_file(temp_path, file_type)
                file_dfs.append(df)
                filenames.append(file.filename)
            finally:
                temp_path.unlink()  # Clean up temporary files

        # Create the pipeline runner
        runner = PipelineRunner(job_code)

        # Run the pipeline in the background
        background_tasks.add_task(
            _run_pipeline_task,
            runner=runner,
            settings=SettingsInput(
                experiment_name=experiment_name,
                species=species,
                alpha=alpha,
                network=network,
                network_file=network_file,
                pathway_file=pathway_file,
                pathway_file_upload=pathway_file_upload,
                min_genes_per_pathway=min_genes_per_pathway,
                max_genes_per_pathway=max_genes_per_pathway,
                fdr_threshold=fdr_threshold,
                run_gsea=run_gsea,
                restrict_to_network=restrict_to_network,
                figure_title=plot_title,
            ),
            file_dfs=file_dfs,
            filenames=filenames,
            job_storage=job_storage,
            job_code=job_code
        )

        return {"job_code": job_code}

    except Exception as e:
        logger.error(f"Error processing pipeline request: {e}")
        raise HTTPException(status_code=422, detail=str(e))


async def _run_pipeline_task(
    runner: PipelineRunner,
    settings: SettingsInput,
    file_dfs: List[pd.DataFrame],
    filenames: List[str],
    job_storage: Dict[str, Dict[str, Any]],
    job_code: str,
):
    """
    Task to run the pipeline in the background.
    """
    try:
        # Execute the pipeline
        results = await runner.run(settings_input=settings, file_dfs=file_dfs, filenames=filenames)

        # Extract the output file from the results
        output_file = results.get("output_file")

        # Update job storage with the download link
        download_url = f"/peanut/download/{job_code}"
        job_storage[job_code].update({
            "status": "Completed",
            "output_file": str(output_file),
            "download_url": download_url,
            "error": None
        })
        logger.info(f"Pipeline completed successfully for job {job_code}. Results available at {download_url}")


    except Exception as e:
        # Update job storage with failure
        job_storage[job_code].update({
            "status": "Failed",
            "output_file": None,
            "download_url": None,
            "error": str(e)
        })
        logger.error(f"Pipeline execution failed for job {job_code}: {e}", exc_info=True)


@router.get("/check-status/{job_code}")
async def check_status(job_code: str) -> JobStatus:
    """
    Check the status of a job by its code.
    """
    if job_code not in job_storage:
        logger.error(f"Job code not found: {job_code}")
        raise HTTPException(status_code=404, detail="Job not found")
        
    job = job_storage[job_code]
    
    # Prepare download URL if job is completed
    download_url = None
    if job["status"] == "Completed" and job["output_file"]:
        download_url = f"/peanut/download/{job_code}"
        logger.debug(f"Download URL prepared for job {job_code}: {download_url}")
    elif job["status"] == "Failed":
        logger.debug(f"Job {job_code} failed with error: {job['error']}")
    else:
        logger.debug(f"Job {job_code} is still in progress: {job['status']}")

    return JobStatus(
        status=job["status"],
        download_url=download_url,
        error=job["error"]
    )


@router.get("/download/{job_code}")
async def download_results(job_code: str):
    """
    Endpoint to download the results of a completed job.
    """
    if job_code not in job_storage or job_storage[job_code]["status"] != "Completed":
        logger.error(f"Download requested for invalid or incomplete job: {job_code}")
        raise HTTPException(status_code=404, detail="Job not found or not completed")

    output_file = job_storage[job_code]["output_file"]
    if not output_file or not Path(output_file).exists():
        logger.error(f"Output file not found for job: {job_code}")
        raise HTTPException(status_code=404, detail="Output file not found")

    return FileResponse(output_file, media_type="application/zip", filename=Path(output_file).name)


@router.get("/", response_class=HTMLResponse)
async def read_index(request: Request):
    """Render the index page."""
    current_date = datetime.now().strftime("%b %d, %Y")
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "current_date": current_date}
    )
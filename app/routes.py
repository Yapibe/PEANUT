import logging
import uuid
from datetime import datetime, timedelta
from typing import List
from pathlib import Path

from fastapi import (
    APIRouter,
    BackgroundTasks,
    File,
    Form,
    HTTPException,
    UploadFile,
)
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .models import JobStatus, PipelineOutput, SettingsInput
from .pipeline_runner import run_pipeline

router = APIRouter()

logger = logging.getLogger(__name__)

static_dir = Path(__file__).resolve().parent / "static"

router.mount(
    "/static",
    StaticFiles(directory=str(static_dir), html=True),
    name="static",
)

# In-memory storage for job status (replace with database in production)
job_storage = {}


@router.get("/")
async def read_index():
    """Read the index page."""
    index_path = Path(__file__).resolve().parent / "templates" / "index.html"
    return FileResponse(str(index_path))


@router.post("/run-pipeline", response_model=PipelineOutput)
async def execute_pipeline(
    background_tasks: BackgroundTasks,
    files: List[UploadFile] = File(...),
    test_name: str = Form(...),
    species: str = Form(...),
    alpha: float = Form(...),
    network: str = Form(...),
    pathway_file: str = Form(...),
    min_genes_per_pathway: int = Form(None),
    max_genes_per_pathway: int = Form(None),
    fdr_threshold: float = Form(...),
    JAC_threshold: float = Form(...),
):
    """Execute the pipeline with the given parameters."""
    logger.info("Received request to run pipeline")
    logger.info(
        f"Files: {[file.filename for file in files]}, "
        f"Test Name: {test_name}, Species: {species}, Alpha: {alpha}, "
        f"Network: {network}, Pathway File: {pathway_file}, "
        f"Min Genes per Pathway: {min_genes_per_pathway}, "
        f"Max Genes per Pathway: {max_genes_per_pathway}, "
        f"FDR Threshold: {fdr_threshold}, JAC Threshold: {JAC_threshold}"
    )

    try:
        settings = SettingsInput(
            test_name=test_name,
            species=species,
            alpha=alpha,
            network=network,
            pathway_file=pathway_file,
            min_genes_per_pathway=min_genes_per_pathway,
            max_genes_per_pathway=max_genes_per_pathway,
            fdr_threshold=fdr_threshold,
            jac_threshold=JAC_threshold,
        )

        job_code = str(uuid.uuid4())[:8]
        job_storage[job_code] = {
            "status": "Processing",
            "output_file": None,
            "expiry": datetime.now() + timedelta(days=7),
        }

        # Read file contents
        file_contents = []
        filenames = []
        for file in files:
            content = await file.read()
            file_contents.append(content)
            filenames.append(file.filename)

        # Run pipeline in the background
        background_tasks.add_task(
            run_pipeline,
            file_contents,
            filenames,
            settings,
            job_storage,
            job_code,
        )

        return PipelineOutput(
            result="Job submitted successfully", job_code=job_code
        )

    except Exception as e:
        logger.error(f"Error occurred while running pipeline: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/check-status/{job_code}", response_model=JobStatus)
async def check_job_status(job_code: str):
    """Check the status of a job."""
    if job_code not in job_storage:
        raise HTTPException(status_code=404, detail="Job not found")

    job = job_storage[job_code]
    if datetime.now() > job["expiry"]:
        del job_storage[job_code]
        raise HTTPException(status_code=404, detail="Job expired")

    return JobStatus(
        status=job["status"],
        download_url=(
            f"/download-file/{job['output_file']}"
            if job["status"] == "Finished" and job["output_file"]
            else None
        ),
    )


@router.get("/download-file/{file_path:path}")
async def download_file(file_path: str):
    """Download a file from the server."""
    logger.info(f"User requested file: {file_path}")
    try:
        full_path = Path(file_path).resolve()
        # Security check to prevent path traversal
        if not full_path.is_file() or not str(full_path).startswith(
            str(Path.cwd())
        ):
            logger.error(f"Invalid file path: {full_path}")
            raise HTTPException(
                status_code=404, detail="File not found"
            )

        logger.info(f"Serving file: {full_path}")
        return FileResponse(str(full_path), filename=full_path.name)
    except Exception as e:
        logger.error(f"Error occurred while sending file: {e}")
        raise HTTPException(status_code=500, detail=str(e))
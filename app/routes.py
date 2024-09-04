from fastapi import APIRouter, HTTPException, UploadFile, File, Form, BackgroundTasks
from .pipeline_runner import run_pipeline
from .models import SettingsInput, PipelineOutput, JobStatus
from typing import Optional
import logging
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
import os
import uuid
from datetime import datetime, timedelta
from starlette.staticfiles import StaticFiles

router = APIRouter()

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

router.mount("/static", StaticFiles(directory=os.path.join(os.getcwd(), "app", "static"), html=True), name="static")

# In-memory storage for job status (replace with database in production)
job_storage = {}

static_dir = os.path.join(os.getcwd(), "app", "static")

@router.get("/")
async def read_index():
    return FileResponse("templates/index.html")


@router.post("/run-pipeline", response_model=PipelineOutput)
async def execute_pipeline(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    test_name: str = Form(...),
    species: str = Form(...),
    alpha: float = Form(...),
    network: str = Form(...),
    pathway_file: str = Form(...),
    min_genes_per_pathway: int = Form(None),
    max_genes_per_pathway: int = Form(None),
    fdr_threshold: float = Form(...),
    JAC_threshold: float = Form(...)
):
    logger.info("Received request to run pipeline")
    logger.info(f"File: {file.filename}, Test Name: {test_name}, Species: {species}, Alpha: {alpha}, Network: {network}, Pathway File: {pathway_file}, Min Genes per Pathway: {min_genes_per_pathway}, Max Genes per Pathway: {max_genes_per_pathway}, FDR Threshold: {fdr_threshold}")
    
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
            JAC_threshold=JAC_threshold
        )

        job_code = str(uuid.uuid4())[:8]  # Generate a unique 8-character job code
        job_storage[job_code] = {
            "status": "Processing",
            "output_file": None,
            "expiry": datetime.now() + timedelta(days=7)  # Set expiry to 7 days from now
        }

        # Read file content
        file_content = await file.read()

        # Run pipeline in the background
        background_tasks.add_task(run_pipeline, file_content, file.filename, settings, job_storage, job_code)

        return PipelineOutput(
            result="Job submitted successfully",
            job_code=job_code
        )

    except Exception as e:
        logger.error(f"Error occurred while running pipeline: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/check-status/{job_code}", response_model=JobStatus)
async def check_job_status(job_code: str):
    if job_code not in job_storage:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = job_storage[job_code]
    if datetime.now() > job["expiry"]:
        del job_storage[job_code]
        raise HTTPException(status_code=404, detail="Job expired")
    
    return JobStatus(
        status=job["status"],
        download_url=f"/download-file/{job['output_file']}" if job["status"] == "Finished" and job["output_file"] else None
    )

@router.get("/download-file/{file_path:path}")
async def download_file(file_path: str):
    logger.info(f"User requested file: {file_path}")
    try:
        full_path = os.path.join(os.getcwd(), file_path)
        if not os.path.exists(full_path):
            logger.error(f"File not found: {full_path}")
            raise HTTPException(status_code=404, detail="File not found")
        
        logger.info(f"Serving file: {full_path}")
        return FileResponse(full_path, filename=os.path.basename(file_path))
    except Exception as e:
        logger.error(f"Error occurred while sending file: {e}")
        raise HTTPException(status_code=500, detail=str(e))
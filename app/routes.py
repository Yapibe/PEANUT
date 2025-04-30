from pathlib import Path
from typing import Dict, List, Optional, Any
from fastapi import APIRouter, UploadFile, HTTPException, BackgroundTasks, Request, Form, File
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.templating import Jinja2Templates
import pandas as pd
import uuid
from datetime import datetime
import traceback
import asyncio
from concurrent.futures import ThreadPoolExecutor
from .models import SettingsInput, JobStatus
from .utils import validate_input_file, setup_logging
from .pipeline_runner import PipelineRunner

logger = setup_logging()
router = APIRouter()

# Set up templates
templates = Jinja2Templates(directory="app/templates")

# In-memory job storage with cleanup
job_storage: Dict[str, Dict] = {}
JOB_EXPIRY_HOURS = 4  # Jobs expire after 4 hours

# Sample data file paths
SAMPLE_DATA_DIR = Path("app/static/sample_data")
SAMPLE_EXPRESSION_FILE = SAMPLE_DATA_DIR / "GSE5281VCX_KEGG_ALZHEIMERS_DISEASE.rnk"
SAMPLE_KEGG_FILE = SAMPLE_DATA_DIR / "kegg.gmt"

# Thread pool for file operations
thread_pool = ThreadPoolExecutor(max_workers=4)

async def cleanup_expired_jobs():
    """Clean up expired jobs from storage."""
    current_time = datetime.now()
    expired_jobs = [
        job_id for job_id, job in job_storage.items()
        if (current_time - datetime.fromisoformat(job['created_at'])).total_seconds() > JOB_EXPIRY_HOURS * 3600
    ]
    for job_id in expired_jobs:
        logger.info(f"Cleaning up expired job: {job_id}")
        del job_storage[job_id]

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
    """Run the PEANUT pipeline with the provided settings and input files."""
    try:
        # Clean up expired jobs
        await cleanup_expired_jobs()
        
        # Generate job code
        job_code = str(uuid.uuid4())
        logger.info(f"Generated job code: {job_code}")
        
        # Initialize job storage
        job_storage[job_code] = {
            "status": "Processing",
            "output_file": None,
            "error": None,
            "created_at": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat()
        }
        logger.debug(f"Initialized job storage for {job_code}")

        # Validate and process input files
        file_dfs, filenames = [], []
        for file in files:
            file_type = file.filename.split('.')[-1].lower()
            if file_type not in ['rnk', 'xlsx']:
                raise HTTPException(
                    status_code=400,
                    detail=f"Unsupported file type: {file_type}. Only .rnk and .xlsx files are supported."
                )

            # Read file content
            content = await file.read()
            temp_path = Path(f"/tmp/{file.filename}")
            
            try:
                # Write to temporary file
                temp_path.write_bytes(content)
                
                # Validate file
                df = await asyncio.get_event_loop().run_in_executor(
                    thread_pool,
                    lambda: validate_input_file(temp_path, file_type)
                )
                file_dfs.append(df)
                filenames.append(file.filename)
            finally:
                # Clean up temporary file
                if temp_path.exists():
                    temp_path.unlink()

        # Create pipeline runner
        runner = PipelineRunner(job_code)
        
        # Add pipeline task to background tasks
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
        
        logger.info(f"Pipeline request processed successfully. Returning job code: {job_code}")
        return {"job_code": job_code}

    except Exception as e:
        logger.error(f"Pipeline request failed with error: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        if isinstance(e, HTTPException):
            raise
        raise HTTPException(status_code=422, detail=str(e))

async def _run_pipeline_task(
    runner: PipelineRunner,
    settings: SettingsInput,
    file_dfs: List[pd.DataFrame],
    filenames: List[str],
    job_storage: Dict[str, Dict[str, Any]],
    job_code: str,
):
    """Run the pipeline in the background."""
    logger.info(f"Starting background pipeline task for job {job_code}")
    
    try:
        # Execute pipeline
        results = await runner.run(settings_input=settings, file_dfs=file_dfs, filenames=filenames)
        
        output_file = results.get("output_file")
        logger.info(f"Pipeline execution completed. Output file: {output_file}")
        
        download_url = f"/peanut/download/{job_code}"
        logger.debug(f"Generated download URL: {download_url}")
        
        # Update job storage
        job_storage[job_code].update({
            "status": "Completed",
            "output_file": str(output_file),
            "download_url": download_url,
            "error": None,
            "last_updated": datetime.now().isoformat()
        })
        logger.info(f"Pipeline completed successfully for job {job_code}. Results available at {download_url}")

    except Exception as e:
        logger.error(f"Pipeline execution failed for job {job_code}")
        logger.error(f"Error details: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        
        job_storage[job_code].update({
            "status": "Failed",
            "output_file": None,
            "download_url": None,
            "error": str(e),
            "last_updated": datetime.now().isoformat()
        })

@router.get("/check-status/{job_code}")
async def check_status(job_code: str) -> JobStatus:
    """Check the status of a job by its code."""
    logger.info(f"Checking status for job {job_code}")
    
    # Clean up expired jobs
    await cleanup_expired_jobs()
    
    if job_code not in job_storage:
        logger.error(f"Job code {job_code} not found in storage")
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = job_storage[job_code]
    logger.debug(f"Retrieved job data: {job}")
    
    download_url = None
    if job["status"] == "Completed" and job["output_file"]:
        download_url = f"/peanut/download/{job_code}"
        logger.info(f"Job completed. Download URL: {download_url}")
    
    return JobStatus(
        status=job["status"],
        download_url=download_url,
        error=job["error"]
    )

@router.get("/download/{job_code}")
async def download_results(job_code: str):
    """Download the results of a completed job."""
    # Clean up expired jobs
    await cleanup_expired_jobs()
    
    if job_code not in job_storage or job_storage[job_code]["status"] != "Completed":
        logger.error(f"Download requested for invalid or incomplete job: {job_code}")
        raise HTTPException(status_code=404, detail="Job not found or not completed")

    output_file = job_storage[job_code]["output_file"]
    if not output_file or not Path(output_file).exists():
        logger.error(f"Output file not found for job: {job_code}")
        raise HTTPException(status_code=404, detail="Output file not found")

    return FileResponse(
        output_file,
        media_type="application/zip",
        filename=Path(output_file).name,
        headers={
            "Content-Disposition": f"attachment; filename={Path(output_file).name}"
        }
    )

@router.get("/", response_class=HTMLResponse)
async def read_index(request: Request):
    """Render the index page."""
    current_date = datetime.now().strftime("%b %d, %Y")
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "current_date": current_date}
    )

@router.get("/sample-data/expression")
async def download_sample_expression():
    """Download sample expression data."""
    if not SAMPLE_EXPRESSION_FILE.exists():
        logger.error("Sample expression file not found")
        raise HTTPException(status_code=404, detail="Sample expression file not found")
    
    return FileResponse(
        SAMPLE_EXPRESSION_FILE,
        media_type="text/plain",
        filename=SAMPLE_EXPRESSION_FILE.name,
        headers={
            "Content-Disposition": f"attachment; filename={SAMPLE_EXPRESSION_FILE.name}"
        }
    )

@router.get("/sample-data/kegg")
async def download_sample_kegg():
    """Download sample KEGG pathway data."""
    if not SAMPLE_KEGG_FILE.exists():
        logger.error("Sample KEGG file not found")
        raise HTTPException(status_code=404, detail="Sample KEGG file not found")
    
    return FileResponse(
        SAMPLE_KEGG_FILE,
        media_type="text/plain",
        filename=SAMPLE_KEGG_FILE.name,
        headers={
            "Content-Disposition": f"attachment; filename={SAMPLE_KEGG_FILE.name}"
        }
    )
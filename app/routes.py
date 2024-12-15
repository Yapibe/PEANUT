import logging
import uuid
from datetime import datetime, timedelta
from typing import List, Optional
from pathlib import Path
from fastapi import (
    APIRouter,
    BackgroundTasks,
    File,
    Form,
    HTTPException,
    UploadFile,
    Request,
)
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.templating import Jinja2Templates
from .utils import validate_rnk_file, validate_xlsx_file, validate_network_file
from .models import JobStatus, PipelineOutput, SettingsInput
from .pipeline_runner import run_pipeline

router = APIRouter()

logger = logging.getLogger(__name__)

templates_dir = Path(__file__).resolve().parent / "templates"
templates = Jinja2Templates(directory=str(templates_dir))
logger.debug(f"Templates directory: {templates_dir}")

# In-memory storage for job status (replace with database in production)
job_storage = {}

@router.get("/", response_class=HTMLResponse)
async def read_index(request: Request):
    """Render the index page."""
    logger.debug(f"Root Path: {request.scope.get('root_path')}")
    current_date = datetime.now().strftime("%b %d, %Y")
    return templates.TemplateResponse("index.html", {"request": request, "url_for": request.url_for, "current_date": current_date})

@router.post("/run-pipeline", response_model=PipelineOutput)
async def execute_pipeline(
    background_tasks: BackgroundTasks,
    files: List[UploadFile] = File(...),
    network_type: str = Form(...),
    network_file: Optional[UploadFile] = None,
    network: Optional[str] = None,
    test_name: str = Form(...),
    species: str = Form(...),
    alpha: float = Form(...),
    pathway_file: str = Form(...),
    min_genes_per_pathway: int = Form(None),
    max_genes_per_pathway: int = Form(None),
    fdr_threshold: float = Form(...),
    run_gsea: bool = Form(False),
    restrict_to_network: bool = Form(False)
):
    """Execute the pipeline with the given parameters."""
    logger.info("Received request to run pipeline")
    
    try:
        # Validate network selection
        if network_type not in ['default', 'custom']:
            raise HTTPException(
                status_code=400,
                detail="Invalid network type. Must be either 'default' or 'custom'."
            )
        
        # Handle network validation
        custom_network_df = None
        if network_type == 'custom':
            if not network_file:
                raise HTTPException(
                    status_code=400,
                    detail="Network file is required when using custom network."
                )
            network_content = await network_file.read()
            custom_network_df = await validate_network_file(network_content, network_file.filename)
        elif network_type == 'default' and not network:
            raise HTTPException(
                status_code=400,
                detail="Network selection is required when using default network."
            )

        # Process input files
        file_dfs = []
        filenames = []
        for file in files:
            content = await file.read()
            if file.filename.endswith('.rnk'):
                df = await validate_rnk_file(content, file.filename)
            else:
                df = await validate_xlsx_file(content, file.filename)
            file_dfs.append(df)
            filenames.append(file.filename)

        # Create settings object
        settings = SettingsInput(
            test_name=test_name,
            species=species,
            alpha=alpha,
            network_type=network_type,
            network=network if network_type == 'default' else None,
            custom_network_df=custom_network_df,
            pathway_file=pathway_file,
            min_genes_per_pathway=min_genes_per_pathway,
            max_genes_per_pathway=max_genes_per_pathway,
            fdr_threshold=fdr_threshold,
            run_gsea=run_gsea,
            restrict_to_network=restrict_to_network
        )

        # Generate job code and set up storage
        job_code = str(uuid.uuid4())
        job_storage[job_code] = {"status": "Running"}

        # Run pipeline in background
        background_tasks.add_task(
            run_pipeline,
            file_dfs,
            filenames,
            settings,
            job_storage,
            job_code,
        )

        return PipelineOutput(
            result="Pipeline execution started",
            job_code=job_code,
        )

    except Exception as e:
        logger.error(f"Error in execute_pipeline: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/check-status/{job_code}", response_model=JobStatus)
async def check_job_status(request: Request, job_code: str):
    """Check the status of a job."""
    try:
        logger.debug(f"Checking status for job_code: {job_code}")
        if job_code not in job_storage:
            logger.warning(f"Job code {job_code} not found in storage.")
            raise HTTPException(status_code=404, detail="Job not found")

        job = job_storage[job_code]
        logger.debug(f"Job data: {job}")

        if datetime.now() > job["expiry"]:
            logger.info(f"Job {job_code} has expired.")
            del job_storage[job_code]
            raise HTTPException(status_code=404, detail="Job expired")

        download_url = None
        if job["status"] == "Finished" and job["output_file"]:
            # Remove leading slash to prevent double slash in URL
            clean_file_path = job['output_file'].lstrip('/')
            # Construct relative download URL based on root_path
            download_url = f"{request.scope.get('root_path')}/download-file/{clean_file_path}"
            logger.debug(f"Download URL for job {job_code}: {download_url}")
        
        return JobStatus(
            status=job["status"],
            download_url=download_url,
        )
    except Exception as e:
        logger.error(f"Error in check_job_status for job_code {job_code}: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail="Internal Server Error")


@router.get("/download-file/{file_path:path}", name='download_file')  # Added 'name' parameter
async def download_file(request: Request, file_path: str):
    """Download a file from the server."""
    logger.info(f"User requested file: {file_path}")
    try:
        full_path = Path(file_path).resolve()
        # Security check to prevent path traversal
        if not full_path.is_file() or not str(full_path).startswith(str(Path.cwd())):
            logger.error(f"Invalid file path: {full_path}")
            raise HTTPException(status_code=404, detail="File not found")

        logger.info(f"Serving file: {full_path}")
        return FileResponse(str(full_path), filename=full_path.name)
    except Exception as e:
        logger.error(f"Error occurred while sending file: {e}")
        raise HTTPException(status_code=500, detail=str(e))
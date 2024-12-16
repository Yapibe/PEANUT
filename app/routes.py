import logging
from pathlib import Path
from typing import Dict, List, Optional
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
    files: List[UploadFile] = File(...)
):
    """
    Endpoint to run the pipeline with the provided settings and input files.
    """
    try:
        # Generate unique job code
        job_code = str(uuid.uuid4())
        logger.info(f"New job request received. Job code: {job_code}")
        
        # Initialize job storage
        job_storage[job_code] = {
            "status": "Processing",
            "output_file": None,
            "error": None
        }
        
        # Handle network file if provided
        if network_file and network == 'custom':
            logger.debug("Reading custom network file")
            network_content = await network_file.read()
            # Create a new UploadFile with the content
            new_network_file = UploadFile(
                filename=network_file.filename,
                file=pd.io.common.BytesIO(network_content)
            )
            network_file = new_network_file
        
        # Determine if a custom network is being used
        is_custom_network = network == 'custom' and network_file is not None

        # Handle custom pathway file if provided
        if pathway_file == 'custom' and pathway_file_upload:
            logger.debug("Processing custom pathway file")
            
            # Read the GMT file content
            content = await pathway_file_upload.read()
            
            # Create temporary directory for processing
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_path = Path(temp_dir)
                
                # Save GMT file temporarily
                gmt_path = temp_path / "custom.gmt"
                gmt_path.write_bytes(content)
                
                try:
                    # Create custom pathways directory
                    custom_pathway_dir = Path("pipeline/Data") / species / "pathways" / "custom"
                    output_file_path = custom_pathway_dir / f"custom_{experiment_name}.tsv"
                    
                    # Convert GMT file to pathway format
                    convert_gmt_to_pathway_format(gmt_path, output_file_path)
                    
                    # Update pathway_file to use the custom file
                    pathway_file = f"custom/custom_{experiment_name}"
                    
                except ValueError as e:
                    logger.error(f"Error processing GMT file: {str(e)}")
                    raise HTTPException(
                        status_code=400,
                        detail=str(e)
                    )

        # Create settings object
        settings = SettingsInput(
            experiment_name=experiment_name,
            species=species,
            alpha=alpha,
            network=network,
            network_file=network_file,
            pathway_file=pathway_file,
            min_genes_per_pathway=min_genes_per_pathway,
            max_genes_per_pathway=max_genes_per_pathway,
            fdr_threshold=fdr_threshold,
            run_gsea=run_gsea,
            restrict_to_network=restrict_to_network,
            create_similarity_matrix=is_custom_network,
            figure_title=plot_title
        )
        
        # Validate input files
        file_dfs = []
        filenames = []
        for file in files:
            file_type = file.filename.split('.')[-1].lower()
            if file_type not in ['rnk', 'xlsx']:
                raise HTTPException(
                    status_code=400,
                    detail=f"Unsupported file type: {file_type}"
                )
            
            # Read and validate file content
            content = await file.read()
            temp_path = Path(f"/tmp/{file.filename}")
            temp_path.write_bytes(content)
            
            try:
                df = validate_input_file(temp_path, file_type)
                file_dfs.append(df)
                filenames.append(file.filename)
            finally:
                temp_path.unlink()  # Clean up temporary file
        
        # Create pipeline runner
        runner = PipelineRunner(job_code)
        
        # Run pipeline in background
        background_tasks.add_task(
            _run_pipeline_task,
            runner=runner,
            settings=settings,
            file_dfs=file_dfs,
            filenames=filenames,
            job_storage=job_storage,
            job_code=job_code
        )
        
        return {"job_code": job_code}
        
    except Exception as e:
        logger.error(f"Error processing pipeline request: {str(e)}")
        raise HTTPException(status_code=422, detail=str(e))

async def _run_pipeline_task(
    runner: PipelineRunner,
    settings: SettingsInput,
    file_dfs: List[pd.DataFrame],
    filenames: List[str],
    job_storage: Dict,
    job_code: str
):
    """Background task to run the pipeline."""
    try:
        # Run pipeline
        results = await runner.run(settings, file_dfs, filenames)
        
        # Update job storage with results
        job_storage[job_code].update({
            "status": "Finished",
            "output_file": str(results),  # results is now a Path object
            "error": None
        })
        
    except Exception as e:
        logger.error(f"Pipeline execution failed for job {job_code}: {str(e)}")
        job_storage[job_code].update({
            "status": "Failed",
            "error": str(e)
        })

@router.get("/check-status/{job_code}")
async def check_status(job_code: str) -> JobStatus:
    """
    Check the status of a job by its code.
    """
    if job_code not in job_storage:
        raise HTTPException(status_code=404, detail="Job not found")
        
    job = job_storage[job_code]
    
    # Prepare download URL if job is finished
    download_url = None
    if job["status"] == "Finished" and job["output_file"]:
        download_url = f"/download/{job_code}"
    
    return JobStatus(
        status=job["status"],
        download_url=download_url,
        error=job["error"]
    )

@router.get("/download/{job_code}")
async def download_results(job_code: str):
    """
    Download results for a completed job.
    """
    if job_code not in job_storage:
        raise HTTPException(status_code=404, detail="Job not found")
        
    job = job_storage[job_code]
    if job["status"] != "Finished":
        raise HTTPException(status_code=400, detail="Job not finished")
        
    if not job["output_file"]:
        raise HTTPException(status_code=404, detail="Output file not found")
    
    # Get the experiment name from the output file path
    output_path = Path(job["output_file"])
    experiment_name = output_path.stem.split('_')[0]  # Get the first part before underscore
        
    return FileResponse(
        job["output_file"],
        filename=f"{experiment_name}_results.zip",
        media_type="application/zip"
    )

@router.get("/", response_class=HTMLResponse)
async def read_index(request: Request):
    """Render the index page."""
    current_date = datetime.now().strftime("%b %d, %Y")
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "current_date": current_date}
    )
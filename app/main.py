from fastapi import FastAPI, File, UploadFile, Form
from pipeline_runner import run_pipeline
from .models import SettingsInput

app = FastAPI()

@app.post("/run-pipeline")
async def run_pipeline_endpoint(
    test_name: str = Form(...),
    species: str = Form(...),
    alpha: float = Form(...),
    network: str = Form(...),
    pathway_file: str = Form(...),
    min_genes_per_pathway: int = Form(15),
    max_genes_per_pathway: int = Form(500),
    fdr_threshold: float = Form(0.05),
    JAC_threshold: float = Form(0.2),
    files: List[UploadFile] = File(...),
):
    job_code = generate_job_code()  # Implement job code generation
    job_storage = get_job_storage()  # Implement job storage retrieval

    # Read file contents
    file_bytes = [await file.read() for file in files]
    filenames = [file.filename for file in files]

    settings = SettingsInput(
        propagation_alpha=alpha,
        fdr_threshold=fdr_threshold,
        species=species,
        network=network,
        pathway_file=pathway_file,
        min_genes_per_pathway=min_genes_per_pathway,
        max_genes_per_pathway=max_genes_per_pathway,
        JAC_threshold=JAC_threshold,
    )

    await run_pipeline(
        files=file_bytes,
        filenames=filenames,
        settings=settings,
        job_storage=job_storage,
        job_code=job_code,
    )

    return {"job_code": job_code}
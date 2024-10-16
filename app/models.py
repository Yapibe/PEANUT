from typing import Optional

from pydantic import BaseModel


class SettingsInput(BaseModel):
    test_name: str
    species: str
    alpha: float
    network: str
    pathway_file: str
    min_genes_per_pathway: Optional[int] = None
    max_genes_per_pathway: Optional[int] = None
    fdr_threshold: float
    jac_threshold: float


class PipelineOutput(BaseModel):
    result: str
    job_code: str



class JobStatus(BaseModel):
    status: str
    download_url: Optional[str] = None  # Make download_url optional


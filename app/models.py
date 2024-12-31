from typing import Optional
from fastapi import UploadFile
from pydantic import BaseModel, Field


class SettingsInput(BaseModel):
    """Settings for the pipeline run."""
    experiment_name: str
    species: str
    alpha: float
    network: str  # Can be 'Anat', 'HumanNet', or 'custom'
    network_file: Optional[UploadFile] = None  # For custom network file upload
    pathway_file: str  # Can be 'KEGG', 'Reactome', 'custom', etc.
    pathway_file_upload: Optional[UploadFile] = None  # For custom pathway file upload
    min_genes_per_pathway: Optional[int] = Field(default=15, ge=1, le=1000)
    max_genes_per_pathway: Optional[int] = Field(default=500, ge=1, le=5000)
    fdr_threshold: float = Field(default=0.05, ge=0, le=1)
    run_gsea: bool = False
    restrict_to_network: bool = False
    create_similarity_matrix: bool = False
    figure_title: str = "Pathway Enrichment"

    def is_custom_network(self) -> bool:
        """Check if using a custom network."""
        return self.network == 'custom' and self.network_file is not None

    def is_custom_pathway(self) -> bool:
        """Check if using a custom pathway file."""
        return self.pathway_file == 'custom' and self.pathway_file_upload is not None

    model_config = {
        "arbitrary_types_allowed": True  # Allow UploadFile and other complex types
    }



class PipelineOutput(BaseModel):
    result: str
    job_code: str


class JobStatus(BaseModel):
    status: str
    download_url: Optional[str] = None
    error: Optional[str] = None

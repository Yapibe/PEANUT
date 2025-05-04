from typing import Optional, Literal, Annotated
from fastapi import UploadFile
from pydantic import BaseModel, Field, validator, StringConstraints
import re


class SettingsInput(BaseModel):
    """Settings for the PEANUT pipeline run.
    
    Attributes:
        experiment_name: Name of the experiment
        species: Species identifier (e.g., 'human', 'mouse')
        alpha: Network propagation parameter (0-1)
        network: Network type ('Anat', 'HumanNet', or 'custom')
        network_file: Optional custom network file
        network_file_path: Optional path to saved network file
        pathway_file: Pathway database ('KEGG', 'Reactome', or 'custom')
        pathway_file_upload: Optional custom pathway file
        pathway_file_path: Optional path to saved pathway file
        min_genes_per_pathway: Minimum genes per pathway
        max_genes_per_pathway: Maximum genes per pathway
        fdr_threshold: False discovery rate threshold
        run_gsea: Whether to run GSEA analysis
        restrict_to_network: Whether to restrict analysis to network genes
        create_similarity_matrix: Whether to create similarity matrix
        figure_title: Title for the output figure
    """
    experiment_name: Annotated[str, StringConstraints(min_length=1, max_length=100)] = Field(
        ..., description="Name of the experiment"
    )
    species: Literal['H_sapiens','mouse'] = Field(..., description="Species for the analysis")
    alpha: float = Field(
        ..., ge=0.0, le=1.0, description="Network propagation parameter"
    )
    network: Literal['Anat', 'HumanNet', 'custom'] = Field(
        ..., description="Network type to use"
    )
    network_file: Optional[UploadFile] = Field(
        None, description="Custom network file (required if network='custom')"
    )
    network_file_path: Optional[str] = Field(
        None, description="Path to saved network file"
    )
    pathway_file: Literal['kegg', 'msig_c2_canonical', 'custom'] = Field(
        ..., description="Pathway database to use"
    )
    pathway_file_upload: Optional[UploadFile] = Field(
        None, description="Custom pathway file (required if pathway_file='custom')"
    )
    pathway_file_path: Optional[str] = Field(
        None, description="Path to saved pathway file"
    )
    min_genes_per_pathway: int = Field(
        default=15, ge=1, le=1000, description="Minimum genes per pathway"
    )
    max_genes_per_pathway: int = Field(
        default=500, ge=1, le=5000, description="Maximum genes per pathway"
    )
    fdr_threshold: float = Field(
        default=0.05, ge=0.0, le=1.0, description="False discovery rate threshold"
    )
    run_gsea: bool = Field(
        default=False, description="Whether to run GSEA analysis"
    )
    restrict_to_network: bool = Field(
        default=False, description="Whether to restrict analysis to network genes"
    )
    create_similarity_matrix: bool = Field(
        default=False, description="Whether to create similarity matrix"
    )
    figure_title: Annotated[str, StringConstraints(min_length=1, max_length=200)] = Field(
        default="Pathway Enrichment", description="Title for the output figure"
    )

    @validator('experiment_name')
    def validate_experiment_name(cls, v):
        """Validate experiment name format."""
        if not re.match(r'^[a-zA-Z0-9_\- ]+$', v):
            raise ValueError('Experiment name can only contain letters, numbers, spaces, underscores, and hyphens')
        return v

    def is_custom_network(self) -> bool:
        """Check if using a custom network."""
        return self.network == 'custom' and (self.network_file is not None or self.network_file_path is not None)

    def is_custom_pathway(self) -> bool:
        """Check if using a custom pathway file."""
        return self.pathway_file == 'custom' and (self.pathway_file_upload is not None or self.pathway_file_path is not None)
        
    def needs_similarity_matrix(self) -> bool:
        """Check if computation needs a similarity matrix.
        
        Returns:
            bool: True if computation requires a similarity matrix
        """
        # No matrix needed for alpha=1 (GSEA equivalent)
        if self.alpha == 1.0:
            return False
            
        return True
        
    def needs_matrix_generation(self) -> bool:
        """Determine if a new similarity matrix needs to be generated.
        
        Returns:
            bool: True if a new matrix needs to be generated
        """
        # Custom network always needs matrix generation
        if self.is_custom_network():
            return True
            
        # No matrix needed for alpha=1 (GSEA equivalent)
        if self.alpha == 1.0:
            return False
            
        # For default networks, only alpha 0.1 and 0.2 have precomputed matrices
        if self.network in ['Anat', 'HumanNet'] and self.alpha not in [0.1, 0.2]:
            return True
            
        return False

    class Config:
        arbitrary_types_allowed = True
        json_schema_extra = {
            "example": {
                "experiment_name": "my_experiment",
                "species": "human",
                "alpha": 0.5,
                "network": "HumanNet",
                "pathway_file": "KEGG",
                "min_genes_per_pathway": 15,
                "max_genes_per_pathway": 500,
                "fdr_threshold": 0.05,
                "run_gsea": False,
                "restrict_to_network": False,
                "create_similarity_matrix": False,
                "figure_title": "Pathway Enrichment Analysis"
            }
        }


class PipelineOutput(BaseModel):
    """Output from the PEANUT pipeline.
    
    Attributes:
        result: Path to the result file
        job_code: Unique identifier for the job
    """
    result: str = Field(..., description="Path to the result file")
    job_code: str = Field(..., description="Unique identifier for the job")


class JobStatus(BaseModel):
    """Status of a PEANUT pipeline job.
    
    Attributes:
        status: Current status of the job
        download_url: URL to download results (if completed)
        error: Error message (if failed)
        job_code: Unique identifier for the job
    """
    status: Literal['Processing', 'Completed', 'Failed'] = Field(
        ..., description="Current status of the job"
    )
    download_url: Optional[str] = Field(
        None, description="URL to download results (if completed)"
    )
    error: Optional[str] = Field(
        None, description="Error message (if failed)"
    )
    job_code: str = Field(
        ..., description="Unique identifier for the job"
    )

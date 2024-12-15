from typing import Optional

from pydantic import BaseModel
import pandas as pd


class SettingsInput(BaseModel):
    """Settings for the pipeline run."""
    test_name: str
    species: str
    create_similarity_matrix: bool
    # network_type: str  # 'default' or 'custom'
    network: Optional[str] = None  # Used when network_type is 'default'
    # custom_network_df: Optional[pd.DataFrame] = None  # Used when network_type is 'custom'
    alpha: float
    pathway_file: str
    min_genes_per_pathway: Optional[int] = None
    max_genes_per_pathway: Optional[int] = None
    fdr_threshold: float
    run_gsea: bool = False
    restrict_to_network: bool = False

    model_config = {
        "arbitrary_types_allowed": True  # Allow pandas DataFrame and other complex types
    }

class PipelineOutput(BaseModel):
    result: str
    job_code: str



class JobStatus(BaseModel):
    status: str
    download_url: Optional[str] = None  # Make download_url optional


import logging
import pandas as pd
from fastapi import HTTPException
from pathlib import Path
from typing import Tuple, List
import yaml
import logging.config

def setup_logging():
    """Load logging configuration from YAML file."""
    config_path = Path(__file__).parent.parent / 'config' / 'log_config.yaml'
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    logging.config.dictConfig(config)
    return logging.getLogger(__name__)

logger = setup_logging()

def validate_columns(df: pd.DataFrame, required_columns: List[str], filename: str) -> None:
    """Validates that a DataFrame contains required columns."""
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in {filename}: {', '.join(missing_cols)}")

def validate_numeric_values(df: pd.DataFrame, numeric_columns: List[str], filename: str) -> None:
    """Validates that specified columns contain numeric values."""
    for col in numeric_columns:
        numeric_series = pd.to_numeric(df[col], errors='coerce')
        non_numeric = df[numeric_series.isna()]
        if not non_numeric.empty:
            raise ValueError(f"Non-numeric values found in {col} column of {filename}")

def validate_network_file(file_path: Path) -> Tuple[pd.DataFrame, bool]:
    """Validates a custom network file."""
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        
        if df.shape[1] < 2:
            raise ValueError("Network file must have at least 2 columns (Source, Target)")
        
        if df.shape[1] == 2:
            df.columns = ['Source', 'Target']
            is_weighted = False
        else:
            df.columns = ['Source', 'Target', 'Weight']
            is_weighted = True
        
        numeric_cols = ['Source', 'Target']
        if is_weighted:
            numeric_cols.append('Weight')
            validate_numeric_values(df, numeric_cols, file_path.name)
            if not ((df['Weight'] >= 0) & (df['Weight'] <= 1)).all():
                raise ValueError("Weights must be between 0 and 1")
        
        return df, is_weighted
        
    except Exception as e:
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Invalid network file: {str(e)}")

def validate_input_file(file_path: Path, file_type: str) -> pd.DataFrame:
    """Validates input files (RNK or XLSX)."""
    try:
        if file_type == 'rnk':
            df = pd.read_csv(file_path, sep='\t', header=None)
            df.columns = ['GeneID', 'Score']
        else:  # xlsx
            df = pd.read_excel(file_path)
        
        validate_columns(df, ['GeneID', 'Score'], file_path.name)
        validate_numeric_values(df, ['Score'], file_path.name)
        
        if df['GeneID'].duplicated().any():
            df = df.drop_duplicates(subset=['GeneID'], keep='first')
        
        return df
        
    except Exception as e:
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Invalid {file_type} file: {str(e)}")

def convert_gmt_to_pathway_format(gmt_path: Path, output_path: Path) -> None:
    """Converts a GMT file to the pathway format used by the pipeline."""
    try:
        new_data = []
        with open(gmt_path, 'r') as file:
            for line_num, line in enumerate(file, 1):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    raise ValueError(f"Invalid GMT format at line {line_num}: each line must have at least 3 columns")
                pathway_name = parts[0]
                genes = parts[2:]  # Skip the description
                new_data.append([pathway_name] + genes)
        
        new_pathways_df = pd.DataFrame(new_data)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        new_pathways_df.to_csv(output_path, index=False, header=False, sep='\t')
        
    except Exception as e:
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Failed to convert GMT file: {str(e)}")

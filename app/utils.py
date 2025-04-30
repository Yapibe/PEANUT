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
    """
    Validates that a DataFrame contains required columns.
    
    Args:
        df: DataFrame to validate
        required_columns: List of required column names
        filename: Name of file being validated for error messages
    
    Raises:
        ValueError: If required columns are missing
    """
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        msg = f"Missing required columns in {filename}: {', '.join(missing_cols)}"
        logger.error(msg)
        raise ValueError(msg)

def validate_numeric_values(df: pd.DataFrame, numeric_columns: List[str], filename: str) -> None:
    """
    Validates that specified columns contain numeric values.
    
    Args:
        df: DataFrame to validate
        numeric_columns: List of columns that should be numeric
        filename: Name of file being validated for error messages
    
    Raises:
        ValueError: If non-numeric values are found
    """
    for col in numeric_columns:
        try:
            # Try to convert to numeric, coercing errors to NaN
            numeric_series = pd.to_numeric(df[col], errors='coerce')
            non_numeric = df[numeric_series.isna()]
            
            if not non_numeric.empty:
                msg = f"Non-numeric values found in {col} column of {filename}"
                logger.error(f"{msg}. First few invalid rows: {non_numeric.head()}")
                raise ValueError(msg)
        except Exception as e:
            logger.error(f"Error validating numeric values in {filename}: {str(e)}")
            raise

def validate_network_file(file_path: Path) -> Tuple[pd.DataFrame, bool]:
    """
    Validates a custom network file.
    
    Args:
        file_path: Path to the network file
    
    Returns:
        Tuple of (validated DataFrame, is_weighted)
    
    Raises:
        ValueError: If validation fails
        HTTPException: If file cannot be read
    """
    try:
        # Read the file with error handling
        try:
            df = pd.read_csv(file_path, sep='\t', header=None)
        except Exception as e:
            raise HTTPException(
                status_code=400,
                detail=f"Could not read network file: {str(e)}"
            )
        
        logger.info(f"Successfully read network file: {file_path}")
        
        # Check minimum columns
        if df.shape[1] < 2:
            raise ValueError("Network file must have at least 2 columns (Source, Target)")
        
        # Rename columns for consistency
        if df.shape[1] == 2:
            df.columns = ['Source', 'Target']
            is_weighted = False
        else:
            df.columns = ['Source', 'Target', 'Weight']
            is_weighted = True
        
        # Validate numeric values
        numeric_cols = ['Source', 'Target']
        if is_weighted:
            numeric_cols.append('Weight')
        
        validate_numeric_values(df, numeric_cols, file_path.name)
        
        if is_weighted:
            # Validate weight range
            if not ((df['Weight'] >= 0) & (df['Weight'] <= 1)).all():
                raise ValueError("Weights must be between 0 and 1")
        
        # Check for duplicate edges
        if df.duplicated(subset=['Source', 'Target']).any():
            logger.warning(f"Duplicate edges found in network file: {file_path}")
        
        logger.info(f"Network file validation successful: {file_path}")
        return df, is_weighted
        
    except Exception as e:
        logger.error(f"Error validating network file {file_path}: {str(e)}")
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Invalid network file: {str(e)}")

def validate_input_file(file_path: Path, file_type: str) -> pd.DataFrame:
    """
    Validates input files (RNK or XLSX).
    
    Args:
        file_path: Path to the input file
        file_type: Type of file ('rnk' or 'xlsx')
    
    Returns:
        Validated DataFrame
    
    Raises:
        ValueError: If validation fails
        HTTPException: If file cannot be read
    """
    try:
        # Read the file with error handling
        try:
            if file_type == 'rnk':
                df = pd.read_csv(file_path, sep='\t', header=None)
                df.columns = ['GeneID', 'Score']
            else:  # xlsx
                df = pd.read_excel(file_path)
        except Exception as e:
            raise HTTPException(
                status_code=400,
                detail=f"Could not read {file_type} file: {str(e)}"
            )
        
        logger.info(f"Successfully read {file_type} file: {file_path}")
        
        # Validate required columns
        validate_columns(df, ['GeneID', 'Score'], file_path.name)
        
        # Validate Score column is numeric
        validate_numeric_values(df, ['Score'], file_path.name)
        
        # Check for duplicate genes
        if df['GeneID'].duplicated().any():
            logger.warning(f"Duplicate genes found in {file_path}")
            # Keep only the first occurrence of each gene
            df = df.drop_duplicates(subset=['GeneID'], keep='first')
        
        logger.info(f"Input file validation successful: {file_path}")
        return df
        
    except Exception as e:
        logger.error(f"Error validating {file_type} file {file_path}: {str(e)}")
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Invalid {file_type} file: {str(e)}")

def convert_gmt_to_pathway_format(gmt_path: Path, output_path: Path) -> None:
    """
    Converts a GMT file to the pathway format used by the pipeline.
    
    Args:
        gmt_path: Path to the input GMT file
        output_path: Path where the converted file should be saved
        
    Raises:
        ValueError: If the GMT file format is invalid
        HTTPException: If file cannot be read or written
    """
    try:
        new_data = []
        try:
            with open(gmt_path, 'r') as file:
                for line_num, line in enumerate(file, 1):
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        raise ValueError(f"Invalid GMT format at line {line_num}: each line must have at least 3 columns")
                    pathway_name = parts[0]
                    genes = parts[2:]  # Skip the description
                    new_data.append([pathway_name] + genes)
        except Exception as e:
            raise HTTPException(
                status_code=400,
                detail=f"Could not read GMT file: {str(e)}"
            )
        
        # Convert to DataFrame and save
        try:
            new_pathways_df = pd.DataFrame(new_data)
            
            # Ensure output directory exists
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Save to TSV format
            new_pathways_df.to_csv(
                output_path,
                index=False,
                header=False,
                sep='\t'
            )
            logger.info(f"Successfully converted GMT file to pathway format: {output_path}")
            
        except Exception as e:
            raise HTTPException(
                status_code=500,
                detail=f"Could not write pathway file: {str(e)}"
            )
        
    except Exception as e:
        logger.error(f"Error converting GMT file: {str(e)}")
        if isinstance(e, HTTPException):
            raise
        raise ValueError(f"Failed to convert GMT file: {str(e)}")

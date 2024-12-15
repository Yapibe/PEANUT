import io
import logging
import pandas as pd
from fastapi import HTTPException
import numpy as np
from pipeline.utils import read_prior_set

logger = logging.getLogger(__name__)

async def validate_rnk_file(file_content: bytes, filename: str) -> pd.DataFrame:
    """Validate the uploaded RNK file."""
    try:
        # Load the RNK file into a DataFrame
        rnk_data = pd.read_csv(io.BytesIO(file_content), sep="\t", header=None)
        
        # Ensure there are exactly two columns
        if rnk_data.shape[1] != 2:
            raise HTTPException(
                status_code=400,
                detail="Invalid RNK file format. Expected two columns (gene identifier and score)."
            )
        
        # Rename columns for consistency
        rnk_data.columns = ["GeneID", "Score"]
        
        # Convert GeneID to string and ensure Score is numeric and float32
        rnk_data["GeneID"] = rnk_data["GeneID"].astype(str)
        rnk_data["Score"] = pd.to_numeric(rnk_data["Score"], errors="coerce", downcast="float").astype(np.float32)
        
        # Check for invalid (non-numeric) values in Score
        if rnk_data["Score"].isna().any():
            raise HTTPException(
                status_code=400,
                detail="The 'Score' column in the RNK file must contain numeric values."
            )
        
        logger.info("RNK file validated successfully")
        return rnk_data
    except Exception as e:
        logger.error(f"Error processing RNK file: {str(e)}")
        raise HTTPException(
            status_code=400,
            detail=f"Error processing RNK file: {str(e)}"
        )


async def validate_xlsx_file(file_content: bytes, filename: str) -> pd.DataFrame:
    """Validate the uploaded Excel file."""
    logger.info("Validating uploaded Excel file")

    # Ensure the file is an Excel file
    if not filename.endswith((".xlsx", ".xls")):
        logger.error("Invalid file type: %s", filename)
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Please upload an Excel file.",
        )

    try:
        # Wrap the file contents in a BytesIO object and read into a DataFrame
        with io.BytesIO(file_content) as excel_data:
            df = read_prior_set(excel_data, is_bytes=True)
        
        logger.info("Excel file preprocessed successfully")

        # Expected columns in the file
        required_columns = ["GeneID", "Score"]
        
        # Check if the required columns are present
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"Missing required columns: {', '.join(missing_columns)}. "
                    f"Expected columns are: {', '.join(required_columns)}"
                ),
            )
        
        # Convert GeneID to string and ensure Score is numeric and float32
        df["GeneID"] = df["GeneID"].astype(str)
        df["Score"] = pd.to_numeric(df["Score"], errors="coerce", downcast="float").astype(np.float32)
        
        # Check for invalid (non-numeric) values in Score
        if df["Score"].isna().any():
            raise HTTPException(
                status_code=400,
                detail="The 'Score' column in the Excel file must contain numeric values.",
            )
        
        logger.info("Excel file validated successfully")
        return df
    except Exception as e:
        logger.error(f"Error processing Excel file: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Error processing Excel file: {str(e)}"
        )


async def validate_network_file(file_content: bytes, filename: str) -> pd.DataFrame:
    """
    Validate the uploaded network file.
    
    Parameters:
    - file_content: The content of the uploaded file in bytes
    - filename: Name of the uploaded file
    
    Returns:
    - pd.DataFrame: A validated DataFrame containing the network edges
    
    Raises:
    - HTTPException: If the file format is invalid or validation fails
    """
    try:
        # Read the network file into a DataFrame
        network_df = pd.read_csv(
            io.BytesIO(file_content),
            sep="\t",
            header=None,
            names=["Source", "Target"]
        )
        
        # Basic validation checks
        # Check number of columns
        if network_df.shape[1] != 2:
            raise HTTPException(
                status_code=400,
                detail="Invalid network file format. Expected exactly two columns (source and target genes)."
            )
        
        # Convert gene IDs to strings
        network_df["Source"] = network_df["Source"].astype(str)
        network_df["Target"] = network_df["Target"].astype(str)
        
        # Remove any self-loops
        network_df = network_df[network_df["Source"] != network_df["Target"]]
        
        # Remove duplicated edges (considering undirected network)
        network_df_reversed = network_df.copy()
        network_df_reversed.columns = ["Target", "Source"]
        network_df = pd.concat([network_df, network_df_reversed])
        network_df = network_df.drop_duplicates()
        
        # Check if network is empty after cleaning
        if network_df.empty:
            raise HTTPException(
                status_code=400,
                detail="Network file is empty or contains only invalid edges."
            )
        
        # Check for missing values
        if network_df.isna().any().any():
            raise HTTPException(
                status_code=400,
                detail="Network file contains missing values."
            )
        
        logger.info(f"Network file validated successfully: {len(network_df)} edges")
        return network_df
        
    except pd.errors.EmptyDataError:
        raise HTTPException(
            status_code=400,
            detail="The uploaded network file is empty."
        )
    except pd.errors.ParserError as e:
        raise HTTPException(
            status_code=400,
            detail=f"Error parsing network file: {str(e)}"
        )
    except Exception as e:
        logger.error(f"Error processing network file: {str(e)}")
        raise HTTPException(
            status_code=500,
            detail=f"Error processing network file: {str(e)}"
        )


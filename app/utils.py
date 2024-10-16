import io
import logging
from pathlib import Path

import pandas as pd
from fastapi import HTTPException

from pipeline.utils import read_prior_set

logger = logging.getLogger(__name__)


async def validate_file(
    file_content: bytes, filename: str
) -> pd.DataFrame:
    """Validate the uploaded file."""
    logger.info("Validating uploaded file")

    # Ensure the file is an Excel file
    if not filename.endswith((".xlsx", ".xls")):
        logger.error("Invalid file type: %s", filename)
        raise HTTPException(
            status_code=400,
            detail="Invalid file type. Please upload an Excel file.",
        )

    try:
        # Wrap the file contents in a BytesIO object
        with io.BytesIO(file_content) as excel_data:
            # Apply preprocessing to the DataFrame
            df = read_prior_set(excel_data, is_bytes=True)
        logger.info("File preprocessed successfully")

        # Expected columns in the file
        required_columns = ["GeneID", "Symbol", "Score", "P-value"]

        # Check if the required columns are present
        missing_columns = [
            col for col in required_columns if col not in df.columns
        ]
        if missing_columns:
            logger.error(
                "Missing required columns: %s",
                ", ".join(missing_columns),
            )
            raise HTTPException(
                status_code=400,
                detail=(
                    f"Missing required columns: {', '.join(missing_columns)}. "
                    f"Expected columns are: {', '.join(required_columns)}"
                ),
            )

        # Validate data types
        if not pd.api.types.is_numeric_dtype(df["Score"]):
            logger.error("The 'Score' column must be numeric")
            raise HTTPException(
                status_code=400,
                detail="The 'Score' column must contain numeric values.",
            )
        if not pd.api.types.is_numeric_dtype(df["P-value"]):
            logger.error("The 'P-value' column must be numeric")
            raise HTTPException(
                status_code=400,
                detail="The 'P-value' column must contain numeric values.",
            )

        logger.info("File validation successful")
        return df

    except Exception as e:
        logger.error(f"Error processing file: {str(e)}")
        raise HTTPException(
            status_code=500, detail=f"Error processing file: {str(e)}"
        )

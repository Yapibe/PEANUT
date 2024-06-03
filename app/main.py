from fastapi import FastAPI, UploadFile, File
from pydantic import BaseModel
import pandas as pd
import os
from pipeline.utils import load_data, process_data

app = FastAPI()

class PathwayAnalysisInput(BaseModel):
    experiment_data: str  # This could be a file path or actual data
    alpha: float = 1.0
    min_gene_per_pathway: int = 5
    max_gene_per_pathway: int = 500
    fdr_threshold: float = 0.05
    jac_threshold: float = 0.3
    p_value_threshold: float = 0.05

@app.post("/analyze/")
async def analyze_pathway(input_data: PathwayAnalysisInput):
    # Assuming load_data and process_data are functions in your utils.py
    experiment_df = load_data(input_data.experiment_data)
    result = process_data(
        experiment_df,
        alpha=input_data.alpha,
        min_genes=input_data.min_gene_per_pathway,
        max_genes=input_data.max_gene_per_pathway,
        fdr_thresh=input_data.fdr_threshold,
        jac_thresh=input_data.jac_threshold,
        p_val_thresh=input_data.p_value_threshold
    )
    return result

@app.post("/uploadfile/")
async def upload_file(file: UploadFile = File(...)):
    file_location = f"inputs/{file.filename}"
    with open(file_location, "wb+") as file_object:
        file_object.write(file.file.read())
    return {"info": "file saved at " + file_location}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse
from app.utils import run_pipeline

router = APIRouter()

@router.post("/run/")
async def run_analysis():
    try:
        result_path = run_pipeline()
        return FileResponse(result_path, filename="results.zip")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

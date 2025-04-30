import logging
import logging.config
import sys
import time
from pathlib import Path
from fastapi import FastAPI, Request, Response, status
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi.staticfiles import StaticFiles
from uvicorn.middleware.proxy_headers import ProxyHeadersMiddleware
from .routes import router

# Configure both application and access logging
logging.basicConfig(
    level=logging.INFO,
    handlers=[
        logging.FileHandler("app.log"),
        logging.StreamHandler(sys.stdout)
    ],
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

# Explicitly configure uvicorn access logger
uvicorn_access_logger = logging.getLogger("uvicorn.access")
uvicorn_access_logger.setLevel(logging.INFO)
if not uvicorn_access_logger.handlers:
    uvicorn_access_logger.addHandler(logging.StreamHandler(sys.stdout))
    uvicorn_access_logger.addHandler(logging.FileHandler("access.log"))

logger = logging.getLogger(__name__)
logger.debug("Logging has been configured with basicConfig")

# FastAPI app lifespan
async def lifespan(app: FastAPI):
    """Manage the lifespan of the FastAPI app."""
    logger.info("Application is starting up")
    yield
    logger.info("Application is shutting down")

# Middleware
class AddCacheControlHeadersMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        response = await call_next(request)
        if request.url.path.startswith("/static/"):
            response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        return response

class RequestLoggingMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request: Request, call_next):
        start_time = time.time()
        try:
            response = await call_next(request)
            process_time = time.time() - start_time
            logger.info(
                f"Request: {request.method} {request.url.path} - Status: {response.status_code} - Time: {process_time:.4f}s"
            )
            return response
        except Exception as e:
            process_time = time.time() - start_time
            logger.error(
                f"Request: {request.method} {request.url.path} - Error: {str(e)} - Time: {process_time:.4f}s"
            )
            return Response(
                content={"detail": "Internal server error"},
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                media_type="application/json"
            )

middleware = [
    Middleware(ProxyHeadersMiddleware),
    Middleware(AddCacheControlHeadersMiddleware),
    Middleware(RequestLoggingMiddleware),
]

# FastAPI app initialization
app = FastAPI(
    lifespan=lifespan,
    middleware=middleware,
    root_path="/peanut",
)

# Health check endpoint
@app.get("/health")
async def health_check():
    """Simple health check endpoint to verify service is running"""
    return {"status": "healthy", "service": "PEANUT"}

# Mount static files directly on the app
static_dir = Path(__file__).resolve().parent / "static"

app.mount(
    "/static",
    StaticFiles(directory=str(static_dir)),
    name="static",
)
logger.info(f"Static files will be served from: {static_dir}")
logger.info(f"Contents of static directory: {list(static_dir.iterdir())}")

# Include routes
app.include_router(router)
import logging
import logging.config
import sys
from pathlib import Path
from fastapi import FastAPI
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi.staticfiles import StaticFiles
from uvicorn.middleware.proxy_headers import ProxyHeadersMiddleware  # Import from uvicorn
from .routes import router

# Configure both application and access logging
logging.basicConfig(
    level=logging.DEBUG,
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

middleware = [
    Middleware(ProxyHeadersMiddleware),  # Use Uvicorn's ProxyHeadersMiddleware
    Middleware(AddCacheControlHeadersMiddleware),
]

# FastAPI app initialization
app = FastAPI(
    lifespan=lifespan,
    middleware=middleware,
    root_path="/peanut",
)

# Mount static files directly on the app
static_dir = Path(__file__).resolve().parent / "static"

app.mount(
    "/static",
    StaticFiles(directory=str(static_dir)),
    name="static",
)

# Include routes
app.include_router(router)
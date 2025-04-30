from pathlib import Path
from fastapi import FastAPI, Request, Response, status
from fastapi.middleware.cors import CORSMiddleware
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi.staticfiles import StaticFiles
from uvicorn.middleware.proxy_headers import ProxyHeadersMiddleware
from .routes import router
from .utils import setup_logging

# Setup logging
logger = setup_logging()

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

class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["X-XSS-Protection"] = "1; mode=block"
        response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"
        return response

middleware = [
    Middleware(ProxyHeadersMiddleware),
    Middleware(AddCacheControlHeadersMiddleware),
    Middleware(RequestLoggingMiddleware),
    Middleware(SecurityHeadersMiddleware),
]

# FastAPI app initialization
app = FastAPI(
    lifespan=lifespan,
    middleware=middleware,
    root_path="/peanut",
    title="PEANUT",
    description="Pathway Enrichment Analysis Using Network Propagation",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, replace with specific origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health check endpoint
@app.get("/health")
async def health_check():
    """Simple health check endpoint to verify service is running"""
    return {"status": "healthy", "service": "PEANUT"}

# Mount static files directly on the app
static_dir = Path(__file__).resolve().parent / "static"
if not static_dir.exists():
    logger.warning(f"Static directory not found: {static_dir}")
    static_dir.mkdir(parents=True, exist_ok=True)

app.mount(
    "/static",
    StaticFiles(directory=str(static_dir)),
    name="static",
)
logger.info(f"Static files will be served from: {static_dir}")
logger.info(f"Contents of static directory: {list(static_dir.iterdir())}")

# Include routes
app.include_router(router)
"""
Main FastAPI application entry point for PEANUT.

This module initializes the FastAPI application with all routes and middleware.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any

from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from uvicorn.middleware.proxy_headers import ProxyHeadersMiddleware
import traceback

from .routes import router
from .utils import setup_logging
from .middleware import (
    AddCacheControlHeadersMiddleware,
    RequestLoggingMiddleware,
    SecurityHeadersMiddleware,
)

# Load configuration from YAML
def load_config() -> Dict[str, Any]:
    """Load configuration from YAML file with environment variable overrides."""
    config_path = Path(__file__).resolve().parent.parent / "config" / "app_config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    # Add computed paths (not in YAML)
    base_dir = Path(__file__).resolve().parent
    config["base_dir"] = base_dir
    config["static_dir"] = base_dir / "static"
    config["sample_data_dir"] = base_dir / "static" / "sample_data"
    
    return config

# Load configuration
config = load_config()

# Initialize app without logging - will be set up in startup event
app = FastAPI(
    title=config["app_name"],
    description=config["app_description"],
    version=config["app_version"],
    docs_url=config["docs_url"],
    redoc_url=None,  # Disable ReDoc UI
    root_path=config["root_path"],
)

# Store logger and config as app state attributes
app.state.logger = None
app.state.config = config


@app.on_event("startup")
async def startup_event() -> None:
    """Initialize application resources on startup."""
    # Set up logging
    app.state.logger = setup_logging()
    app.state.logger.info(f"Starting {config['app_name']} application")
    
    # Ensure static directories exist
    config["static_dir"].mkdir(parents=True, exist_ok=True)
    config["sample_data_dir"].mkdir(parents=True, exist_ok=True)
    app.state.logger.info(f"Static directory configured: {config['static_dir']}")


@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Handle all unhandled exceptions globally."""
    if app.state.logger:
        app.state.logger.error(
            f"Unhandled exception in {request.method} {request.url.path}: {str(exc)}"
        )
        app.state.logger.error(traceback.format_exc())
    
    # Don't expose internal details in production
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={"detail": "Internal server error"},
    )


# Register middleware (order matters)
app.add_middleware(ProxyHeadersMiddleware)
app.add_middleware(AddCacheControlHeadersMiddleware)
app.add_middleware(RequestLoggingMiddleware)
app.add_middleware(SecurityHeadersMiddleware)

# Add CORS middleware with configurable origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=config["cors_origins"],
    allow_credentials=config["allow_credentials"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount(
    config["static_url_prefix"],
    StaticFiles(directory=str(config["static_dir"])),
    name="static",
)

# Simple health check endpoint
@app.get("/health")
async def health_check() -> Dict[str, str]:
    """Simple health check endpoint."""
    return {"status": "healthy", "service": config["app_name"]}

# Include routes
app.include_router(router)
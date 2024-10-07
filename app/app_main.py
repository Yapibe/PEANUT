import logging
import logging.config
import yaml
from fastapi import FastAPI
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from .routes import router

try:
    # Load logging configuration from the config directory
    with open("config/log_config.yaml", "r") as f:
        log_config = yaml.safe_load(f.read())
        logging.config.dictConfig(log_config)
except Exception as e:
    print(f"Failed to load log configuration: {e}")
    # Fallback to basic logging
    logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)


async def lifespan(app: FastAPI):
    """Manage the lifespan of the FastAPI app."""
    logger.info("Application is starting up")
    yield
    logger.info("Application is shutting down")


class AddCacheControlHeadersMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        response = await call_next(request)
        if request.url.path.startswith("/static/"):
            response.headers["Cache-Control"] = (
                "no-cache, no-store, must-revalidate"
            )
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        return response


middleware = [Middleware(AddCacheControlHeadersMiddleware)]

app = FastAPI(
    lifespan=lifespan,
    middleware=middleware,
    root_path="/peanut",  # Added root_path to handle the proxy prefix
)

# Include the routes from routes.py
app.include_router(router)
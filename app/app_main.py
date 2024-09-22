import logging
from fastapi import FastAPI
from .routes import router
from contextlib import asynccontextmanager
from starlette.middleware import Middleware
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("Application is starting up")
    yield
    logger.info("Application is shutting down")


class AddCacheControlHeadersMiddleware(BaseHTTPMiddleware):
    async def dispatch(self, request, call_next):
        response = await call_next(request)
        if request.url.path.startswith("/static/"):
            response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        return response


middleware = [Middleware(AddCacheControlHeadersMiddleware)]

app = FastAPI(lifespan=lifespan, middleware=middleware)

# Include the routes from routes.py
app.include_router(router)

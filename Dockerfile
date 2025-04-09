# Use the full Python runtime image
FROM python:3.10

# Set the working directory in the container
WORKDIR /app

# Install curl for healthcheck
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy the requirements file and install dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . /app

# Copy the config directory
COPY config /app/config

# Create necessary directories and set permissions
RUN mkdir -p /app/pipeline/Outputs/logs && \
    touch /app/pipeline/Outputs/logs/app.log && \
    touch /app/pipeline/Outputs/logs/access.log && \
    chmod -R 775 /app/pipeline/Outputs/logs && \
    mkdir -p /app/app/static/sample_data && \
    chown -R root:root /app

# Expose the port your app runs on
EXPOSE 8000

# Set default root path (can be overridden at runtime)
ENV ROOT_PATH=/peanut

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PORT=8000
ENV WORKERS=4
ENV TIMEOUT=300

# Add health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:${PORT}/peanut/health || exit 1

# Command to run the application
CMD uvicorn app.app_main:app --host=0.0.0.0 --port=${PORT} --workers=${WORKERS} --timeout-keep-alive=${TIMEOUT}

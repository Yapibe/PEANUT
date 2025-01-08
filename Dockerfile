# Use the full Python runtime image
FROM python:3.10

# Set the working directory in the container
WORKDIR /app

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
    chown -R root:root /app

# Expose the port your app runs on
EXPOSE 8000

# Set default root path (can be overridden at runtime)
ENV ROOT_PATH=/peanut

# Command to run your application
CMD ["sh", "-c", "uvicorn app.app_main:app --host 0.0.0.0 --port 8000 --root-path ${ROOT_PATH} --log-config /app/config/log_config.yaml"]

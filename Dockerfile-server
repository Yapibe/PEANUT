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

# Create the Outputs/logs directory
RUN mkdir -p /app/pipeline/Outputs/logs

# Expose the port your app runs on
EXPOSE 8000

# Command to run your application with proxy support
CMD ["uvicorn", "app.app_main:app", "--host", "0.0.0.0", "--port", "8000", "--root-path", "/peanut", "--log-config", "/app/config/log_config.yaml"]

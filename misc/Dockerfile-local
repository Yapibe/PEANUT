# Use the full Python runtime image
FROM python:3.10

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file and install dependencies
COPY requirements.txt ./ 
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . /app

# Include the entire pathways and networks directories
COPY pipeline/Data/H_sapiens/pathways /app/pipeline/Data/H_sapiens/pathways
COPY pipeline/Data/H_sapiens/network /app/pipeline/Data/H_sapiens/network

# Include only the specific file from the matrix directory
COPY pipeline/Data/H_sapiens/matrix/Anat_0.1.npz /app/pipeline/Data/H_sapiens/matrix/

# Create the Outputs/logs directory
RUN mkdir -p /app/pipeline/Outputs/logs

# Expose the port your app runs on
EXPOSE 8000

# Command to run your application locally with reload support
CMD ["uvicorn", "app.app_main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]

# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements file and install dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . /app

# Expose the port your app runs on
EXPOSE 8000

# Command to run your application
CMD ["uvicorn", "app.app_main:app", "--host", "0.0.0.0", "--port", "8000"]

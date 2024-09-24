# PEANUT: Pathway Enrichment Analysis through Network UTilization

## Overview

PEANUT is a web-based tool designed to perform pathway enrichment analysis on RNA-seq data using network propagation techniques. By leveraging protein-protein interaction (PPI) networks, PEANUT helps identify significant biological pathways associated with a given experiment. The web interface allows easy customization of parameters and processing of multiple conditions, while providing downloadable results and visualizations.

## Features
- **Web Interface for Ease of Use**: Upload gene sets and customize parameters directly through a web interface.
- **Multiple Condition Support**: Analyze multiple gene sets in one run and compare the results across different conditions.
- **Network Propagation**: Propagate gene scores through a PPI network to enhance biological data analysis.
- **Statistical Enrichment**: Perform statistical tests (Kolmogorov-Smirnov, Mann-Whitney) for pathway enrichment, with False Discovery Rate (FDR) correction.
- **Downloadable Results**: Results and visualizations are packaged and downloadable as a zip file.
- **Pathway Comparison**: View and compare pathway trends across conditions in generated plots.

## Installation

1. **Clone the repository**:
   ```sh
   git clone https://github.com/Yapibe/PEANUT.git
   ```
2. **Install required packages**:
   Ensure you have Python 3.8 or higher. Install dependencies using pip:
   ```sh
   pip install -r requirements.txt
   ```
3. **Run the application**:
   Use FastAPI to start the web server
   ```sh
   uvicorn app.app_main:app --reload
   ```

4. **Access the web interface**:
    Open a browser and go to `http://127.0.0.1:8000` to access the PEANUT interface.

## Web Interface Usage

1. **Start a New Job**:
    - Enter a test name and upload one or more preranked gene set files (Excel format).
    - Select species, network type, propagation alpha, and gene set (e.g., KEGG, Reactome).
    - (Optional) Customize advanced settings, including minimum/maximum genes per pathway, FDR threshold, and Jaccard index threshold.
    - Submit the form to start the analysis.
2. **Monitor Job Status**:
    - After submitting a job, you will receive a job code to check the status.
    - Use the provided job code to monitor progress. Once the job is completed, you can download the results as a zip file.
3. **Download Results**:
    - The results will include significant pathways and gene information for each condition, along with comparative visualizations in a PNG file.

## Backend Usage (Optional)

1. **Configure the Program**: Modify config.yaml to specify input/output directories and statistical parameters for advanced usage outside the web interface.

2. **Run the Backend Pipeline**: For those who prefer command-line execution:
   ```sh
   python pipeline/pipeline_main.py
   ```

## Pipeline Overview
### Web Pipeline
- **Gene Set Upload**: Upload multiple preranked tables for analysis.
- **Network Propagation**: Propagate gene scores through the PPI network to identify enriched pathways.
- **Pathway Enrichment**: Perform statistical tests (Kolmogorov-Smirnov, Mann-Whitney) to assess significance.
- **Visualization**: Generate plots comparing pathway significance across conditions.


### Backend Pipeline (Detailed)
For users preferring command-line execution, the pipeline consists of:
1. **Data Loading**: Load experimental data from CSV files with columns for GeneID, Symbol, Score, and P-value.
2. **Network Propagation**: Propagate gene scores through a PPI network and filter based on significance.
3. **Statistical Enrichment**: Perform statistical tests to identify significant pathways.
4. **Result Compilation**: Compile and save results and plots for each condition.


### Constants and parameters
- **Thresholds:**
  - `Alpha`: Hyper parameter for the network propagation. Default Alpha = 1 means no propagation.
  - `MIN_GENE_PER_PATHWAY`: Minimum number of genes per pathway to consider.
  - `MAX_GENE_PER_PATHWAY`: Maximum number of genes per pathway to consider.
  - `FDR_THRESHOLD`: FDR threshold for significance.
  - `JAC_THRESHOLD`: Jaccard index threshold for comparing sets.
  - `P_VALUE_THRESHOLD`: P-value threshold for statistical significance.

### Directory Structure
- **YAIR_PROPAGATION**
  - **app/**: Contains the FastAPI web server and associated routes.
    - `static/`: Static assets for the web interface.
    - `app_main.py`: FastAPI entry point.
    - `models.py`: Defines the data models.
    - `pipeline_runner.py`: Executes the backend pipeline.
    - `routes.py`: Handles file uploads, job status, and result download routes.
    - `utils.py`: Utility functions for data processing and file handling.
  - **pipeline/**: Backend code for pathway analysis.
    - **Data/**
      - `Human`: contains gene names, pathways, network, and matrix data
    - **Inputs/**
      - `experiments_data`: contains experimental data for each condition
    - **Outputs/**: Generates plots for pathway comparison.
    - `pipeline_main.py`: Main script for pathway analysis.
    - `settings.py`: Contains settings for the pipeline.
    - `propagation_routines.py`: Handles network propagation.
    - `pathway_enrichment.py`: Contains statistical analysis methods.
    - `visualization_tools.py`: Generates plots for pathway comparison.
    - `utils.py`: Utility functions for data processing and file handling.
    - `statistical_methods.py`: Contains statistical analysis methods.
    
### Statistical Methods
- **Hypergeometric Test:** <br>
  Used to identify statistically significant pathways by comparing the observed number of genes in a pathway to the expected number under a null model.

- **Kolmogorov-Smirnov Test:** <br>
    Used to score each pathway by assessing if the expression changes of its genes deviate significantly from all other genes.

- **Mann-Whitney Test:** <br>
    Non-parametric test to compare differences between two independent groups.

- **FDR Correction:** <br>
    Adjusts p-values to account for multiple testing, controlling the false discovery rate.

## Plot Explanation
The output plots display pathway trends across conditions:
- **Colored Bars**: Show the mean score for each pathway in a given condition.
- **Solid Bars**:  Indicate significant pathways based on the FDR threshold.
- **Text File**: Includes full details on p-values, trends, and significant genes.

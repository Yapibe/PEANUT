# PEANUT: Pathway Enrichment Analysis through Network UTilization

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
- [Usage](#usage)
  - [Web Interface Usage](#web-interface-usage)
  - [Backend Usage](#backend-usage)
- [Features](#features)
- [Troubleshooting](#troubleshooting)
- [Citing PEANUT](#citing-peanut)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

## Overview
PEANUT is a web-based tool designed to perform pathway enrichment analysis on RNA-seq data using network propagation techniques. By leveraging protein-protein interaction (PPI) networks, PEANUT helps identify significant biological pathways associated with a given experiment. The tool integrates network propagation with a robust statistical framework, including permutation tests and multiple comparison corrections, to enhance the detection of relevant pathways.

The PEANUT tool is available at [GitHub Repository](https://github.com/Yapibe/PEANUT) and has been assigned a Zenodo DOI: [10.5281/zenodo.15184862](https://doi.org/10.5281/zenodo.15184862).

## Installation

### Prerequisites
- Python 3.8 or higher
- Required libraries (listed in `requirements.txt`)

1. **Clone the repository**:
   ```sh
   git clone https://github.com/Yapibe/PEANUT.git
   ```
2. **Install required packages**:
   Ensure you have Python 3.8 or higher. Install dependencies using pip:
   ```sh
   pip install -r requirements.txt
   ```
   The following packages are required:
   - `fastapi==0.115.0`
   - `gseapy==1.1.3`
   - `matplotlib==3.8.2`
   - `networkx==3.2.1`
   - `numpy>=1.21,<2`
   - `pandas==2.2.3`
   - `pydantic==2.9.2`
   - `PyYAML==6.0.2`
   - `scipy==1.14.1`
   - `starlette>=0.37.2,<0.39.0`
   - `statsmodels==0.14.2`
   - `uvicorn==0.25.0`
   - `python-multipart==0.0.6`
   - `openpyxl==3.1.2`
   - `jinja2==3.1.2`
3. **Run the application**:
   ```sh
   uvicorn app.app_main:app --reload
   ```
4. **Access the web interface**:
   Open a browser and go to `http://127.0.0.1:8000` to access the PEANUT interface.

## Usage

### Web Interface Usage
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

### Backend Usage (Optional)
1. **Configure the Program**: Modify `config.yaml` to specify input/output directories and statistical parameters for advanced usage outside the web interface.
2. **Run the Backend Pipeline**:
   ```sh
   python pipeline/pipeline_main.py
   ```

## Features
- **Web Interface for Ease of Use**: Upload gene sets and customize parameters directly through a web interface.
- **Multiple Condition Support**: Analyze multiple gene sets in one run and compare the results across different conditions.
- **Network Propagation**: Propagate gene scores through a PPI network to enhance biological data analysis.
- **Statistical Enrichment**: Perform statistical tests (Kolmogorov-Smirnov, Mann-Whitney, permutation tests) for pathway enrichment, with False Discovery Rate (FDR) correction using the Benjamini-Hochberg method.
- **Downloadable Results**: Results and visualizations are packaged and downloadable as a zip file.
- **Pathway Comparison**: View and compare pathway trends across conditions in generated plots.

## Troubleshooting
- **Common Issues**: List common installation or usage issues and their solutions.

## Citing PEANUT
If you use PEANUT in your research, please cite it as follows:
- DOI: [10.5281/zenodo.15184862](https://doi.org/10.5281/zenodo.15184862)

## Contributing
We welcome contributions! Please see our [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For support or inquiries, please contact [your_email@example.com].

## Acknowledgments
- Acknowledge any funding sources, collaborators, or institutions that contributed to the development of the tool.

## Pipeline Overview

### Web Pipeline
- **Gene Set Upload**: Upload multiple preranked tables for analysis.
- **Network Propagation**: Propagate gene scores through the PPI network to identify enriched pathways.
- **Pathway Enrichment**: Perform statistical tests (Kolmogorov-Smirnov, Mann-Whitney, permutation tests) to assess significance.
- **Visualization**: Generate plots comparing pathway significance across conditions.

### Backend Pipeline (Detailed)
For users preferring command-line execution, the pipeline consists of:
1. **Data Loading**: Load experimental data from CSV files with columns for GeneID, Symbol, Score, and P-value.
2. **Network Propagation**: Propagate gene scores through a PPI network and filter based on significance.
3. **Statistical Enrichment**: Perform statistical tests to identify significant pathways.
4. **Result Compilation**: Compile and save results and plots for each condition.

### Constants and Parameters
- **Thresholds:**
  - `Alpha`: Hyperparameter for the network propagation. Default Alpha = 0.2.
  - `MIN_GENE_PER_PATHWAY`: Minimum number of genes per pathway to consider.
  - `MAX_GENE_PER_PATHWAY`: Maximum number of genes per pathway to consider.
  - `FDR_THRESHOLD`: FDR threshold for significance.
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
      - `H_sapiens`: Contains gene names, pathways, network, and matrix data.
    - **Inputs/**
      - `experiments_data`: Contains experimental data for each condition.
    - **Outputs/**: Generates plots for pathway comparison.
    - `pipeline_main.py`: Main script for pathway analysis.
    - `settings.py`: Contains settings for the pipeline.
    - `propagation_routines.py`: Handles network propagation.
    - `pathway_enrichment.py`: Contains statistical analysis methods.
    - `visualization_tools.py`: Generates plots for pathway comparison.
    - `utils.py`: Utility functions for data processing and file handling.
    - `statistical_methods.py`: Contains statistical analysis methods.

### Statistical Methods
- **Kolmogorov-Smirnov Test:** Used to score each pathway by assessing if the expression changes of its genes deviate significantly from all other genes.
- **Mann-Whitney Test:** Non-parametric test to compare differences between two independent groups.
- **Permutation Test:** Empirically evaluates the significance of observed pathway scores by generating a null distribution.
- **FDR Correction:** Adjusts p-values to account for multiple testing, controlling the false discovery rate.

## Plot Explanation
The output plots display pathway trends across conditions:
- **Colored Bars**: Show the mean score for each pathway in a given condition.
- **Solid Bars**: Indicate significant pathways based on the FDR threshold.
- **Text File**: Includes full details on p-values, trends, and significant genes.

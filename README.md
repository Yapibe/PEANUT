# Network Propagation to Enhance Biological Data

## Overview

This program is designed for pathway analysis in RNA-seq data, focusing on filtering and statistical enrichment analysis within a biological protein-protein interaction (PPI) network. It leverages various statistical methods to identify significant pathways associated with a given experiment, particularly tailored for Parkinson's disease data analysis.

## Features

- **Pathway Data Loading**: Loads and parses pathway data.
- **Gene Filtering**: Applies thresholds to filter genes based on pathway size.
- **Statistical Analysis**: Performs hypergeometric and Mann-Whitney tests to identify significant pathways.
- **Multiple Testing Correction**: Adjusts p-values using the False Discovery Rate (FDR) method.
- **Output Generation**: Saves results in specified directories for further analysis.

## Installation

1. **Clone the repository**:
   ```sh
   git clone https://github.com/Yapibe/Yair_propagation.git
    ```
2. **Install required packages**:
   Ensure you have Python 3.8 or higher. Install dependencies using pip:
   ```sh
    pip install -r requirements.txt
    ```
   requirements.txt should include:
   ```plaintext
    numpy==1.21.2
    pandas==1.3.3
    scipy==1.7.3
    statsmodels==0.12.2
    networkx==2.6.3
    matplotlib==3.4.3
    seaborn==0.11.2
    ```
   
## Usage

1. **Configure the program**:
   Modify the configuration file `config.yaml` to specify input and output directories, as well as filtering thresholds and statistical parameters.
2. **Run the program**:
   Execute the main script `main.py` to perform pathway analysis on the specified data.
   ```sh
   python main.py
   ```
3. **View the results**:
   Results will be saved in the Outputs directory, including filtered pathways and statistical analysis results.

## Detailed Explanation

### Pipeline Stages
#### Pathway Filtering
Initially, pathways are filtered based on their overlap with differentially expressed genes (DEGs) to focus on those
most relevant to the study's conditions using a hypergeometric test.

#### Statistical Tests for Enrichment
Each pathway is scored using the Kolmogorov-Smirnov test to assess if the expression changes of its genes deviate
significantly from all other genes. Those pathways which are found significant after an FDR correction
(p<0.05) are maintained.

#### Additional Filtering
To control for potential differences between genes in annotated pathways vs. other genes, we perform a 
final Mann-Whitney test with the remaining pathways to compare the distribution of their expression changes to
that of other annotated genes. Pathways with FDR-corrected p-values smaller than 0.05 are reported.

### Constants and parameters
- **Thresholds:**
  - `MIN_GENE_PER_PATHWAY`: Minimum number of genes per pathway to consider.
  - `MAX_GENE_PER_PATHWAY`: Maximum number of genes per pathway to consider.
  - `FDR_THRESHOLD`: FDR threshold for significance.
  - `JAC_THRESHOLD`: Jaccard index threshold for comparing sets.
  - `P_VALUE_THRESHOLD`: P-value threshold for statistical significance.

### Directory Structure
- **Data:**
  - `H_sapiens`:
    - `gene_names`: contains dictionaries of gene names and their corresponding IDs.
    - `pathways`: contains lists of pathways and their corresponding genes.
    - `network`: contains the PPI network data.
- **Inputs:**
  - `experiments_data/`:
    - `Parkinson`: contains the Parkinson's disease data.
- **Outputs:**
  - `Plots`: contains plots generated during the analysis.
  - `Text`: contains text files with the results of the analysis.

### Statisical Methods
- **Hypergeometric Test:** <br>
  Used to identify statistically significant pathways by comparing the observed number of genes in a pathway to the expected number under a null model.

- **Kolmogorov-Smirnov Test:** <br>
    Used to score each pathway by assessing if the expression changes of its genes deviate significantly from all other genes.

- **Mann-Whitney Test:** <br>
    Non-parametric test to compare differences between two independent groups.

- **FDR Correction:** <br>
    Adjusts p-values to account for multiple testing, controlling the false discovery rate.

## Plot Explanation

For each pathway in the generated plots, you can see its trends for each condition. If the bar is solid, it means that 
it's significantly changing in that condition. In the accompanying text file, you can find full information on
pathway p-values, trends, and scores of significant genes that are part of the pathways.
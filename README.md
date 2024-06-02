# Network Propagation to Enhance Biological Data

## Overview

This program is designed for pathway analysis in RNA-seq data, focusing on filtering and statistical enrichment analysis within a biological protein-protein interaction (PPI) network. It leverages various statistical methods and network propagation to identify significant pathways associated with a given experiment.

## Features
- **Experimental Data Loading**: Reads and processes RNA-seq data from CSV files.
- **Network Propagation**: Propagates gene scores through a PPI network to enhance biological data analysis.
- **Pathway Data Loading and filtering**: Loads and parses pathway data and filters based on various criteria.
- **Statistical Analysis**: Performs a series of statistical tests to identify significant pathways.
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
####  Data Loading
Load the experimental data from the specified input directory. The data should be in the form of a CSV file with columns
for GeneID, Symbol, Score and P-value for each condition.

#### Propagation 
Propagate gene scores through a protein-protein interaction (PPI) network to obtain updated scores which are saved for
further analysis. The default network is the Human PPI network which is loaded from the data directory.

#### Pathway Filtering
Initially, pathways are filtered based on their size and overlap with differentially expressed genes (DEGs) to focus on those
most relevant to the study's conditions using a hypergeometric test.

#### Statistical Tests for Enrichment
Each pathway is scored using the Kolmogorov-Smirnov test to assess if the expression changes of its genes deviate
significantly from all other genes. Those pathways which are found significant after an FDR correction
(p<0.05) are maintained.

#### Additional Filtering
To control for potential differences between genes in annotated pathways vs. other genes, we perform a 
final Mann-Whitney test with the remaining pathways to compare the distribution of their expression changes to
that of other annotated genes. Pathways with FDR-corrected p-values smaller than 0.05 are reported.

#### Propagation Steps
1. **Load Prior Data:** Load the initial gene scores and P-values from the experimental data.
2. **Filter Network:** Filter the PPI network to include only genes present in the prior data.
3. **Create or Load Similarity Matrix:** Generate or load the similarity matrix for the network.
4. **Propagate Scores:** Propagate the gene scores through the network.
5. **Normalize Scores:** Normalize the propagated scores using ones vector.
6. **Save Results:** Save the propagated scores and the updated gene data.

#### Enrichment Steps
1. **Identify Significant Genes:** Identify genes with P-values below the significance threshold.
2. **Filter Pathways:** Filter pathways based on gene count criteria.
3. **Hypergeometric Test:** Calculate hypergeometric P-values for each pathway with enough genes.
4. **Kolmogorov-Smirnov Test:** Perform the KS test to compare distributions of scores.
5. **Mann-Whitney Test:** Perform the Mann-Whitney test for additional filtering.
6. **Adjust P-values:** Apply Benjamini-Hochberg correction to the P-values.
7. **Save Results:** Save the significant pathways and their associated genes.

### Constants and parameters
- **Thresholds:**
  - `Alpha`: Hyper parameter for the network propagation. Default Alpha = 1 means no propagation.
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

For each pathway in the generated plots, you can see its trends for each condition. If the bar is solid, it means that 
it's significantly changing in that condition. In the accompanying text file, you can find full information on
pathway p-values, trends, and scores of significant genes that are part of the pathways.
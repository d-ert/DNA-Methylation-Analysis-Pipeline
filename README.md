# Methylation Array Analysis for Differential Methylation

## Project Overview

This repository contains an R-based analysis of Illumina HumanMethylation450k array data to identify differentially methylated probes (DMPs) between a disease and a control group. The project's goal is to uncover epigenetic variations that could serve as potential biomarkers for the studied condition.

### Key Features:

*   **Data Quality Control:** Rigorous quality checks to ensure data integrity.
*   **Normalization:** Application of the `preprocessNoob` method to correct for background noise and dye-bias.
*   **Statistical Analysis:** Identification of DMPs using the `dmpFinder` function from the `minfi` package.
*   **Visualization:** Generation of volcano plots, Manhattan plots, and heatmaps to interpret the results.


## Getting Started

To replicate this analysis, clone the repository and follow the steps below.

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your_username/your_repository.git
    ```
2.  **Install R packages:**
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install(c("minfi", "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19"))

    install.packages(c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel", "knitr", "tinytex"))
    ```

### Execution

1.  Place your raw data files (`.idat`) in the `Input_Data/` directory.
2.  Open the `methylation_analysis.Rmd` file in RStudio.
3.  Run the code chunks sequentially to perform the analysis.

# ComeBack
Pipeline to compute Methylation Risk Scores (MRS) using co-methylation pruning (ComeBack) and compare it to standard methods. Includes CpG selection, score generation from EWAS summary stats, predictive evaluation, and visualization.

# ComeBack-MRS: Pipeline for Generating Co-Methylation-Based Methylation Risk Scores

This repository contains a complete and reproducible pipeline to compute **Methylation Risk Scores (MRS)** using the **ComeBack method** (co-methylation pruning), and to compare it with conventional methods like **Methylscore**.

The workflow includes:
- Preprocessing of DNA methylation beta matrices  
- Identification of Co-Methylated Regions (CMRs)  
- Generation of MRS using EWAS summary statistics  
- Evaluation of predictive performance (R²)  
- Comparison between methods (ComeBack vs Methylscore)  
- Sensitivity analysis across p-value thresholds

## Requirements

- R version ≥ 4.0  
- R packages:  
  `comeback`, `robustHD`, `methylscore`, `data.table`, `dplyr`, `readxl`, `tibble`, `ggplot2`, `pROC`, `scales`, etc.

You will also need:
- A matrix of DNA methylation beta values (rows = CpGs, columns = samples)  
- EWAS summary statistics including `Marker`, `BETA`, `SE`, and `Pvalue` columns  

## Pipeline Overview
The ComeBack R package can be downloaded and manually installed from: https://bitbucket.org/flopflip/comeback/src/master/comeback_1.0.1.tar.gz

As explained in the MRS_func.R script on GitHub (https://github.com/jche453/Pruning-Thresholding-MRS/blob/main/MRS_func.R), this function provides the base logic to generate MRS using optional co-methylation pruning. 

### 🔹 Step 1: Load Data and Define Co-Methylated Regions (CMRs)

To define Co-Methylated Regions (CMRs), we use the cmr() function from the ComeBack R package, which performs correlation-based pruning of CpGs. This groups CpGs into regions based on both their methylation similarity and genomic proximity.

- **Input data required**
You must first load your beta matrix (methylation values, with CpGs as columns and samples as rows).

- **Running cmr**
  - You can customize two key parameters:
      -corlo: the minimum correlation threshold between CpGs to be grouped together. Example values: 0.3 (less strict, more CpGs grouped)         or 0.5 (more strict, fewer CpGs grouped). Higher values make the selection more conservative, requiring stronger co-methylation.
      - maxprbdst: the maximum allowed physical distance (in base pairs) between two CpGs to be considered part of the same CMR. Default          is often set to 2000 (i.e., CpGs within 2kb).
        
- **List of Co-Methylated Regions** (using GenCoMeRegion)
This allows you to check the number of regions identified and how many CpGs are in each. 

### 🔹 Step 2: Generate MRS
This step integrates EWAS results with individual-level methylation data (beta matrix).

- **Required Input**: You need an EWAS summary statistics file with four columns named exactly as follows:
  - Marker: CpG identifier (e.g., cg00000029)
  - BETA: effect size
  - SE: standard error
  - Pvalue: p-value of the association
    
- **Applying the GenMRS Function**
The function allows you to generate a MRS_ComeBack:
- Apply a p-value threshold to select CpGs (Pthred, e.g., 2:minpvalue = 5e-2)
- Optionally apply co-methylation pruning using CoMeRegion
- Calculate the MRS using the selected CpGs and their effect sizes

### 🔹 Step 3: Evaluate Predictive Performance

- Association of MRS with phenotype (e.g., binary outcome)
- Calculation of pseudo R² or AUC
- Visualization of performance across p-value thresholds

### 🔹 Step 4: Compare ComeBack to Methylscore

The pipeline supports head-to-head comparison of two approaches for Methylation Risk Score (MRS) generation:
- `MRS_Methylscore`: standard Methylscore (no pruning) -> (https://github.com/agonse/methylscore)
- `MRS_ComeBack`: pruning by co-methylation

## Output

The pipeline produces:
- `.csv` files with MRS and metadata per sample
- `.RData` files with CMRs and selected CpGs
- `.png` figures showing:
  - Predictive performance (R²)
  - Sensitivity to p-value thresholds
  - CpG count by method

## Contact

For questions or suggestions, please open an issue or contact the author (noguasch@recerca.clinic.cat).


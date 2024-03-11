# CSE284_Final_Project

## Personal Genomics Project: Python Implementation of GWAS Linear Regression
#### Team Members: Andrew Ghafari (A59020215) and Morgan Kelley (A59019598)

# Introduction

Genome-Wide Association Studies (GWAS) analyze the association between genetic variants (usually Single Nucleotide Polymorphisms, or SNPs) across different individuals to discover genetic associations with observable traits. This project aims to replicate the functionality of the Plink command-line tool in Python, specifically focusing on the use of linear regression to assess the association between SNPs and a quantitative trait while controlling for confounding factors such as ancestry. Our Python script processes VCF files to conduct a GWAS using the --linear option, similar to the command 
```
plink --vcf [file] --linear --maf 0.05 --pheno [file] --allow-no-sex --out ps3_gwas
```
to identify common SNPs (Minor Allele Frequency or MAF > 0.05) associated with the trait of interest.

In this project, we directlt compare our results with those obtained from the Plink command-line tool. Our focus is on estimating the association between genetic variants and studied traits using a compressed vcf file as input.

The output of our script will provide statistical analysis including coefficients and p-values, to highlight potential genetic associations with the traits of interest. This endeavor not only aims to replicate the efficiency and accuracy of Plink's native capabilities but also to explore the performance of our Python-based approach versus other programming languages. 


For this analysis, we utilize the 1000 Genomes dataset, covering all 22 chromosomes available from the Public Datahub directory. This project stands as an excellent opportunity for us to dive into optimization techniques within Python, leveraging libraries like pandas for data manipulation and statsmodels for statistical analysis.

This initiative is conducted under Option #1 of the project guidelines for CSE 284 - Personal Genomics, offered by Professor Gymrek.

# Install Instructions

Before running our GWAS analysis script, ensure that your environment has the necessary Python libraries. Our script primarily utilizes pandas for data manipulation, numpy for numerical operations, and statsmodels for conducting linear regression. These can be installed via pip if not already available:

```
pip install pandas numpy statsmodels
```
This is the exhaustive list of libraries required for our model:
```
import argparse
import pandas as pd
import numpy as np
import gzip
from io import StringIO
from tqdm import tqdm
import statsmodels.api as sm
from multiprocessing import Pool
```

# Algorithm

#### Data Preprocessing
* VCF Parsing: The algorithm starts by parsing the VCF file to extract SNP information. This process involves reading the VCF, filtering out header lines, and then converting the relevant SNP data into a pandas DataFrame for easier manipulation.
* Allele Frequency Calculation: For each SNP, we calculate allele counts by parsing genotype data. We then compute the minor allele frequency (MAF) by dividing the count of the minor allele by the total number of alleles observed.
* Filtering: SNPs with a MAF less than 0.05 are filtered out to focus on common variants, reducing the complexity and improving the efficiency of the analysis.

#### Genotype and Phenotype Integration
* Phenotype Data Processing: The phenotype data, provided in a separate file, is read and integrated with the genotype data. This step ensures that each sample's genotype and phenotype information is aligned for the regression analysis.
* Data Transformation: The genotype data is transformed to numeric format, accommodating the regression analysis requirements. This involves converting genotype strings to integers that represent allele counts.

#### Linear Regression
* Regression Model: With the prepared dataset, we perform linear regression for each SNP, using its allele counts as the independent variable and the phenotype value as the dependent variable. The regression analysis is conducted using the statsmodels library, which provides detailed statistics for each SNP, including beta coefficients and p-values.
* Result Filtering and Sorting: The results of the regression analysis are filtered and sorted based on the p-values to identify the SNPs with the most significant associations with the phenotype.

#### Metrics Evaluation

* Comparison with Plink Results: To validate our algorithm's effectiveness, we compare our linear regression results with those obtained from Plink. This involves calculating the mean absolute error (MAE), mean squared error (MSE), and root mean squared error (RMSE) for both p-values and beta coefficients, ensuring our Python implementation's accuracy aligns closely with Plink's established results.

Our Python-based GWAS analysis algorithm is designed to be efficient, accurate, and user-friendly, providing a viable alternative to traditional command-line tools like Plink with the added flexibility and accessibility of Python for genomic research.


# Results

The output of our script includes a detailed summary of the linear regression analysis, listing each SNP analyzed, its beta coefficient (indicating the effect size), and the p-value (indicating the significance of the association). This allows for the identification of genetic variants that may have significant associations with the trait of interest.

In our analysis, we calculate the metrics:

* MAE for p-values: 0.0001235241545070243
* MSE for p-values: 5.163739526631637e-05
* RMSE for p-values: 0.0071859164527787525
* MAE for beta coefficients: 0.11027235180768107
* MSE for beta coefficients: 0.059025047364823356
* RMSE for beta coefficients: 0.24295070974340321

# Datasets

The datasets that we used and produced during our implementation can be found through our Google Drive link:

https://drive.google.com/drive/folders/1xLiMV_VE6ZGjsPoFl-XqzJvUCEiqqqkg?usp=sharing


# Steps to Run the Code. 

To run the GWAS analysis using our Python script:

**1.** Navigate to the Google Drive link (in section Datasets) and download the files **ps3_gwas.vcf, and ps3_gwas.assoc.linear**.

* ps3_gwas.vcf contains the unformatted, unprocessed input vcf data.
* ps3_gwas.assoc.linear contains the output results from Plink that we compare the results to.

**2.** Download and save **(RUN_ME)GWAS_Analysis.sh, FinalProjPlinkPYTHON1.py, and FinalProjPlinkPYTHON2.py** from the files on GitHub.

* (RUN_ME)GWAS_Analysis.sh contains the entire shell script to run both of the python codes
* FinalProjPlinkPYTHON1.py contains the initial code to preprocess the data
* FinalProjPlinkPYTHON2.py contains the code to merge the data and perform linear regression evaluation.

**Note that we had to break our .py code into two sections because of breaking from memory.** If you are on a system with more memory, feel free to combine the two .py codes into one file and adjust the shell script accordingly. If your system breaks after the 



**2. (optional)** 
We already ran the Plink command for you and saved the results in **ps3_gwas.assoc.linear**, which works if you are using the exact same data set we provided, **ps3_gwas.vcf**. If you are using different .vcf data, you **must** run the original Plink command for linear regression to compare with our output results.


**3.** Navigate to the file you downloaded the data to.
The shell script you will be running (found in the Google Drive link) is titled **(RUN_ME)GWAS_Analysis.sh**. 

* In order to run, you must edit the above script to update **"VCF_PATH", "PHENO_PATH","PYTHON_SCRIPT_PATH_PREPROCESSING", "PYTHON_SCRIPT_PATH_1", "PYTHON_SCRIPT_PATH_2", and "PLINK_RESULTS_PATH"**. 

**4.** You must make the shell script executable by navigating to the directory containing the script and running: `chmod +x (RUN_ME)GWAS_Analysis.sh`.
  
**5.** Once the paths are set and the script is executable, you can run the script from the terminal with `./(RUN_ME)GWAS_Analysis.sh`.
  
* Note that this script will use both of the .py files we have provided, so you must save both in the same directory you are running the script in.
* At the end of the first py script, you should get four dataframe outputs **ps3_gwas_output_df1.csv, ps3_gwas_output_df2.csv, ps3_gwas_output_df3.csv, ps3_gwas_output_df4.csv**. The second .py code with then upload these 4 dataframes again (in order to not break on memory). If you have lots of memory, feel free to delete the download portion from **FinalProjPlinkPYTHON1.py** and the upload portion from **FinalProjPlinkPYTHON2.py**.
* If your code breaks during the second .py script but you do have the 4 datasets downloaded already, you should restart your server and continue by running just the following code in your terminal:
  ```
  #!/bin/bash
  VCF_PATH="$HOME/1_284FinalProject/ps3_gwas.vcf.gz"
  PHENO_PATH="$HOME/1_284FinalProject/ps3_gwas.phen"
  OUTPUT_PREFIX="ps3_gwas_output"
  PYTHON_SCRIPT_PATH_1="$HOME/1_284FinalProject/FinalProjPlinkPYTHON1.py"
  PYTHON_SCRIPT_PATH_2="$HOME/1_284FinalProject/FinalProjPlinkPYTHON2.py"
  PLINK_RESULTS_PATH="$HOME/1_284FinalProject/ps3_gwas.assoc.linear"
  echo "Data preprocessing complete. Intermediate results saved."
  echo "Starting processing and analysis with intermediate results..."
  python $PYTHON_SCRIPT_PATH_2 --out $OUTPUT_PREFIX --plink_results $PLINK_RESULTS_PATH --pheno $PHENO_PATH
  echo "GWAS analysis complete."
  ```
* This code can take up to an hour to run so please be patient and please let us know if you face any difficulties, we are happy to assist. 


# Contributors

This project was developed by Andrew Ghafari and Morgan Kelley, under the guidance of Professor Melissa Gymrek for the CSE 284 - Personal Genomics course at UC San Diego. Our goal was not only to replicate the functionality of existing GWAS tools but also to explore the potential for Python in genomics research.


# Future Work:

Future directions for this project include extending the analysis to other programming languages, such as C++ or Java, for performance comparison. Additionally, we aim to incorporate more complex statistical models and machine learning algorithms to enhance the analysis of genetic data.

# Acknowledgements

We extend our gratitude to Professor Melissa Gymrek and the teaching staff of CSE 284 for their support and guidance throughout this project. We also thank the contributors of the 1000 Genomes Project for providing a rich dataset for our analysis.

# References
* Purcell, S., et al. (2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. American Journal of Human Genetics, 81(3), 559-575.
* Python Software Foundation. Python Language Reference, version 3.x. https://python.org



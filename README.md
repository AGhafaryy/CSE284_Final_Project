# CSE284_Final_Project

## Personal Genomics Project: Python Implementation of GWAS Linear Regression
#### Team Members: Andrew Ghafari (A59020215) and Morgan Kelley (A59019598)

# Introduction

Genome-Wide Association Studies (GWAS) analyze the association between genetic variants (usually Single Nucleotide Polymorphisms, or SNPs) across different individuals to discover genetic associations with observable traits. This project aims to replicate the functionality of the Plink command-line tool in Python, specifically focusing on the use of linear regression to assess the association between SNPs and a quantitative trait while controlling for confounding factors such as ancestry. Our Python script processes VCF files to conduct a GWAS using the --linear option, similar to the command plink --vcf [file] --linear --maf 0.05 --pheno [file] --allow-no-sex --out ps3_gwas, to identify common SNPs (Minor Allele Frequency or MAF > 0.05) associated with the trait of interest.

In this project, we directlt compare our results with those obtained from the Plink command-line tool. Our focus is on estimating the association between genetic variants and studied traits using a compressed vcf file as input.

The output of our script will provide statistical analysis including coefficients and p-values, to highlight potential genetic associations with the traits of interest. This endeavor not only aims to replicate the efficiency and accuracy of Plink's native capabilities but also to explore the performance of our Python-based approach versus other programming languages. 

Future work includes replicating these scripts in other languages such as C++ or Java and comparing the results in terms of execution speed.

For this analysis, we utilize the 1000 Genomes dataset, covering all 22 chromosomes available from the Public Datahub directory. This project stands as an excellent opportunity for us to dive into optimization techniques within Python, leveraging libraries like pandas for data manipulation and statsmodels for statistical analysis.

This initiative is conducted under Option #1 of the project guidelines for CSE 284 - Personal Genomics, offered by Professor Gymrek.

# Install Instructions

Before running our GWAS analysis script, ensure that your environment has the necessary Python libraries. Our script primarily utilizes pandas for data manipulation, numpy for numerical operations, and statsmodels for conducting linear regression. These can be installed via pip if not already available:

```
pip install pandas numpy statsmodels
```

# Basic Usage

To run the GWAS analysis using our Python script, navigate to the directory containing the script and execute it with the appropriate command-line arguments. The basic command structure is as follows:

```
python gwas_analysis.py --vcf input_file.vcf.gz --pheno phenotype.phen --out output_prefix
```
    --vcf: Path to the VCF file containing SNP data.
    --pheno: Path to the file containing phenotype data.
    --out: Prefix for the output files.
    
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

# Datasets

The datasets that we used and produced during our implementation can be found through our Google Drive link:

https://drive.google.com/drive/folders/1qRffK7Sr2Z7onwI9QyoyMFd-UoVbxQjF?usp=drive_link

We describe in our .ipynb file when each dataset was produced. Many of these were created as checkpoints to ease the implementation.



# Steps to run the code. 

Running the code is pretty straightforward, the first part of it is running the  plink --vcf --linear --maf 0.05 --pheno --allow-no-sex command (the one we used for PA3), which will generate a .linear file. This should be pretty easy, since we have already done this in class before. 

For the python file, everything you need is present inside the notebook, you just need to run everything sequentially and make sure you have the proper libraries installed before. 

(You can also check out the pdf for an already compiled version)

Please let us know if you face any difficulties, we are happy to assist. 

# Results

The output of our script includes a detailed summary of the linear regression analysis, listing each SNP analyzed, its beta coefficient (indicating the effect size), and the p-value (indicating the significance of the association). This allows for the identification of genetic variants that may have significant associations with the trait of interest.

In our analysis, we calculate the metrics:

* MAE for p-values: 0.0001235241545070243
* MSE for p-values: 5.163739526631637e-05
* RMSE for p-values: 0.0071859164527787525
* MAE for beta coefficients: 0.11027235180768107
* MSE for beta coefficients: 0.059025047364823356
* RMSE for beta coefficients: 0.24295070974340321

# Contributors

This project was developed by Andrew Ghafari and Morgan Kelley, under the guidance of Professor Melissa Gymrek for the CSE 284 - Personal Genomics course at UC San Diego. Our goal was not only to replicate the functionality of existing GWAS tools but also to explore the potential for Python in genomics research.


# Future Work:

Future directions for this project include extending the analysis to other programming languages, such as C++ or Java, for performance comparison. Additionally, we aim to incorporate more complex statistical models and machine learning algorithms to enhance the analysis of genetic data.

# Acknowledgements

We extend our gratitude to Professor Melissa Gymrek and the teaching staff of CSE 284 for their support and guidance throughout this project. We also thank the contributors of the 1000 Genomes Project for providing a rich dataset for our analysis.

# References
* Purcell, S., et al. (2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. American Journal of Human Genetics, 81(3), 559-575.
* Python Software Foundation. Python Language Reference, version 3.x. https://python.org



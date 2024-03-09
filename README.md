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
    

# Datasets
ps3_gwas.vcf: https://drive.google.com/file/d/17Lr1YpH88fVytZj7FYFau4S0zLo9vxLJ/view?usp=drive_link

vcfdfout.csv: https://drive.google.com/file/d/1_Ll7ZcPNY5UCpEOb4ndbxSE9VqzXrURT/view?usp=drive_link

ps3_gwas.assoc.linear: https://drive.google.com/file/d/1XQZSAZjCgWIm8pi_rOTNqGdh2Uc1vmhy/view?usp=drive_link

numeric_df_first_quarter.csv:
https://drive.google.com/file/d/1w7eOc5qI6_sxKxn2lzYfzbv8vsISvjVN/view?usp=drive_link

numeric_df_second_quarter.csv
https://drive.google.com/file/d/1LQ8RhJOD_BC7dzqEuZVurPZZJyg0ygBS/view?usp=drive_link

numeric_second_half_part1_df.csv:
https://drive.google.com/file/d/1KRy5SDyvoFGuqpsHCw-_sj9ydPIH39zu/view?usp=drive_link

numeric_second_half_part2_df.csv:
https://drive.google.com/file/d/1BwNeRY9lwiZ6J2Ia67dArEEeG2Jupz-0/view?usp=drive_link

filtered_df.csv:
https://drive.google.com/file/d/1Jg8NsQVbtrJPLgMFz3cxR1y_1J0vGV0F/view?usp=drive_link

concat_w_phen.csv:
https://drive.google.com/file/d/1leGlbsmFMIfjceoUgjO5rnapMFuiDi_R/view?usp=drive_link

linRegResults.csv:
https://drive.google.com/file/d/1BBqy3v9A2ROZh351Lx-gMj-9lFStlGgK/view?usp=drive_link

concatenated_df_final.csv:
https://drive.google.com/file/d/1qZrUtbmKuY-p5wX59QP54sgMpuYBtRNP/view?usp=drive_link

ps3_gwas.phen:
https://drive.google.com/file/d/1kWi1_M5T-0f84XX_DNvAlGVEifTWSrwk/view?usp=drive_link





# Steps to run the code. 

Running the code is pretty straightforward, the first part of it is running the  plink --vcf --linear --maf 0.05 --pheno --allow-no-sex command (the one we used for PA3), which will generate a .linear file. This should be pretty easy, since we have already done this in class before. 

For the python file, everything you need is present inside the notebook, you just need to run everything sequentially and make sure you have the proper libraries installed before. 

(You can also check out the pdf for an already compiled version)

Please let us know if you face any difficulties, we are happy to assist. 

# Results

The output of our script includes a detailed summary of the linear regression analysis, listing each SNP analyzed, its beta coefficient (indicating the effect size), and the p-value (indicating the significance of the association). This allows for the identification of genetic variants that may have significant associations with the trait of interest.

In our analysis, we calculate the metrics:

MAE for p-values: 0.0001235241545070243
MSE for p-values: 5.163739526631637e-05
RMSE for p-values: 0.0071859164527787525
MAE for beta coefficients: 0.11027235180768107
MSE for beta coefficients: 0.059025047364823356
RMSE for beta coefficients: 0.24295070974340321

# Contributors

This project was developed by Andrew Ghafari and Morgan Kelley, under the guidance of Professor Melissa Gymrek for the CSE 284 - Personal Genomics course at UC San Diego. Our goal was not only to replicate the functionality of existing GWAS tools but also to explore the potential for Python in genomics research.


# Future Work:

Future directions for this project include extending the analysis to other programming languages, such as C++ or Java, for performance comparison. Additionally, we aim to incorporate more complex statistical models and machine learning algorithms to enhance the analysis of genetic data.

# Acknowledgements

We extend our gratitude to Professor Melissa Gymrek and the teaching staff of CSE 284 for their support and guidance throughout this project. We also thank the contributors of the 1000 Genomes Project for providing a rich dataset for our analysis.

# References

    Purcell, S., et al. (2007). PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. American Journal of Human Genetics, 81(3), 559-575.
    Python Software Foundation. Python Language Reference, version 3.x. https://python.org



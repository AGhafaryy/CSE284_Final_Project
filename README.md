# CSE284_Final_Project

## Personal Genomics Project: Python Implementation of Linear Regression for SNP Analysis
#### Team Members: Andrew Ghafari (A59020215) and Morgan Kelley (A59019598)

In this project, we aim to create a Python script that implements linear regression with multiple covariates to analyze SNP (Single Nucleotide Polymorphisms) data, directly comparing our results with those obtained from the Plink command-line tool. Our focus is on estimating the association between genetic variants and studied traits using a compressed Variant Call Format (VCF) file as input.

The output of our script will provide statistical analysis including coefficients and p-values, to highlight potential genetic associations with the traits of interest. This endeavor not only aims to replicate the efficiency and accuracy of Plink's native capabilities but also to explore the performance of our Python-based approach versus other programming languages. 
Future work includes replicating these scripts in other languages such as C++ or Java and comparing the results in terms of execution speed.

For this analysis, we utilize the 1000 Genomes dataset, covering all 22 chromosomes available from the Public Datahub directory. This project stands as an excellent opportunity for us to dive into optimization techniques within Python, leveraging libraries like pandas for data manipulation and statsmodels for statistical analysis.

This initiative is conducted under Option #1 of the project guidelines for CSE 284 - Personal Genomics, offered by Professor Gymrek.


# Steps to run the code. 

Running the code is pretty straightforward, the first part of it is running the  plink --vcf --linear --maf 0.05 --pheno --allow-no-sex command (the one we used for PA3), which will generate a .linear file. This should be pretty easy, since we have already done this in class before. 

For the python file, everything you need is present inside the notebook, you just need to run everything sequentially and make sure you have the proper libraries installed before. 

Please let us know if you face any difficulties, we are happy to assist. 


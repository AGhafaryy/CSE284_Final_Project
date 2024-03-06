# CSE284_Final_Project

## Personal Genomics Project: Python Implementation of Linear Regression for SNP Analysis
#### Team Members: Andrew Ghafari (A59020215) and Morgan Kelley (A59019598)

In this project, we aim to create a Python script that implements linear regression with multiple covariates to analyze SNP (Single Nucleotide Polymorphisms) data, directly comparing our results with those obtained from the Plink command-line tool. Our focus is on estimating the association between genetic variants and studied traits using a compressed Variant Call Format (VCF) file as input.

The output of our script will provide statistical analysis including coefficients and p-values, to highlight potential genetic associations with the traits of interest. This endeavor not only aims to replicate the efficiency and accuracy of Plink's native capabilities but also to explore the performance of our Python-based approach versus other programming languages. 
Future work includes replicating these scripts in other languages such as C++ or Java and comparing the results in terms of execution speed.

For this analysis, we utilize the 1000 Genomes dataset, covering all 22 chromosomes available from the Public Datahub directory. This project stands as an excellent opportunity for us to dive into optimization techniques within Python, leveraging libraries like pandas for data manipulation and statsmodels for statistical analysis.

This initiative is conducted under Option #1 of the project guidelines for CSE 284 - Personal Genomics, offered by Professor Gymrek.

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

# Future Work:

We are still working on the repo, and will be cleaning the code and adding more analysis before we deliver the final project. Thank you for your patience.





#!/bin/bash


VCF_PATH="$HOME/1_284FinalProject/ps3_gwas.vcf.gz"
PHENO_PATH="$HOME/1_284FinalProject/ps3_gwas.phen"
OUTPUT_PREFIX="ps3_gwas_output"
PYTHON_SCRIPT_PATH_1="$HOME/1_284FinalProject/FinalProjPlinkPYTHON1.py"
PYTHON_SCRIPT_PATH_2="$HOME/1_284FinalProject/FinalProjPlinkPYTHON2.py"
PLINK_RESULTS_PATH="$HOME/1_284FinalProject/ps3_gwas.assoc.linear"



echo "Starting data preprocessing and saving intermediate results..."
python $PYTHON_SCRIPT_1 --vcf $VCF_PATH --pheno $PHENO_PATH --out $OUTPUT_PREFIX

echo "Data preprocessing complete. Intermediate results saved."

echo "Starting processing and analysis with intermediate results..."
python $PYTHON_SCRIPT_PATH_2 --out $OUTPUT_PREFIX --plink_results $PLINK_RESULTS_PATH --pheno $PHENO_PATH


echo "GWAS analysis complete."


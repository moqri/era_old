#!/bin/sh

# Generate EWAS data table
echo "Collecting DNA methylation data..."
python ewas_prep.py

# Run EWAS analysis using GLINT
echo "Running EWAS..."
cd GLINT_1
python glint.py --datafile data.txt --covarfile covariates.txt --phenofile phenotypes.txt --gsave
jwakim$ python glint.py --datafile data.glint --ewas --pheno treatment --covar --minstd 0.02
python glint.py --plot --manhattan --results results.glint.linreg.txt


import pandas as pd

# Download data from GEO
rep = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142439/matrix//GSE142439_series_matrix.txt.gz'

# Extract DNA methylation data (starting from line 71, excluding last line)
df = pd.read_csv(rep, sep='\t', skiprows=71, index_col=0)
df = df.drop('!series_matrix_table_end')
df.to_csv("GLINT_1/data.txt")

# Generate the covariate file
meta=pd.read_csv(rep,sep='\t',skiprows=40,index_col=0,nrows=3) # meta data starts from line 40
meta.columns=df.columns
meta.index=['age','cell','treatment']
meta=meta.T
meta['age']=meta.age.str.strip('age: ')
meta['cell']=meta.cell.str.strip('cell type: ')
meta['treatment']=meta.treatment.str.strip('treatmentstatus: ')

# Encode covariates
coding = {"cell": {"skin fibroblasts" : 0, "vein endothelial cells": 1},
        "treatment": {"Normal": 0, "Treated": 1}}
meta = meta.replace(coding)
covars = meta.drop(columns=['treatment'])
covars.to_csv("GLINT_1/covariates.txt")
pheno = meta.drop(columns=['age', 'cell'])
pheno.to_csv("GLINT_1/phenotypes.txt")

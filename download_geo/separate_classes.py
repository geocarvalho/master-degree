import pandas as pd
import numpy as np
import os

pheno = "GSE51032/GSE51032_all_phenotype.csv"
pheno_df = pd.read_csv(pheno)
direc = os.path.dirname(os.path.abspath(pheno))

pheno_df["cancer type (icd-10):ch1"] = pheno_df["cancer type (icd-10):ch1"].replace(
    np.nan, "normal", regex=True)
# pheno_df["cancer type (icd-10):ch1"].value_counts()
pheno_df.to_csv(
    os.path.join(direc, "GSE51032_all_phenotype_normals.csv"), index=False)

pheno_df[pheno_df["cancer type (icd-10):ch1"]=="normal"]["Unnamed: 0"].to_csv(
    os.path.join(direc, "GSE51032_all_phenotype_onlynormals.csv"), index=False, header=False)
# breast cancer
pheno_df[pheno_df["cancer type (icd-10):ch1"]=="C50"]["Unnamed: 0"].to_csv(
    os.path.join(direc, "GSE51032_all_phenotype_onlymama.csv"), index=False, header=False)
# colon cancer
# pheno_df[pheno_df["cancer type (icd-10):ch1"]=="C18"]["Unnamed: 0"].to_csv(
#     os.path.join(direc, "GSE51032_all_phenotype_onlycolon.csv"), index=False, header=False)
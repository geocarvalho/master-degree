import numpy as np
import pandas as pd

def threshold_bvalue(value):
    if float(value) <= 0.2:
        return "low"
    elif float(value) >=0.8:
        return "high"
    else:
        return "middle"

def prepare_df(gse_df, pheno):
    """ Merge GSE and phenotype files into a DataFrame
    """
    # GSE
    gse_df["geo_accession"] = gse_df["samples"].apply(lambda lst: lst.split("_")[0])

    # Transform b-values in methylation rate groups
    cols = [col for col in gse_df.columns if col.startswith("cg")]
    gse_df[cols] = gse_df[cols].applymap(threshold_bvalue)
    
    # phenotype
    pheno_df = pd.read_csv(pheno)
    pheno_df = pheno_df.rename({
        "Unnamed: 0": "samples", "subtype_ihc:ch1": "subtype"
        }, axis=1)  

    l_pheno_df = pheno_df[["geo_accession", "subtype"]]
    gse_pheno_df = gse_df.merge(l_pheno_df, on="geo_accession")
    gse_pheno_df.drop(["geo_accession", "samples"], axis=1, inplace=True)

    return gse_pheno_df

# Open datasets
gse1 = "GSE72245/GSE72245_filtered_bvalues.csv"
pheno1 = "GSE72245/GSE72245_all_phenotype.csv"

gse2 = "GSE72251/GSE72251_filtered_bvalues.csv"
pheno2 = "GSE72251/GSE72251_all_phenotype.csv"
# gse3 = "GSE72254/GSE72254_filtered_bvalues.csv"

gse1_df = pd.read_csv(gse1, index_col=0, header=None).transpose()
gse2_df = pd.read_csv(gse2, index_col=0, header=None).transpose()
# gse3_df = pd.read_csv(gse3)

gse1_df.rename(columns={np.NaN: "samples"}, inplace=True)
gse2_df.rename(columns={np.NaN: "samples"}, inplace=True)

comm_cpgs = ["samples"] + list(set(list(gse1_df_trans.columns)).intersection(list(gse2_df_trans.columns)))

gse1_df_filtered = gse1_df[comm_cpgs]
gse_pheno_df1 = prepare_df(gse1_df_filtered, pheno1)

gse2_df_filtered = gse2_df[comm_cpgs]
gse_pheno_df2 = prepare_df(gse2_df_filtered, pheno2)

# Concatenate DFs
gse_pheno_df = pd.concat([gse_pheno_df1, gse_pheno_df2], sort=False).reset_index(drop=True) 

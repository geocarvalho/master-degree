import pandas as pd
import numpy as np
import math

def calculate_classes_res(series):
    """ Calculate the number of classes and the threshold"""
    n = len(series.unique())
    ## Check number of classes 
    for c in range(1, 100):
        if 2**c >= n:
            break
    # Calculate the threshold
    min = series.min()
    max = series.max()
    minus = max - min
    ## Always round up the result
    res = math.ceil(minus / c)
    # res = round(minus/c, 2)
    print("classes, threshold: ", c, res)
    return c, res

def create_classes(value, inter, classes, thres):
    """ Create classes based on the minimun value, number of classes and inteval
    between classes """
    if value < inter:
        return "<%s" % inter
    for c in range(1,classes+1):
        if inter <= value < inter+thres:
            return ">=%s,<%s" % (str(inter), str(inter+thres))
        inter+=thres

def col_classes(df, col):
    """ Create a new column with the classes for the given column """
    print("Column: ", col)
    new_col = col + "_classes"
    ## Number of data points
    classes, res = calculate_classes_res(df[col])
    ## Create new column with the classe
    df[new_col] = df.apply(lambda row: create_classes(row[col], df[col].min(), classes, res), axis=1)
    return df.copy()

dataset = "../data/GSE51032_all_phenotype.csv"
pheno_df = pd.read_csv(dataset)
cell_count_data = "GSE51032/GSE51032_cell_estimation.csv"
cell_count_df = pd.read_csv(cell_count_data)
epismoker = "GSE51032/epismoker_SSt_GSE109381.csv"
epismoker_df = pd.read_csv(epismoker)
age = "GSE51032/age_prediction_GSE109381.csv"
age_df = pd.read_csv(age)


# Replace negative values in cell count by 0
cell_count_df["CD4T"] = np.where(cell_count_df["CD4T"] < 0.0, 0.0, cell_count_df["CD4T"])

# Create sample_id
cell_count_df = cell_count_df.rename(columns={"Unnamed: 0": "sample_id"})
pheno_df["sample_id"] = pheno_df["geo_accession"] + "_" + pheno_df["title"]

# Select just samples from breast cancer and normal
pheno_df["cancer type (icd-10):ch1"] = pheno_df["cancer type (icd-10):ch1"].replace(
    np.nan, "normal", regex=True)
df_breast = pheno_df.loc[pheno_df["cancer type (icd-10):ch1"].isin(["C50", "normal"])].copy()
df_breast = df_breast[df_breast["gender:ch1"]=="F"]

# Merge df_breast and cell_count_df
result = pd.merge(df_breast, cell_count_df, on="sample_id", how="left")
result = pd.merge(result, epismoker_df, left_on="sample_id", right_on="SampleName", how="left")
result = pd.merge(result, age_df, left_on="sample_id", right_on="SampleID", how="left")

# Rename columns
result = result.rename(columns={
    "age:ch1": "age", "time to diagnosis:ch1": "time_to_diagnosis",
    "age at menarche:ch1": "age_at_menarche", "cancer type (icd-10):ch1": "cancer_type"})

# Change cell count to percentage
for cell in ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"]:
  result[cell] = result[cell] * 100

# Create classes for the columns
cols = ["age", "time_to_diagnosis", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "mAge_Hovath", "mAge_Hannum", "PhenoAge"]
for col in cols:
    result = col_classes(result, col)

# Replace nan as a class
result["time_to_diagnosis_classes"].fillna(-1, inplace=True)
result["age_at_menarche"].fillna(-1, inplace=True)
result["age_at_menarche_int"] = result["age_at_menarche"].astype(int)

# Select important columns
output = result[[
    "sample_id", "age_at_menarche_int", "time_to_diagnosis", "time_to_diagnosis_classes", "age", "age_classes",
    "CD8T", "CD8T_classes", "CD4T", "CD4T_classes", "NK", "NK_classes", "Bcell", "Bcell_classes", "Mono", "Mono_classes",
    "Gran", "Gran_classes", "cancer_type", "PredictedSmokingStatus", "mAge_Hovath_classes", "mAge_Hannum_classes",
    "PhenoAge_classes"]]
output.to_csv("../download_geo/GSE51032/GSE51032_classes_design.csv", index=False)
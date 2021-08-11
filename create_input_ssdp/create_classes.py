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
    print("Number of classes: ", c)
    # Calculate the threshold
    min = series.min()
    max = series.max()
    minus = max - min
    ## Always round up the result
    res = math.ceil(minus / c)
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

dataset = "../data/GSE51032_all_phenotype.csv"
df = pd.read_csv(dataset)
# Select just samples from breast cancer and normal
df["cancer type (icd-10):ch1"] = df["cancer type (icd-10):ch1"].replace(
    np.nan, "normal", regex=True)
df_breast = df.loc[df["cancer type (icd-10):ch1"].isin(["C50", "normal"])]
df_breast = df_breast[df_breast["gender:ch1"]=="F"]
print("Just female, normals and C50 (breast cancer)")
print(df_breast["gender:ch1"].value_counts())
print(df_breast["cancer type (icd-10):ch1"].value_counts())

# Create classes for age
## Number of data points
classes, res = calculate_classes_res(df["age:ch1"])
print("Classes n e threshold: ", classes, res)
## Create new column with the classe
df_breast["age:ch1"] = df_breast["age:ch1"].round(2)
df_breast["age"] = df_breast.apply(lambda row: create_classes(row["age:ch1"], df_breast["age:ch1"].min(), classes, res), axis=1)
print("\nAge:")
print(df_breast["age"].value_counts())

# Create classes for time to diagnosis
## Number of data points
classes, res = calculate_classes_res(df["time to diagnosis:ch1"])
## Create new column with the classe
df_breast["time to diagnosis:ch1"] = df_breast["time to diagnosis:ch1"].round(2)
df_breast["time_to_diagnosis"] = df_breast.apply(lambda row: create_classes(row["time to diagnosis:ch1"], df_breast["time to diagnosis:ch1"].min(), classes, res), axis=1)
# Replace nan as a class
df_breast["time_to_diagnosis"].fillna(-1, inplace=True)
print("\nTime of diagnosis:")
print(df_breast["time_to_diagnosis"].value_counts())
df_breast["age at menarche:ch1"].fillna(-1, inplace=True)
df_breast["age at menarche:ch1"] = df_breast["age at menarche:ch1"].astype(int)
print("\nAge of menarche:")
print(df_breast["age at menarche:ch1"].value_counts())

# Select important columns
output = df_breast[["title", "geo_accession", "age at menarche:ch1", "time to diagnosis:ch1", "time_to_diagnosis", "age:ch1", "age", "cancer type (icd-10):ch1"]]
output = output.rename(columns={"age at menarche:ch1": "age_at_menarche",
    "cancer type (icd-10):ch1": "cancer_type"})
output.to_csv("../download_geo/GSE51032/GSE51032_classes_design.csv", index=False)
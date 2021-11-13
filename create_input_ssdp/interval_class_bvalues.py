from sklearn import preprocessing
import concurrent.futures
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
    for c in range(1,classes+2):
        if inter <= value < inter+thres:
            return ">=%s,<%s" % (str(inter), str(inter+thres))
        inter+=thres

def col_classes(series):
    """ Create a new column with the classes for the given column """
    minimal = series.min()
    ## Number of data points
    classes, res = calculate_classes_res(series)
    ## Create new column with the classe
    new_series = series.apply(lambda row: create_classes(row, minimal, classes, res))
    return new_series

def factorize_df(df):
    """ Factorize column classes into integers """
    new_df = pd.DataFrame()
    le = preprocessing.LabelEncoder()
    cols = [col for col in df.columns if col != "target"]
    for col in cols:
        print(col)
        le.fit(df[col])
        new_df[col] = le.transform(df[col])
    new_df.index = df.index
    if "target" in df.columns:
        new_df["target"] = df["target"].copy()
    return new_df

data = "../download_geo/GSE51032/GSE51032_bvalues_filtered_cn.csv"
pheno = "../download_geo/GSE51032/GSE51032_classes_design.csv"
output = data.replace(".csv", "_pheno.csv")
pheno_df = pd.read_csv(pheno)
pheno_df_filtered = pheno_df[[
    "sample_id", "age_at_menarche_int", "age_classes", "CD8T_classes",
    "CD4T_classes", "NK_classes", "Bcell_classes", "Mono_classes", "Gran_classes", "PredictedSmokingStatus",
    "mAge_Hannum_classes", "mAge_Hovath_classes", "PhenoAge_classes", "cancer_type"]] # time_to_diagnosis_classes
pheno_df_filtered.set_index("sample_id", inplace=True)
pheno_df_filtered = pheno_df_filtered.rename(columns={"cancer_type": "target"})
pheno_df_filtered["target"] = pheno_df_filtered["target"].map({"normal": "0", "C50": "1"}) 
df = pd.read_csv(data)
df.set_index("probes", inplace=True)
df = df.astype(float)
# df = df.T

cols = df.columns.tolist()
series_lst = []
with concurrent.futures.ProcessPoolExecutor() as executer:
    results = [executer.submit(col_classes, df[col]*100) for col in cols]
    for f in concurrent.futures.as_completed(results):
        series_lst.append(f.result())
print("Process' results done.")
class_df = pd.concat(series_lst, axis=1)
# Check for empty cells
for col in class_df.columns:
    if class_df[col].isnull().any():
        print("Empty column: ", col)
# Factorize columns https://stackoverflow.com/questions/47595268/convert-classes-to-numeric-in-a-pandas-dataframe
new_class_df = factorize_df(class_df)
new_pheno_df_filtered = factorize_df(pheno_df_filtered)
# Transpose
class_df = class_df.T
new_class_df = new_class_df.T
# merge with phenotipes
merge = pd.concat([class_df, pheno_df_filtered], axis=1)
new_merge = pd.concat([new_class_df, new_pheno_df_filtered], axis=1)
merge.to_csv(output, index=False)
new_merge.to_csv(output.replace(".csv", "_factor.csv"), index=False)

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
    # res = math.ceil(minus / c)
    res = round(minus/c, 2)
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

def col_classes(series):
    """ Create a new column with the classes for the given column """
    minimal = series.min()
    ## Number of data points
    classes, res = calculate_classes_res(series)
    ## Create new column with the classe
    new_series = series.apply(lambda row: create_classes(row, minimal, classes, res), axis=1)
    return new_series

data = "./GSE51032/GSE51032_bvalues_filtered_cn.csv"
output = data.replace(".csv", "_pheno.csv")
df = pd.read_csv(data)
df.set_index("probes", inplace=True)
df = df.astype(float)

cols = df.columns.tolist()
series_lst = []
with concurrent.futures.ProcessPoolExecutor() as executer:
    results = [executer.submit(col_classes, df[col]) for col in cols]
    for f in concurrent.futures.as_completed(results):
        series_lst.append(f.result())

new_df = pd.concat(series_lst, axis=1)
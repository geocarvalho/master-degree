import numpy as np
import pandas as pd

# see example in website
# http://htmlpreview.github.io/?https://github.com/sailalithabollepalli/EpiSmokEr/blob/master/vignettes/epismoker.html
# open pheno data
ss = pd.read_csv("GSE51032_all_phenotype.csv")
# Just select samples from breast cancer and normal
ss["cancer type (icd-10):ch1"] = ss["cancer type (icd-10):ch1"].replace(
    np.nan, "normal", regex=True)
ss_filter = ss.loc[ss["cancer type (icd-10):ch1"].isin(["C50", "normal"])]
# select only necessary columns
ss_filter = ss_filter[["geo_accession", "gender:ch1", "age:ch1", "title"]].copy()
# create empty columns as in example
ss_filter[""]= ss_filter["geo_accession"] + "_" + ss_filter["title"]
# create age column with int
ss_filter["Age"] = ss_filter["age:ch1"].astype(int)
# create gender and sex columns
ss_filter["Gender"] = ss_filter["gender:ch1"].replace({"M": "m", "F": "f"})
ss_filter["sex"] = ss_filter["gender:ch1"].replace({"M": 1, "F": 2})
# create slide and array columns
ss_filter["Slide"] = ss_filter.apply(lambda sample: sample["title"].split("_")[0], axis=1)
ss_filter["Array"] = ss_filter.apply(lambda sample: sample["title"].split("_")[1], axis=1)
# select only the columns created and create output csv file
output = ss_filter[["", "Age", "Gender", "sex", "Slide", "Array"]]
output.to_csv("idat/normal_mama/samplesheet_GSE109381_simple.csv", index=False)

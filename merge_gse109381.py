import numpy as np
import pandas as pd
from dask import dataframe as dd

def normalize_values(df):
    "Classify the values in low, middle and high methylated"
    df[df <= 0.2] = 0
    df[df >= 0.8] = 2
    df[(df > 0.2) & (df < 0.8)] = 1
    return df


# csv1 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part1.csv"
# csv2 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part2.csv"
# csv3 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part3.csv"
# csv4 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part4.csv"
# pheno = "/home/watson/george/master-degree/GSE109381/GSE109381_all_phenotype_filtered_tclasses.csv"
# to_transpose = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_formatted.csv"
# Create transposed file
# csv_cols_list = {}
# for csv in [csv1, csv2, csv3, csv4]:
#     mylst = []
#     for chunk in pd.read_csv(csv, index_col=0, header=None, low_memory=False, chunksize=100000):
#         mylst.append(chunk)
#     gse_df = pd.concat(mylst, axis=0)
#     del mylst
#     gse_df = gse_df.transpose()
#     gse_df.rename(columns={np.NaN: "samples"}, inplace=True)
#     csv_cols_list[csv] = gse_df.columns
#     csv_output_name = csv.replace(".csv", "_t.csv")
#     gse_df.to_csv(csv_output_name, index=False)

# Select common columns
# print("Opening datasets")
# gse1_df = dd.read_csv(csv1, blocksize=64000000)
# gse1_df = gse1_df.set_index("probes")
# gse1_df = normalize_values(gse1_df)
# 
# gse2_df = dd.read_csv(csv2, blocksize=64000000)
# gse2_df = gse2_df.set_index("probes")
# gse2_df = normalize_values(gse2_df)
# 
# gse3_df = dd.read_csv(csv3, blocksize=64000000)
# gse3_df = gse3_df.set_index("probes")
# gse3_df = normalize_values(gse3_df)
# 
# gse4_df = dd.read_csv(csv4, blocksize=64000000)
# gse4_df = gse4_df.set_index("probes")
# gse4_df = normalize_values(gse4_df)
# 
# gse_merge = gse1_df.merge(gse2_df, left_index=True, right_index=True)
# gse_merge = gse_merge.merge(gse3_df, left_index=True, right_index=True)
# gse_merge = gse_merge.merge(gse4_df, left_index=True, right_index=True)
# 
# gse_cols = [col.split("_")[0] for col in gse_merge.columns]
# gse_merge.columns = gse_cols
# 
# pheno_df = dd.read_csv(pheno)
# pheno_df["sample_name"] = "classes"
# pheno_df = pheno_df.set_index("sample_name")
# gse_out = dd.concat([gse_merge, pheno_df])
# gse_out.to_csv(to_transpose, single_file=True)
# print(gse_out.columns)

# Create a merged csv
# from glob import glob
# path = /home/watson/george/master-degree/GSE109381/GSE109381_bvalues_formatted.csv 
# filenames = glob(path) 
# with open(to_transpose, 'w') as out: 
#   for fn in filenames: 
#       with open(fn) as f: 
#           out.write(f.read())

# Create transposed file
# mylst = []
# for chunk in pd.read_csv(to_transpose, index_col=0, header=None, low_memory=False, chunksize=100000):
#     mylst.append(chunk)
#     gse_df = pd.concat(mylst, axis=0)
# del mylst
# gse_df = gse_df.transpose()
# csv_output_name = to_transpose.replace(".csv", "_t.csv")
# gse_df.to_csv(csv_output_name, index=False)

# Work with the most common class
gse = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_gbm.csv"
csv = "/home/watson/george/master-degree/GSE109381/GBM_class_sample_2cols.csv"
# csv
csv_df = pd.read_csv(csv, names=["samples", "class"])
csv_df["class"] = pd.Categorical(csv_df["class"])
csv_df["class_num"] = csv_df["class"].cat.codes
csv_df.to_csv("/home/watson/george/master-degree/GSE109381/GBM_class_sample_categories.csv", index=False)
classes = csv_df[["samples", "class_num"]]

# gse
gse_df = pd.read_csv(gse, index_col=0, header=None).transpose()
gse_df.rename(columns={np.NaN: "samples"}, inplace=True)
gse_df["samples"] = [sample.split("_")[0] for sample in gse_df["samples"]]
# the majority class is "GBM, RTK II" that now is 5

# merge
merged_df = gse_df.merge(classes, on="sample")
merged_df = merged_df.set_index("samples")
cols = gse_df.columns.to_list()
gse_df[cols] = gse_df[cols].astype(np.float32)
gse_df = normalize_values(gse_df)

mylst = []
for chunk in pd.read_csv(csv, low_memory=False, chunksize=100000):
    cols = [col for col in chunk.columns if col != "class_num"]
    chunk[cols] = chunk[cols].astype(float)
    chunk[chunk[cols] <= 0.2] = 0
    chunk[chunk[cols] >= 0.8] = 2
    chunk[(chunk[cols] > 0.2) & (chunk[cols] < 0.8)] = 1
    mylst.append(chunk)
gse_df = pd.concat(mylst, axis=0)
gse_df.to_csv("GSE109381_bvalues_gbm_merged_classes_formatted.csv", index=False)

# Remove float points to have just integers
# cat GSE109381_bvalues_gbm_merged_classes_formatted.csv | sed 's/\.0//g' > GSE109381_bvalues_gbm_merged_classes_formatted_nopoints.csv

# Rodar o SSDP+
# java -jar -Xmx50G -XX:+UseConcMarkSweepGC ssdp_plus/out/artifacts/ssdp_plus_jar/ssdp_plus.jar master-degree/GSE109381_bvalues_gbm_merged_classes_formatted_nopoints.csv 5 0.9 5 | tee -a  master-degree/GSE109381_log.txt
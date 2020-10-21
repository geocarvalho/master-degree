import numpy as np
import pandas as pd

csv1 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part1.csv"
csv2 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part2.csv"
csv3 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part3.csv"
csv4 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part4.csv"

# Create transposed file
csv_cols_list = {}
for csv in [csv2, csv3, csv4]:
    mylst = []
    for chunk in pd.read_csv(csv, index_col=0, header=None, low_memory=False, chunksize=100000):
        mylst.append(chunk)
    gse_df = pd.concat(mylst, axis=0)
    del mylst
    gse_df = gse_df.transpose()
    gse_df.rename(columns={np.NaN: "samples"}, inplace=True)
    csv_cols_list[csv] = gse_df.columns
    csv_output_name = csv.replace(".csv", "_t.csv")
    gse_df.to_csv(csv_output_name, index=False)

print(csv_cols_list)
# print("Take common features") 
# comm_cpgs1 = list(set(list(gse1_df.columns)).intersection(list(gse2_df.columns)))
# comm_cpgs2 = list(set(list(gse3_df.columns)).intersection(list(gse4_df.columns)))
# comm_cpgs = comm_cpgs1.intersection(comm_cpgs2)

# print(comm_cpgs)
# print(len(comm_cpgs2))
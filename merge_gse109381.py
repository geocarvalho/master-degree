import numpy as np
import pandas as pd

csv1 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part1.csv"
csv2 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part2.csv"
csv3 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part3.csv"
csv4 = "/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_part4.csv"

gse1_df = pd.read_csv(csv1, index_col=0, header=None).transpose()
gse2_df = pd.read_csv(csv2, index_col=0, header=None).transpose()
gse3_df = pd.read_csv(csv3, index_col=0, header=None).transpose()
gse4_df = pd.read_csv(csv4, index_col=0, header=None).transpose()

gse1_df.rename(columns={np.NaN: "samples"}, inplace=True)
gse2_df.rename(columns={np.NaN: "samples"}, inplace=True)
gse3_df.rename(columns={np.NaN: "samples"}, inplace=True)
gse4_df.rename(columns={np.NaN: "samples"}, inplace=True)

comm_cpgs1 = list(set(list(gse1_df.columns)).intersection(list(gse2_df.columns)))
comm_cpgs2 = list(set(list(gse3_df.columns)).intersection(list(gse4_df.columns)))
comm_cpgs = comm_cpgs1.intersection(comm_cpgs2)

print(comm_cpgs)
print(len(comm_cpgs2))
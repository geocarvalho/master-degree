import pandas as pd

diff_probes = "/home/watson/george/master-degree/download_geo/GSE51032/diff_probes.csv"
dp_df = pd.read_csv(diff_probes, header=0, names=["probes", "C50-normal"])
filtered = dp_df[(dp_df["C50-normal"]==-1) | (dp_df["C50-normal"]==1)]
probes = filtered["probes"]
probes.to_csv("diff_probes_list.csv", header=False, index=False)
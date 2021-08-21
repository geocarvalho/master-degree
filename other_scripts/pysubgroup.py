import pysubgroup as ps
import pandas as pd
import numpy as np

data = "../download_geo/GSE51032/GSE51032_bvalues_filtered_cn_pheno_factor.csv"
df = pd.read_csv(data)

target = ps.BinaryTarget ('target', True)
searchspace = ps.create_selectors(df, ignore=['target'])
task = ps.SubgroupDiscoveryTask (
    df, 
    target, 
    searchspace, 
    result_set_size=5, 
    depth=2, 
    qf=ps.WRAccQF())
result = ps.BeamSearch().execute(task)
result_df = result.to_dataframe()
restult_df.to_csv("../download_geo/GSE51032/GSE51032_beamsearch_result.csv", index=False)
import pysubgroup as ps
import pandas as pd
import numpy as np

data = "../download_geo/GSE51032/GSE51032_bvalues_filtered_cn_pheno_factor.csv"
df = pd.read_csv(data)
cols = [col for col in df.columns if col != "time_to_diagnosis_classes"]
df = df[cols]

target = ps.BinaryTarget('target', True)
searchspace = ps.create_selectors(df, ignore=['target'])
task = ps.SubgroupDiscoveryTask (
    df, 
    target, 
    searchspace, 
    result_set_size=10, 
    depth=50, 
    qf=ps.WRAccQF())
result = ps.BeamSearch().execute(task)
result_df = result.to_dataframe()
result_df.to_csv("../download_geo/GSE51032/GSE51032_beamsearch_result_dp50_sz10_nodiagnosis.csv", index=False)
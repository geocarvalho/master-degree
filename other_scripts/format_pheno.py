import pandas as pd

data = "GSE109381/GSE109381_all_phenotype.csv"
df = pd.read_csv(data)
df = df.rename(
    {"Unnamed: 0": "sample", 
    "methylation class:ch1": "methylation_class"}, axis=1)
df_filtered = df[["sample", "methylation_class"]]
df_filtered["class"] = pd.factorize(df_filtered["methylation_class"])[0]+1
df_out = df_filtered.set_index("sample").transpose()
df_out.to_csv("GSE109381/GSE109381_all_phenotype_filtered_t.csv", index=False)
# just classes
df_filtered = df_filtered.drop(columns=["methylation_class"]) 
df_out = df_filtered.set_index("sample").transpose()
df_out.to_csv("GSE109381/GSE109381_all_phenotype_filtered_tclasses.csv", index=False)

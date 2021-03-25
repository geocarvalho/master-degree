import pandas as pd

data = "/home/watson/george/master-degree/download_geo/GSE51032/GSE51032_bvalues_filtered.csv"
output = data.replace(".csv", "_pheno.csv")
df = pd.read_csv(data)
df.set_index("probes", inplace=True)
df = df.astype(float)
cols = [col.split("_")[0] for col in df.columns]
df.columns = cols
df[df<=0.2] = 0
df[df>=0.8] = 2
df[(df>0.2) & (df<0.8)] = 1
df = df.astype(int)

pheno = "/home/watson/george/master-degree/download_geo/GSE51032/GSE51032_all_phenotype_normals.csv"
pheno_df = pd.read_csv(pheno)
# cols = ['Unnamed: 0', 'age at menarche:ch1', 'age:ch1', 'gender:ch1', 'time to diagnosis:ch1', 'cancer type (icd-10):ch1']
cols = ['Unnamed: 0', 'gender:ch1', 'cancer type (icd-10):ch1']
pheno_df = pheno_df[cols]
# cols = ['probes', 'age_menarche', 'age', 'gender', 'time_diagnosis', 'cancer_type']
cols = ['probes', 'gender', 'cancer_type']
pheno_df.columns = cols
pheno_df.set_index("probes", inplace=True)
pheno_df = pheno_df[(pheno_df["cancer_type"] == "normal") | (pheno_df["cancer_type"] == "C50")]
pheno_df = pheno_df.T

merge = pd.concat([df, pheno_df])
merge = merge.T
merge.to_csv(output)
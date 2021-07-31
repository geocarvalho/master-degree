import pandas as pd

def class_age(row):
    """ Create class for age <54 and >=54 """
    if row["age"] < 54:
        return 0
    else:
        return 1
    
def class_diagnosis(row):
    """ Create class for time_diagnosis <5 and >=5 """
    if row["time_diagnosis"] < 5:
        return 0
    else:
        return 1

def class_menarche(row):
    """ Create class for age_menarche <13 and >=13 """
    if row["age_menarche"] < 13:
        return 0
    else:
        return 1


data = "../download_geo/GSE51032/GSE51032_bvalues_filtered.csv"
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

pheno = "../download_geo/GSE51032/GSE51032_all_phenotype_normals.csv"
pheno_df = pd.read_csv(pheno)
# cols = ['Unnamed: 0', 'age:ch1', 'gender:ch1', 'time to diagnosis:ch1', 'cancer type (icd-10):ch1']
cols = ['Unnamed: 0', 'age at menarche:ch1', 'age:ch1', 'gender:ch1', 'time to diagnosis:ch1', 'cancer type (icd-10):ch1']
# cols = ['Unnamed: 0', 'gender:ch1', 'cancer type (icd-10):ch1']
pheno_df = pheno_df[cols]
# cols = ['probes', 'age', 'gender', 'time_diagnosis', 'cancer_type']
cols = ['probes', 'age_menarche', 'age', 'gender', 'time_diagnosis', 'cancer_type']
# cols = ['probes', 'gender', 'cancer_type']
pheno_df.columns = cols
pheno_df.set_index("probes", inplace=True)
pheno_df = pheno_df[(pheno_df["cancer_type"] == "normal") | (pheno_df["cancer_type"] == "C50")]
# Criar binário da classe alvo
pheno_df["target"] = pheno_df["cancer_type"].map(
    {"normal": 0, "C50": 1}
    # {"normal": 0, "C18": 1}
)
pheno_df["bin_gender"] = pheno_df["gender"].map(
    {"F": 0, "M": 1}
)
pheno_df["bin_age"] = pheno_df.apply(class_age, axis=1)
pheno_df["bin_diagnosis"] = pheno_df.apply(class_diagnosis, axis=1)
pheno_df["bin_menarche"] = pheno_df.apply(class_menarche, axis=1)
pheno_df = pheno_df[["bin_gender", "bin_age", "bin_diagnosis", "bin_menarche", "target"]]
pheno_df = pheno_df.T

merge = pd.concat([df, pheno_df])
merge = merge.T
merge.to_csv(output)

# Rodar o SSDP+
# java -jar -Xmx50G -XX:+UseConcMarkSweepGC ssdp_plus/out/artifacts/ssdp_plus_jar/ssdp_plus.jar /home/watson/george/master-degree/download_geo/GSE51032/GSE51032_bvalues_filtered_pheno.csv 5 0.5 1 | tee a  /home/watson/george/master-degree/create_input_ssdp/log_files/GSE51032_50_log.txt
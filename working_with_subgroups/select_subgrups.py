import pandas as pd

file = "GSE109381_bvalues_gbm_merged_classes_formatted_nopoints.csv"
csv = "/home/watson/george/master-degree/GSE109381/GBM_class_sample_categories.csv"

# similiridade de 50
## subgrupo 1
## cg23847017 = 0,cg11669516 = 0,cg17146473 = 2,cg08206623 = 1
cols = ["cg23847017", "cg11669516", "cg17146473", "cg08206623", "class_num"]
df = pd.read_csv(file, usecols=cols)
csv_df = pd.read_csv(csv, skiprows=1, names=["samples", "class", "class_num_2"])

df_merge = pd.concat([df, csv_df], axis=1)
df_subgroup = df_merge.query("cg23847017 == 0 & cg11669516 == 0 & cg17146473 == 2 & cg08206623 == 1 & class_num == 5")

samples = df_merge["samples"].to_list()
with open('samples_50_sg1.txt', 'w') as f:
    f.write("%s\n" % samples)
import pandas as pd


ss = pd.read_csv("samplesheet_GSE109381.csv")
ss["sample"] = ss.apply(lambda link: link["supplementary_file"].split("/")[-1].split(".")[0], axis=1)
ss["Slide"] = ss.apply(lambda sample: sample["sample"].split("_")[1], axis=1)
ss["Array"] = ss.apply(lambda sample: sample["sample"].split("_")[2], axis=1)
ss["sample"] = ss.apply(lambda sample: sample["sample"].replace("_Grn", ""), axis=1)

output = ss[["sample", "gender", "sex", "Slide", "Array"]]
output.to_csv("samplesheet_GSE109381_simple.csv", index=False)
import pandas as pd

list_of_probes = "./GSE51032/GSE51032_sigCpGs_cn.csv"
b_values = "./GSE51032/GSE51032_bvalues.csv"
output = b_values.replace(".csv", "_filtered_cn.csv")

with open(list_of_probes, "r") as buffer:
    probes_list = [probe.strip("\n") for probe in buffer.readlines()]
    buffer.close()

bvalues_df = pd.read_csv(b_values)
bvalues_df.rename(columns={"Unnamed: 0": "probes"}, inplace=True)
bvalues_df.set_index("probes", inplace=True)
bvalues_filtered_df  = bvalues_df[bvalues_df.index.isin(probes_list)]
bvalues_filtered_df.to_csv(output)
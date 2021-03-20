import pandas as pd

list_of_probes = "/home/watson/george/master-degree/create_input_ssdp/diff_probes_list.csv"
b_values = "/home/watson/george/master-degree/download_geo/GSE51032/GSE51032_bvalues.csv"

bvalues_df = pd.read_csv(b_values)
bvalues_df.rename(columns={"Unnamed: 0": "probes"}, inplace=True)
bvalues_df.set_index("probes", inplace=True)
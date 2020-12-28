import pandas as pd
import json


def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

json_file = "subgroups_1.json"
with open(json_file, "r") as read_file:
    data = json.load(read_file)

cols = list(data["10"].keys())
df_10 = pd.DataFrame(index=cols, columns=cols)
df_50 = df_10.copy()
df_90 = df_10.copy()
dataframes = {"10": df_10, "50": df_50, "90": df_90}

for similarity in data.keys():
    df = dataframes[similarity]
    for subgroup, subgroup_lst in data[similarity].items():
        for second in range(int(subgroup)+1, len(data[similarity].keys())+1):
            jac = jaccard_similarity(subgroup_lst, data[similarity][str(second)])
            df[str(second)][subgroup] = jac

writer = pd.ExcelWriter("jaccard_index_subgoups.xlsx",engine='xlsxwriter')   
for dataframe, sheet in zip(list(dataframes.values()), list(dataframes.keys())):
    dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0)   
writer.save()

# Valdidation:
# print("Validation for 8 and 10")
# val_90_9 = data["90"]["8"]
# val_90_10 = data["90"]["10"]
# jac_val = jaccard_similarity(val_90_9, val_90_10)
# print(jac_val)
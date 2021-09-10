import pandas as pd

data = "../download_geo/GSE51032/GSE51032_bvalues_filtered_cn_pheno_factor.csv"
df = pd.read_csv(data)

# pysubgroup common cpgs 50-10
df[(df['cg27544799']==0) & (df['cg04551440']==0) & (df['cg00475161']==0) & (df['cg24154336']==16) & (df['ch.9.2042397F']==0) & (df['cg21279459']==9) & (df['cg20972117']==0) & (df['cg07610406']==0) & (df['cg15717853']==0) & (df['cg15859496']==14) & (df['cg03801179']==0) & (df['cg10437571']==0) & (df['cg02647825']==15) & (df['cg06721712']==15) & (df['cg14651446']==0) & (df['cg26821539']==0) & (df['cg01892567']==16) & (df['cg06932756']==0) & (df['cg00015930']==16) & (df['cg09449449']==15) & (df['cg01622006']==0) & (df['cg07859799']==15) & (df['cg08905629']==0) & (df['cg25373297']==0) & (df['cg01946760']==0) & (df['cg07062711']==0) & (df['cg03161767']==0) & (df['cg03255749']==15) & (df['cg00505073']==15) & (df['cg00262415']==0) & (df['cg10318121']==0) & (df['cg01962937']==0) & (df['cg01016662']==0) & (df['cg05046026']==0) & (df['cg06770554']==0) & (df['cg05006947']==9) & (df['cg12071806']==0) & (df['cg12589307']==0) & (df['cg01031032']==0) & (df['cg12510286']==0)]
# nothing interesting about the phenotypes
# 'age_at_menarche_int', 'time_to_diagnosis_classes', 'age_classes', 'CD8T_classes', 'CD4T_classes', 'NK_classes', 'Bcell_classes', 'Mono_classes', 'Gran_classes', 'PredictedSmokingStatus', 'target']
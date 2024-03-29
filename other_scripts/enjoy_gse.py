import pandas as pd
import numpy as np
# import seaborn as sns
import xgboost as xgb
# import matplotlib.pyplot as plt
from sklearn.svm import LinearSVC
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, train_test_split

def threshold_bvalue(value):
    if float(value) <= 0.2:
        return "low"
    elif float(value) >=0.8:
        return "high"
    else:
        return "middle"

def prepare_df(gse, pheno):
    """ Merge GSE and phenotype files into a DataFrame
    """
    # GSE
    gse_df = pd.read_csv(gse, index_col=0, header=None).T
    gse_df = gse_df.rename(columns={np.NaN: "samples"})
    gse_df["geo_accession"] = gse_df["samples"].apply(lambda lst: lst.split("_")[0])

    # Transform b-values in methylation rate groups
    cols = [col for col in gse_df.columns if col.startswith("cg")]
    gse_df[cols] = gse_df[cols].applymap(threshold_bvalue)
    
    # phenotype
    pheno_df = pd.read_csv(pheno)
    pheno_df = pheno_df.rename({
        "Unnamed: 0": "samples", "subtype_ihc:ch1": "subtype"
        }, axis=1)

    l_pheno_df = pheno_df[["geo_accession", "subtype"]]
    gse_pheno_df = gse_df.merge(l_pheno_df, on="geo_accession")
    gse_pheno_df.drop(["geo_accession", "samples"], axis=1, inplace=True)

    return gse_pheno_df

# Organize dataset
gse_file1 = "/home/genomika/george/master-degree/GSE72245/GSE72245_bvalues.csv"
pheno1 = "/home/genomika/george/master-degree/GSE72245/GSE72245_all_phenotype.csv"
gse_pheno_df1 = prepare_df(gse_file1, pheno1)

gse_file2 = "/home/genomika/george/master-degree/GSE72251/GSE72251_bvalues.csv"
pheno2 = "/home/genomika/george/master-degree/GSE72251/GSE72251_all_phenotype.csv"
gse_pheno_df2 = prepare_df(gse_file2, pheno2)

# gse_file3 = "/home/genomika/george/master-degree/GSE72254/GSE72254_bvalues.csv"
# pheno3 = "/home/genomika/george/master-degree/GSE72254/GSE72254_all_phenotype.csv"
# gse_pheno_df3 = prepare_df(gse_file3, pheno3)

# Concatenate DFs
gse_pheno_df = pd.concat([gse_pheno_df1, gse_pheno_df2], sort=False).reset_index(drop=True) 

gse_pheno_df["subtype"] = gse_pheno_df["subtype"].map({
    "LumB":1, "Basal":2, "HER2":3, "LumA":4
    })


subtype = gse_pheno_df["subtype"]
gse_pheno_df.drop("subtype", axis=1, inplace=True)
subtype.drop(subtype.index[[141, 216]], inplace=True)
gse_pheno_df.drop(gse_pheno_df.index[[141, 216]], inplace=True)

gse_pheno_df = gse_pheno_df.astype(float)
subtype = subtype.astype(float)

# TODO: PROBLEM WITH NAN values
subtype
# Plot subtype distribution
# plt.figure(figsize=(12,5))
# sns.countplot(x=subtype, color='mediumseagreen')
# plt.title('ECancer subtype class distribution', fontsize=16)
# plt.ylabel('Class Counts', fontsize=16)
# plt.xlabel('Class Label', fontsize=16)
# plt.xticks(rotation='vertical')
 
# ML algorithms based on:
# https://www.freecodecamp.org/news/multi-class-classification-with-sci-kit-learn-xgboost-a-case-study-using-brainwave-data-363d7fca5f69/
# Run RF
# %%time
pl_random_forest = Pipeline(steps=[('random_forest', RandomForestClassifier())])
scores = cross_val_score(pl_random_forest, gse_pheno_df, subtype, cv=10,scoring='accuracy')
print('Accuracy for RandomForest : ', scores.mean())

# Dataset 1
# Accuracy for RandomForest :  0.5681818181818181
# CPU times: user 1min 2s, sys: 3.61 s, total: 1min 6s
# Wall time: 1min 6s

# Both datasets
# Accuracy for RandomForest :  0.6855072463768116

# Run LR
# %%time
pl_log_reg = Pipeline(steps=[('scaler',StandardScaler()),
                             ('log_reg', LogisticRegression(multi_class='multinomial', solver='saga', max_iter=200))])
scores = cross_val_score(pl_log_reg, gse_pheno_df, subtype, cv=10,scoring='accuracy')
print('Accuracy for Logistic Regression: ', scores.mean())

# Dataset1
# Accuracy for Logistic Regression:  0.5840909090909091
# CPU times: user 55min 14s, sys: 1min 7s, total: 56min 21s
# Wall time: 55min 1s

# Both datasets
# Accuracy for Logistic Regression:  0.6942028985507247

scaler = StandardScaler()
scaled_df = scaler.fit_transform(gse_pheno_df)
pca = PCA(n_components = 20)
pca_vectors = pca.fit_transform(scaled_df)
for index, var in enumerate(pca.explained_variance_ratio_):
    print("Explained Variance ratio by Principal Component ", (index+1), " : ", var)
# Scatter plot from the 20 components
# plt.figure(figsize=(25,8))
# sns.scatterplot(x=pca_vectors[:, 0], y=pca_vectors[:, 1], hue=subtype)
# plt.title('Principal Components vs Class distribution', fontsize=16)
# plt.ylabel('Principal Component 2', fontsize=16)
# plt.xlabel('Principal Component 1', fontsize=16)
# plt.xticks(rotation='vertical')

# %%time
# pl_log_reg_pca = Pipeline(steps=[('scaler',StandardScaler()),
#                              ('pca', PCA(n_components = 2)),
#                              ('log_reg', LogisticRegression(multi_class='multinomial', solver='saga', max_iter=200))])
# scores = cross_val_score(pl_log_reg_pca, gse_pheno_df, subtype, cv=10,scoring='accuracy')
# print('Accuracy for Logistic Regression with 2 Principal Components: ', scores.mean())

# Accuracy for Logistic Regression with 2 Principal Components:  0.40681818181818186

# %%time
# pl_log_reg_pca_10 = Pipeline(steps=[('scaler',StandardScaler()),
#                              ('pca', PCA(n_components = 10)),
#                              ('log_reg', LogisticRegression(multi_class='multinomial', solver='saga', max_iter=200))])
# scores = cross_val_score(pl_log_reg_pca_10, gse_pheno_df, subtype, cv=10,scoring='accuracy')
# print('Accuracy for Logistic Regression with 10 Principal Components: ', scores.mean())

# Accuracy for Logistic Regression with 10 Principal Components:  0.5515151515151515
# CPU times: user 29min 11s, sys: 5min 22s, total: 34min 34s
# Wall time: 3min 17s

# %%time
pl_mlp = Pipeline(steps=[('scaler',StandardScaler()),
                             ('mlp_ann', MLPClassifier(hidden_layer_sizes=(1275, 637)))])
scores = cross_val_score(pl_mlp, gse_pheno_df, subtype, cv=10,scoring='accuracy')
print('Accuracy for ANN : ', scores.mean())

# Accuracy for ANN :  0.5613636363636364
# CPU times: user 4h 47min 35s, sys: 7h 21min 40s, total: 12h 9min 15s
# Wall time: 6h 3min 10s

# %%time
pl_svm = Pipeline(steps=[('scaler',StandardScaler()),
                             ('pl_svm', LinearSVC())])
scores = cross_val_score(pl_svm, gse_pheno_df, subtype, cv=10,scoring='accuracy')
print('Accuracy for Linear SVM : ', scores.mean())

# Dataset 1
# Accuracy for Linear SVM :  0.6098484848484848
# CPU times: user 6min 6s, sys: 16.9 s, total: 6min 23s
# Wall time: 4min 55s

# %%time
pl_xgb = Pipeline(steps=
                  [('xgboost', xgb.XGBClassifier(objective='multi:softmax'))])
scores = cross_val_score(pl_xgb, gse_pheno_df, subtype, cv=10)
print('Accuracy for XGBoost Classifier : ', scores.mean())

# Accuracy for XGBoost Classifier :  0.6787878787878789
# CPU times: user 15h 57min 11s, sys: 1h 15min 57s, total: 17h 13min 8s
# Wall time: 23min 59s

# Both datasets
# Accuracy for XG Boost Classifier :  0.7416666666666667
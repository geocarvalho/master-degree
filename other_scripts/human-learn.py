# pip install human-learn
# https://calmcode.io/model-mining/introduction.html
import numpy as np
import pandas as pd
from hulearn.classification import FunctionClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import precision_score, recall_score, accuracy_score, make_scorer

def subgroup_1(dataf):
    """ Rule with first subgroup """
    return np.array(df["cg00015930"]==16) & (df["cg00262415"]==0) & (df["cg00475161"]==0) & (df["cg00505073"]==15) & \
    (df["cg01016662"]==0) & (df["cg01031032"]==0) & (df["cg01622006"]==0) & (df["cg01892567"]==16) & (df["cg01946760"]==0) & \
    (df["cg01962937"]==0) & (df["cg02647825"]==15) & (df["cg03161767"]==0) & (df["cg03255749"]==15) & (df["cg03801179"]==0) & \
    (df["cg04104256"]==0) & (df["cg04551440"]==0) & (df["cg05006947"]==9) & (df["cg05046026"]==0) & (df["cg05664296"]==15) & \
    (df["cg06721712"]==15) & (df["cg06770554"]==0) & (df["cg06932756"]==0) & (df["cg07062711"]==0) & (df["cg07610406"]==0) & \
    (df["cg07859799"]==15) & (df["cg08905629"]==0) & (df["cg09449449"]==15) & (df["cg10318121"]==0) & (df["cg10437571"]==0) & \
    (df["cg12071806"]==0) & (df["cg12510286"]==0) & (df["cg12589307"]==0) & (df["cg14651446"]==0) & (df["cg15717853"]==0) & \
    (df["cg15859496"]==14) & (df["cg17386812"]==15) & (df["cg20972117"]==0) & (df["cg21279459"]==9) & (df["cg24154336"]==16) & \
    (df["cg25373297"]==0) & (df["cg26821539"]==0) & (df["cg27544799"]==0) & (df["ch.9.2042397F"]==0).astype(int)

def subgroup_2(df):
    """ Rule with secondsubgroup """
    return np.array(df[(df["cg00015930"]==16) & (df["cg00262415"]==0) & (df["cg00475161"]==0) & (df["cg00505073"]==15) & (df["cg01016662"]==0) & \
    (df["cg01031032"]==0) & (df["cg01622006"]==0) & (df["cg01892567"]==16) & (df["cg01946760"]==0) & (df["cg01962937"]==0) & \
    (df["cg02647825"]==15) & (df["cg03161767"]==0) & (df["cg03255749"]==15) & (df["cg03801179"]==0) & (df["cg04104256"]==0) & \
    (df["cg04551440"]==0) & (df["cg05006947"]==9) & (df["cg05046026"]==0) & (df["cg06721712"]==15) & (df["cg06770554"]==0) & \
    (df["cg06932756"]==0) & (df["cg07062711"]==0) & (df["cg07610406"]==0) & (df["cg07859799"]==15) & (df["cg08905629"]==0) & \
    (df["cg09449449"]==15) & (df["cg10318121"]==0) & (df["cg10437571"]==0) & (df["cg12071806"]==0) & (df["cg12510286"]==0) & \
    (df["cg12589307"]==0) & (df["cg14651446"]==0) & (df["cg15717853"]==0) & (df["cg15859496"]==14) & (df["cg16661579"]==0) & \
    (df["cg17386812"]==15) & (df["cg20972117"]==0) & (df["cg21279459"]==9) & (df["cg24154336"]==16) & (df["cg25373297"]==0) & \
    (df["cg26821539"]==0) & (df["cg27544799"]==0) & (df["ch.9.2042397F"]==0)]).astype(int)

data = "../download_geo/GSE51032/GSE51032_bvalues_filtered_cn_pheno_factor.csv"
df = pd.read_csv(data)
X, y = df.drop(columns=['target']), df['target']
mod = FunctionClassifier(subgroup_1)
mod.fit(X, y).predict(X)

grid = GridSearchCV(mod,
                    cv=2,
                    scoring={'accuracy': make_scorer(accuracy_score),
                            'precision': make_scorer(precision_score),
                            'recall': make_scorer(recall_score)},
                    refit='accuracy')
grid.fit(X, y)

grid_res = (pd.DataFrame(grid.cv_results_)
  [['param_threshold', 'mean_test_recall', 'mean_test_accuracy', 'mean_test_precision']]
  .set_index('param_threshold'))
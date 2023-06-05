# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from autogluon.tabular import TabularDataset, TabularPredictor
import json

if not os.path.exists("./2-output"):
        os.makedirs("./2-output")

#### Prepare the data, feature extraction (dimension reduction)
# read in, transpose and fit PCA to the training set
beta_train= pd.read_feather('./1-output/beta_train.feather').transpose() # read in data frame and transpose, rownames are lost, but beta_train and _test are in same order
components = beta_train.shape[0] -1 # number of samples -1
pca = PCA(n_components = components) # pca object
train_data = pca.fit_transform(beta_train)

# read in, transpose and project test set to PCA fit test set
beta_test = pd.read_feather('./1-output/beta_test.feather').transpose() # reads in dictionary where keys are the name of objects 
test_data = pca.transform(beta_test)

# convert np.arrays to pd.dataframes
train_data = pd.DataFrame(train_data)
train_data['basename'] = list(beta_train.index)
test_data = pd.DataFrame(test_data)
test_data['basename'] = list(beta_test.index)

# add pheno labels and age/ic_retrain as additional features
pheno_train = pd.read_feather('./1-output/pheno_train.feather')
train_data = pd.merge(pheno_train, train_data)

pheno_test = pd.read_feather('./1-output/pheno_test.feather')
test_data = pd.merge(pheno_test, test_data)

#### train models to predict breast,ovarian or endometrial cancer
# remove basename from training data, this is not a feature to train on
train_data =  train_data.drop(columns=['basename'])

#arguments for training
label = 'type' # variable to train on
print("Summary of type variable: \n", train_data[label].describe())
time_limit = 10800 #time limit in seconds, exagerated, but care little about runtime
save_path = "./2-output/agModels"

#train
predictor = TabularPredictor(label=label,path=save_path).fit(train_data,time_limit=time_limit,presets='best_quality')

#### test (best) model(s) on external dataset
# load trained predictor
#label = 'type' # variable to train on
#save_path = "./2-output/agModels"
#predictor = TabularPredictor.load(save_path)
predictor.get_model_best() #Best model 'WeightedEnsemble_L3'

# prepare test data
y_test = test_data[label]
test_data_nolab =  test_data.drop(columns=[label])
# check
test_data_nolab.head()

# predict using best model
y_pred = predictor.predict(test_data_nolab)
print("Predictions: \n", y_pred)
y_pred_proba = predictor.predict_proba(test_data_nolab)
print("Probabilities: \n", y_pred_proba)
perf = predictor.evaluate_predictions(y_true=y_test, y_pred=y_pred, auxiliary_metrics=True)
with open(os.path.join('./2-output', 'performance.txt'), 'w') as performance:
     performance.write(json.dumps(perf))

out1 = predictor.leaderboard(test_data,silent=True)
out1.to_csv("./2-output/leaderboard.csv")

# save predicitons and probabilities for 3C external validation
pred = pd.DataFrame({"basename": test_data["basename"],
                    "WID-3C_PC_autogluon_prediction": y_pred})# make a df with basename, y_pred
pred = pred.join(y_pred_proba) # add probabilities
pred.to_csv("./2-output/predictions.csv")

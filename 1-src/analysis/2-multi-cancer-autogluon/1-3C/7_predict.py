# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

import os
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor
import json

if not os.path.exists("./7-output"):
        os.makedirs("./7-output")

# load trained predictor
label = 'type' # variable to train on
save_path = "./6-output/agModels"
predictor = TabularPredictor.load(save_path)
predictor.get_model_best() #Best model: 'WeightedEnsemble_L2'

# load internal validation set
test_data =  TabularDataset('./5-output/test3C.csv')
y_test = test_data[label]
test_data_nolab =  test_data.drop(columns=[label])
# check
test_data_nolab.head()

## deal with missing features, will have to do for external datasets
train_data = TabularDataset('./5-output/train3Cdisc.csv') # note that age and estimated ic content are input features
features_train = list(train_data.columns.values)
features_test = list(test_data.columns.values)
missing = set(sorted(features_train)).difference(sorted(features_test))
empty_missing = pd.DataFrame(np.nan, index= list(test_data.index.values), columns=list(missing))
test = test_data_nolab.join(empty_missing)

# predict using best model
y_pred = predictor.predict(test)
print("Predictions: \n", y_pred)
y_pred_proba = predictor.predict_proba(test)
print("Probabilityes: \n", y_pred_proba)
perf = predictor.evaluate_predictions(y_true=y_test, y_pred=y_pred, auxiliary_metrics=True)
with open(os.path.join('./7-output', 'performance.txt'), 'w') as performance:
     performance.write(json.dumps(perf))

test2 = test_data.join(empty_missing)
out1 = predictor.leaderboard(test2,silent=True) 
out1.to_csv("./7-output/leaderboard.csv") 

# save predicitons and probabilities for 3C external validation
pred = pd.DataFrame({"basename": test["basename"], 
                    "WID-3C_prediction": y_pred})# make a df with basename, y_pred
pred = pred.join(y_pred_proba)# add probabilities
pred.to_csv("./7-output/predictions.csv")

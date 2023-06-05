# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

import os
import pandas as pd
import numpy as np
from autogluon.tabular import TabularDataset, TabularPredictor

if not os.path.exists("./4-output"):
        os.makedirs("./4-output")

# load trained predictor
label = 'type' # variable to train on
save_path = "./3-output/agModels"
predictor = TabularPredictor.load(save_path) #Best model: "WeightedEnsemble_L2"

# load internal validation set
test_data =  TabularDataset('./2-output/validate3Cdisc.csv')
y_test = test_data[label]
test_data_nolab =  test_data.drop(columns=[label])
# check
test_data_nolab.head()

# predict using best model
y_pred = predictor.predict(test_data_nolab)
print("Predictions: \n", y_pred)
perf = predictor.evaluate_predictions(y_true=y_test, y_pred=y_pred, auxiliary_metrics=True)

out1 = predictor.leaderboard(test_data,silent=True) 
out1.to_csv("./4-output/leaderboard_internalval.csv") 

### save importance features for pruning and pruned training set next round
importance = predictor.feature_importance(test_data)
importance.to_csv("./4-output/importance.csv")
#importance = pd.read_csv ("./4-output/importance.csv")
importance[["importance"]].describe()

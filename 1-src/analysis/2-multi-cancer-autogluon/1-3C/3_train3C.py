# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

import os
from autogluon.tabular import TabularDataset, TabularPredictor

if not os.path.exists("./3-output"):
        os.makedirs("./3-output")

### train models to predict control, breast, ovarian or endometrial cancer

# initial features selection was done by selecting top 5000 immune and top 5000 epithelial db for each cancer ~17000 CpGs
# also age and estimated immune cell composition will be given as features
# internal validation set is excluded here, and will be used for checking feature importance

#read in training data saved as .csv
train_data = TabularDataset('./2-output/train3Cdisc.csv')

#arguments for training
label = 'type' # variable to train on
print("Summary of type variable: \n", train_data[label].describe())
time_limit = 3600 #time limit in seconds
save_path = "./3-output/agModels"
#metric = 'roc_auc' # does not work for multi-class

#train
predictor = TabularPredictor(label=label,path=save_path).fit(train_data,time_limit=time_limit,presets='best_quality')


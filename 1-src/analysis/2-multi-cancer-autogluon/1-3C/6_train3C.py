# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

import os
from autogluon.tabular import TabularDataset, TabularPredictor

if not os.path.exists("./6-output"):
        os.makedirs("./6-output")

### train models to predict control, breast, ovarian or endometrial cancer, second round

# full set of samples (training + internal validation) used
# only 28 features pruned compared to first training round

#read in training data saved as .csv
train_data = TabularDataset('./5-output/train3Cdisc.csv')

#arguments for training
label = 'type' # variable to train on
print("Summary of type variable: \n", train_data[label].describe())
time_limit = 10800 #time limit in seconds, little desperate
save_path = "./6-output/agModels"

#train
predictor = TabularPredictor(label=label,path=save_path).fit(train_data,time_limit=time_limit,presets='best_quality')


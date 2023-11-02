# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Check performance with increasing training set size
# penalized regression, 10K cross-validation for each method
# full FHS set used for training (2639 peripheral blood samples, GSE55763)
# Hannum subset (479 samples withinsame age range as training, GSE40279) for testing accuracy
# Technical reps from FHS for testing repeatability (2x36 reps)
# for training size sample input choose three levels of granularity (95 sample sizes/data points checked):
# fine: n min = 34, n max = 54, n step = 5 
# medium: n min = 59, n max = 1619, n step = 20 
# coarse: n min = 1639, n max = 2639, n step = 100 

### final output 1. dataframe:
## training parameters kept:
# method: PC.clock, ElNet (no PC)
# input.CpGs: all
# train.sample.size: nr training samples
## model parameters kept:
# n.features : nr CpGs, or components in trained model
## test results kept:
# val.slope.regress: slope regression age ~ predicted age (accuracy)
# val.intercept.regress: intercept regression age ~ predicted age
# val.MSE: mean squared error age ~ predicted age
# rep.ICC: ICC rep set smears (precision)
# rep.lbound.ICC
# rep.ubound.ICC
### final output 2. list 
## training parameters kept:
# age distribution training samples

# Same experiment is repeated three times, results are then collated results

library(tidyverse)
library(ggplot2)

dir.create("./7-output")

### Source scripts for training+predicting age in validation sets for testing accuracy/precision ###----

source("../../functions/regression-train-validate.R")
source("../../functions/PC-train-validate.R")

### Prepare data ###----
## Beta
beta_tr_full <- readRDS("../../../0-data/beta/beta_BloodFHS_450K_imputed.Rds") # Training set - FHS, missing values imputed
beta_val <- readRDS("../../../0-data/beta/beta_Hannum_ExtVal.Rds") # Validation set accuracy - Hannum subset
beta_rep <- readRDS("../../../0-data/beta/beta_BloodRep_450K_imputed.Rds") # validation set precision - technical reps FHS, missing values imputed

## arrange pheno files
#training
pheno_tr_full <- readRDS("../../../0-data/pheno/pheno_BloodFHS_450K.Rds")
pheno_tr_full <- pheno_tr_full %>% 
  select(basename,age) %>% 
  filter(pheno_tr_full$basename %in% colnames(beta_tr_full)) %>% 
  droplevels()
pheno_tr_full <- pheno_tr_full[match(colnames(beta_tr_full),pheno_tr_full$basename),]
identical(colnames(beta_tr_full),pheno_tr_full$basename)

# validation accuracy
pheno_val <- readRDS("../../../0-data/pheno/pheno_Hannum_ExtVal.Rds") %>% dplyr::rename(age=Age,basename=sampleID)
pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(colnames(beta_val),pheno_val$basename)
age_val <- pheno_val %>% pull(age) %>% as.numeric()

# validation precision
pheno_rep <- readRDS("../../../0-data/pheno/pheno_BloodRep_450K.Rds")
pheno_rep <- pheno_rep %>%
  dplyr::rename(patient_ID=sample) %>%
  select(basename,patient_ID,rep) %>% #only keep necessary stuff
  droplevels()
pheno_rep <- pheno_rep[match(colnames(beta_rep),pheno_rep$basename),]
identical(colnames(beta_rep),pheno_rep$basename)

### Define list of samples to loop over and save age distribution ###----
nr <- c(seq(from=34,to=54,by=5),seq(from=59,to=1619,by=20),seq(from=1639,to=2639,by=100))

sampleList <- list()
for (i in 1:length(nr)){
  set.seed(12+i)
  sampleList[[i]] <- sample(pheno_tr_full$basename,nr[i])
}
names(sampleList) <- nr
saveRDS(sampleList, file="./7-output/sampleListloop.Rds") # fix list to repeat

# out.list - age distribution for different input sample sets - length (nmax-nmin/n) +1:
out.list <- list()
for (j in 1:length(sampleList)){
  out.list[[j]] <- pheno_tr_full %>% filter(basename %in% sampleList[[j]]) %>% select(age,basename)
}
names(out.list) <- names(sampleList)
saveRDS(out.list, file="./7-output/sampleListloopAges.Rds")

### Train and test different methods ###----
# initialize out.df
out.df <- data.frame(train.sample.size=numeric(0),
                     method=character(0),
                     n.features=numeric(0),
                     sample.type=character(0),
                     val.slope.regress=numeric(0),
                     val.intercept.regress=numeric(0),
                     val.MSE=numeric(0),
                     rep.ICC=numeric(0),
                     rep.lbound.ICC=numeric(0),
                     rep.ubound.ICC=numeric(0))

# PC/Elastic Net
for (j in 1:length(sampleList)){
  
  # take sample subset
  beta_tr <- beta_tr_full %>% select(all_of(sampleList[[j]]))
  pheno_tr <- pheno_tr_full %>% filter(basename %in% all_of(sampleList[[j]]))
  pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
  identical(colnames(beta_tr),pheno_tr$basename)
  age_tr <- pheno_tr %>% pull(age) %>% as.numeric()
  
  ## Stick to approach Higgins to subset all common CpGs from train/test/sets
  # note this is actually a suboptimal machine learning approach, as data is leaking from the test to the train set
  CpGs <- intersect(rownames(beta_tr),intersect(rownames(beta_val),rownames(beta_rep)))
  beta_tr <- beta_tr %>% filter(rownames(beta_tr) %in% CpGs)
  
  # Elastic Net
  result <- regressionTrainvalidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,0.5)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "ElNet"
  result$sample.type <- "blood"
  out.df <- full_join(out.df,result)
  
  # Elastic net on PCA-derived feature space
  result <- PCTrainValidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,0.5)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "PC"
  result$sample.type <- "blood"
  out.df <- full_join(out.df,result)
  
  saveRDS(out.df, file="./7-output/results_training_size.Rds") # overwrite results at end of each loop as a means of saving temporary result
}

#!/usr/bin/env Rscript
# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(dplyr)
library(tidyverse)
library(data.table)

dir.create("./5-output")

## training set, beta, pheno, selected features
#read in
beta <- readRDS("../../../../0-data/beta/beta_3CDisc.Rds")
pheno <- readRDS("./0-output/pheno_3CDisc_ic.Rds")

# prune features based on first training round
features <- readRDS("./1-output/selected_features_3C.Rds") #selected 17863 features first round
importance <- read.csv("./4-output/importance.csv") #feature importance calculated for internal validation set
keep <- importance %>% filter(importance >= 0) %>% pull(X) 
length(keep) # 17837, removing only 28, not very promising
features <- intersect(features,keep)

#save Training set (full set of samples)
datTrain <- beta %>% filter(rownames(beta) %in% features) %>% select(c(pheno$basename)) %>% droplevels()
datTrain <- t(as.matrix(datTrain)) #transpose
datTrain <- as.data.frame(datTrain) %>% rownames_to_column(var="basename")
pheno <- pheno %>% select(basename,age,type, ic_retrain) %>% droplevels()
datTrain <- full_join(pheno,datTrain) %>% drop_na() %>% column_to_rownames(var="basename")
#fwrite(datTrain,file="./5-output/train3Cdisc.csv", row.names=FALSE)

# test set last stage
beta <- readRDS("../../../0-preprocess/0-output/beta_3CExtVal.Rds")
pheno <- readRDS("../../../3-remove-unreliable-probes/2-output/pheno_3CExtVal_rem_EPICV2_lowMI.Rds")
dat <- t(as.matrix(beta))
dat <- as.data.frame(dat) %>% rownames_to_column(var="basename")
pheno <- pheno %>% select(basename,age, type, ic_retrain) %>% droplevels()

#save space, might be still features missing in test set
features_common <- intersect(features,colnames(dat))
dat <- full_join(pheno,dat) %>% drop_na() %>% select(c(features_common,c(basename,type,age,ic_retrain))) 
fwrite(dat,file="./5-output/test3C.csv",row.names=TRUE)

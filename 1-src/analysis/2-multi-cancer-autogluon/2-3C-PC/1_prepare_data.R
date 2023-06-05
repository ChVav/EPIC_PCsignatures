# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

##1_prepare_data so that can read it in Python
# remove already unreliable probes
# replace missing features (CpGs) in test_dataset with max(density(beta of feature in training set))

#!/usr/bin/env Rscript

library(dplyr)
library(tidyverse)
library(stats)
library(feather)

dir.create("./1-output")

## training set pheno
pheno <- readRDS("../1-3C/0-output/pheno_3CDisc_ic.Rds")
pheno <- pheno %>% select(basename,type,age,ic_retrain) %>% droplevels() # use also age and ic_retrain1 as training features
write_feather(pheno,"./1-output/pheno_train.feather")

## training set beta
beta <- readRDS("../../../../0-data/beta/beta_3CDisc.Rds")
removedEPICv2 <- read.csv("../../../../../0-data/flagged_probes/notonEPICv2_names.csv", header=FALSE) %>% pull(V1)
lowMI <- read.csv("../../../../../0-data/flagged_probes/lowMI_names.csv", header=FALSE) %>% pull(V1)
toremove <- unique(c(removedEPICv2,lowMI))

# remove unreliable probes
beta_train <- beta %>% filter(!rownames(beta) %in% toremove) ; rm(beta) ; gc()

## test set pheno
beta <- readRDS("../../../0-preprocess/0-output/beta_3CExtVal.Rds")
pheno <- pheno %>% select(basename,type,age,ic_retrain) %>% droplevels()
write_feather(pheno,"./1-output/pheno_test.feather")

## test set beta
# read in
beta <- readRDS("../../../../0-data/beta/beta_3CExtVal.Rds")

# only keep features also in training set
beta<- beta %>% filter(rownames(beta) %in% rownames(beta_train))

# deal with missing features
features_common <- intersect(rownames(beta),rownames(beta_train))
features_missing <- setdiff(sort(rownames(beta_train)),sort(features_common))
df <- as.data.frame(matrix(ncol=length(colnames(beta)),nrow=length(features_missing)))
rownames(df) <- features_missing
colnames(df) <- colnames(beta)
for (i in 1:length(features_missing)){ 
  d <- beta_train %>% filter(rownames(beta_train) == features_missing[i]) %>% as.vector() %>% as.numeric() %>% density()
  df[i,] <- rep(max(d$x), times=length(colnames(df)))
}
identical(colnames(beta), colnames(df))
beta <- bind_rows(beta,df)

# fix order
beta <- beta[match(rownames(beta_train), rownames(beta)),]

# double check that CpGs (features) are in the same order, because rownames will be lost
identical(rownames(beta), rownames(beta_train))

#save betas
write_feather(beta_train,"./1-output/beta_train.feather")
write_feather(beta,"./1-output/beta_test.feather")
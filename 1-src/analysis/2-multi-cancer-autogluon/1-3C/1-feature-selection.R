# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Description: 
# Feature selection of those that are top 5000 db immune/epithelial for breast, ovarian and endometrial cancer separately
# Note that this is done for the training set only (internal validation set excluded)
# Unreliable probes were excluded

library(dplyr)
library(tidyverse)

# tried top 15000 and this did not work well for training models max 3600 sec

dir.create("./1-output")

## Breast
# read in top 5000 beta subsets and save CpG names
beta_tr_epi <- readRDS(file = './0-output/BC/beta_tr_epi.Rds')
beta_tr_imm <- readRDS(file = './0-output/BC/beta_tr_imm.Rds')
features_breast <- unique(c(rownames(beta_tr_epi)[1:5000],rownames(beta_tr_imm)[1:5000])) #7357 unique CpGs

## Ovarian
# read in top 5000 beta subsets and save CpG names
beta_tr_epi <- readRDS(file = './0-output/OC/beta_tr_epi.Rds')
beta_tr_imm <- readRDS(file = './0-output/OC/beta_tr_imm.Rds')
features_ovarian <- unique(c(rownames(beta_tr_epi)[1:5000],rownames(beta_tr_imm)[1:5000])) #8292 unique CpGs

## Endometrial
# read in top 5000 beta subsets and save CpG names
beta_tr_epi <- readRDS(file = './0-output/EC/beta_tr_epi.Rds')
beta_tr_imm <- readRDS(file = './0-output/EC/beta_tr_imm.Rds')
features_endometrial <- unique(c(rownames(beta_tr_epi)[1:5000],rownames(beta_tr_imm)[1:5000])) #8292 unique CpGs

features <- unique(c(features_breast,unique(c(features_ovarian,features_endometrial)))) #17863 final features
  
# double check additionally flagged probes are removed
removedEPICv2 <- read.csv("../../../../../0-data/flagged_probes/notonEPICv2_names.csv", header=FALSE) %>% pull(V1)
lowMI <- read.csv("../../../../../0-data/flagged_probes/lowMI_names.csv", header=FALSE) %>% pull(V1)
toremove <- unique(c(removedEPICv2,lowMI))
intersect(toremove,features) #all good :)

saveRDS(features, file="./1-output/selected_features_3C.Rds")

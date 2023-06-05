#!/usr/bin/env Rscript

# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(dplyr)
library(tidyverse)
library(data.table)

#dir.create("./2-output")

## training set, beta, pheno, selected features
#read in
beta <- readRDS("../../../../0-data/beta/beta_3CDisc.Rds")
pheno <- readRDS("./0-output/pheno_3CDisc_ic.Rds")
features <- readRDS("./1-output/selected_features_3C.Rds")

# split train/test
pheno1 <- pheno %>% filter(split=="train") %>% droplevels() # training set for this round
pheno2 <- pheno %>% filter(split=="test") %>% droplevels() # internal validation set for testing feature importance and subsequent pruning

#save Training set
datTrain <- beta %>% filter(rownames(beta) %in% features) %>% select(c(pheno1$basename)) %>% droplevels()
datTrain <- t(as.matrix(datTrain)) #transpose
datTrain <- as.data.frame(datTrain) %>% rownames_to_column(var="basename")
pheno1 <- pheno1 %>% select(basename,age,type,ic_retrain) %>% droplevels()
datTrain <- full_join(pheno1,datTrain) %>% drop_na() %>% column_to_rownames(var="basename")
fwrite(datTrain,file="./2-output/train3Cdisc.csv", row.names=FALSE)

## "internal" validation set for pruning features
datVal <- beta %>% filter(rownames(beta) %in% features) %>% select(c(pheno2$basename)) %>% droplevels()
datVal <- t(as.matrix(datVal)) #transpose
datVal <- as.data.frame(datVal) %>% rownames_to_column(var="basename")
pheno2 <- pheno2 %>% select(basename,age,type,ic_retrain) %>% droplevels()
datVal <- full_join(pheno2,datVal) %>% drop_na() %>% column_to_rownames(var="basename")
fwrite(datVal,file="./2-output/validate3Cdisc.csv", row.names=FALSE)

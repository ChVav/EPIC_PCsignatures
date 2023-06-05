# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Check performance with increasing training set size
# penalized regression, 10K cross-validation for each method
# full 3C discovery set used for training (no train-test split, cancers included)
# 3C external validation/smears repeatability used for testing performance
# for sample input choose, n min = 57 (random), n max = 1657, n step =20 (random) -> 81 sample sizes checked
# later added more fine-grained steps for smaller sizes n =32,37,...,52
### final output 1. dataframe:
## training parameters kept:
# method: PC.clock, ElNet, Lasso
# input.CpGs: all, unreliable.removed
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

# note that since this is computationally intensive, the final execution was performed with Rscript on a linux workstation, running steps from the loop in parallel
# This is just an example script

# Same experiment is repeated three times, results are then collated results

library(dplyr)
library(tidyverse)
library(ggplot2)

dir.create("./1-output")

### Source scripts for training+predicting age in validation sets for testing accuracy/precision ###----

source("../../functions/regression-train-validate.R")
source("../../functions/PC-train-validate.R")

### Prepare data ###----
## Training set - 3C Discovery
# Beta
beta_tr_full <- readRDS("../../../0-data/beta/beta_3CDisc.Rds")

# pheno
pheno_tr_full <- readRDS("../../../0-data/pheno/pheno_3CDisc.Rds")
pheno_tr_full <- pheno_tr_full %>% select(basename,age)

# arrange age for training 
pheno_tr_full <- pheno_tr_full %>% filter(pheno_tr_full$basename %in% colnames(beta_tr_full)) %>% droplevels()
pheno_tr_full <- pheno_tr_full[match(colnames(beta_tr_full),pheno_tr_full$basename),]
identical(colnames(beta_tr_full),pheno_tr_full$basename)

# set probes to remove
removedEPICv2 <- read.csv("../../../0-data/flagged_probes/notonEPICv2_names.csv", header=FALSE) %>% pull(V1)
lowMI <- read.csv("../../../0-data/flagged_probes/lowMI_names.csv", header=FALSE) %>% pull(V1)
toremove <- unique(c(removedEPICv2,lowMI))

## Validation set accuracy - 3C External validation
beta_val<- readRDS("../../../0-data/beta/beta_3CExtVal.Rds") #778511 CpGs, 449 samples

#pheno
pheno_val <- readRDS("./../../0-data/pheno/pheno_3CDisc.Rds")
pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(colnames(beta_val),pheno_val$basename)
age_val <- pheno_val %>% pull(age) %>% as.numeric()

## Validation set precision - Cervical smears Repeatability
beta_rep <- readRDS("../../../0-data/beta/beta_repeatability.Rds")
pheno_rep <- readRDS("../../../0-data/pheno/pheno_repeatability.Rds")
pheno_rep <- pheno %>%
  filter(sample_type=="Cervical_smear") %>% # kick out blood
  select(basename,patient_ID,rep) %>% #only keep necessary stuff
  droplevels()
beta_rep <- beta_rep %>% select(c(pheno_rep$basename))  #keep smears only
pheno_rep <- pheno_rep[match(colnames(beta_rep),pheno_rep$basename),]
identical(colnames(beta_rep),pheno_rep$basename)

### Define list of samples to loop over and save age distribution ###----
nr <- seq(from=57,to=1657,by=20)
# gradually increase sample size, for each step new set of random samples will be selected
sampleList <- list()
for (i in 1:length(nr)){
  set.seed(27+i)
  sampleList[[i]] <- sample(pheno_tr_full$basename,nr[i])
}
names(sampleList) <- nr
saveRDS(sampleList, file="./1-output/sampleListloop.Rds") # fix list to repeat

# add more fine-grained steps at lower end
nr <- seq(from=32,to=52,by=5)
sampleList2 <- list()
for (i in 1:length(nr)){
  set.seed(27+i)
  sampleList2[[i]] <- sample(pheno_tr_full$basename,nr[i])
}
names(sampleList2) <- nr
saveRDS(sampleList2, file="./1-output/sampleListloop2.Rds") # fix list to repeat

# out.list - age distribution for different input sample sets - length (nmax-nmin/n) +1:
sampleList <- c(sampleList2,sampleList)
out.list <- list()
for (j in 1:length(sampleList)){
  out.list[[j]] <- pheno_tr_full %>% filter(basename %in% sampleList[[j]]) %>% select(age,basename)
}
names(out.list) <- names(sampleList)
saveRDS(out.list, file="./1-output/sampleListloopAges.Rds")

### Train and test different methods ###----
# initialize out.df
out.df <- data.frame(train.sample.size=numeric(0),
                     method=character(0),
                     input.CpGs=character(0),
                     n.features=numeric(0),
                     val.slope.regress=numeric(0),
                     val.intercept.regress=numeric(0),
                     val.MSE=numeric(0),
                     rep.ICC=numeric(0),
                     rep.lbound.ICC=numeric(0),
                     rep.ubound.ICC=numeric(0))

# Lasso/Elastic Net
for (j in 1:length(sampleList)){
  
  # take sample subset
  beta_tr <- beta_tr_full %>% select(all_of(sampleList[[j]]))
  pheno_tr <- pheno_tr_full %>% filter(basename %in% all_of(sampleList[[j]]))
  pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
  identical(colnames(beta_tr),pheno_tr$basename)
  age_tr <- pheno_tr %>% pull(age) %>% as.numeric()
  
  # Lasso - All
  result <- regressionTrainvalidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,1.0)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "Lasso"
  result$input.CpGs <- "all"
  out.df <- full_join(out.df,result)
  
  # Elastic Net -All
  result <- regressionTrainvalidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,0.5)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "ElNet"
  result$input.CpGs <- "all"
  out.df <- full_join(out.df,result)
  
  # PC - All
  result <- PCTrainValidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "PC"
  result$input.CpGs <- "all"
  out.df <- full_join(out.df,result)
  
  # take CpG subset from sample subset
  beta_tr <- beta_tr %>% filter(!rownames(beta_tr) %in% toremove) #634898 CpGs
  identical(colnames(beta_tr),pheno_tr$basename)
  
  # Lasso - unreliable probes removed
  result <- regressionTrainvalidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,1.0)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "Lasso"
  result$input.CpGs <- "unreliable.removed"
  out.df <- full_join(out.df,result)
  
  # Elastic Net - Unreliable probes removed
  result <- regressionTrainvalidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep,0.5)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "ElNet"
  result$input.CpGs <- "unreliable.removed"
  out.df <- full_join(out.df,result)
  
  # PC - unreliable probes removed
  result <- PCTrainValidate(beta_tr,age_tr,beta_val,age_val,beta_rep,pheno_rep)
  result$train.sample.size <- names(sampleList[j]) %>% as.numeric()
  result$method <- "PC"
  result$input.CpGs <- "unreliable.removed"
  out.df <- full_join(out.df,result)
  
  saveRDS(out.df, file="./1-output/results_training_size.Rds") # overwrite results at end of each loop # save temporary results for lasso/elnet
}

### Check output one experiment ##----

# accuracy/precision/features

results_training_size %>%
  ggplot(aes(x=train.sample.size,y=val.slope.regress,color=input.CpGs)) +
  geom_point() +
  geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95) +
  facet_wrap(~method)

results_training_size %>%
  ggplot(aes(x=train.sample.size,y=val.MSE,color=input.CpGs)) +
  geom_point() +
  geom_smooth(method="loess", se=FALSE) +
  facet_wrap(~method)

results_training_size %>%
  ggplot(aes(x=train.sample.size,y=val.MSE,color=method)) +
  geom_point() +
  geom_smooth(method="loess", se=FALSE) +
  facet_wrap(~input.CpGs)

results_training_size %>%
  ggplot(aes(x=train.sample.size,y=rep.ICC,color=input.CpGs)) +
  geom_point(stat="identity") +
  geom_smooth(method="loess", se=FALSE) +
  #geom_errorbar(aes(ymin=rep.lbound.ICC, ymax=rep.ubound.ICC), width=.2)+
  facet_wrap(~method)

results_training_size %>%
  ggplot(aes(x=train.sample.size,y=n.features,color=input.CpGs)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~method)

# plot age density and mean for a given subset
mu <- mean(out.list[[81]]$age)
ggplot(data=out.list[[81]], aes(x=age)) +
  geom_density() +
  geom_vline(xintercept=mu,linetype="dashed") +
  xlab(paste0("age in training set with ",names(out.list[81])," samples")) +
  theme_minimal()

mu <- mean(out.list[[10]]$age)
ggplot(data=out.list[[10]], aes(x=age)) +
  geom_density() +
  geom_vline(xintercept=mu,linetype="dashed") +
  xlab(paste0("age in training set with ",names(out.list[10])," samples")) +
  theme_minimal()



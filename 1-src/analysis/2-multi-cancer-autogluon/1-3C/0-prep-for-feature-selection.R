# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Save CpGs with top db for each cancer

library(dplyr)
library(tidyverse)

# source helper functions
source("../../../functions/estimateDeltaBeta.R")
source("../../../functions/saveBetaSubsets.R") 

# create output folder
dir.create("./0-output")

## Calculate IC after removing unreliable probes ##----

beta <- readRDS("../../../../0-data/beta/beta_3CDisc.Rds")
pheno <- readRDS(file = "../../../../0-data/pheno/pheno_3CDisc.Rds")

# remove flagged probes
removedEPICv2 <- read.csv("../../../../0-data/flagged_probes/notonEPICv2_names.csv", header=FALSE) %>% pull(V1)
lowMI <- read.csv("../../../../0-data/flagged_probes/lowMI_names.csv", header=FALSE) %>% pull(V1)
toremove <- unique(c(removedEPICv2,lowMI))

beta <- beta %>% as.data.frame() %>% filter(!rownames(beta) %in% probes_to_remove) %>% as.matrix()

# calculate IC proportion
out.l <- epidish(beta.m = beta, ref.m = centEpiFibIC.m, method = "RPC")

# save pheno with recalculated ic
ic <- out.l$estF %>% 
  as.data.frame() %>%
  select(IC) %>%
  rownames_to_column(var="basename")

pheno <- full_join(pheno,ic) %>%
  drop_na() %>%
  dplyr::rename(ic_retrain=IC)

saveRDS(pheno, file='./0-output/pheno_3CDisc_ic.Rds') 

## Endometrial db ##----
dir.create("./0-output/EC")

path_to_output <- "./0-output/EC"

# define beta and pheno, match samples
pheno_tr <- pheno %>% 
  filter(type %in% c("Control","Endometrial")) %>% 
  filter(split == "train") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()
pheno_val <- pheno %>% 
  filter(type %in% c("Control","Endometrial")) %>% 
  filter(split == "test") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()

beta_tr <- beta %>% select(c(pheno_tr$basename)) %>% droplevels()
beta_val <- beta %>% select(c(pheno_val$basename)) %>% droplevels()

identical(rownames(beta_tr),rownames(beta_val))

# arrange age and ic for training and validation
pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_tr),pheno_tr$basename)

pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_val),pheno_val$basename)

# delta-beta estimates
db <- estimateDeltaBeta(beta_tr,pheno_tr,toremove)
saveRDS(db, file=paste0(path_to_output,"/delta-beta.Rds"))

## Breast ##----
dir.create("./0-output/BC")

path_to_output <- "./0-output/BC"

# define beta and pheno, match samples
pheno_tr <- pheno %>% 
  filter(type %in% c("Control","Breast")) %>% 
  filter(split == "train") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()
pheno_val <- pheno %>% 
  filter(type %in% c("Control","Breast")) %>% 
  filter(split == "test") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()

beta_tr <- beta %>% select(c(pheno_tr$basename)) %>% droplevels()
beta_val <- beta %>% select(c(pheno_val$basename)) %>% droplevels()

identical(rownames(beta_tr),rownames(beta_val))

# arrange age and ic for training and validation
pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_tr),pheno_tr$basename)

pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_val),pheno_val$basename)

# delta-beta estimates
db <- estimateDeltaBeta(beta_tr,pheno_tr,toremove)
saveRDS(db, file=paste0(path_to_output,"/delta-beta.Rds"))


## Ovarian ##----
dir.create("./0-output/OC")

path_to_output <- "./0-output/OC"

# define beta and pheno, match samples
pheno_tr <- pheno %>% 
  filter(type %in% c("Control","Ovarian")) %>% 
  filter(split == "train") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()
pheno_val <- pheno %>% 
  filter(type %in% c("Control","Ovarian")) %>% 
  filter(split == "test") %>% 
  column_to_rownames(var="basename") %>%
  select(type,age,ic_retrain) %>%
  droplevels()

beta_tr <- beta %>% select(c(pheno_tr$basename)) %>% droplevels()
beta_val <- beta %>% select(c(pheno_val$basename)) %>% droplevels()

identical(rownames(beta_tr),rownames(beta_val))

# arrange age and ic for training and validation
pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_tr),pheno_tr$basename)

pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(rownames(beta_tr),rownames(beta_val))
identical(colnames(beta_val),pheno_val$basename)

# delta-beta estimates
db <- estimateDeltaBeta(beta_tr,pheno_tr,toremove)
saveRDS(db, file=paste0(path_to_output,"/delta-beta.Rds"))
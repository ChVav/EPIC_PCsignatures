
### Retrain hannum PC clock on BloodHannum_450K (GSE40279)
# test on BloodFull_450K/BloodRep_450K (GSE55763)
# 1. prep data
# 2. subset common CpGs shared in training and test sets, train
# 3. calcluate PC clock predictions for BloodFull_450K
# 4. calculate PC clock predictions for BloodRep_450K (control, this time using more CpGs)

#### 1. prep training set, subset common CpGs shared in training and test sets, train ----

library(tidyverse)
library(glmnet)

dir.create("./3-output")

### 1. Prepare data ###----
## Beta
beta_tr <- readRDS("../../../0-data/beta/beta_Hannum.Rds") # Training set
#beta_tr <- readRDS("~/Desktop/Datasets-downloaded/GSE40279/beta_Hannum.Rds")
beta_val <- readRDS("../../../0-data/beta/beta_BloodFull_450K_imputed.Rds") # Validation set accuracy - full GSE55763, missing data imputed
#beta_val <- readRDS("~/Desktop/Datasets-downloaded/GSE55763/beta_BloodFull_450K_imputed.Rds")
beta_rep <- readRDS("../../../0-data/beta/beta_BloodRep_450K_imputed.Rds") # validation set precision - technical reps GSE55763, missing data imputed
#beta_rep <- readRDS("~/Desktop/Datasets-downloaded/GSE55763/beta_BloodRep_450K_imputed.Rds")

## arrange pheno files
#training
pheno_tr <- readRDS("../../../0-data/pheno/pheno_Hannum.Rds") %>% 
  dplyr::rename(age=Age,basename=sampleID)
pheno_tr <- pheno_tr %>% 
  select(basename,age) %>% 
  filter(pheno_tr$basename %in% colnames(beta_tr)) %>% 
  droplevels()
pheno_tr <- pheno_tr[match(colnames(beta_tr),pheno_tr$basename),]
identical(colnames(beta_tr),pheno_tr$basename)
pheno_tr <- pheno_tr %>% pull(age) %>% as.numeric()

# validation accuracy
pheno_val <- readRDS("./1-output/pheno_BloodFull_450K_Hannum.Rds") 
pheno_val <- pheno_val[match(colnames(beta_val),pheno_val$basename),]
identical(colnames(beta_val),pheno_val$basename)
age_val <- pheno_val %>% pull(age) %>% as.numeric()

# validation precision
pheno_rep <- readRDS("./1-output/pheno_BloodRep_450K_Hannum.Rds")
pheno_rep <- pheno_rep[match(colnames(beta_rep),pheno_rep$basename),]
identical(colnames(beta_rep),pheno_rep$basename)

## Stick to approach Higgins to subset all CpGs from train/test/sets
CpGs = intersect(rownames(beta_tr),intersect(rownames(beta_val),rownames(beta_rep))) #473,034 CpGs
saveRDS(CpGs, file="./3-output/commonCpgs.Rds")
beta_tr = beta_tr %>% filter(rownames(beta_tr) %in% CpGs)
beta_val = beta_val %>% filter(rownames(beta_val) %in% CpGs)
beta_val = beta_val[match(CpGs, rownames(beta_val)),]
beta_rep = beta_rep %>% filter(rownames(beta_rep) %in% CpGs)
beta_rep = beta_rep[match(CpGs, rownames(beta_rep)),]

# check order CpGs in all datasets
if (identical(rownames(beta_tr),rownames(beta_val)) & identical(rownames(beta_tr),rownames(beta_rep))){
  ## Perform PCA and projections. Remove last PC.
  PCA = prcomp(t(as.matrix(beta_tr)),scale.=F) #data centered, but not scaled
  TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
  
  ## train with 10-fold cross-validation (default settings Lambda)
  TrainAge = as.numeric(pheno_tr)
  
  #Train PC clock, elastic net
  nFolds = 10
  foldid = sample(rep(seq(nFolds), length.out = nrow(TrainPCData)))
  cv = cv.glmnet(TrainPCData, TrainAge, nfolds=nFolds, foldid=foldid, alpha=0.5, family="gaussian") 
  
  CalcPCAge = vector(mode = "list",length = 0)
  temp = as.matrix(coef(cv,s = cv$lambda.min))
  CalcPCAge$model = temp[temp!=0,][-1]
  CalcPCAge$intercept = temp[1,1]
  CalcPCAge$center = PCA$center
  CalcPCAge$rotation = PCA$rotation[,names(CalcPCAge$model)]
  
  ## predict age in validation sets
  predAge <- function(beta,CalcPCAge,CpGs){ # helper function for predicting age using trained model
    # calculate projections
    PCAge = sweep(t(as.matrix(beta)),2,CalcPCAge$center) %*% CalcPCAge$rotation %*% CalcPCAge$model + CalcPCAge$intercept
    return(PCAge)
  }
  
  ageVal = predAge(beta_val,CalcPCAge,CpGs) # validation set 1 ~ accuracy
  ageRep = predAge(beta_rep,CalcPCAge,CpGs) # validation set 2 ~ precision
  
  pheno_val$Hannum.PC2 <- ageVal[,1]
  pheno_rep$Hannum.PC2 <- ageRep[,1]
  
  saveRDS(pheno_val, "./3-output/pheno_BloodFull_450K_Hannum_PC.Rds")
  saveRDS(pheno_rep, "./3-output/pheno_BloodRep_450K_Hannum_PC.Rds")
  
}


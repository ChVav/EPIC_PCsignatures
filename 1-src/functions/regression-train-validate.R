# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Train Lasso regression + predict age in validation sets + calculate precision/accuracy

regressionTrainvalidate <- function(beta_tr,pheno_tr,beta_val,pheno_val,beta_rep,pheno_rep,alpha) {
  require(dplyr)
  require(tidyverse)
  require(glmnet)
  require(irr)
  
  # set alpha=0.0 for Ridge, 0.5 for Elastic net, 1.0 for Lasso
  
  ## train with 10-fold cross-validation (default settings Lambda)
  nFolds = 10
  foldid = sample(rep(seq(nFolds), length.out = ncol(beta_tr))) #beta_tr or nrow for TrainPCData
  cv = cv.glmnet(t(beta_tr), as.numeric(pheno_tr), nfolds=nFolds, foldid=foldid, alpha=alpha, family="gaussian",type.measure="mse")
  index_coef = coef(cv,s = cv$lambda.min)
  rnames = rownames(index_coef)
  index_coef = as.numeric(index_coef)
  names(index_coef) = rnames
  index_coef = index_coef[index_coef!=0] # keep only non-zero coefficients
  intercept = index_coef[1]
  w = index_coef[2:length(index_coef)]
  
  ## predict age in validation sets
  predAge <- function(beta,w,intercept){ # helper function for predicting age using trained model
    ind = match(names(w),rownames(beta))
    B = beta[ind[!is.na(ind)],]
    w2 = w[!is.na(ind)]
    B1 = B*w2
    age = intercept + apply(B1, MARGIN = 2, FUN = 'sum')
    return(age)
  }
  
  ageVal = predAge(beta_val,w,intercept) # validation set 1 ~ accuracy
  ageRep = predAge(beta_rep,w,intercept) # validation set 2 ~ precision
  
  ## regression validation accuracy
  modelVal <- lm(ageVal ~ pheno_val)
  
  ## ICC repeatablity set precision
  dat = pheno_rep
  dat$predAge = as.numeric(ageRep)
  dat = dat %>% pivot_wider(id_cols=c("patient_ID"),names_from = "rep",values_from = "predAge") %>% select(-patient_ID)
  out = icc(dat, model="twoway",type="agreement")
  
  ### final output 
  df <- data.frame(n.features=length(w),
                   val.slope.regress=summary(modelVal)$coefficients[2],
                   val.intercept.regress=summary(modelVal)$coefficients[1],
                   val.MSE=mean((ageVal-pheno_val)^2),
                   rep.ICC=out$value,
                   rep.lbound.ICC=out$lbound,
                   rep.ubound.ICC=out$ubound
                   )
  
  return(df)
}



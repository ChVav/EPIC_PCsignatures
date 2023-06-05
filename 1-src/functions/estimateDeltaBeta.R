# This code is available for research use only. Commercial use of the code or any data related to it is prohibited.
# author: James E. Barrett, edited Charlotte Vavourakis

# Description: infer epithelial and immune delta-beta estimates

estimateDeltaBeta <- function(beta, pheno, probes_to_remove) {
  require(dplyr)
  require(tidyverse)
  
  cat('Beginning delta-beta estimation script...\n\n')
  
  cat('Removing unreliable probes\n')
  
  beta <- beta %>% as.data.frame() %>% filter(!rownames(beta) %in% probes_to_remove) %>% as.matrix()
  ind <- match(colnames(beta), rownames(pheno))
  pheno <- pheno[ind,]
  identical(rownames(pheno),colnames(beta))
  
  cat('Statistical adjustment will be made for:\n')
  for(i in colnames(pheno)){cat(i,'\n')}
  cat('\n')
  
  # Set factor levels and data types correctly
  pheno$type <- as.factor(as.character(pheno$type))
  pheno$ic_retrain <- as.numeric(pheno$ic_retrain)
  pheno$age <- as.numeric(pheno$age)
  
  #==============================================================================#
  # fit linear models
  
  cat('Estimating delta-betas...\n')
  db <- matrix(NA, nrow=nrow(beta),ncol=2+ncol(pheno))
  rownames(db) <- rownames(beta)
  colnames(db) <- c('db_epithelial','db_immune',colnames(pheno))
  
  ic.control <- pheno$ic_retrain[which(pheno$type=='Control')]
  ic.case <- pheno$ic_retrain[which(pheno$type!='Control')]
  
  pB <- txtProgressBar(min=1,max=nrow(beta), width =50L, style = 3)
  for (i in 1:nrow(beta)){
    setTxtProgressBar(pB, i)
    
    ldat <- data.frame(beta=as.numeric(beta[i,]),
                       pheno)
    
    lfit <- lm(beta ~ ., data=ldat)
    db[i,3:ncol(db)] <- summary(lfit)$coefficients[2:(ncol(pheno)+1),4]
    
    beta.control <- as.numeric(beta[i,which(pheno$type=='Control')])
    beta.case <- as.numeric(beta[i,which(pheno$type!='Control')])
    
    dat.control <- data.frame(beta=beta.control, ic=ic.control)
    dat.case <- data.frame(beta=beta.case, ic=ic.case)
    
    fit.control <- lm(beta ~ ic, data=dat.control)
    fit.case <- lm(beta ~ ic, data=dat.case)
    
    delta_beta_epithelial <- round(fit.case$coefficients[1]-fit.control$coefficients[1],digits=3)
    delta_beta_immune <- round(-fit.control$coefficients[1] - fit.control$coefficients[2]
                               +fit.case$coefficients[1] + fit.case$coefficients[2],digits=3)
    db[i,1] <- delta_beta_epithelial
    db[i,2] <- delta_beta_immune
  }
  close(pB)
  cat('done\n\n')
  return(db)
}
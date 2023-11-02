# edited code copied from ../0-data/PaperCode_PC-Clocks-main/PC-Clocks-main/TrainPCClocks.R
# training either on age 

PCTrainValidate <- function(beta_tr, pheno_tr,beta_val,pheno_val,beta_rep,pheno_rep) {
  
  require(tidyverse)
  require(glmnet)
  require(irr)
  
  ## Stick to approach Higgins to subset all CpGs from train/test/sets
  CpGs = intersect(rownames(beta_tr),intersect(rownames(beta_val),rownames(beta_rep)))
  beta_tr = beta_tr %>% filter(rownames(beta_tr) %in% CpGs)
  beta_val = beta_val %>% filter(rownames(beta_val) %in% CpGs)
  beta_val = beta_val[match(CpGs, rownames(beta_val)),]
  beta_rep = beta_rep %>% filter(rownames(beta_rep) %in% CpGs)
  
  # check order CpGs in all datasets
  if (identical(rownames(beta_tr),rownames(beta_val)) & identical(rownames(beta_tr),rownames(beta_rep))){
    ## Perform PCA and projections. Remove last PC.
    PCA = prcomp(t(as.matrix(beta_tr)),scale.=F) #data centered, but not scaled
    TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
    
    ## train with 10-fold cross-validation (default settings Lambda)
    TrainAge = as.numeric(pheno_tr)
    
    #Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
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
    
    ## regression validation accuracy
    modelVal <- lm(ageVal ~ pheno_val)
    
    ## ICC repeatablity set precision
    dat = pheno_rep
    dat$predAge = as.numeric(ageRep)
    dat = dat %>% pivot_wider(id_cols=c("patient_ID"),names_from = "rep",values_from = "predAge") %>% select(-patient_ID)
    out = icc(dat, model="twoway",type="agreement")
    
    ### final output 
    df <- data.frame(n.features=length(CalcPCAge$model),
                     val.slope.regress=summary(modelVal)$coefficients[2],
                     val.intercept.regress=summary(modelVal)$coefficients[1],
                     val.MSE=mean((ageVal-pheno_val)^2),
                     rep.ICC=out$value,
                     rep.lbound.ICC=out$lbound,
                     rep.ubound.ICC=out$ubound
    )
  } else {
    stop("CpGs in train and validation sets are not in the correct order")}

  return(df)
}

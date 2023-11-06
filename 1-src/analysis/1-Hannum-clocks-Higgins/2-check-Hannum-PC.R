
### Train hannum PC and PC clocks on BloodHannum_450K (GSE40279), test on BloodRep_450K (GSE55763) and compare to results calculation in paper Higgins-Chen, Nature Aging 2022
# 1. prep training set, subset only 78,464 CpGs that were shared in all training sets used for the paper - GSE40279, calculate Hannum
# 2. prep test set (transpose), calculate Hannum
# 3. (impute +) calculate PC clock using actual age
# 4. (impute +) train Hannum PC proxy, elastic net regression should select 390 of 655 Hannum PCs? (last PC ignored)
# code copied from ../0-data/PaperCode_PC-Clocks-main/PC-Clocks-main/TrainPCClocks.R
# 5. calculate PC clock and Hannum PC clock proxy for blood replication test set ("project test set on PCA space")

#### 1. prep training set ----

library(tidyverse)
library(openxlsx)
library(glmnet)
library(methylclock)

# get data from GEO
beta <- readRDS("../../../0-data/beta/beta_Hannum.Rds")
pheno <- readRDS("../../../0-data/pheno/pheno_Hannum.Rds")

# beta (datMeth), rownames should be sample names, colnames CpGs, so will have to transpose
beta <- beta %>% column_to_rownames(var="ID_REF")# move CpGs to rownames
datMeth <- t(as.matrix(beta)) #656 samples, 473034 CpGs

# check order samples
datMeth <- datMeth[match(pheno$sampleID,rownames(datMeth)),]
identical(pheno$sampleID, rownames(datMeth)) 

# subset 78464 CpGs used in paper, prep pheno, check order and names to proceed with step 3.
# download 43587_2022_248_MOESM3_ESM.xlsx from Supplementary Table 1
subset <- read.xlsx("<path-to-download>/43587_2022_248_MOESM3_ESM.xlsx", sheet=8)
subset <- subset[,1]
head(subset)
n <- length(subset)
subset <- subset[2:n]

datMethTrain <- as.data.frame(datMeth)
datMethTrain <- datMethTrain %>% select(all_of(subset))
rm(datMeth);gc()

datPhenoTrain <- pheno %>% select(Age)
rownames(datPhenoTrain) <- pheno$sampleID
rm(trainset_pheno);gc()

# calculate original Hannum clock
predicted.age <- DNAmAge(beta, clocks="Hannum")
predicted.age <- predicted.age[match(rownames(datPhenoTrain),predicted.age$id),]
identical(predicted.age$id,rownames(datPhenoTrain))
datPhenoTrain$Hannum <- predicted.age$Hannum
rm(beta);gc()

identical(rownames(datPhenoTrain),rownames(datMethTrain))

#### 2. prep test set ----

#beta, transpose, subset 78464 CpGs used in paper, prep pheno, check order and names to proceed with step 3.
beta <- readRDS("../../../0-data/beta/beta_BloodRep_450K.Rds")
datMethTest <- t(as.matrix(beta)) #72 samples,473864 CpGs
datMethTest <- as.data.frame(datMethTest)
datMethTest <- datMethTest %>% select(all_of(subset))

#pheno
datPhenoTest <- readRDS("../../../0-data/pheno/pheno_BloodRep_450K.Rds")
rownames(datPhenoTest) <- datPhenoTest$basename
datPhenoTest <- datPhenoTest %>% select(age) %>% dplyr::rename(Age=age)

# calculate original Hannum clock
beta <- beta %>% rownames_to_column(var="CpG")
predicted.age <- DNAmAge(beta, clocks="Hannum")
predicted.age <- predicted.age[match(rownames(datPhenoTest),predicted.age$id),]
identical(predicted.age$id,rownames(datPhenoTest))
datPhenoTest$Hannum <- predicted.age$Hannum
rm(beta);gc()

identical(rownames(datMethTest),rownames(datPhenoTest))

#### 3. Train PC clock (predicting age) ----

CpGs <- subset #78464 CpGs already subsetted

identical(colnames(datMethTrain),colnames(datMethTest))

#Impute missing values 
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain <- apply(datMethTrain,2,meanimpute)

#check order of data is correct
if(all(colnames(datMethTrain) == CpGs)){
  message("CpGs are all in order")
} else(message(paste("Only",sum(colnames(datMethTrain) == CpGs),"CpGs are in order")))

if(all(rownames(datMethTrain) == rownames(datPhenoTrain))){ #may need to change rownames(datPhenoTrain) to a column name
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(datMethTrain) == rownames(datPhenoTrain)),"Samples are in order")))

#Perform PCA and projections. Remove last PC.
PCA = prcomp(datMethTrain,scale.=F) #data centered, but not scaled
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TestPCData = predict(PCA,datMethTest)[,1:(dim(PCA$x)[2]-1)]

#Select phenotype to be predicted. For example, here we are predicting age.
TrainAge = as.numeric(datPhenoTrain$Age)
TestAge = as.numeric(datPhenoTest$Age)

#Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainAge, nfolds=10,alpha=0.5, family="gaussian") 
fit = glmnet(TrainPCData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)

#Examine full model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.min),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.min))

#Examine sparse model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.1se),xlab = "Age",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.1se))

#Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCAge <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCAge$model = temp[temp!=0,][-1]
CalcPCAge$intercept = temp[1,1]
CalcPCAge$center = PCA$center
CalcPCAge$rotation = PCA$rotation[,names(CalcPCAge$model)]

length(CalcPCAge$model)

#### 4. Train Hannum PC clock proxy (predicting Hannum) ----

CpGs <- subset #78464 CpGs already subsetted

identical(colnames(datMethTrain),colnames(datMethTest))

#Impute missing values .
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain <- apply(datMethTrain,2,meanimpute)

#check order of data is correct
if(all(colnames(datMethTrain) == CpGs)){
  message("CpGs are all in order")
} else(message(paste("Only",sum(colnames(datMethTrain) == CpGs),"CpGs are in order")))

if(all(rownames(datMethTrain) == rownames(datPhenoTrain))){ #may need to change rownames(datPhenoTrain) to a column name
  message("Samples are all in order")
} else(message(paste("Only",sum(rownames(datMethTrain) == rownames(datPhenoTrain)),"Samples are in order")))

#Perform PCA and projections. Remove last PC.
PCA = prcomp(datMethTrain,scale.=F) #data centered, but not scaled
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TestPCData = predict(PCA,datMethTest)[,1:(dim(PCA$x)[2]-1)]

#Select phenotype to be predicted. For example, here we are predicting Hannum Age
TrainAge = as.numeric(datPhenoTrain$Hannum)
TestAge = as.numeric(datPhenoTest$Hannum)

#Train PC clock. Can test different models using different alpha and lambda parameters (see glmnet documentation)
cv = cv.glmnet(TrainPCData, TrainAge, nfolds=10,alpha=0.5, family="gaussian") 
fit = glmnet(TrainPCData, TrainAge, family="gaussian", alpha=0.5, nlambda=100)
plot(cv)

#Examine full model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min),xlab = "Hannum",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.min))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.min),xlab = "Hannum",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.min))

#Examine sparse model
plot(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se),xlab = "Hannum",ylab = "Predicted Age", main = "Training")
cor(TrainAge,predict(fit,TrainPCData,s = cv$lambda.1se))
plot(TestAge,predict(fit,TestPCData,s = cv$lambda.1se),xlab = "Hannum",ylab = "Predicted Age", main = "Testing")
cor(TestAge,predict(fit,TestPCData,s = cv$lambda.1se))

#Most likely your final model will only use a small subset of PCs. Thus you can compress your model:
CalcPCHannum <- vector(mode = "list",length = 0)
temp = as.matrix(coef(cv,s = cv$lambda.min))
CalcPCHannum$model = temp[temp!=0,][-1]
CalcPCHannum$intercept = temp[1,1]
CalcPCHannum$center = PCA$center
CalcPCHannum$rotation = PCA$rotation[,names(CalcPCHannum$model)]

length(CalcPCHannum$model)
# I found 369 PCs, while the paper described "elastic net regression should select 390 of 655 Hannum PCs (655 out of 666, last PC ignored)
# every time you do k-fold cross validation (cv.glmnet), folds or picked randomly, so results will differ
# as such impossible to replicate figure 3c exactly

#### 5. calculate PC Hannum PC proxy clocks for blood replication test set ----

PCAge <- sweep(as.matrix(datMethTest),2,CalcPCAge$center) %*% CalcPCAge$rotation %*% CalcPCAge$model + CalcPCAge$intercept
PCAge <- data.frame(basename = rownames(PCAge), PCAge.EUTOPS = c(PCAge[,1]))

# PC proxy
PCHannum <- sweep(as.matrix(datMethTest),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept
PCHannum <- data.frame(basename = rownames(PCHannum), PCHannum.EUTOPS = c(PCHannum[,1]))

# add to pheno
pheno <- readRDS("./1-output/pheno_BloodRep_450K_Hannum.Rds")
pheno <- full_join(pheno,PCAge)
pheno <- full_join(pheno,PCHannum)

#rename
pheno <- pheno %>% dplyr::rename(Hannum_PC.EUTOPS=PCAge.EUTOPS, Hannum_PCproxy.EUTOPS=PCHannum.EUTOPS)

# 5. compare to data paper, they only give the proxy results
# download 43587_2022_248_MOESM5_ESM.xlsx Source Data Fig.3
data_higgins <- read.xlsx("<path-to-download>/43587_2022_248_MOESM5_ESM.xlsx", sheet=1)
data_higgins <- data_higgins %>% 
  filter(Clock== "PCHannum") %>% gather(key="rep", value="PCHannum", "1":"2") %>% 
  select(-Clock) %>% arrange(sample,rep) %>% droplevels()

identical(pheno$sample,data_higgins$sample)
identical(pheno$rep,data_higgins$rep)

pheno$Hannum_PCproxy.higgins <- data_higgins$PCHannum

identical(pheno$Hannum_PCproxy.EUTOPS, pheno$Hannum_PCproxy.higgins) #nope does not correspond, but almost impossible because of cv.glmnet anyways

saveRDS(pheno, "./2-output/pheno_BloodRep_450K_Hannum_PC.Rds")



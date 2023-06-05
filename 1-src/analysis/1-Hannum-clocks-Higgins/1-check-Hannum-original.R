
### calculate hannum original clock for BloodRep_450K and compare to calculation in paper Higgins-Chen, Nature Aging 2022
#Note, the calculation for Horvath differed from what was published, 
#and we failed to reproduce also the results for Horvath using the online calculator 
#https://dnamage.genetics.ucla.edu/ (various settings tested)

library(tidyverse)
library(dplyr)
library(methylclock)
library(openxlsx)

## load pheno/beta blood_rep set
pheno <- readRDS("../../../0-data/pheno/pheno_BloodRep_450K.Rds")
beta <- readRDS("../../../0-data/beta/beta_BloodRep_450K.Rds")

## calculate clocks and add to pheno
#cpgs.missing <- checkClocks(beta) # test for missing CpGs (80% of CpGs in clock should be there)
beta <- beta %>% rownames_to_column(var="ProbeID")
Hannum.EUTOPS <- DNAmAge(beta, clocks="Hannum") # missing data is imputed, normalization = set to FALSE by default (assumed data was already normalized)
identical(Hannum.EUTOPS$id,pheno$basename) #sanity-check, should be TRUE
pheno$Hannum.EUTOPS <- Hannum.EUTOPS$Hannum

## compare with data presented in paper, download 43587_2022_248_MOESM4_ESM.xlsx from Source Data Fig. 1
data_higgins <- read.xlsx("<path-to-download-folder>/43587_2022_248_MOESM4_ESM.xlsx", sheet=3)

# make long format and add basename to data_paper
data_higgins <- data_higgins %>% filter(Clock== "Hannum") %>% gather(key="rep", value="Hannum", "1":"2") %>% select(-Clock) %>% arrange(sample,rep) %>% droplevels()

# basename not given, sort and check whether age and gender corresponds to pheno download from GEO, then add Hannum ages to pheno
data_higgins <- data_higgins  %>% arrange(sample,rep) 
pheno <- pheno %>% arrange(sample,rep)

identical(pheno$sample,data_higgins$sample)
identical(pheno$age,data_higgins$age)
identical(pheno$gender,data_higgins$gender)

pheno$Hannum.higgins <- data_higgins$Hannum

# compare
identical(round(pheno$Hannum.EUTOPS, digits=4), round(pheno$Hannum.higgins,digits=4)) #TRUE

# update pheno
saveRDS(pheno, file="./1-output/pheno_BloodRep_450K_Hannum.Rds")



















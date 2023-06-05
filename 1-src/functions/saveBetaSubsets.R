# This code is available for research use only. Commercial use of the code or any data related to it is prohibited.
# author: James E. Barrett, edited Charlotte Vavourakis

# Description: script to save the top 30,000 CpGs ranked by epithelial, immune delta-betas or p-values.
# The smaller beta matrices require less memory to train classifiers in parallel.

saveBetaSubsets <- function(beta_tr,beta_val,db,path_to_output) {
 
  cat('Beginning beta subset script...\n\n')
  
  ord <- order(abs(db[,1]),decreasing = TRUE)[1:30000]
  cpgs <- rownames(db[ord,])
  beta_tr_epi <- beta_tr[match(cpgs,rownames(beta_tr)),]
  beta_val_epi <- beta_val[match(cpgs,rownames(beta_val)),]
  saveRDS(beta_tr_epi, file=paste0(path_to_output,"/beta_tr_epi.Rds"))
  saveRDS(beta_val_epi, file=paste0(path_to_output,"/beta_val_epi.Rds"))
  
  ord <- order(abs(db[,2]),decreasing = TRUE)[1:30000]
  cpgs <- rownames(db[ord,])
  beta_tr_imm <- beta_tr[match(cpgs,rownames(beta_tr)),]
  beta_val_imm <- beta_val[match(cpgs,rownames(beta_val)),]
  saveRDS(beta_tr_imm, file=paste0(path_to_output,"/beta_tr_imm.Rds"))
  saveRDS(beta_val_imm, file=paste0(path_to_output,"/beta_val_imm.Rds"))
  
  ord <- order(db[,3],decreasing = FALSE)[1:30000]
  cpgs <- rownames(db[ord,])
  beta_tr_pv <- beta_tr[match(cpgs,rownames(beta_tr)),]
  beta_val_pv <- beta_val[match(cpgs,rownames(beta_val)),]
  saveRDS(beta_tr_pv, file=paste0(path_to_output,"/beta_tr_pv.Rds"))
  saveRDS(beta_val_pv, file=paste0(path_to_output,"/beta_val_pv.Rds"))
  
  cat(' done\n\n')
}
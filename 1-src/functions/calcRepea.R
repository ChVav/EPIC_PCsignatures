# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

#helper function for plotting coefficient of variation

calcRepea <- function(data,indices){
  
  require(irr) # this package is designed for interrater variablity?
  require(dplyr)
  require(tidyverse)
  
  ## Sample specific estimates
  #reformat data
  
  N = length(unique(data$rep))
  
  comp_long <- data %>%
    select(c(sample_ID, rep, sample_type, all_of(indices))) %>%
    pivot_longer(all_of(indices),
                 names_to = "index",
                 values_to = "value")
  
  comp_wide <- comp_long %>%
    pivot_wider(id_cols = c("sample_ID", "index", "sample_type"),
                names_from = "rep",
                values_from = "value")
  
  # calculate mean, SD, SEM and dispersion for each sample
  df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("sample_ID", "index", "Mean.index", "SD", "median","MAD")) #init
  
  for (i in unique(comp_long$sample_ID)) {
    for (j in unique(comp_long$index)) {
      values <- comp_long %>% filter(sample_ID == i, index == j) %>% pull()
      mean <- values %>% mean() 
      SD <- values %>% sd()
      median <- values %>% median() 
      MAD <- values %>% mad() 
      df[nrow(df) + 1,] = c(i,j,mean,SD,median,MAD)
    }
  }
  
  df[,3:6] <- sapply(df[,3:6],as.numeric)
  
  df$SEM <- as.numeric(df$SD/sqrt(N))
  df$disp <- as.numeric(df$MAD/abs(df$median))
  
  ## Sample type specific estimates
  # calculate mean(SD)
  # calculate icc ~ Higgins et al. single rater, absolute-agreements, two way random-effects
  df2 <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("index","sample_type", "ICC","CI0.95_lower","CI0.95_higher")) #init
  
  ind <- unique(comp_long$index)
  sample_types <- unique(comp_long$sample_type)
  
  for (i in ind) {
    for (j in sample_types){
      dat <- comp_wide %>% filter(index == i & sample_type==j) %>% select(-c(sample_ID,index,sample_type))
      out <- icc(dat, model="twoway",type="agreement")
      df2[nrow(df2) + 1,] = c(i,j,out$value,out$lbound,out$ubound)
    }
  }
  
  df2$ICC <- as.numeric(df2$ICC)
  df2$CI0.95_lower <- as.numeric(df2$CI0.95_lower)
  df2$CI0.95_higher <- as.numeric(df2$CI0.95_higher)
  
  list <- list(df,df2)
  
  return(list)
  
}
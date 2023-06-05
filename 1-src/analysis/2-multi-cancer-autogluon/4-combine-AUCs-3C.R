# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

# Calculate retrained autogluon clock results and update pheno for all data sets

library(dplyr)
library(tidyverse)

dir.create("./7-output")

results <- list.files("./6-output", pattern = ".csv")

df <- read.csv(file= paste0("./6-output/",results[1])) %>% select(-X)

for (f in results[-1]){
  df2 <-read.csv(file= paste0("./6-output/",f)) %>% select(-X)
  df <- full_join(df,df2)
}

saveRDS(df, file = "./7-output/results_multi3C_AUC.Rds")

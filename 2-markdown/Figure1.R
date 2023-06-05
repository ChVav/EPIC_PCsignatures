# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(rcartocolor)

palette = carto_pal(8,"Safe")

# source helper functions
source("../1-src/functions/calcRepea.R") # helper function for estimating precision/repeatability
source("../1-src/functions/regresAge.R") # helper function for calculating root mean squared error

### a) performance Hannum original, Hannum_PC, Hannum PC_proxy in 450K blood data ##----
result_BloodRep <- readRDS("../1-src/analysis/1-Hannum-clocks-Higgins/2-output/pheno_BloodRep_450K_Hannum_PC.Rds")

# create data frame for plotting
df <- as.data.frame(matrix(data=NA, nrow=2,ncol=2))
colnames(df) <- c("clock","RMSE")
clocks <- c("Hannum.EUTOPS","Hannum_PC.EUTOPS")

# accuracy for each clock, regress and save RMSE
for(i in 1:length(clocks)){
  df[i,1] <- clocks[i]
  df[i,2] <- sqrt(mean(regresAge(result_BloodRep,clocks[i],"", mode="fitting")$residuals^2)) 
}

# precision for each clock, calculate ICC
result_BloodRep$sample_ID <- result_BloodRep$sample
result_BloodRep$sample_type <- "Blood"
precision <- calcRepea(result_BloodRep,clocks)
precision <- precision[[2]] %>% dplyr::rename(clock=index) %>% select(-sample_type)

df <- full_join(df, precision)

#plot
cols <- c(palette[1],palette[5])
df$clock <- factor(df$clock, levels = clocks)
labels <- c("Hannum", "Hannum_PC")

plota <- df %>%
  ggplot(aes(x=RMSE,y=ICC,color=clock, shape=clock)) +
  geom_point(size=4) +
  labs(color  = "", shape = "") +
  geom_errorbar(mapping=aes(ymin=CI0.95_lower, ymax=CI0.95_higher), width=0.2, size=0.2) + 
  xlim(3,4.5) +
  ylim(0.95,1) +
  scale_color_manual(labels=labels, values=cols)+
  scale_shape_manual(labels=labels, values=c(15,16))+
  theme_minimal()

# Original Hannum vs chronological age
fit <- lm(round(Hannum.EUTOPS, digits=4) ~ round(age, digits=4), data = result_BloodRep)
anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
               "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))

#CI for slope
#fit <- lm(Hannum.EUTOPS ~ age, data = result_BloodRep)
#confint(fit,'age',level=0.95)

plotb <- result_BloodRep %>%
  ggplot(aes(x = round(age, digits=4), y = round(Hannum.EUTOPS, digits=4))) +
  geom_point(color=cols[1], shape=15) +
  theme(panel.background = element_blank()) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              size = 0.5,
              se = FALSE,
              alpha = 0.3,
              colour = "gray40") +
  xlab("Age") +
  ylab("Hannum predicted age")  +
  xlim(35,85) +
  ylim(35,85) +
  annotate("text",
           x=50,
           y = 40,
           label = anno,
           hjust = 0,
           size = 3) 

# PC Age vs chronological age
fit <- lm(round(Hannum_PC.EUTOPS, digits=4) ~ round(age, digits=4), data = result_BloodRep)
anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
               "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))

#CI for slope
#fit <- lm(Hannum_PC.EUTOPS ~ age, data = result_BloodRep)
#confint(fit,'age',level=0.95)
#fit <- lm(Hannum_PCproxy.EUTOPS ~ age, data = result_BloodRep)
#confint(fit,'age',level=0.95)

plotc <- result_BloodRep %>%
  ggplot(aes(x = round(age, digits=4), y = round(Hannum_PC.EUTOPS, digits=4))) +
  geom_point(color=cols[2], shape=16) +
  theme(panel.background = element_blank(),plot.margin = margin(0, 0, 0, 0, "pt")) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              size = 0.5,
              se = FALSE,
              alpha = 0.3,
              colour = "gray40") +
  xlab("Age") +
  ylab("Hannum_PC predicted age")  +
  xlim(35,85) +
  ylim(35,85) +
  annotate("text",
           x=50,
           y = 40,
           label = anno,
           hjust = 0,
           size = 3) 


## combined plot ##----

des <- "
ABBCC
ABBCC"

plot <- plota +  plotb + plotc + plot_layout(design = des) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave(filename='Figure1.pdf', plot, width=9, height=2.5, dpi = 600)
ggsave(filename= 'Figure1.jpeg', plot, width=180, height=60, units= "mm", dpi=600)
